import os
import re
import sys
import glob
import json
import gzip
import zlib
import heapq
import queue
import struct
import locale
import tempfile
import threading
import functools
import contextlib
import webbrowser
from collections import namedtuple
from pyliftover import LiftOver

try:
    import tkinter as tk
    import tkinter.font as tkfont
    from tkinter import filedialog, messagebox, ttk
    from tkinterdnd2 import DND_FILES, TkinterDnD
    HAS_TK = True
except ImportError:
    HAS_TK = False

try:
    from pyfaidx import Fasta
    HAS_PYFAIDX = True
except ImportError:
    HAS_PYFAIDX = False


CHAIN_DOWNLOAD_URL = "https://hgdownload.soe.ucsc.edu/goldenPath"

__version__ = "1.0.2"


def detect_system_language():
    if sys.platform.startswith('win'):
        try:
            import ctypes
            lcid = ctypes.windll.kernel32.GetUserDefaultUILanguage()
            loc_name = locale.windows_locale.get(lcid)
            if loc_name:
                code = loc_name.split('_')[0].lower()
                if len(code) == 2:
                    return code
        except Exception:
            pass

    for env_var in ('LANGUAGE', 'LC_ALL', 'LC_MESSAGES', 'LANG'):
        val = os.environ.get(env_var)
        if val:
            first_entry = val.split(':')[0].split('.')[0]
            code = first_entry.split('_')[0].lower()
            if len(code) == 2:
                return code

    try:
        loc = locale.getlocale()[0]
        if loc:
            code = loc.split('_')[0].lower()
            if len(code) == 2:
                return code
    except Exception:
        pass

    return "en"

@functools.lru_cache(maxsize=None)
def chrom_sort_key(chrom):
    c = chrom.lower()
    if c.startswith('chr'):
        c = c[3:]
    if c.isdigit():
        return (0, int(c), "")
    special = {'x': 23, 'y': 24, 'm': 25, 'mt': 25}
    if c in special:
        return (1, special[c], "")
    return (2, 0, c)

REV_COMP_TRANS = str.maketrans("ACGTRYKMBDHVNacgtrykmbdhvn", "TGCAYRMKVHDBNtgcayrmkvhdbn")

_M_CONTIG_NAMES = ('chrM', 'CHRM', 'ChrM', 'M', 'm')
_MT_CONTIG_NAMES = ('chrMT', 'CHRMT', 'ChrMT', 'MT', 'mt')
M_NAMES_UPPER = ('MT', 'CHRM', 'M', 'CHRMT')

def _mt_kind(name):
    value = str(name).lower()
    if value.startswith('chr'):
        value = value[3:]
    return value if value in ('m', 'mt') else None

def _mt_alias_candidates(name, known_sizes=None):
    kind = _mt_kind(name)
    if kind is None:
        return []

    m_names = tuple(_M_CONTIG_NAMES)
    mt_names = tuple(_MT_CONTIG_NAMES)
    primary = m_names if kind == 'm' else mt_names
    other = mt_names if kind == 'm' else m_names

    def group_names(sizes, group):
        found = []
        for raw_name in (sizes or {}):
            if _mt_kind(raw_name) == group:
                found.append(raw_name)
        return found

    if known_sizes is None:
        explicit_name = str(name).lower().startswith('chr')
        ordered = list(primary) if explicit_name else list(primary) + list(other)
    else:
        present_m = group_names(known_sizes, 'm')
        present_mt = group_names(known_sizes, 'mt')
        explicit_chr_name = str(name).lower().startswith('chr')
        if present_m and present_mt:
            ordered = list(primary)
        else:
            ordered = list(primary)
            other_names = group_names(known_sizes, 'mt' if kind == 'm' else 'm')
            if explicit_chr_name:
                other_names = []
            if other_names:
                other_sizes = {int(known_sizes[n]) for n in other_names if known_sizes.get(n) is not None}
                if other_sizes and other_sizes == {16571} and kind == 'mt':
                    other_names = []
                if other_names:
                    ordered.extend(other)

    result = [name]
    result.extend(x for x in ordered if x.lower() != str(name).lower())
    return list(dict.fromkeys(result))


def _mt_target_candidates(name, sizes):
    """Кандидаты имён митохондрии в целевой FASTA / таблице длин.
    Пары chrM <-> chrMT (оба с префиксом chr) по-прежнему не смешиваются (legacy hg19 chrM = 16 571 п.н.),
    но контиги без префикса (MT / M — стиль Ensembl/NCBI, это всегда rCRS) допустимы для любого имени с chr."""
    candidates = _mt_alias_candidates(name, sizes)
    if str(name).lower().startswith('chr'):
        for alias in _M_CONTIG_NAMES + _MT_CONTIG_NAMES:
            if not alias.lower().startswith('chr') and alias not in candidates:
                candidates.append(alias)
    return candidates


def _fasta_contig_sizes(fasta_obj):
    sizes = {}
    try:
        for name in fasta_obj.keys():
            try:
                sizes[str(name)] = len(fasta_obj[name])
            except Exception:
                continue
    except Exception:
        return None
    return sizes


def _chain_source_sizes(lo):
    chain_file = getattr(lo, 'chain_file', None)
    chains = getattr(chain_file, 'chains', None)
    if chains is None:
        return None
    sizes = {}
    try:
        if isinstance(chains, dict):
            for chrom, chain_list in chains.items():
                if chain_list:
                    c = chain_list[0]
                    sizes.setdefault(str(getattr(c, 'source_name', chrom)), int(c.source_size))
        else:
            for chain in chains:
                sizes.setdefault(str(getattr(chain, 'source_name', '')), int(chain.source_size))
    except Exception:
        pass
    return sizes if sizes else None

def is_symbolic(allele):
    return allele.startswith('<') or allele in ('*', '.')

def is_breakend(allele):
    return '[' in allele or ']' in allele

def rev_comp_allele(allele):
    if is_symbolic(allele) or allele == '-':
        return allele
    return allele.translate(REV_COMP_TRANS)[::-1]


def get_fasta_seq(fasta_obj, chrom, start_0, end_0, target_lengths=None):
    try:
        if start_0 < 0 or end_0 < start_0:
            return None
        if _mt_kind(chrom):
            sizes = target_lengths if target_lengths is not None else _fasta_contig_sizes(fasta_obj)
            candidates = _mt_target_candidates(chrom, sizes)
        else:
            candidates = [chrom]
            chrom_lower = chrom.lower()
            if chrom_lower.startswith("chr"):
                suffix = chrom[3:]
                candidates.extend((suffix, "chr" + suffix, "CHR" + suffix, "Chr" + suffix))
            else:
                candidates.extend(("chr" + chrom, "CHR" + chrom, "Chr" + chrom))
        for candidate in dict.fromkeys(candidates):
            if candidate in fasta_obj:
                return str(fasta_obj[candidate][start_0:end_0]).upper()
    except (KeyError, IndexError, TypeError, ValueError):
        pass
    except Exception:
        pass
    return None

def add_info(info, tags):
    extra = ';'.join(tags)
    return extra if info in ('.', '') else info + ';' + extra

def _fmt_num(x):
    return format(x, '.6g')

def swap_genotypes_and_format(format_str, sample_str):
    if not format_str or not sample_str or sample_str == '.':
        return sample_str
    keys = format_str.split(':')
    vals = sample_str.split(':')

    def idx(name):
        i = keys.index(name) if name in keys else -1
        return i if 0 <= i < len(vals) else -1

    i = idx('GT')
    if i != -1:
        tokens = re.split(r'([/|])', vals[i])
        vals[i] = ''.join(('1' if t == '0' else '0' if t == '1' else t) for t in tokens)

    for name in ('PL', 'GL'):
        i = idx(name)
        if i != -1:
            items = vals[i].split(',')
            if len(items) in (2, 3):
                vals[i] = ','.join(reversed(items))

    for name in ('AD', 'ADF', 'ADR'):
        i = idx(name)
        if i != -1:
            items = vals[i].split(',')
            if len(items) == 2:
                vals[i] = ','.join(reversed(items))

    i = idx('AF')
    if i != -1 and ',' not in vals[i] and vals[i] not in ('.', ''):
        try:
            vals[i] = _fmt_num(1.0 - float(vals[i]))
        except ValueError:
            pass

    return ':'.join(vals)

def swap_info_fields(info):
    if info in ('.', ''):
        return info
    items = info.split(';')
    an = None
    for it in items:
        if it.startswith('AN='):
            try:
                an = int(it[3:])
            except ValueError:
                pass
    for n, it in enumerate(items):
        key, eq, val = it.partition('=')
        if not eq or ',' in val:
            continue
        try:
            if key in ('AF', 'MLEAF'):
                items[n] = f"{key}={_fmt_num(1.0 - float(val))}"
            elif key in ('AC', 'MLEAC') and an is not None:
                items[n] = f"{key}={an - int(val)}"
        except ValueError:
            pass
    return ';'.join(items)

_BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")

class BgzfWriter:
    BLOCK = 0xff00

    def __init__(self, path):
        self._f = open(path, 'wb')
        self._buf = bytearray()
        self._closed = False

    def _write_block(self, data):
        comp = zlib.compressobj(6, zlib.DEFLATED, -15)
        cdata = comp.compress(data) + comp.flush()
        bsize = len(cdata) + 25
        header = (b'\x1f\x8b\x08\x04' + b'\x00\x00\x00\x00' + b'\x00\xff' +
                  b'\x06\x00' + b'BC' + b'\x02\x00' + struct.pack('<H', bsize))
        footer = struct.pack('<II', zlib.crc32(data) & 0xffffffff, len(data) & 0xffffffff)
        self._f.write(header + cdata + footer)

    def write(self, text):
        if self._closed:
            raise ValueError("write to closed BGZF stream")
        if not isinstance(text, str):
            raise TypeError("BGZF text stream accepts str, not bytes")
        self._buf += text.encode('utf-8')
        while len(self._buf) >= self.BLOCK:
            self._write_block(bytes(self._buf[:self.BLOCK]))
            del self._buf[:self.BLOCK]

    def writelines(self, lines):
        for line in lines:
            self.write(line)

    def flush(self):
        if self._closed:
            raise ValueError("flush of closed BGZF stream")
        if self._buf:
            self._write_block(bytes(self._buf))
            self._buf.clear()
        self._f.flush()

    def close(self):
        if self._closed:
            return
        self._closed = True
        try:
            if self._buf:
                self._write_block(bytes(self._buf))
                self._buf.clear()
            self._f.write(_BGZF_EOF)
        finally:
            self._f.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

def _is_gzip_file(filepath):
    try:
        with open(filepath, 'rb') as f:
            return f.read(2) == b'\x1f\x8b'
    except OSError:
        return False

def open_file_read(filepath):
    if _is_gzip_file(filepath):
        return gzip.open(filepath, 'rt', encoding='utf-8-sig', errors='replace')
    return open(filepath, 'rt', encoding='utf-8-sig', errors='replace')

def open_file_write(filepath):
    if filepath.lower().endswith(('.gz', '.bgz')):
        return BgzfWriter(filepath)
    return open(filepath, 'w', encoding='utf-8', newline='\n')

def read_chain_target_sizes(chain_path):
    sizes = {}
    try:
        opener = gzip.open if _is_gzip_file(chain_path) else open
        with opener(chain_path, 'rb') as f:
            for line in f:
                if line.startswith(b'chain '):
                    p = line.split()
                    if len(p) >= 9:
                        sizes[p[7].decode('utf-8', 'replace')] = int(p[8])
    except Exception:
        pass
    return sizes

def lookup_length(lengths, name):
    if _mt_kind(name):
        candidates = _mt_target_candidates(name, lengths)
    else:
        candidates = [name, name[3:] if name.lower().startswith('chr') else 'chr' + name]
    for candidate in candidates:
        if candidate in lengths:
            return lengths[candidate]
        candidate_lower = candidate.lower()
        for existing, length in lengths.items():
            if existing.lower() == candidate_lower:
                return length
    return None

def detect_format(filename):
    name = filename.lower()
    for suf in ('.gz', '.bgz'):
        if name.endswith(suf):
            name = name[:-len(suf)]
            break
    if name.endswith(('.g.vcf', '.gvcf')): return 'gvcf'
    if name.endswith('.vcf'): return 'vcf'
    if name.endswith(('.bed', '.bedgraph', '.bdg')): return 'bed'
    if name.endswith(('.gff', '.gff3', '.gtf')): return 'gff'
    if name.endswith('.wig'): return 'wig'
    if name.endswith('.maf'): return 'maf'
    return None

def suggest_output_path(filename):
    name, comp = filename, ''
    for suf in ('.gz', '.bgz'):
        if name.lower().endswith(suf):
            comp = name[-len(suf):]
            name = name[:-len(suf)]
            break
    base, ext = os.path.splitext(name)
    return f"{base}_lifted{ext}{comp}"

def parse_dnd_path(root, event_data):
    try:
        items = root.tk.splitlist(event_data)
    except Exception:
        items = [event_data.strip()]
    return items[0] if items else ""

def find_preset_chain(filename):
    search_dirs = []
    if hasattr(sys, '_MEIPASS'):
        search_dirs.extend([os.path.join(sys._MEIPASS, "chains"), sys._MEIPASS])
    base_dir = os.path.dirname(os.path.abspath(sys.executable if getattr(sys, 'frozen', False) else __file__))
    search_dirs.extend([os.path.join(base_dir, "chains"), base_dir])
    for d in search_dirs:
        candidate = os.path.join(d, filename)
        if os.path.isfile(candidate):
            return candidate
    return None

class ExternalSorter:
    def __init__(self, max_bytes=96 * 1024 * 1024):
        self.max_bytes = max_bytes
        self.buf = []
        self.size = 0
        self.chunks = []

    @staticmethod
    def _key(r):
        return (chrom_sort_key(r[0]), r[1], r[2])

    def add(self, chrom, pos, end, payload):
        self.buf.append((chrom, pos, end, payload))
        self.size += len(payload) + len(chrom) + 80
        if self.size >= self.max_bytes:
            self._spill()

    def _spill(self):
        if not self.buf:
            return
        self.buf.sort(key=self._key)
        f = tempfile.TemporaryFile('w+', encoding='utf-8', newline='\n')
        f.writelines(f"{c}\t{p}\t{e}\t{pl}" for c, p, e, pl in self.buf)
        f.flush()
        f.seek(0)
        self.chunks.append(f)
        self.buf = []
        self.size = 0

    def items(self):
        if not self.chunks:
            self.buf.sort(key=self._key)
            yield from self.buf
            return
        self._spill()
        def gen(fh):
            for line in fh:
                c, p, e, pl = line.split('\t', 3)
                yield (c, int(p), int(e), pl)
        yield from heapq.merge(*[gen(f) for f in self.chunks], key=self._key)

    def close(self):
        for f in self.chunks:
            try: f.close()
            except Exception: pass
        self.chunks = []
        self.buf = []
        self.size = 0

class I18nManager:
    def __init__(self, default_lang=None):
        self.languages = {}
        self.translations = {}
        self.fallback_translations = {}
        self._load_available_languages()

        if "en" in self.languages:
            try:
                with open(self.languages["en"]["path"], 'r', encoding='utf-8-sig') as f:
                    self.fallback_translations = json.load(f)
            except Exception:
                pass

        if default_lang is None:
            sys_lang = detect_system_language()
            default_lang = sys_lang if sys_lang in self.languages else "en"

        self.current_lang = default_lang if default_lang in self.languages else "en"
        self.set_language(self.current_lang)

    def _get_locales_dirs(self):
        dirs = []
        if hasattr(sys, '_MEIPASS'):
            dirs.append(os.path.join(sys._MEIPASS, "locales"))
        base_dir = os.path.dirname(os.path.abspath(sys.executable if getattr(sys, 'frozen', False) else __file__))
        dirs.append(os.path.join(base_dir, "locales"))
        result, seen = [], set()
        for directory in dirs:
            normalized = os.path.normcase(os.path.abspath(directory))
            if normalized not in seen:
                seen.add(normalized)
                result.append(directory)
        return result

    def _load_available_languages(self):
        for loc_dir in reversed(self._get_locales_dirs()):
            if os.path.exists(loc_dir):
                for fpath in sorted(glob.glob(os.path.join(loc_dir, "*.json"))):
                    code = os.path.splitext(os.path.basename(fpath))[0].lower()
                    try:
                        with open(fpath, 'r', encoding='utf-8-sig') as f:
                            data = json.load(f)
                            if code not in self.languages:
                                self.languages[code] = {"name": data.get("lang_name", code.upper()), "path": fpath}
                    except Exception:
                        pass
        self.languages = dict(sorted(self.languages.items()))

    def set_language(self, code):
        if code in self.languages:
            self.current_lang = code
            try:
                with open(self.languages[code]["path"], 'r', encoding='utf-8-sig') as f:
                    self.translations = json.load(f)
            except Exception:
                self.translations = {}

    def t(self, key, **kwargs):
        text = self.translations.get(key, self.fallback_translations.get(key, key))
        if kwargs:
            try: text = text.format(**kwargs)
            except (KeyError, IndexError, ValueError): pass
        return text

Lifted = namedtuple('Lifted', 'chrom start0 end0 strand exact')
END_RE = re.compile(r'(?:^|;)END=(\d+)(?=;|$)')
END_SUB = re.compile(r'((?:^|;)END=)\d+(?=;|$)')
MAF_ALLELE_COLUMNS = ('Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
                      'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2')

class LiftOverEngine:
    def __init__(self, lo, fasta_obj, cmap, auto_chr, save_unmapped, log, t,
                 chain_file="", target_lengths=None, protected_paths=None):
        self.lo = lo
        self.chain_source_sizes = _chain_source_sizes(lo)
        self.fasta = fasta_obj
        self.strict = fasta_obj is not None
        self.cmap = cmap or {}
        self.auto_chr = auto_chr
        self.save_unm = save_unmapped
        self.log = log
        self.t = t
        self.chain_file = chain_file
        self.target_lengths = target_lengths
        self.protected_paths = list(protected_paths or [])
        self.f_unm = None
        self.f_rej = None
        self._reset_stats()

    def _reset_stats(self):
        self.total = self.mapped = self.unmapped = 0
        self.swapped = self.rejected = 0
        self.inverted = self.inverted_skipped = 0

    @property
    def stats(self):
        return dict(total=self.total, mapped=self.mapped, unmapped=self.unmapped, swapped=self.swapped,
                    rejected=self.rejected, inverted=self.inverted, inverted_skipped=self.inverted_skipped)

    def _lift_coord(self, chrom, pos_0):
        if pos_0 < 0: return None
        source_chrom = self.cmap.get(chrom, chrom)
        candidates = [source_chrom]
        if self.auto_chr:
            if _mt_kind(source_chrom):
                candidates.extend(_mt_alias_candidates(source_chrom, self.chain_source_sizes))
            else:
                source_lower = source_chrom.lower()
                if source_lower.startswith('chr'):
                    suffix = source_chrom[3:]
                    candidates.extend((suffix, 'chr' + suffix, 'CHR' + suffix, 'Chr' + suffix))
                else:
                    candidates.extend(('chr' + source_chrom, 'CHR' + source_chrom, 'Chr' + source_chrom))
        for candidate in dict.fromkeys(candidates):
            res = self.lo.convert_coordinate(candidate, pos_0)
            if res: return res[0]
        return None

    def _lift_interval(self, chrom, start0, end0):
        if start0 < 0 or end0 < start0: return None
        cs = self._lift_coord(chrom, start0)
        if not cs: return None
        ce = cs if end0 == start0 else self._lift_coord(chrom, end0)
        if not ce or cs[0] != ce[0] or cs[2] != ce[2]: return None
        s, e = (ce[1], cs[1]) if cs[2] == '-' else (cs[1], ce[1])
        if e < s: return None
        return Lifted(cs[0], s, e, cs[2], (e - s) == (end0 - start0))

    @staticmethod
    def _format_output_chrom(orig_chrom, target_chrom, explicit=False):
        if explicit: return target_chrom
        if _mt_kind(orig_chrom) or _mt_kind(target_chrom):
            if orig_chrom.lower().startswith('chr'):
                return target_chrom          # UCSC-стиль на входе: фактическое имя контига целевой сборки
            if _mt_kind(orig_chrom):
                return orig_chrom            # Ensembl/NCBI-стиль (MT, M): как и у остальных хромосом, имя входа сохраняется
            return target_chrom[3:] if target_chrom.lower().startswith('chr') else target_chrom
        orig_has_chr = orig_chrom.lower().startswith('chr')
        target_has_chr = target_chrom.lower().startswith('chr')
        if not orig_has_chr and target_has_chr:
            return target_chrom[3:]
        elif orig_has_chr and not target_has_chr:
            return 'chr' + target_chrom
        return target_chrom

    def _tick(self):
        self.total += 1
        if self.total % 25000 == 0:
            self.log(self.t("log_progress", total=self.total, mapped=self.mapped, unmapped=self.unmapped))

    def _unmapped(self, line):
        self.unmapped += 1
        if self.f_unm:
            self.f_unm.write(line if line.endswith('\n') else line + '\n')
        return None

    def _reject(self, parts, reason):
        self.rejected += 1
        if self.f_rej:
            p = list(parts)
            p[7] = add_info(p[7], ["LIFTOVER_REJECT_REASON=" + reason])
            self.f_rej.write('\t'.join(p) + '\n')
        return None

    def _write_sorted(self, f_out, sorter):
        if self.mapped:
            self.log(self.t("log_sorting", count=self.mapped))
        for _c, _p, _e, payload in sorter.items():
            f_out.write(payload)

    def run(self, cmd_type, in_file, out_file):
        self._reset_stats()
        protected = [self.chain_file, *self.protected_paths]
        if self.fasta is not None:
            fasta_path = getattr(self.fasta, "filename", None)
            if isinstance(fasta_path, (str, os.PathLike)):
                protected.append(os.fspath(fasta_path))
        path_error = check_paths(
            in_file, out_file, cmd_type=cmd_type,
            save_unmapped=self.save_unm, strict=self.strict,
            protected_paths=protected,
        )
        if path_error:
            raise ValueError(self.t(path_error))

        if self.target_lengths is None:
            if self.strict:
                try: self.target_lengths = {str(k): len(self.fasta[k]) for k in self.fasta.keys()}
                except Exception: self.target_lengths = {}
            else:
                self.target_lengths = read_chain_target_sizes(self.chain_file) if self.chain_file else {}

        unm_file = out_file + ".unmapped" if self.save_unm else None
        rej_file = (out_file + ".rejected.vcf") if cmd_type in ('vcf', 'gvcf') else (out_file + ".rejected")
        if not self.strict: rej_file = None

        handlers = {'vcf': self._run_vcf, 'gvcf': self._run_vcf, 'bed': self._run_bed,
                    'gff': self._run_gff, 'wig': self._run_wig, 'maf': self._run_maf}
        opened_outputs = []
        try:
            with contextlib.ExitStack() as stack:
                f_in = stack.enter_context(open_file_read(in_file))
                f_out = stack.enter_context(open_file_write(out_file))
                opened_outputs.append(out_file)
                if unm_file:
                    self.f_unm = stack.enter_context(open_file_write(unm_file))
                    opened_outputs.append(unm_file)
                if rej_file:
                    self.f_rej = stack.enter_context(open_file_write(rej_file))
                    opened_outputs.append(rej_file)
                handlers[cmd_type](f_in, f_out)
        except BaseException:
            for path in opened_outputs:
                if os.path.exists(path):
                    try: os.remove(path)
                    except OSError: pass
            raise
        finally:
            self.f_unm = self.f_rej = None

        if unm_file and self.unmapped == 0 and os.path.exists(unm_file): os.remove(unm_file)
        if rej_file and self.rejected == 0 and os.path.exists(rej_file): os.remove(rej_file)
        return self.stats

    def _run_vcf(self, f_in, f_out):
        sorter = ExternalSorter()
        header, chrom_line, seen = [], None, set()
        try:
            for line in f_in:
                if line.startswith('##'):
                    header.append(line if line.endswith('\n') else line + '\n')
                    if self.f_rej: self.f_rej.write(header[-1])
                    continue
                if line.startswith('#'):
                    chrom_line = line if line.endswith('\n') else line + '\n'
                    if self.f_rej:
                        self.f_rej.write('##INFO=<ID=LIFTOVER_REJECT_REASON,Number=1,Type=String,Description="Reason for liftover rejection">\n')
                        self.f_rej.write(chrom_line)
                    continue
                if not line.strip(): continue
                self._tick()
                rec = self._convert_vcf_record(line)
                if rec:
                    chrom, pos, out_line = rec
                    seen.add(chrom)
                    sorter.add(chrom, pos, pos, out_line)
                    self.mapped += 1

            out_header = [
                h for h in header
                if h.split('=', 1)[0] not in ('##contig', '##reference', '##assembly')
            ]
            existing = set(re.findall(r'^##INFO=<ID=([^,>]+)', ''.join(out_header), re.M))
            new_info = []
            if self.strict and 'LIFTOVER_SWAPPED' not in existing:
                new_info.append('##INFO=<ID=LIFTOVER_SWAPPED,Number=0,Type=Flag,Description="Reference and alternate alleles were swapped during liftover to match target genome reference">\n')
            if 'LIFTOVER_REVERSE_COMPLEMENTED' not in existing:
                new_info.append('##INFO=<ID=LIFTOVER_REVERSE_COMPLEMENTED,Number=0,Type=Flag,Description="Alleles were reverse-complemented during liftover due to negative strand alignment">\n')
            out_header.extend(new_info)
            chain_name = os.path.basename(self.chain_file) if self.chain_file else "unknown"
            out_header.append(f"##LiftOverZ_chain={chain_name}\n")
            for c in sorted(seen, key=chrom_sort_key):
                length = lookup_length(self.target_lengths or {}, c)
                out_header.append(f"##contig=<ID={c},length={length}>\n" if length else f"##contig=<ID={c}>\n")
            for h in out_header: f_out.write(h)
            if chrom_line: f_out.write(chrom_line)
            self._write_sorted(f_out, sorter)
        finally:
            sorter.close()

    def _normalize_left_anchor(self, chrom, pos1, ref, alts):
        als = [ref] + alts
        while all(len(a) > 0 for a in als) and len({a[-1] for a in als}) == 1:
            als = [a[:-1] for a in als]
        if any(len(a) == 0 for a in als):
            if pos1 <= 1: return None
            base = get_fasta_seq(self.fasta, chrom, pos1 - 2, pos1 - 1, self.target_lengths)
            if not base: return None
            als = [base + a for a in als]
            pos1 -= 1
        return pos1, als[0], als[1:]

    def _convert_vcf_record(self, line):
        parts = line.rstrip('\r\n').split('\t')
        if len(parts) < 5: return self._unmapped(line)
        try: pos = int(parts[1])
        except ValueError: return self._unmapped(line)
        while len(parts) < 8: parts.append('.')
        orig = list(parts)
        chrom = parts[0]
        ref = parts[3].upper()
        alts = parts[4].upper().split(',')

        m = END_RE.search(parts[7])
        end1 = int(m.group(1)) if m else None
        has_end = end1 is not None and end1 >= pos
        span_end1 = end1 if has_end else pos + len(ref) - 1

        lifted = self._lift_interval(chrom, pos - 1, span_end1 - 1)
        if lifted is None: return self._unmapped(line)
        if not has_end and not lifted.exact: return self._unmapped(line)

        neg = lifted.strand == '-'
        target_c = self._format_output_chrom(chrom, lifted.chrom, chrom in self.cmap)
        new_pos = lifted.start0 + 1
        new_end = lifted.end0 + 1
        tags = []
        same_len = True

        if neg:
            if any(is_breakend(a) for a in alts): return self._unmapped(line)
            cur_alts = [rev_comp_allele(a) for a in alts]
            if has_end:
                cur_ref = None
            else:
                same_len = all(is_symbolic(a) or len(a) == len(ref) for a in alts)
                if not same_len and not self.strict:
                    self.inverted_skipped += 1
                    return self._unmapped(line)
                cur_ref = rev_comp_allele(ref)
            tags.append('LIFTOVER_REVERSE_COMPLEMENTED')
        else:
            cur_ref, cur_alts = ref, alts

        swapped_now = False

        if not self.strict:
            if neg:
                parts[3] = cur_ref if cur_ref is not None else 'N'
                parts[4] = ','.join(cur_alts)
        else:
            if has_end and neg:
                base = get_fasta_seq(self.fasta, lifted.chrom, lifted.start0, lifted.start0 + 1, self.target_lengths)
                if base is None: return self._reject(orig, "ContigNotFoundInFasta")
                cur_ref = base
            else:
                fasta_ref = get_fasta_seq(self.fasta, lifted.chrom, lifted.start0, lifted.start0 + len(cur_ref), self.target_lengths)
                if fasta_ref is None: return self._reject(orig, "ContigNotFoundInFasta")
                if cur_ref == fasta_ref:
                    if neg and not same_len and not any(is_symbolic(a) for a in cur_alts):
                        norm = self._normalize_left_anchor(lifted.chrom, new_pos, cur_ref, cur_alts)
                        if norm is None: return self._reject(orig, "CannotAnchorInvertedIndel")
                        new_pos, cur_ref, cur_alts = norm
                elif len(cur_alts) == 1 and cur_alts[0] == fasta_ref and len(cur_ref) == len(fasta_ref):
                    cur_ref, cur_alts = cur_alts[0], [cur_ref]
                    swapped_now = True
                    tags.insert(0, 'LIFTOVER_SWAPPED')
                else:
                    return self._reject(orig, "MismatchedRefAllele")
            parts[3] = cur_ref
            parts[4] = ','.join(cur_alts)

        if swapped_now:
            if len(parts) >= 10:
                fmt = parts[8]
                for s_idx in range(9, len(parts)):
                    parts[s_idx] = swap_genotypes_and_format(fmt, parts[s_idx])
            parts[7] = swap_info_fields(parts[7])
            self.swapped += 1

        if has_end: parts[7] = END_SUB.sub(lambda mm: mm.group(1) + str(new_end), parts[7])
        if tags: parts[7] = add_info(parts[7], tags)
        if neg: self.inverted += 1

        parts[0] = target_c
        parts[1] = str(new_pos)
        return (target_c, new_pos, '\t'.join(parts) + '\n')

    def _lift_thick(self, chrom, ts, te, start, end, new_s, new_e, neg, exact):
        if ts == te: return new_s, new_s
        if exact:
            if neg: return new_s + (end - te), new_s + (end - ts)
            return new_s + (ts - start), new_s + (te - start)
        a = self._lift_interval(chrom, ts, te - 1)
        if a: return a.start0, a.end0 + 1
        return new_s, new_e

    def _run_bed(self, f_in, f_out):
        sorter, header = ExternalSorter(), []
        try:
            for line in f_in:
                if line.startswith(('#', 'track', 'browser')):
                    header.append(line if line.endswith('\n') else line + '\n')
                    continue
                if not line.strip(): continue
                self._tick()
                rec = self._convert_bed_line(line)
                if rec:
                    c, s, e, out = rec
                    sorter.add(c, s, e, out)
                    self.mapped += 1
            f_out.writelines(header)
            self._write_sorted(f_out, sorter)
        finally:
            sorter.close()

    def _convert_bed_line(self, line):
        parts = line.rstrip('\r\n').split('\t')
        if len(parts) < 3: return self._unmapped(line)
        try:
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if end < start: return self._unmapped(line)
            zero = end == start
            lifted = self._lift_interval(chrom, start, start if zero else end - 1)
            if lifted is None: return self._unmapped(line)
            new_s = lifted.start0
            new_e = new_s if zero else lifted.end0 + 1
            neg = lifted.strand == '-'
            exact = (new_e - new_s) == (end - start)
            target_c = self._format_output_chrom(chrom, lifted.chrom, chrom in self.cmap)

            if neg and len(parts) >= 6 and parts[5] in ('+', '-'):
                parts[5] = '-' if parts[5] == '+' else '+'
            if len(parts) >= 8:
                ts, te = self._lift_thick(chrom, int(parts[6]), int(parts[7]), start, end, new_s, new_e, neg, exact)
                parts[6], parts[7] = str(ts), str(te)
            if len(parts) >= 12:
                if not exact: return self._unmapped(line)
                trailing = ',' if parts[10].endswith(',') else ''
                sizes = [int(x) for x in parts[10].strip(',').split(',')]
                starts = [int(x) for x in parts[11].strip(',').split(',')]
                blocks = list(zip(starts, sizes))
                if neg: blocks = sorted(((end - start) - (st + sz), sz) for st, sz in blocks)
                parts[10] = ','.join(str(sz) for _, sz in blocks) + trailing
                parts[11] = ','.join(str(st) for st, _ in blocks) + trailing
        except ValueError:
            return self._unmapped(line)

        if neg: self.inverted += 1
        parts[0], parts[1], parts[2] = target_c, str(new_s), str(new_e)
        return (target_c, new_s, new_e, '\t'.join(parts) + '\n')

    def _run_gff(self, f_in, f_out):
        sorter, header, trailer, in_fasta = ExternalSorter(), [], [], False
        try:
            for line in f_in:
                if in_fasta:
                    trailer.append(line if line.endswith('\n') else line + '\n')
                    continue
                if line.startswith('##FASTA'):
                    in_fasta = True
                    trailer.append(line if line.endswith('\n') else line + '\n')
                    continue
                if line.startswith('#'):
                    header.append(line if line.endswith('\n') else line + '\n')
                    continue
                if not line.strip(): continue
                self._tick()
                rec = self._convert_gff_line(line)
                if rec:
                    c, s, e, out = rec
                    sorter.add(c, s, e, out)
                    self.mapped += 1
            f_out.writelines(header)
            self._write_sorted(f_out, sorter)
            f_out.writelines(trailer)
        finally:
            sorter.close()

    def _convert_gff_line(self, line):
        parts = line.rstrip('\r\n').split('\t')
        if len(parts) < 5: return self._unmapped(line)
        try: chrom, s1, e1 = parts[0], int(parts[3]), int(parts[4])
        except ValueError: return self._unmapped(line)
        lifted = self._lift_interval(chrom, s1 - 1, e1 - 1)
        if lifted is None: return self._unmapped(line)
        target_c = self._format_output_chrom(chrom, lifted.chrom, chrom in self.cmap)
        new_s, new_e = lifted.start0 + 1, lifted.end0 + 1
        if lifted.strand == '-':
            self.inverted += 1
            if len(parts) > 6 and parts[6] in ('+', '-'):
                parts[6] = '-' if parts[6] == '+' else '+'
        parts[0], parts[3], parts[4] = target_c, str(new_s), str(new_e)
        return (target_c, new_s, new_e, '\t'.join(parts) + '\n')

    def _run_wig(self, f_in, f_out):
        state = {'sorter': ExternalSorter()}
        chrom = mode = None
        start = step = span = 1
        idx = 0

        def flush():
            last = None
            for c, p, _e, payload in state['sorter'].items():
                sp, val = payload.rstrip('\n').split('\t', 1)
                if (c, sp) != last:
                    f_out.write(f"variableStep chrom={c}" + (f" span={sp}" if sp != '1' else '') + "\n")
                    last = (c, sp)
                f_out.write(f"{p + 1}\t{val}\n")
            state['sorter'].close()
            state['sorter'] = ExternalSorter()

        try:
            for line in f_in:
                sline = line.strip()
                if not sline: continue
                if sline.startswith(('#', 'browser')):
                    f_out.write(sline + '\n')
                    continue
                if sline.startswith('track'):
                    flush()
                    f_out.write(sline + '\n')
                    chrom = mode = None
                    continue
                if sline.startswith(('variableStep', 'fixedStep')):
                    tokens = sline.split()
                    mode = 'variable' if tokens[0] == 'variableStep' else 'fixed'
                    kv = dict(t.split('=', 1) for t in tokens[1:] if '=' in t)
                    chrom = kv.get('chrom')
                    try:
                        span, start, step = int(kv.get('span', 1)), int(kv.get('start', 1)), int(kv.get('step', 1))
                    except ValueError:
                        span, start, step = 1, 1, 1
                    idx = 0
                    continue

                self._tick()
                parts = sline.split()
                try:
                    if mode == 'fixed':
                        pos1, val = start + idx * step, parts[0]
                        idx += 1
                    else:
                        pos1, val = int(parts[0]), parts[1]
                    if chrom is None or mode is None: raise ValueError
                except (ValueError, IndexError):
                    self._unmapped(f"{chrom}\t{sline}\n")
                    continue

                if pos1 < 1 or span < 1:
                    self._unmapped(f"{chrom}\t{pos1}\t{val}\n")
                    continue
                lifted = self._lift_interval(chrom, pos1 - 1, pos1 - 1 + span - 1)
                if lifted is None or not lifted.exact or (lifted.end0 - lifted.start0 + 1) != span:
                    self._unmapped(f"{chrom}\t{pos1}\t{val}\n")
                    continue
                target_c = self._format_output_chrom(chrom, lifted.chrom, chrom in self.cmap)
                if lifted.strand == '-': self.inverted += 1
                state['sorter'].add(target_c, lifted.start0, lifted.start0, f"{span}\t{val}\n")
                self.mapped += 1
            if self.mapped:
                self.log(self.t("log_sorting", count=self.mapped))
            flush()
        finally:
            state['sorter'].close()

    def _run_maf(self, f_in, f_out):
        sorter, header = ExternalSorter(), []
        chr_idx, start_idx, end_idx = 4, 5, 6
        allele_idx = []
        header_line = None
        try:
            for line in f_in:
                if line.startswith('#'):
                    header.append(line if line.endswith('\n') else line + '\n')
                    continue
                if not line.strip(): continue
                parts = line.rstrip('\r\n').split('\t')
                if header_line is None:
                    header_line = line if line.endswith('\n') else line + '\n'
                    if 'Chromosome' in parts: chr_idx = parts.index('Chromosome')
                    if 'Start_Position' in parts: start_idx = parts.index('Start_Position')
                    if 'End_Position' in parts: end_idx = parts.index('End_Position')
                    allele_idx = [parts.index(n) for n in MAF_ALLELE_COLUMNS if n in parts]
                    continue
                self._tick()
                rec = self._convert_maf_line(line, parts, chr_idx, start_idx, end_idx, allele_idx)
                if rec:
                    c, s, e, out = rec
                    sorter.add(c, s, e, out)
                    self.mapped += 1
            f_out.writelines(header)
            if header_line: f_out.write(header_line)
            self._write_sorted(f_out, sorter)
        finally:
            sorter.close()

    def _convert_maf_line(self, line, parts, chr_idx, start_idx, end_idx, allele_idx):
        if len(parts) <= max(chr_idx, start_idx, end_idx): return self._unmapped(line)
        try:
            chrom = parts[chr_idx]
            s1, e1 = int(parts[start_idx]), int(parts[end_idx])
        except ValueError: return self._unmapped(line)
        if e1 < s1: return self._unmapped(line)
        lifted = self._lift_interval(chrom, s1 - 1, e1 - 1)
        if lifted is None: return self._unmapped(line)
        target_c = self._format_output_chrom(chrom, lifted.chrom, chrom in self.cmap)
        new_s, new_e = lifted.start0 + 1, lifted.end0 + 1
        if lifted.strand == '-':
            self.inverted += 1
            for i in allele_idx:
                if i < len(parts): parts[i] = rev_comp_allele(parts[i])
        parts[chr_idx], parts[start_idx], parts[end_idx] = target_c, str(new_s), str(new_e)
        return (target_c, new_s, new_e, '\t'.join(parts) + '\n')


def load_chrom_map(path):
    cmap = {}
    if not path:
        return cmap
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Chromosome mapping file does not exist or is not a file: {path}")
    with open(path, 'r', encoding='utf-8-sig', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                cmap[parts[0]] = parts[1]
    return cmap


def _same_path(path_a, path_b):
    if not path_a or not path_b:
        return False
    a = os.path.normcase(os.path.realpath(os.path.abspath(os.fspath(path_a))))
    b = os.path.normcase(os.path.realpath(os.path.abspath(os.fspath(path_b))))
    if a == b:
        return True
    if os.path.exists(path_a) and os.path.exists(path_b):
        try:
            return os.path.samefile(path_a, path_b)
        except OSError:
            pass
    return False


def _output_artifact_paths(out, cmd_type=None, save_unmapped=False, strict=False):
    paths = [out]
    if save_unmapped:
        paths.append(out + ".unmapped")
    if strict:
        paths.append(out + (".rejected.vcf" if cmd_type in ("vcf", "gvcf") else ".rejected"))
    return paths


def check_paths(inp, out, cmd_type=None, save_unmapped=False, strict=False, protected_paths=None):
    if not inp or not os.path.isfile(inp):
        return "err_no_input"
    if not out:
        return "err_no_output"
    if not os.path.isdir(os.path.dirname(os.path.abspath(out))):
        return "err_bad_output_dir"

    output_paths = _output_artifact_paths(out, cmd_type, save_unmapped, strict)
    if _same_path(inp, out):
        return "err_same_file"

    protected = [inp]
    if protected_paths:
        protected.extend(p for p in protected_paths if p)
    for output_path in output_paths:
        for protected_path in protected:
            if _same_path(output_path, protected_path):
                return "err_output_conflict"
    for i, first in enumerate(output_paths):
        if any(_same_path(first, later) for later in output_paths[i + 1:]):
            return "err_output_conflict"
    return None


class LiftOverZApp:
    FORMAT_MAPPING = {
        "VCF (.vcf, .vcf.gz)": "vcf",
        "gVCF (.gvcf, .gvcf.gz)": "gvcf",
        "BED / BEDGraph (.bed, .bedgraph)": "bed",
        "GFF / GTF (.gff, .gtf, .gff3)": "gff",
        "Wiggle (.wig)": "wig",
        "MAF (.maf)": "maf"
    }

    CHAIN_PRESETS = [
        ("preset_custom", None),
        ("preset_hg19_to_hg38", "hg19ToHg38.over.chain.gz"),
        ("preset_hg38_to_hg19", "hg38ToHg19.over.chain.gz"),
        ("preset_hg38_to_hs1", "hg38ToHs1.over.chain.gz"),
        ("preset_hs1_to_hg38", "hs1ToHg38.over.chain.gz"),
        ("preset_hg19_to_hs1", "hg19ToHs1.over.chain.gz"),
        ("preset_hs1_to_hg19", "hs1ToHg19.over.chain.gz"),
        ("preset_mm10_to_mm39", "mm10ToMm39.over.chain.gz"),
        ("preset_mm39_to_mm10", "mm39ToMm10.over.chain.gz"),
        ("preset_mm9_to_mm10", "mm9ToMm10.over.chain.gz")
    ]

    GENOMIC_PATTERNS = (
        "*.vcf *.vcf.gz *.vcf.bgz *.g.vcf *.g.vcf.gz *.gvcf *.gvcf.gz "
        "*.bed *.bed.gz *.bedgraph *.bedgraph.gz *.gff *.gff.gz *.gff3 *.gff3.gz "
        "*.gtf *.gtf.gz *.wig *.wig.gz *.maf *.maf.gz"
    )

    def __init__(self, root):
        self.root = root
        self.i18n = I18nManager()
        self.ui_queue = queue.Queue()

        self.root.geometry("780x760")
        self.root.minsize(720, 680)

        self.input_path = tk.StringVar()
        self.chain_path = tk.StringVar()
        self.target_fasta_path = tk.StringVar()
        self.out_path = tk.StringVar()
        self.chrom_map_path = tk.StringVar()
        self.selected_format = tk.StringVar(value=list(self.FORMAT_MAPPING.keys())[0])

        self.auto_chr = tk.BooleanVar(value=True)
        self.save_unmapped = tk.BooleanVar(value=True)

        self._build_ui()
        self._update_ui_texts()
        self.cb_presets.current(1)
        self._on_preset_selected()
        self.root.after(100, self._poll_ui_queue)

    def _build_ui(self):
        self.main_frame = ttk.Frame(self.root, padding=12)
        self.main_frame.pack(fill=tk.BOTH, expand=True)

        top_bar = ttk.Frame(self.main_frame)
        top_bar.pack(fill=tk.X, pady=(0, 6))

        self.lbl_lang = ttk.Label(top_bar)
        self.lbl_lang.pack(side=tk.LEFT, padx=(0, 5))

        lang_options = [f"{info['name']} ({code})" for code, info in self.i18n.languages.items()]
        self.lang_cb = ttk.Combobox(top_bar, values=lang_options, state="readonly", width=18)
        current_name = self.i18n.languages.get(self.i18n.current_lang, {}).get("name", "English")
        self.lang_cb.set(f"{current_name} ({self.i18n.current_lang})")
        self.lang_cb.pack(side=tk.LEFT)
        self.lang_cb.bind("<<ComboboxSelected>>", self._on_language_changed)

        self.fmt_frame = ttk.LabelFrame(self.main_frame, padding=8)
        self.fmt_frame.pack(fill=tk.X, pady=(0, 8))

        self.cb_fmt = ttk.Combobox(
            self.fmt_frame,
            textvariable=self.selected_format,
            values=list(self.FORMAT_MAPPING.keys()),
            state="readonly"
        )
        self.cb_fmt.pack(fill=tk.X, expand=True)

        self.files_frame = ttk.LabelFrame(self.main_frame, padding=8)
        self.files_frame.pack(fill=tk.X, pady=(0, 8))
        self.files_frame.columnconfigure(1, weight=1)

        self.lbl_input = ttk.Label(self.files_frame)
        self.lbl_input.grid(row=0, column=0, sticky="w", pady=3)
        self.ent_inp = ttk.Entry(self.files_frame, textvariable=self.input_path)
        self.ent_inp.grid(row=0, column=1, sticky="ew", padx=5, pady=3)
        self.btn_browse_inp = ttk.Button(self.files_frame, command=self._browse_input)
        self.btn_browse_inp.grid(row=0, column=2, pady=3)
        self.ent_inp.drop_target_register(DND_FILES)
        self.ent_inp.dnd_bind('<<Drop>>', self._on_drop_input)

        self.lbl_chain_preset = ttk.Label(self.files_frame)
        self.lbl_chain_preset.grid(row=1, column=0, sticky="w", pady=3)
        self.cb_presets = ttk.Combobox(self.files_frame, state="readonly")
        self.cb_presets.grid(row=1, column=1, columnspan=2, sticky="ew", padx=5, pady=3)
        self.cb_presets.bind("<<ComboboxSelected>>", self._on_preset_selected)

        self.lbl_chain = ttk.Label(self.files_frame)
        self.lbl_chain.grid(row=2, column=0, sticky="w", pady=3)
        self.ent_chain = ttk.Entry(self.files_frame, textvariable=self.chain_path)
        self.ent_chain.grid(row=2, column=1, sticky="ew", padx=5, pady=3)
        self.btn_browse_chain = ttk.Button(self.files_frame, command=self._browse_chain)
        self.btn_browse_chain.grid(row=2, column=2, pady=3)
        self.ent_chain.drop_target_register(DND_FILES)
        self.ent_chain.dnd_bind('<<Drop>>', self._on_drop_chain)

        link_font = tkfont.nametofont("TkDefaultFont").copy()
        link_font.configure(underline=True)
        self._link_font = link_font
        self.lnk_chain_dl = ttk.Label(self.files_frame, foreground="#1a73e8", cursor="hand2", font=link_font)
        self.lnk_chain_dl.grid(row=3, column=1, columnspan=2, sticky="w", padx=5, pady=(0, 3))
        self.lnk_chain_dl.bind("<Button-1>", self._open_chain_download_page)

        self.lbl_target_fasta = ttk.Label(self.files_frame)
        self.lbl_target_fasta.grid(row=4, column=0, sticky="w", pady=3)
        self.ent_fasta = ttk.Entry(self.files_frame, textvariable=self.target_fasta_path)
        self.ent_fasta.grid(row=4, column=1, sticky="ew", padx=5, pady=3)
        self.btn_browse_fasta = ttk.Button(self.files_frame, command=self._browse_fasta)
        self.btn_browse_fasta.grid(row=4, column=2, pady=3)
        self.ent_fasta.drop_target_register(DND_FILES)
        self.ent_fasta.dnd_bind('<<Drop>>', lambda e: self.target_fasta_path.set(parse_dnd_path(self.root, e.data)))

        self.lbl_output = ttk.Label(self.files_frame)
        self.lbl_output.grid(row=5, column=0, sticky="w", pady=3)
        self.ent_out = ttk.Entry(self.files_frame, textvariable=self.out_path)
        self.ent_out.grid(row=5, column=1, sticky="ew", padx=5, pady=3)
        self.btn_save_out = ttk.Button(self.files_frame, command=self._save_out_dialog)
        self.btn_save_out.grid(row=5, column=2, pady=3)

        self.opts_frame = ttk.LabelFrame(self.main_frame, padding=8)
        self.opts_frame.pack(fill=tk.X, pady=(0, 8))
        self.opts_frame.columnconfigure(1, weight=1)

        self.lbl_chrom_map = ttk.Label(self.opts_frame)
        self.lbl_chrom_map.grid(row=0, column=0, sticky="w", pady=3)
        self.ent_cmap = ttk.Entry(self.opts_frame, textvariable=self.chrom_map_path)
        self.ent_cmap.grid(row=0, column=1, sticky="ew", padx=5, pady=3)
        self.btn_browse_map = ttk.Button(self.opts_frame, command=self._browse_chrom_map)
        self.btn_browse_map.grid(row=0, column=2, pady=3)
        self.ent_cmap.drop_target_register(DND_FILES)
        self.ent_cmap.dnd_bind('<<Drop>>', lambda e: self.chrom_map_path.set(parse_dnd_path(self.root, e.data)))

        self.chk_auto_chr_w = ttk.Checkbutton(self.opts_frame, variable=self.auto_chr)
        self.chk_auto_chr_w.grid(row=1, column=0, columnspan=3, sticky="w", pady=2)

        self.chk_save_unm_w = ttk.Checkbutton(self.opts_frame, variable=self.save_unmapped)
        self.chk_save_unm_w.grid(row=2, column=0, columnspan=3, sticky="w", pady=2)

        btn_action_frame = ttk.Frame(self.main_frame)
        btn_action_frame.pack(fill=tk.X, pady=5)
        btn_action_frame.columnconfigure(0, weight=3)
        btn_action_frame.columnconfigure(1, weight=1)

        self.btn_run = ttk.Button(btn_action_frame, command=self._start_processing)
        self.btn_run.grid(row=0, column=0, sticky="ew", padx=(0, 5), ipady=4)

        self.btn_reset = ttk.Button(btn_action_frame, command=self._reset_all)
        self.btn_reset.grid(row=0, column=1, sticky="ew", ipady=4)

        self.log_frame = ttk.LabelFrame(self.main_frame, padding=5)
        self.log_frame.pack(fill=tk.BOTH, expand=True, pady=(5, 0))

        self.log_area = tk.Text(self.log_frame, height=9, state='disabled', wrap='word')
        self.log_area.pack(fill=tk.BOTH, expand=True)

    def _update_ui_texts(self):
        self.root.title(f'{self.i18n.t("app_title")}  v{__version__}')
        self.lbl_lang.config(text=self.i18n.t("language_label"))
        self.fmt_frame.config(text=self.i18n.t("step_1_format"))
        self.files_frame.config(text=self.i18n.t("step_2_files"))
        self.opts_frame.config(text=self.i18n.t("step_3_options"))
        self.log_frame.config(text=self.i18n.t("log_title"))

        self.lbl_input.config(text=self.i18n.t("lbl_input"))
        self.lbl_chain_preset.config(text=self.i18n.t("lbl_chain_preset"))
        self.lbl_chain.config(text=self.i18n.t("lbl_chain"))
        self.lbl_target_fasta.config(text=self.i18n.t("lbl_target_fasta"))
        self.lnk_chain_dl.config(text=self.i18n.t("lnk_chain_download"))
        self.lbl_output.config(text=self.i18n.t("lbl_output"))
        self.lbl_chrom_map.config(text=self.i18n.t("lbl_chrom_map"))
        self.chk_auto_chr_w.config(text=self.i18n.t("chk_auto_chr"))
        self.chk_save_unm_w.config(text=self.i18n.t("chk_save_unmapped"))

        self.btn_browse_inp.config(text=self.i18n.t("btn_browse"))
        self.btn_browse_chain.config(text=self.i18n.t("btn_browse"))
        self.btn_browse_fasta.config(text=self.i18n.t("btn_browse"))
        self.btn_browse_map.config(text=self.i18n.t("btn_browse"))
        self.btn_save_out.config(text=self.i18n.t("btn_save_as"))
        self.btn_run.config(text=self.i18n.t("btn_run"))
        self.btn_reset.config(text=self.i18n.t("btn_reset"))

        curr_idx = self.cb_presets.current()
        preset_labels = [self.i18n.t(item[0]) for item in self.CHAIN_PRESETS]
        self.cb_presets.config(values=preset_labels)
        if curr_idx >= 0:
            self.cb_presets.current(curr_idx)

    def _reset_all(self):
        self.input_path.set("")
        self.chain_path.set("")
        self.target_fasta_path.set("")
        self.out_path.set("")
        self.chrom_map_path.set("")
        self.selected_format.set(list(self.FORMAT_MAPPING.keys())[0])
        self.cb_presets.current(1)
        self._on_preset_selected()
        self.auto_chr.set(True)
        self.save_unmapped.set(True)

        self.log_area.config(state='normal')
        self.log_area.delete('1.0', tk.END)
        self.log_area.config(state='disabled')

    def _on_language_changed(self, event=None):
        code = self.lang_cb.get().split('(')[-1].rstrip(')').strip()
        self.i18n.set_language(code)
        self._update_ui_texts()

    def _on_preset_selected(self, event=None):
        idx = self.cb_presets.current()
        if idx <= 0:
            self.chain_path.set("")
            return
        _, filename = self.CHAIN_PRESETS[idx]
        actual_path = find_preset_chain(filename)
        if actual_path:
            self.chain_path.set(actual_path)
        else:
            base_dir = os.path.dirname(
                os.path.abspath(sys.executable if getattr(sys, 'frozen', False) else __file__)
            )
            self.chain_path.set(os.path.join(base_dir, "chains", filename))

    def _on_drop_input(self, event):
        path = parse_dnd_path(self.root, event.data)
        if path:
            self.input_path.set(path)
            self._auto_detect_format(path)

    def _on_drop_chain(self, event):
        path = parse_dnd_path(self.root, event.data)
        if path:
            self.chain_path.set(path)
            self.cb_presets.current(0)

    def _open_chain_download_page(self, event=None):
        try:
            opened = webbrowser.open(CHAIN_DOWNLOAD_URL)
        except Exception:
            opened = False
        if not opened:
            self._append_log(CHAIN_DOWNLOAD_URL)

    def _browse_input(self):
        filename = filedialog.askopenfilename(
            title=self.i18n.t("dialog_select_input"),
            filetypes=[(self.i18n.t("ft_genomic"), self.GENOMIC_PATTERNS), (self.i18n.t("ft_all"), "*.*")]
        )
        if filename:
            self.input_path.set(filename)
            self._auto_detect_format(filename)

    def _browse_chain(self):
        filename = filedialog.askopenfilename(
            title=self.i18n.t("dialog_select_chain"),
            filetypes=[(self.i18n.t("ft_chain"), "*.chain *.chain.gz"), (self.i18n.t("ft_all"), "*.*")]
        )
        if filename:
            self.chain_path.set(filename)
            self.cb_presets.current(0)

    def _browse_fasta(self):
        filename = filedialog.askopenfilename(
            title=self.i18n.t("dialog_select_fasta"),
            filetypes=[(self.i18n.t("ft_fasta"), "*.fa *.fasta *.fa.gz *.fasta.gz *.fna"), (self.i18n.t("ft_all"), "*.*")]
        )
        if filename:
            self.target_fasta_path.set(filename)

    def _browse_chrom_map(self):
        filename = filedialog.askopenfilename(
            title=self.i18n.t("dialog_select_map"),
            filetypes=[(self.i18n.t("ft_map"), "*.txt *.map *.tsv"), (self.i18n.t("ft_all"), "*.*")]
        )
        if filename:
            self.chrom_map_path.set(filename)

    def _save_out_dialog(self):
        default_ext = ""
        inp = self.input_path.get()
        if inp:
            _, ext = os.path.splitext(inp)
            default_ext = ext

        filename = filedialog.asksaveasfilename(
            title=self.i18n.t("dialog_save_result"),
            defaultextension=default_ext,
            filetypes=[(self.i18n.t("ft_all"), "*.*")]
        )
        if filename:
            self.out_path.set(filename)

    def _auto_detect_format(self, filename):
        fmt = detect_format(filename)
        if fmt:
            for label, code in self.FORMAT_MAPPING.items():
                if code == fmt:
                    self.selected_format.set(label)
                    break
        self.out_path.set(suggest_output_path(filename))

    def log(self, text):
        self.ui_queue.put(('log', text))

    def _append_log(self, text):
        self.log_area.config(state='normal')
        self.log_area.insert(tk.END, text + "\n")
        self.log_area.see(tk.END)
        self.log_area.config(state='disabled')

    def _poll_ui_queue(self):
        try:
            while True:
                kind, payload = self.ui_queue.get_nowait()
                if kind == 'log':
                    self._append_log(payload)
                elif kind == 'info':
                    messagebox.showinfo(*payload)
                elif kind == 'error':
                    messagebox.showerror(*payload)
                elif kind == 'done':
                    self.btn_run.config(state='normal')
                    self.btn_reset.config(state='normal')
        except queue.Empty:
            pass
        self.root.after(100, self._poll_ui_queue)

    def _load_chrom_map(self):
        return load_chrom_map(self.chrom_map_path.get().strip())

    def _start_processing(self):
        inp = self.input_path.get().strip()
        chain = self.chain_path.get().strip()
        fasta_file = self.target_fasta_path.get().strip()
        out = self.out_path.get().strip()
        map_file = self.chrom_map_path.get().strip()
        cmd_type = self.FORMAT_MAPPING[self.selected_format.get()]
        err = self.i18n.t("err_title")

        if not inp or not os.path.isfile(inp):
            messagebox.showerror(err, self.i18n.t("err_no_input"))
            return
        if not chain or not os.path.isfile(chain):
            messagebox.showerror(err, self.i18n.t("err_no_chain"))
            return
        if fasta_file and not os.path.isfile(fasta_file):
            messagebox.showerror(err, self.i18n.t("err_no_fasta"))
            return
        try:
            cmap = self._load_chrom_map()
        except OSError as e:
            messagebox.showerror(err, str(e))
            return

        path_error = check_paths(
            inp, out, cmd_type=cmd_type,
            save_unmapped=self.save_unmapped.get(), strict=bool(fasta_file),
            protected_paths=(chain, fasta_file, map_file),
        )
        if path_error:
            messagebox.showerror(err, self.i18n.t(path_error))
            return

        auto_chr = self.auto_chr.get()
        save_unm = self.save_unmapped.get()

        self.btn_run.config(state='disabled')
        self.btn_reset.config(state='disabled')
        threading.Thread(
            target=self._run_conversion,
            args=(cmd_type, chain, fasta_file, inp, out, cmap, auto_chr, save_unm, map_file),
            daemon=True
        ).start()

    def _run_conversion(self, cmd_type, chain_file, fasta_file, in_file, out_file, cmap, auto_chr, save_unm, map_file=None):
        t = self.i18n.t
        fasta_obj = None
        try:
            self.log(t("log_started", fmt=cmd_type.upper()))
            self.log(t("log_loading_chain", path=chain_file))
            lo = LiftOver(chain_file)

            if fasta_file:
                self.log(t("log_loading_fasta", path=fasta_file))
                try:
                    if not HAS_PYFAIDX:
                        raise ImportError("pyfaidx is not installed")
                    fasta_obj = Fasta(fasta_file, as_raw=True, sequence_always_upper=True, read_ahead=100000)
                except Exception as e:
                    raise RuntimeError(t("err_fasta_open", error=e))

            engine = LiftOverEngine(
                lo, fasta_obj, cmap, auto_chr, save_unm, self.log, t,
                chain_file=chain_file, protected_paths=(chain_file, fasta_file, map_file),
            )
            self.log(t("log_processing", path=in_file))
            st = engine.run(cmd_type, in_file, out_file)

            if st['inverted']:
                self.log(t("log_inverted", count=st['inverted']))
            if st['inverted_skipped']:
                self.log(t("log_inverted_skipped", count=st['inverted_skipped']))

            if fasta_obj is not None:
                self.log(t("log_summary_fasta", total=st['total'], mapped=st['mapped'], swapped=st['swapped'],
                           rejected=st['rejected'], unmapped=st['unmapped']))
                msg = t("success_message_fasta", mapped=st['mapped'], total=st['total'], swapped=st['swapped'],
                        rejected=st['rejected'], unmapped=st['unmapped'])
            else:
                self.log(t("log_summary", total=st['total'], mapped=st['mapped'], unmapped=st['unmapped']))
                msg = t("success_message", mapped=st['mapped'], total=st['total'], unmapped=st['unmapped'])
            self.log(t("log_finished"))
            self.ui_queue.put(('info', (t("success_title"), msg)))

        except Exception as e:
            self.log(t("log_error", error=e))
            self.ui_queue.put(('error', (t("err_title"), str(e))))
        finally:
            if fasta_obj is not None:
                try:
                    fasta_obj.close()
                except Exception:
                    pass
            self.ui_queue.put(('done', None))

def _prepare_stdio():
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors='replace')
        except Exception:
            pass


def run_cli(args):
    _prepare_stdio()
    i18n = I18nManager(default_lang=args.lang)
    t = i18n.t

    fmt = args.format
    if fmt == 'auto':
        fmt = detect_format(args.input)
        if not fmt:
            print(f"[{t('err_title')}] {t('err_no_input')} (Format could not be auto-detected)", file=sys.stderr)
            sys.exit(1)

    if not args.input or not os.path.isfile(args.input):
        print(f"[{t('err_title')}] {t('err_no_input')}", file=sys.stderr)
        sys.exit(1)

    chain_path = args.chain
    if not os.path.isfile(chain_path):
        preset = find_preset_chain(chain_path)
        if preset and os.path.isfile(preset):
            chain_path = preset
        else:
            print(f"[{t('err_title')}] {t('err_no_chain')}", file=sys.stderr)
            sys.exit(1)

    try:
        cmap = load_chrom_map(args.map)
    except OSError as e:
        print(t("log_error", error=e), file=sys.stderr)
        sys.exit(1)

    path_error = check_paths(
        args.input, args.output, cmd_type=fmt,
        save_unmapped=not args.no_unmapped, strict=bool(args.fasta),
        protected_paths=(chain_path, args.fasta, args.map),
    )
    if path_error:
        print(f"[{t('err_title')}] {t(path_error)}", file=sys.stderr)
        sys.exit(1)

    fasta_obj = None
    if args.fasta:
        if not os.path.isfile(args.fasta):
            print(f"[{t('err_title')}] {t('err_no_fasta')}", file=sys.stderr)
            sys.exit(1)
        if not HAS_PYFAIDX:
            print(f"[{t('err_title')}] pyfaidx library is missing.", file=sys.stderr)
            sys.exit(1)
        try:
            fasta_obj = Fasta(args.fasta, as_raw=True, sequence_always_upper=True, read_ahead=100000)
        except Exception as e:
            print(f"[{t('err_title')}] {t('err_fasta_open', error=e)}", file=sys.stderr)
            sys.exit(1)

    try:
        print(t("log_started", fmt=fmt.upper()))
        lo = LiftOver(chain_path)
        engine = LiftOverEngine(
            lo, fasta_obj, cmap,
            auto_chr=not args.no_auto_chr,
            save_unmapped=not args.no_unmapped,
            log=print, t=t, chain_file=chain_path,
            protected_paths=(chain_path, args.fasta, args.map),
        )
        st = engine.run(fmt, args.input, args.output)
    except Exception as e:
        print(t("log_error", error=e), file=sys.stderr)
        sys.exit(1)
    finally:
        if fasta_obj is not None:
            try:
                fasta_obj.close()
            except Exception:
                pass

    if st['inverted']:
        print(t("log_inverted", count=st['inverted']))
    if st['inverted_skipped']:
        print(t("log_inverted_skipped", count=st['inverted_skipped']))

    if fasta_obj is not None:
        print(t("log_summary_fasta", total=st['total'], mapped=st['mapped'], swapped=st['swapped'],
                   rejected=st['rejected'], unmapped=st['unmapped']))
    else:
        print(t("log_summary", total=st['total'], mapped=st['mapped'], unmapped=st['unmapped']))

    print(t("log_finished"))


if __name__ == "__main__":
    import argparse

    if len(sys.argv) > 1:
        parser = argparse.ArgumentParser(
            description="LiftOverZ - Standalone Genomic Coordinate Converter",
            formatter_class=argparse.RawTextHelpFormatter
        )
        parser.add_argument("-V", "--version", action="version", version=f"LiftOverZ {__version__}")
        parser.add_argument("-i", "--input", required=True, help="Input genomic file path")
        parser.add_argument("-c", "--chain", required=True, help="Chain file path OR built-in preset name (e.g. hg19ToHg38.over.chain.gz)")
        parser.add_argument("-o", "--output", required=True, help="Output file path")
        parser.add_argument("-f", "--format", choices=['vcf', 'gvcf', 'bed', 'gff', 'wig', 'maf', 'auto'], default='auto', help="Format of the input data (default: auto)")
        parser.add_argument("--fasta", help="Target reference FASTA for strict biological validation (optional)")
        parser.add_argument("--map", help="Chromosome mapping file (e.g., NC_000001.10 -> chr1) (optional)")
        parser.add_argument("--no-auto-chr", action="store_true", help="Disable auto-matching of 'chr' and 'MT' prefixes")
        parser.add_argument("--no-unmapped", action="store_true", help="Do not save unmapped records to a .unmapped file")
        parser.add_argument("--lang", default="en", help="Language code for logs (en, ru, zh, hi, es, ar, etc.) (default: en)")
        
        args = parser.parse_args()
        run_cli(args)
    elif not HAS_TK:
        print("Tkinter / tkinterdnd2 are not available: the graphical mode cannot start. "
              "Use the command-line mode (run with --help).", file=sys.stderr)
        sys.exit(2)
    else:
        root = TkinterDnD.Tk()
        app = LiftOverZApp(root)
        root.mainloop()
