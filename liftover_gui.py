import os
import sys
import glob
import json
import gzip
import locale
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinterdnd2 import DND_FILES, TkinterDnD
from pyliftover import LiftOver
from pyfaidx import Fasta


def detect_system_language():
    """Определяет двухбуквенный код языка интерфейса операционной системы"""
    # 1. Windows: системный вызов WinAPI
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

    # 2. Linux / macOS: переменные окружения
    for env_var in ('LANGUAGE', 'LC_ALL', 'LC_MESSAGES', 'LANG'):
        val = os.environ.get(env_var)
        if val:
            first_entry = val.split(':')[0].split('.')[0]
            code = first_entry.split('_')[0].lower()
            if len(code) == 2:
                return code

    # 3. Стандартная локаль Python
    try:
        loc = locale.getlocale()[0]
        if loc:
            code = loc.split('_')[0].lower()
            if len(code) == 2:
                return code
    except Exception:
        pass

    return "en"


def chrom_sort_key(chrom):
    """Канонический порядок сортировки хромосом: 1..22, X, Y, M/MT, затем остальные контиги"""
    c = chrom.lower()
    if c.startswith('chr'):
        c = c[3:]
    if c.isdigit():
        return (0, int(c), "")
    special = {'x': 23, 'y': 24, 'm': 25, 'mt': 25}
    if c in special:
        return (1, special[c], "")
    return (2, 0, c)


REV_COMP_TRANS = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rev_comp_allele(allele):
    """Обратный комплемент нуклеотидов с пропуском спец-символов VCF"""
    if allele.startswith('<') or allele == '*' or allele == '.':
        return allele
    return allele.translate(REV_COMP_TRANS)[::-1]


def get_fasta_seq(fasta_obj, chrom, start_0, end_0):
    """Извлекает срез нуклеотидов из индекса FASTA с авто-сопоставлением имен хромосом"""
    try:
        if chrom in fasta_obj:
            return str(fasta_obj[chrom][start_0:end_0]).upper()
        if chrom.startswith("chr") and chrom[3:] in fasta_obj:
            return str(fasta_obj[chrom[3:]][start_0:end_0]).upper()
        if not chrom.startswith("chr") and ("chr" + chrom) in fasta_obj:
            return str(fasta_obj["chr" + chrom][start_0:end_0]).upper()
        if chrom.upper() in ("MT", "CHRM", "M", "CHRMT"):
            for m_alias in ("chrM", "MT", "M", "chrMT"):
                if m_alias in fasta_obj:
                    return str(fasta_obj[m_alias][start_0:end_0]).upper()
    except Exception:
        pass
    return None


def swap_genotypes_and_format(format_str, sample_str):
    """Инвертирует GT (0<->1), а также пересчитывает PL и AD при замене REF/ALT"""
    if not format_str or not sample_str or sample_str == '.':
        return sample_str
    fmt_keys = format_str.split(':')
    sm_vals = sample_str.split(':')

    gt_idx = fmt_keys.index('GT') if 'GT' in fmt_keys else -1
    pl_idx = fmt_keys.index('PL') if 'PL' in fmt_keys else -1
    ad_idx = fmt_keys.index('AD') if 'AD' in fmt_keys else -1

    if gt_idx != -1 and gt_idx < len(sm_vals):
        gt_val = sm_vals[gt_idx]
        sep = '/' if '/' in gt_val else ('|' if '|' in gt_val else None)
        if sep:
            alleles = gt_val.split(sep)
            new_alleles = []
            for a in alleles:
                if a == '0':
                    new_alleles.append('1')
                elif a == '1':
                    new_alleles.append('0')
                else:
                    new_alleles.append(a)
            sm_vals[gt_idx] = sep.join(new_alleles)

    if pl_idx != -1 and pl_idx < len(sm_vals):
        pls = sm_vals[pl_idx].split(',')
        if len(pls) == 3:
            sm_vals[pl_idx] = f"{pls[2]},{pls[1]},{pls[0]}"

    if ad_idx != -1 and ad_idx < len(sm_vals):
        ads = sm_vals[ad_idx].split(',')
        if len(ads) == 2:
            sm_vals[ad_idx] = f"{ads[1]},{ads[0]}"

    return ':'.join(sm_vals)


class I18nManager:
    """Менеджер локализации без захардкоженных языковых строк с fallback на en.json"""
    def __init__(self, default_lang=None):
        self.languages = {}
        self.translations = {}
        self.fallback_translations = {}
        self._load_available_languages()

        if "en" in self.languages:
            try:
                with open(self.languages["en"]["path"], 'r', encoding='utf-8') as f:
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
        return dirs

    def _load_available_languages(self):
        for loc_dir in self._get_locales_dirs():
            if os.path.exists(loc_dir):
                for fpath in glob.glob(os.path.join(loc_dir, "*.json")):
                    code = os.path.splitext(os.path.basename(fpath))[0].lower()
                    try:
                        with open(fpath, 'r', encoding='utf-8') as f:
                            data = json.load(f)
                            self.languages[code] = {"name": data.get("lang_name", code.upper()), "path": fpath}
                    except Exception:
                        pass

    def set_language(self, code):
        if code in self.languages:
            self.current_lang = code
            try:
                with open(self.languages[code]["path"], 'r', encoding='utf-8') as f:
                    self.translations = json.load(f)
            except Exception:
                self.translations = {}

    def t(self, key, **kwargs):
        text = self.translations.get(key, self.fallback_translations.get(key, key))
        if kwargs:
            try:
                text = text.format(**kwargs)
            except KeyError:
                pass
        return text


def parse_dnd_path(event_data):
    path = event_data.strip()
    if path.startswith('{') and path.endswith('}'):
        path = path[1:-1]
    return path


def find_preset_chain(filename):
    """Поиск встроенного chain-файла в папке chains внутри сборки или рядом с ней"""
    search_dirs = []
    if hasattr(sys, '_MEIPASS'):
        search_dirs.extend([
            os.path.join(sys._MEIPASS, "chains"),
            sys._MEIPASS
        ])
    base_dir = os.path.dirname(os.path.abspath(sys.executable if getattr(sys, 'frozen', False) else __file__))
    search_dirs.extend([
        os.path.join(base_dir, "chains"),
        base_dir
    ])

    for d in search_dirs:
        candidate = os.path.join(d, filename)
        if os.path.exists(candidate):
            return candidate
    return None


def open_file_read(filepath):
    return gzip.open(filepath, 'rt', encoding='utf-8', errors='replace') if filepath.endswith('.gz') else open(filepath, 'rt', encoding='utf-8', errors='replace')


def open_file_write(filepath):
    return gzip.open(filepath, 'wt', encoding='utf-8') if filepath.endswith('.gz') else open(filepath, 'w', encoding='utf-8')


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

    def __init__(self, root):
        self.root = root
        self.i18n = I18nManager()

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

        # 1. Входной файл
        self.lbl_input = ttk.Label(self.files_frame)
        self.lbl_input.grid(row=0, column=0, sticky="w", pady=3)
        self.ent_inp = ttk.Entry(self.files_frame, textvariable=self.input_path)
        self.ent_inp.grid(row=0, column=1, sticky="ew", padx=5, pady=3)
        self.btn_browse_inp = ttk.Button(self.files_frame, command=self._browse_input)
        self.btn_browse_inp.grid(row=0, column=2, pady=3)
        self.ent_inp.drop_target_register(DND_FILES)
        self.ent_inp.dnd_bind('<<Drop>>', self._on_drop_input)

        # 2. Пресет Chain
        self.lbl_chain_preset = ttk.Label(self.files_frame)
        self.lbl_chain_preset.grid(row=1, column=0, sticky="w", pady=3)
        self.cb_presets = ttk.Combobox(self.files_frame, state="readonly")
        self.cb_presets.grid(row=1, column=1, columnspan=2, sticky="ew", padx=5, pady=3)
        self.cb_presets.bind("<<ComboboxSelected>>", self._on_preset_selected)

        # 3. Активный Chain файл
        self.lbl_chain = ttk.Label(self.files_frame)
        self.lbl_chain.grid(row=2, column=0, sticky="w", pady=3)
        self.ent_chain = ttk.Entry(self.files_frame, textvariable=self.chain_path)
        self.ent_chain.grid(row=2, column=1, sticky="ew", padx=5, pady=3)
        self.btn_browse_chain = ttk.Button(self.files_frame, command=self._browse_chain)
        self.btn_browse_chain.grid(row=2, column=2, pady=3)
        self.ent_chain.drop_target_register(DND_FILES)
        self.ent_chain.dnd_bind('<<Drop>>', self._on_drop_chain)

        # 4. Целевой референс FASTA (Опционально)
        self.lbl_target_fasta = ttk.Label(self.files_frame)
        self.lbl_target_fasta.grid(row=3, column=0, sticky="w", pady=3)
        self.ent_fasta = ttk.Entry(self.files_frame, textvariable=self.target_fasta_path)
        self.ent_fasta.grid(row=3, column=1, sticky="ew", padx=5, pady=3)
        self.btn_browse_fasta = ttk.Button(self.files_frame, command=self._browse_fasta)
        self.btn_browse_fasta.grid(row=3, column=2, pady=3)
        self.ent_fasta.drop_target_register(DND_FILES)
        self.ent_fasta.dnd_bind('<<Drop>>', lambda e: self.target_fasta_path.set(parse_dnd_path(e.data)))

        # 5. Выходной файл
        self.lbl_output = ttk.Label(self.files_frame)
        self.lbl_output.grid(row=4, column=0, sticky="w", pady=3)
        self.ent_out = ttk.Entry(self.files_frame, textvariable=self.out_path)
        self.ent_out.grid(row=4, column=1, sticky="ew", padx=5, pady=3)
        self.btn_save_out = ttk.Button(self.files_frame, command=self._save_out_dialog)
        self.btn_save_out.grid(row=4, column=2, pady=3)

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
        self.ent_cmap.dnd_bind('<<Drop>>', lambda e: self.chrom_map_path.set(parse_dnd_path(e.data)))

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
        self.root.title(self.i18n.t("app_title"))
        self.lbl_lang.config(text=self.i18n.t("language_label"))
        self.fmt_frame.config(text=self.i18n.t("step_1_format"))
        self.files_frame.config(text=self.i18n.t("step_2_files"))
        self.opts_frame.config(text=self.i18n.t("step_3_options"))
        self.log_frame.config(text=self.i18n.t("log_title"))

        self.lbl_input.config(text=self.i18n.t("lbl_input"))
        self.lbl_chain_preset.config(text=self.i18n.t("lbl_chain_preset"))
        self.lbl_chain.config(text=self.i18n.t("lbl_chain"))
        self.lbl_target_fasta.config(text=self.i18n.t("lbl_target_fasta"))
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
        """Полная очистка полей и логов"""
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
            return
        _, filename = self.CHAIN_PRESETS[idx]
        actual_path = find_preset_chain(filename)
        if actual_path:
            self.chain_path.set(actual_path)
        else:
            self.chain_path.set(os.path.join("chains", filename))

    def _on_drop_input(self, event):
        path = parse_dnd_path(event.data)
        self.input_path.set(path)
        self._auto_detect_format(path)

    def _on_drop_chain(self, event):
        path = parse_dnd_path(event.data)
        self.chain_path.set(path)
        self.cb_presets.current(0)

    def _browse_input(self):
        filename = filedialog.askopenfilename(
            title=self.i18n.t("dialog_select_input"),
            filetypes=[("Genomic Files", "*.vcf *.vcf.gz *.gvcf *.gvcf.gz *.bed *.bedgraph *.gff *.gtf *.gff3 *.wig *.maf"), ("All Files", "*.*")]
        )
        if filename:
            self.input_path.set(filename)
            self._auto_detect_format(filename)

    def _browse_chain(self):
        filename = filedialog.askopenfilename(
            title=self.i18n.t("dialog_select_chain"),
            filetypes=[("Chain Files", "*.chain *.chain.gz"), ("All Files", "*.*")]
        )
        if filename:
            self.chain_path.set(filename)
            self.cb_presets.current(0)

    def _browse_fasta(self):
        filename = filedialog.askopenfilename(
            title=self.i18n.t("dialog_select_fasta"),
            filetypes=[("FASTA Files", "*.fa *.fasta *.fa.gz *.fasta.gz"), ("All Files", "*.*")]
        )
        if filename:
            self.target_fasta_path.set(filename)

    def _browse_chrom_map(self):
        filename = filedialog.askopenfilename(
            title=self.i18n.t("dialog_select_map"),
            filetypes=[("Map Files", "*.txt *.map *.tsv"), ("All Files", "*.*")]
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
            filetypes=[("All Files", "*.*")]
        )
        if filename:
            self.out_path.set(filename)

    def _auto_detect_format(self, filename):
        ext = filename.lower()
        if ext.endswith(('.vcf', '.vcf.gz')):
            self.selected_format.set("VCF (.vcf, .vcf.gz)")
        elif ext.endswith(('.gvcf', '.gvcf.gz')):
            self.selected_format.set("gVCF (.gvcf, .gvcf.gz)")
        elif ext.endswith(('.bed', '.bedgraph')):
            self.selected_format.set("BED / BEDGraph (.bed, .bedgraph)")
        elif ext.endswith(('.gff', '.gff3', '.gtf')):
            self.selected_format.set("GFF / GTF (.gff, .gtf, .gff3)")
        elif ext.endswith('.wig'):
            self.selected_format.set("Wiggle (.wig)")
        elif ext.endswith('.maf'):
            self.selected_format.set("MAF (.maf)")

        if filename.endswith('.gz'):
            base, ext_inner = os.path.splitext(filename[:-3])
            ext_name = ext_inner + '.gz'
        else:
            base, ext_name = os.path.splitext(filename)
        self.out_path.set(f"{base}_lifted{ext_name}")

    def log(self, text):
        self.log_area.config(state='normal')
        self.log_area.insert(tk.END, text + "\n")
        self.log_area.see(tk.END)
        self.log_area.config(state='disabled')

    def _start_processing(self):
        inp = self.input_path.get()
        chain = self.chain_path.get()
        fasta_file = self.target_fasta_path.get().strip()
        out = self.out_path.get()
        cmd_type = self.FORMAT_MAPPING[self.selected_format.get()]

        if not inp or not os.path.exists(inp):
            messagebox.showerror(self.i18n.t("err_title"), self.i18n.t("err_no_input"))
            return
        if not chain or not os.path.exists(chain):
            messagebox.showerror(self.i18n.t("err_title"), self.i18n.t("err_no_chain"))
            return
        if fasta_file and not os.path.exists(fasta_file):
            messagebox.showerror(self.i18n.t("err_title"), self.i18n.t("dialog_select_fasta"))
            return
        if not out:
            messagebox.showerror(self.i18n.t("err_title"), self.i18n.t("err_no_output"))
            return

        self.btn_run.config(state='disabled')
        self.btn_reset.config(state='disabled')
        threading.Thread(target=self._run_conversion, args=(cmd_type, chain, fasta_file, inp, out), daemon=True).start()

    def _load_chrom_map(self):
        m_path = self.chrom_map_path.get()
        cmap = {}
        if m_path and os.path.exists(m_path):
            with open(m_path, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        cmap[parts[0]] = parts[1]
        return cmap

    def _lift_coord(self, lo, chrom, pos_0, cmap, auto_chr):
        target_chrom = cmap.get(chrom, chrom)
        res = lo.convert_coordinate(target_chrom, pos_0)
        if res:
            return res[0]

        if auto_chr:
            if target_chrom.upper() in ('MT', 'CHRM', 'M', 'CHRMT'):
                for m_alias in ('chrM', 'MT', 'M', 'chrMT'):
                    if m_alias != target_chrom:
                        res = lo.convert_coordinate(m_alias, pos_0)
                        if res:
                            return res[0]

            alt_chrom = target_chrom[3:] if target_chrom.startswith('chr') else 'chr' + target_chrom
            res = lo.convert_coordinate(alt_chrom, pos_0)
            if res:
                return res[0]

        return None

    @staticmethod
    def _format_output_chrom(orig_chrom, target_chrom):
        if orig_chrom.upper() == 'MT':
            return 'MT'
        if orig_chrom.lower() == 'chrm':
            return 'chrM'
        if orig_chrom.upper() == 'M':
            return 'M'

        orig_has_chr = orig_chrom.lower().startswith('chr')
        target_has_chr = target_chrom.lower().startswith('chr')

        if not orig_has_chr and target_has_chr:
            stripped = target_chrom[3:]
            if stripped.upper() == 'M':
                return 'MT' if orig_chrom.upper() == 'MT' else 'M'
            return stripped
        elif orig_has_chr and not target_has_chr:
            if target_chrom.upper() == 'MT':
                return 'chrM'
            return 'chr' + target_chrom
        return target_chrom

    def _run_conversion(self, cmd_type, chain_file, fasta_file, in_file, out_file):
        try:
            self.log(self.i18n.t("log_started", fmt=cmd_type.upper()))
            self.log(self.i18n.t("log_loading_chain", path=chain_file))
            lo = LiftOver(chain_file)

            fasta_obj = None
            if fasta_file:
                self.log(self.i18n.t("log_loading_fasta", path=fasta_file))
                fasta_obj = Fasta(fasta_file, as_raw=True, sequence_always_upper=True, read_ahead=100000)

            cmap = self._load_chrom_map()
            auto_chr = self.auto_chr.get()
            save_unm = self.save_unmapped.get()

            unm_file = out_file + ".unmapped" if save_unm else None
            f_unm = open_file_write(unm_file) if unm_file else None

            rej_file = (out_file + ".rejected.vcf") if cmd_type in ['vcf', 'gvcf'] else (out_file + ".rejected")
            f_rej = open_file_write(rej_file) if fasta_obj else None

            self.log(self.i18n.t("log_processing", path=in_file))
            total = 0
            mapped = 0
            unmapped = 0
            swapped = 0
            rejected = 0

            converted_records = []

            with open_file_read(in_file) as f_in, open_file_write(out_file) as f_out:
                if cmd_type in ['vcf', 'gvcf']:
                    for line in f_in:
                        if line.startswith('##'):
                            f_out.write(line)
                            if f_rej:
                                f_rej.write(line)
                            continue
                        if line.startswith('#CHROM'):
                            if fasta_obj:
                                f_out.write('##INFO=<ID=LIFTOVER_SWAPPED,Number=0,Type=Flag,Description="Reference and alternate alleles were swapped during liftover to match target genome reference">\n')
                                f_out.write('##INFO=<ID=LIFTOVER_REVERSE_COMPLEMENTED,Number=0,Type=Flag,Description="Alleles were reverse-complemented during liftover due to negative strand alignment">\n')
                                if f_rej:
                                    f_rej.write('##INFO=<ID=LIFTOVER_REJECT_REASON,Number=1,Type=String,Description="Reason for liftover rejection">\n')
                            f_out.write(line)
                            if f_rej:
                                f_rej.write(line)
                            continue

                        total += 1
                        parts = line.rstrip('\r\n').split('\t')
                        if len(parts) < 5:
                            continue

                        chrom = parts[0]
                        pos_1 = int(parts[1])
                        ref = parts[3].upper()
                        alt_str = parts[4].upper()
                        alts = alt_str.split(',')

                        conv = self._lift_coord(lo, chrom, pos_1 - 1, cmap, auto_chr)
                        if not conv:
                            unmapped += 1
                            if f_unm:
                                f_unm.write(line)
                            continue

                        target_c = self._format_output_chrom(chrom, conv[0])
                        new_pos = conv[1] + 1
                        strand = conv[2] if len(conv) > 2 else '+'

                        # 1. Быстрый режим без FASTA
                        if not fasta_obj:
                            parts[0] = target_c
                            parts[1] = str(new_pos)
                            converted_records.append((target_c, new_pos, '\t'.join(parts) + '\n'))
                            mapped += 1
                        else:
                            # 2. Строгий режим с валидацией по FASTA
                            strand_flipped = False
                            if strand == '-':
                                cur_ref = rev_comp_allele(ref)
                                cur_alts = [rev_comp_allele(a) for a in alts]
                                strand_flipped = True
                            else:
                                cur_ref = ref
                                cur_alts = alts

                            fasta_ref = get_fasta_seq(fasta_obj, conv[0], new_pos - 1, new_pos - 1 + len(cur_ref))

                            if fasta_ref is None:
                                rejected += 1
                                if f_rej:
                                    parts_rej = list(parts)
                                    tag = "LIFTOVER_REJECT_REASON=ContigNotFoundInFasta"
                                    parts_rej[7] = tag if parts_rej[7] == '.' else parts_rej[7] + ';' + tag
                                    f_rej.write('\t'.join(parts_rej) + '\n')
                                continue

                            # Случай 1: Референс совпал
                            if cur_ref == fasta_ref:
                                parts[0] = target_c
                                parts[1] = str(new_pos)
                                parts[3] = cur_ref
                                parts[4] = ','.join(cur_alts)
                                if strand_flipped:
                                    parts[7] = "LIFTOVER_REVERSE_COMPLEMENTED" if parts[7] == '.' else parts[7] + ";LIFTOVER_REVERSE_COMPLEMENTED"
                                converted_records.append((target_c, new_pos, '\t'.join(parts) + '\n'))
                                mapped += 1

                            # Случай 2: Референс сменился на мутантный аллель (Swap REF/ALT)
                            elif len(cur_alts) == 1 and cur_alts[0] == fasta_ref and len(cur_ref) == len(fasta_ref):
                                swapped_ref = cur_alts[0]
                                swapped_alt = cur_ref
                                parts[0] = target_c
                                parts[1] = str(new_pos)
                                parts[3] = swapped_ref
                                parts[4] = swapped_alt

                                if len(parts) >= 10:
                                    fmt = parts[8]
                                    for s_idx in range(9, len(parts)):
                                        parts[s_idx] = swap_genotypes_and_format(fmt, parts[s_idx])

                                info_tags = ["LIFTOVER_SWAPPED"]
                                if strand_flipped:
                                    info_tags.append("LIFTOVER_REVERSE_COMPLEMENTED")
                                add_info = ';'.join(info_tags)
                                parts[7] = add_info if parts[7] == '.' else parts[7] + ';' + add_info

                                converted_records.append((target_c, new_pos, '\t'.join(parts) + '\n'))
                                mapped += 1
                                swapped += 1

                            # Случай 3: Несовпадение референса (Mismatch)
                            else:
                                rejected += 1
                                if f_rej:
                                    parts_rej = list(parts)
                                    tag = "LIFTOVER_REJECT_REASON=MismatchedRefAllele"
                                    parts_rej[7] = tag if parts_rej[7] == '.' else parts_rej[7] + ';' + tag
                                    f_rej.write('\t'.join(parts_rej) + '\n')

                        if total % 25000 == 0:
                            self.log(self.i18n.t("log_progress", total=total, mapped=mapped, unmapped=unmapped))

                    if converted_records:
                        self.log(self.i18n.t("log_sorting", count=len(converted_records)))
                        converted_records.sort(key=lambda x: (chrom_sort_key(x[0]), x[1]))
                        for rec in converted_records:
                            f_out.write(rec[2])
                        converted_records.clear()

                elif cmd_type == 'bed':
                    for line in f_in:
                        if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                            f_out.write(line)
                            continue
                        total += 1
                        parts = line.rstrip('\r\n').split('\t')
                        if len(parts) < 3:
                            continue
                        chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                        conv_s = self._lift_coord(lo, chrom, start, cmap, auto_chr)
                        conv_e = self._lift_coord(lo, chrom, end - 1, cmap, auto_chr)
                        if conv_s and conv_e and conv_s[0] == conv_e[0]:
                            target_c = self._format_output_chrom(chrom, conv_s[0])
                            new_s = min(conv_s[1], conv_e[1])
                            new_e = max(conv_s[1], conv_e[1]) + 1
                            parts[0] = target_c
                            parts[1] = str(new_s)
                            parts[2] = str(new_e)
                            converted_records.append((target_c, new_s, new_e, '\t'.join(parts) + '\n'))
                            mapped += 1
                        else:
                            unmapped += 1
                            if f_unm:
                                f_unm.write(line)
                        if total % 25000 == 0:
                            self.log(self.i18n.t("log_progress", total=total, mapped=mapped, unmapped=unmapped))

                    if converted_records:
                        self.log(self.i18n.t("log_sorting", count=len(converted_records)))
                        converted_records.sort(key=lambda x: (chrom_sort_key(x[0]), x[1], x[2]))
                        for rec in converted_records:
                            f_out.write(rec[3])
                        converted_records.clear()

                elif cmd_type == 'gff':
                    for line in f_in:
                        if line.startswith('#'):
                            f_out.write(line)
                            continue
                        total += 1
                        parts = line.rstrip('\r\n').split('\t')
                        if len(parts) < 5:
                            continue
                        chrom, start, end = parts[0], int(parts[3]) - 1, int(parts[4]) - 1
                        conv_s = self._lift_coord(lo, chrom, start, cmap, auto_chr)
                        conv_e = self._lift_coord(lo, chrom, end, cmap, auto_chr)
                        if conv_s and conv_e and conv_s[0] == conv_e[0]:
                            target_c = self._format_output_chrom(chrom, conv_s[0])
                            new_s = min(conv_s[1], conv_e[1]) + 1
                            new_e = max(conv_s[1], conv_e[1]) + 1
                            parts[0] = target_c
                            parts[3] = str(new_s)
                            parts[4] = str(new_e)
                            converted_records.append((target_c, new_s, new_e, '\t'.join(parts) + '\n'))
                            mapped += 1
                        else:
                            unmapped += 1
                            if f_unm:
                                f_unm.write(line)
                        if total % 25000 == 0:
                            self.log(self.i18n.t("log_progress", total=total, mapped=mapped, unmapped=unmapped))

                    if converted_records:
                        self.log(self.i18n.t("log_sorting", count=len(converted_records)))
                        converted_records.sort(key=lambda x: (chrom_sort_key(x[0]), x[1], x[2]))
                        for rec in converted_records:
                            f_out.write(rec[3])
                        converted_records.clear()

                elif cmd_type == 'wig':
                    curr_chrom = None
                    for line in f_in:
                        sline = line.strip()
                        if not sline or sline.startswith('#') or sline.startswith('track'):
                            f_out.write(line)
                            continue
                        if sline.startswith('variableStep') or sline.startswith('fixedStep'):
                            tokens = sline.split()
                            for t in tokens:
                                if t.startswith('chrom='):
                                    curr_chrom = t.split('=')[1]
                            f_out.write(line)
                            continue
                        parts = sline.split()
                        if len(parts) >= 2 and curr_chrom:
                            total += 1
                            pos = int(parts[0]) - 1
                            conv = self._lift_coord(lo, curr_chrom, pos, cmap, auto_chr)
                            if conv:
                                parts[0] = str(conv[1] + 1)
                                f_out.write('\t'.join(parts) + '\n')
                                mapped += 1
                            else:
                                unmapped += 1
                                if f_unm:
                                    f_unm.write(line)
                        else:
                            f_out.write(line)

                elif cmd_type == 'maf':
                    chr_idx = 4
                    start_idx = 5
                    end_idx = 6
                    header_found = False
                    for line in f_in:
                        if line.startswith('#'):
                            f_out.write(line)
                            continue
                        parts = line.rstrip('\r\n').split('\t')
                        if not header_found:
                            header_found = True
                            if 'Chromosome' in parts:
                                chr_idx = parts.index('Chromosome')
                            if 'Start_Position' in parts:
                                start_idx = parts.index('Start_Position')
                            if 'End_Position' in parts:
                                end_idx = parts.index('End_Position')
                            f_out.write(line)
                            continue
                        total += 1
                        if len(parts) > max(chr_idx, start_idx, end_idx):
                            chrom = parts[chr_idx]
                            start = int(parts[start_idx]) - 1
                            end = int(parts[end_idx]) - 1
                            conv_s = self._lift_coord(lo, chrom, start, cmap, auto_chr)
                            conv_e = self._lift_coord(lo, chrom, end, cmap, auto_chr)
                            if conv_s and conv_e and conv_s[0] == conv_e[0]:
                                target_c = self._format_output_chrom(chrom, conv_s[0])
                                new_s = min(conv_s[1], conv_e[1]) + 1
                                new_e = max(conv_s[1], conv_e[1]) + 1
                                parts[chr_idx] = target_c
                                parts[start_idx] = str(new_s)
                                parts[end_idx] = str(new_e)
                                converted_records.append((target_c, new_s, new_e, '\t'.join(parts) + '\n'))
                                mapped += 1
                            else:
                                unmapped += 1
                                if f_unm:
                                    f_unm.write(line)
                        if total % 25000 == 0:
                            self.log(self.i18n.t("log_progress", total=total, mapped=mapped, unmapped=unmapped))

                    if converted_records:
                        self.log(self.i18n.t("log_sorting", count=len(converted_records)))
                        converted_records.sort(key=lambda x: (chrom_sort_key(x[0]), x[1], x[2]))
                        for rec in converted_records:
                            f_out.write(rec[3])
                        converted_records.clear()

            if f_unm:
                f_unm.close()
                if unmapped == 0 and os.path.exists(unm_file):
                    os.remove(unm_file)

            if f_rej:
                f_rej.close()
                if rejected == 0 and os.path.exists(rej_file):
                    os.remove(rej_file)

            if fasta_obj:
                self.log(self.i18n.t("log_summary_fasta", total=total, mapped=mapped, swapped=swapped, rejected=rejected, unmapped=unmapped))
                self.log(self.i18n.t("log_finished"))
                messagebox.showinfo(self.i18n.t("success_title"), self.i18n.t("success_message_fasta", mapped=mapped, total=total, swapped=swapped, rejected=rejected, unmapped=unmapped))
            else:
                self.log(self.i18n.t("log_summary", total=total, mapped=mapped, unmapped=unmapped))
                self.log(self.i18n.t("log_finished"))
                messagebox.showinfo(self.i18n.t("success_title"), self.i18n.t("success_message", mapped=mapped, total=total, unmapped=unmapped))

        except Exception as e:
            self.log(f"ERROR: {str(e)}")
            messagebox.showerror(self.i18n.t("err_title"), str(e))
        finally:
            self.btn_run.config(state='normal')
            self.btn_reset.config(state='normal')


if __name__ == "__main__":
    root = TkinterDnD.Tk()
    app = LiftOverZApp(root)
    root.mainloop()
