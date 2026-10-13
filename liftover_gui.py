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


def detect_system_language():
    """Определяет двухбуквенный код языка интерфейса операционной системы"""
    # 1. Windows: чтение системной локали через WinAPI
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

    # 2. Linux / macOS / POSIX: системные переменные окружения
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


class I18nManager:
    """Универсальный менеджер локализации с каскадным fallback на английский"""
    def __init__(self, default_lang=None):
        self.languages = {}
        self.translations = {}
        self.fallback_translations = {}
        self._load_available_languages()

        # Загружаем базовый английский словарь для fallback
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
        # 1. Текущий язык -> 2. Английский fallback -> 3. Имя самого ключа
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
    """Поиск встроенного chain-файла в ресурсах сборки и локальной директории"""
    search_dirs = []
    if hasattr(sys, '_MEIPASS'):
        search_dirs.extend([
            os.path.join(sys._MEIPASS, "chains"),
            os.path.join(sys._MEIPASS, "changes"),
            sys._MEIPASS
        ])
    base_dir = os.path.dirname(os.path.abspath(sys.executable if getattr(sys, 'frozen', False) else __file__))
    search_dirs.extend([
        os.path.join(base_dir, "chains"),
        os.path.join(base_dir, "changes"),
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

        self.root.geometry("760x730")
        self.root.minsize(700, 650)

        self.input_path = tk.StringVar()
        self.chain_path = tk.StringVar()
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

        self.lbl_output = ttk.Label(self.files_frame)
        self.lbl_output.grid(row=3, column=0, sticky="w", pady=3)
        self.ent_out = ttk.Entry(self.files_frame, textvariable=self.out_path)
        self.ent_out.grid(row=3, column=1, sticky="ew", padx=5, pady=3)
        self.btn_save_out = ttk.Button(self.files_frame, command=self._save_out_dialog)
        self.btn_save_out.grid(row=3, column=2, pady=3)

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

        # Панель кнопок: Запуск и Сброс
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
        self.lbl_output.config(text=self.i18n.t("lbl_output"))
        self.lbl_chrom_map.config(text=self.i18n.t("lbl_chrom_map"))
        self.chk_auto_chr_w.config(text=self.i18n.t("chk_auto_chr"))
        self.chk_save_unm_w.config(text=self.i18n.t("chk_save_unmapped"))

        self.btn_browse_inp.config(text=self.i18n.t("btn_browse"))
        self.btn_browse_chain.config(text=self.i18n.t("btn_browse"))
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
        """Очищает поля, восстанавливает настройки по умолчанию и сбрасывает логи"""
        self.input_path.set("")
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

        # Всегда формируем актуальное имя файла результата при смене входного файла
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
        out = self.out_path.get()
        cmd_type = self.FORMAT_MAPPING[self.selected_format.get()]

        if not inp or not os.path.exists(inp):
            messagebox.showerror(self.i18n.t("err_title"), self.i18n.t("err_no_input"))
            return
        if not chain or not os.path.exists(chain):
            messagebox.showerror(self.i18n.t("err_title"), self.i18n.t("err_no_chain"))
            return
        if not out:
            messagebox.showerror(self.i18n.t("err_title"), self.i18n.t("err_no_output"))
            return

        self.btn_run.config(state='disabled')
        self.btn_reset.config(state='disabled')
        threading.Thread(target=self._run_conversion, args=(cmd_type, chain, inp, out), daemon=True).start()

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
            # 1. Автоматическая обработка митохондрий (MT, chrM, M, chrMT)
            if target_chrom.upper() in ('MT', 'CHRM', 'M', 'CHRMT'):
                for m_alias in ('chrM', 'MT', 'M', 'chrMT'):
                    if m_alias != target_chrom:
                        res = lo.convert_coordinate(m_alias, pos_0)
                        if res:
                            return res[0]

            # 2. Стандартное сопоставление chr <-> non-chr (1 <-> chr1)
            alt_chrom = target_chrom[3:] if target_chrom.startswith('chr') else 'chr' + target_chrom
            res = lo.convert_coordinate(alt_chrom, pos_0)
            if res:
                return res[0]

        return None

    @staticmethod
    def _format_output_chrom(orig_chrom, target_chrom):
        """
        Сохраняет исходный стиль именования хромосом входного файла.
        Если входной файл был без префикса 'chr' (например, Ensembl/NCBI: 1, 2, MT),
        то результат конверсии (UCSC 'chr1', 'chrM') возвращается без 'chr' (1, 2, MT).
        """
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

    def _run_conversion(self, cmd_type, chain_file, in_file, out_file):
        try:
            self.log(self.i18n.t("log_started", fmt=cmd_type.upper()))
            self.log(self.i18n.t("log_loading_chain", path=chain_file))
            lo = LiftOver(chain_file)

            cmap = self._load_chrom_map()
            auto_chr = self.auto_chr.get()
            save_unm = self.save_unmapped.get()

            unm_file = out_file + ".unmapped" if save_unm else None
            f_unm = open_file_write(unm_file) if unm_file else None

            self.log(self.i18n.t("log_processing", path=in_file))
            total = 0
            mapped = 0
            unmapped = 0

            with open_file_read(in_file) as f_in, open_file_write(out_file) as f_out:
                if cmd_type in ['vcf', 'gvcf']:
                    for line in f_in:
                        if line.startswith('#'):
                            f_out.write(line)
                            continue
                        total += 1
                        parts = line.rstrip('\r\n').split('\t')
                        if len(parts) < 2:
                            continue
                        chrom, pos = parts[0], int(parts[1]) - 1
                        conv = self._lift_coord(lo, chrom, pos, cmap, auto_chr)
                        if conv:
                            parts[0] = self._format_output_chrom(chrom, conv[0])
                            parts[1] = str(conv[1] + 1)
                            f_out.write('\t'.join(parts) + '\n')
                            mapped += 1
                        else:
                            unmapped += 1
                            if f_unm:
                                f_unm.write(line)
                        if total % 25000 == 0:
                            self.log(self.i18n.t("log_progress", total=total, mapped=mapped, unmapped=unmapped))

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
                            parts[0] = self._format_output_chrom(chrom, conv_s[0])
                            new_s = min(conv_s[1], conv_e[1])
                            new_e = max(conv_s[1], conv_e[1]) + 1
                            parts[1] = str(new_s)
                            parts[2] = str(new_e)
                            f_out.write('\t'.join(parts) + '\n')
                            mapped += 1
                        else:
                            unmapped += 1
                            if f_unm:
                                f_unm.write(line)
                        if total % 25000 == 0:
                            self.log(self.i18n.t("log_progress", total=total, mapped=mapped, unmapped=unmapped))

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
                            parts[0] = self._format_output_chrom(chrom, conv_s[0])
                            parts[3] = str(min(conv_s[1], conv_e[1]) + 1)
                            parts[4] = str(max(conv_s[1], conv_e[1]) + 1)
                            f_out.write('\t'.join(parts) + '\n')
                            mapped += 1
                        else:
                            unmapped += 1
                            if f_unm:
                                f_unm.write(line)
                        if total % 25000 == 0:
                            self.log(self.i18n.t("log_progress", total=total, mapped=mapped, unmapped=unmapped))

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
                                parts[chr_idx] = self._format_output_chrom(chrom, conv_s[0])
                                parts[start_idx] = str(min(conv_s[1], conv_e[1]) + 1)
                                parts[end_idx] = str(max(conv_s[1], conv_e[1]) + 1)
                                f_out.write('\t'.join(parts) + '\n')
                                mapped += 1
                            else:
                                unmapped += 1
                                if f_unm:
                                    f_unm.write(line)
                        if total % 25000 == 0:
                            self.log(self.i18n.t("log_progress", total=total, mapped=mapped, unmapped=unmapped))

            if f_unm:
                f_unm.close()

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
