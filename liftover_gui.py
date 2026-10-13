import os
import sys
import glob
import json
import shutil
import subprocess
import tempfile
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinterdnd2 import DND_FILES, TkinterDnD

try:
    from cm.main import crossmap_entry
except ImportError:
    try:
        from CrossMap.main import crossmap_entry
    except ImportError:
        crossmap_entry = None

try:
    import pybigtools
except ImportError:
    pybigtools = None


class I18nManager:
    """Менеджер локализации: загружает языки из PyInstaller bundle и внешней папки locales/"""
    def __init__(self, default_lang="en"):
        self.languages = {}
        self.current_lang = default_lang
        self.translations = {}
        self._load_available_languages()
        self.set_language(default_lang)

    def _get_locales_dirs(self):
        dirs = []
        # Внутренняя папка PyInstaller
        if hasattr(sys, '_MEIPASS'):
            dirs.append(os.path.join(sys._MEIPASS, "locales"))
        # Внешняя папка рядом со скриптом / бинарником
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
                            display_name = data.get("lang_name", code.upper())
                            self.languages[code] = {"name": display_name, "path": fpath}
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
        text = self.translations.get(key, key)
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


def get_bundled_executable(name):
    if hasattr(sys, '_MEIPASS'):
        bundle_path = os.path.join(sys._MEIPASS, name)
        if os.path.exists(bundle_path):
            return bundle_path
        if os.path.exists(bundle_path + ".exe"):
            return bundle_path + ".exe"
    return shutil.which(name)


class TextRedirector:
    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, str_val):
        self.text_widget.config(state='normal')
        self.text_widget.insert(tk.END, str_val)
        self.text_widget.see(tk.END)
        self.text_widget.config(state='disabled')

    def flush(self):
        pass


class LiftOverZApp:
    FORMAT_MAPPING = {
        "BED / BEDGraph (.bed, .bedgraph)": "bed",
        "BAM / SAM / CRAM (.bam, .sam, .cram)": "bam",
        "GFF / GTF (.gff, .gtf, .gff3)": "gff",
        "Wiggle (.wig)": "wig",
        "BigWig (.bw, .bigwig)": "bigwig",
        "VCF (.vcf, .vcf.gz)": "vcf",
        "gVCF (.gvcf, .gvcf.gz)": "gvcf",
        "MAF (.maf)": "maf"
    }

    def __init__(self, root):
        self.root = root
        self.i18n = I18nManager(default_lang="ru" if "ru" in I18nManager().languages else "en")

        self.root.geometry("770x760")
        self.root.minsize(720, 680)

        self.input_path = tk.StringVar()
        self.chain_path = tk.StringVar()
        self.fasta_path = tk.StringVar()
        self.out_path = tk.StringVar()
        self.selected_format = tk.StringVar(value=list(self.FORMAT_MAPPING.keys())[0])

        self.chrom_map_path = tk.StringVar()
        self.min_mapq = tk.IntVar(value=0)
        self.min_insert_size = tk.IntVar(value=0)
        self.max_insert_size = tk.IntVar(value=1000)
        self.stdev_insert_size = tk.IntVar(value=300)
        self.append_tags = tk.BooleanVar(value=False)
        self.no_comp_as_chrom = tk.BooleanVar(value=False)

        self._build_ui()
        self._update_ui_texts()

    def _build_ui(self):
        self.main_frame = ttk.Padding(self.root, padding=12)
        self.main_frame.pack(fill=tk.BOTH, expand=True)

        # Верхняя панель: Переключатель языка
        top_bar = ttk.Frame(self.main_frame)
        top_bar.pack(fill=tk.X, pady=(0, 6))

        self.lbl_lang = ttk.Label(top_bar, text=self.i18n.t("language_label"))
        self.lbl_lang.pack(side=tk.LEFT, padx=(0, 5))

        lang_options = [f"{info['name']} ({code})" for code, info in self.i18n.languages.items()]
        self.lang_cb = ttk.Combobox(top_bar, values=lang_options, state="readonly", width=18)
        current_name = self.i18n.languages.get(self.i18n.current_lang, {}).get("name", "English")
        self.lang_cb.set(f"{current_name} ({self.i18n.current_lang})")
        self.lang_cb.pack(side=tk.LEFT)
        self.lang_cb.bind("<<ComboboxSelected>>", self._on_language_changed)

        # 1. Формат
        self.fmt_frame = ttk.LabelFrame(self.main_frame, padding=8)
        self.fmt_frame.pack(fill=tk.X, pady=(0, 8))

        self.cb_fmt = ttk.Combobox(
            self.fmt_frame, 
            textvariable=self.selected_format, 
            values=list(self.FORMAT_MAPPING.keys()), 
            state="readonly"
        )
        self.cb_fmt.pack(fill=tk.X, expand=True)
        self.cb_fmt.bind("<<ComboboxSelected>>", self._on_format_changed)

        # 2. Файлы
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

        self.lbl_chain = ttk.Label(self.files_frame)
        self.lbl_chain.grid(row=1, column=0, sticky="w", pady=3)
        self.ent_chain = ttk.Entry(self.files_frame, textvariable=self.chain_path)
        self.ent_chain.grid(row=1, column=1, sticky="ew", padx=5, pady=3)
        self.btn_browse_chain = ttk.Button(self.files_frame, command=self._browse_chain)
        self.btn_browse_chain.grid(row=1, column=2, pady=3)
        self.ent_chain.drop_target_register(DND_FILES)
        self.ent_chain.dnd_bind('<<Drop>>', lambda e: self.chain_path.set(parse_dnd_path(e.data)))

        self.lbl_output = ttk.Label(self.files_frame)
        self.lbl_output.grid(row=2, column=0, sticky="w", pady=3)
        self.ent_out = ttk.Entry(self.files_frame, textvariable=self.out_path)
        self.ent_out.grid(row=2, column=1, sticky="ew", padx=5, pady=3)
        self.btn_save_out = ttk.Button(self.files_frame, command=self._save_out_dialog)
        self.btn_save_out.grid(row=2, column=2, pady=3)

        # 3. Динамические опции
        self.opts_frame = ttk.LabelFrame(self.main_frame, padding=8)
        self.opts_frame.pack(fill=tk.X, pady=(0, 8))
        self._rebuild_dynamic_options()

        # 4. Запуск
        self.btn_run = ttk.Button(self.main_frame, command=self._start_processing)
        self.btn_run.pack(fill=tk.X, pady=5, ipady=4)

        # 5. Лог
        self.log_frame = ttk.LabelFrame(self.main_frame, padding=5)
        self.log_frame.pack(fill=tk.BOTH, expand=True, pady=(5, 0))

        self.log_area = tk.Text(self.log_frame, height=8, state='disabled', wrap='word')
        self.log_area.pack(fill=tk.BOTH, expand=True)

    def _update_ui_texts(self):
        """Обновляет все заголовки и надписи в интерфейсе"""
        self.root.title(self.i18n.t("app_title"))
        self.lbl_lang.config(text=self.i18n.t("language_label"))
        self.fmt_frame.config(text=self.i18n.t("step_1_format"))
        self.files_frame.config(text=self.i18n.t("step_2_files"))
        self.opts_frame.config(text=self.i18n.t("step_3_options"))
        self.log_frame.config(text=self.i18n.t("log_title"))

        self.lbl_input.config(text=self.i18n.t("lbl_input"))
        self.lbl_chain.config(text=self.i18n.t("lbl_chain"))
        self.lbl_output.config(text=self.i18n.t("lbl_output"))

        self.btn_browse_inp.config(text=self.i18n.t("btn_browse"))
        self.btn_browse_chain.config(text=self.i18n.t("btn_browse"))
        self.btn_save_out.config(text=self.i18n.t("btn_save_as"))
        self.btn_run.config(text=self.i18n.t("btn_run"))

        self._rebuild_dynamic_options()

    def _on_language_changed(self, event=None):
        selected = self.lang_cb.get()
        # Извлекаем код из скобок, напр. "Русский (ru)" -> "ru"
        code = selected.split('(')[-1].rstrip(')').strip()
        self.i18n.set_language(code)
        self._update_ui_texts()

    def _rebuild_dynamic_options(self):
        for child in self.opts_frame.winfo_children():
            child.destroy()

        cmd = self.FORMAT_MAPPING[self.selected_format.get()]
        self.opts_frame.columnconfigure(1, weight=1)
        row_idx = 0

        if cmd in ["vcf", "gvcf", "maf"]:
            ttk.Label(self.opts_frame, text=self.i18n.t("lbl_fasta_req")).grid(row=row_idx, column=0, sticky="w", pady=3)
            ent_fa = ttk.Entry(self.opts_frame, textvariable=self.fasta_path)
            ent_fa.grid(row=row_idx, column=1, sticky="ew", padx=5, pady=3)
            ttk.Button(self.opts_frame, text=self.i18n.t("btn_browse"), command=self._browse_fasta).grid(row=row_idx, column=2, pady=3)
            ent_fa.drop_target_register(DND_FILES)
            ent_fa.dnd_bind('<<Drop>>', lambda e: self.fasta_path.set(parse_dnd_path(e.data)))
            row_idx += 1

        ttk.Label(self.opts_frame, text=self.i18n.t("lbl_chrom_map")).grid(row=row_idx, column=0, sticky="w", pady=3)
        ent_cmap = ttk.Entry(self.opts_frame, textvariable=self.chrom_map_path)
        ent_cmap.grid(row=row_idx, column=1, sticky="ew", padx=5, pady=3)
        ttk.Button(self.opts_frame, text=self.i18n.t("btn_browse"), command=self._browse_chrom_map).grid(row=row_idx, column=2, pady=3)
        ent_cmap.drop_target_register(DND_FILES)
        ent_cmap.dnd_bind('<<Drop>>', lambda e: self.chrom_map_path.set(parse_dnd_path(e.data)))
        row_idx += 1

        if cmd == "bam":
            bam_panel = ttk.Frame(self.opts_frame)
            bam_panel.grid(row=row_idx, column=0, columnspan=3, sticky="ew", pady=5)

            ttk.Label(bam_panel, text=self.i18n.t("lbl_mapq")).grid(row=0, column=0, sticky="w", padx=(0, 2))
            ttk.Spinbox(bam_panel, from_=0, to=255, textvariable=self.min_mapq, width=5).grid(row=0, column=1, padx=(0, 10))

            ttk.Label(bam_panel, text=self.i18n.t("lbl_min_insert")).grid(row=0, column=2, sticky="w", padx=(0, 2))
            ttk.Spinbox(bam_panel, from_=0, to=100000, textvariable=self.min_insert_size, width=7).grid(row=0, column=3, padx=(0, 10))

            ttk.Label(bam_panel, text=self.i18n.t("lbl_max_insert")).grid(row=0, column=4, sticky="w", padx=(0, 2))
            ttk.Spinbox(bam_panel, from_=0, to=100000, textvariable=self.max_insert_size, width=7).grid(row=0, column=5, padx=(0, 10))

            ttk.Label(bam_panel, text=self.i18n.t("lbl_stdev")).grid(row=1, column=0, sticky="w", padx=(0, 2), pady=5)
            ttk.Spinbox(bam_panel, from_=0, to=10000, textvariable=self.stdev_insert_size, width=5).grid(row=1, column=1, padx=(0, 10), pady=5)

            ttk.Checkbutton(bam_panel, text=self.i18n.t("chk_append_tags"), variable=self.append_tags).grid(row=1, column=2, columnspan=4, sticky="w", pady=5)
            row_idx += 1

        if cmd in ["vcf", "gvcf"]:
            ttk.Checkbutton(
                self.opts_frame, 
                text=self.i18n.t("chk_no_comp"), 
                variable=self.no_comp_as_chrom
            ).grid(row=row_idx, column=0, columnspan=3, sticky="w", pady=4)

    def _on_format_changed(self, event=None):
        self._rebuild_dynamic_options()

    def _on_drop_input(self, event):
        path = parse_dnd_path(event.data)
        self.input_path.set(path)
        self._auto_detect_format(path)

    def _browse_input(self):
        filename = filedialog.askopenfilename(
            title=self.i18n.t("dialog_select_input"),
            filetypes=[("Genomic files", "*.vcf *.vcf.gz *.gvcf *.gvcf.gz *.bam *.sam *.cram *.bed *.gff *.gtf *.wig *.bw *.bigwig *.maf"), ("All Files", "*.*")]
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

    def _browse_fasta(self):
        filename = filedialog.askopenfilename(
            title=self.i18n.t("dialog_select_fasta"),
            filetypes=[("FASTA Files", "*.fa *.fasta *.fa.gz *.fasta.gz"), ("All Files", "*.*")]
        )
        if filename:
            self.fasta_path.set(filename)

    def _browse_chrom_map(self):
        filename = filedialog.askopenfilename(
            title=self.i18n.t("dialog_select_map"),
            filetypes=[("Mapping Files", "*.txt *.map *.tsv"), ("All Files", "*.*")]
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
        elif ext.endswith(('.bam', '.sam', '.cram')):
            self.selected_format.set("BAM / SAM / CRAM (.bam, .sam, .cram)")
        elif ext.endswith(('.bed', '.bedgraph')):
            self.selected_format.set("BED / BEDGraph (.bed, .bedgraph)")
        elif ext.endswith(('.gff', '.gff3', '.gtf')):
            self.selected_format.set("GFF / GTF (.gff, .gtf, .gff3)")
        elif ext.endswith('.wig'):
            self.selected_format.set("Wiggle (.wig)")
        elif ext.endswith(('.bw', '.bigwig')):
            self.selected_format.set("BigWig (.bw, .bigwig)")
        elif ext.endswith('.maf'):
            self.selected_format.set("MAF (.maf)")

        if not self.out_path.get():
            base, ext_name = os.path.splitext(filename)
            if filename.endswith('.gz'):
                base = base[:-4]
            self.out_path.set(f"{base}_lifted{ext_name}")

        self._rebuild_dynamic_options()

    def log(self, text):
        self.log_area.config(state='normal')
        self.log_area.insert(tk.END, text + "\n")
        self.log_area.see(tk.END)
        self.log_area.config(state='disabled')

    def _start_processing(self):
        inp = self.input_path.get()
        chain = self.chain_path.get()
        out = self.out_path.get()
        fasta = self.fasta_path.get()
        cmd_type = self.FORMAT_MAPPING[self.selected_format.get()]

        if not inp or not os.path.exists(inp):
            messagebox.showerror(self.i18n.t("err_title"), self.i18n.t("err_no_input"))
            return
        if not chain or not os.path.exists(chain):
            messagebox.showerror(self.i18n.t("err_title"), self.i18n.t("err_no_chain"))
            return
        if cmd_type in ["vcf", "gvcf", "maf"] and (not fasta or not os.path.exists(fasta)):
            messagebox.showerror(self.i18n.t("err_title"), self.i18n.t("err_no_fasta", fmt=cmd_type.upper()))
            return
        if not out:
            messagebox.showerror(self.i18n.t("err_title"), self.i18n.t("err_no_output"))
            return

        self.btn_run.config(state='disabled')
        threading.Thread(target=self._run_pipeline, args=(cmd_type, chain, inp, fasta, out), daemon=True).start()

    def _run_pipeline(self, cmd_type, chain, inp, fasta, out):
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        redirector = TextRedirector(self.log_area)

        try:
            sys.stdout = redirector
            sys.stderr = redirector

            self.log(self.i18n.t("log_started", fmt=cmd_type.upper()))
            self.log(self.i18n.t("log_input", path=inp))
            self.log(self.i18n.t("log_chain", path=chain))
            if cmd_type in ["vcf", "gvcf", "maf"]:
                self.log(self.i18n.t("log_fasta", path=fasta))
            self.log(self.i18n.t("log_output", path=out) + "\n")

            # BigWig через pybigtools
            if cmd_type == 'bigwig':
                self.log(self.i18n.t("log_converting_bigwig"))
                with tempfile.TemporaryDirectory() as tmpdir:
                    tmp_wig = os.path.join(tmpdir, "temp.wig")
                    tmp_lifted_wig = os.path.join(tmpdir, "temp_lifted.wig")
                    
                    reader = pybigtools.open(inp)
                    with open(tmp_wig, 'w') as f_out:
                        for chrom, length in reader.chroms().items():
                            for interval in reader.intervals(chrom):
                                f_out.write(f"{chrom}\t{interval[0]}\t{interval[1]}\t{interval[2]}\n")

                    sys.argv = ["CrossMap", "bed", chain, tmp_wig, tmp_lifted_wig]
                    crossmap_entry()

                    self.log(self.i18n.t("log_packing_bigwig"))
                    pybigtools.write(out, tmp_lifted_wig)

                self.log("\n" + self.i18n.t("log_finished"))
                messagebox.showinfo(self.i18n.t("success_title"), self.i18n.t("success_message", fmt="BigWig"))
                return

            # BAM/CRAM на Windows через samtools
            if cmd_type == 'bam' and sys.platform.startswith('win'):
                samtools_bin = get_bundled_executable("samtools")
                if samtools_bin and os.path.exists(samtools_bin):
                    with tempfile.TemporaryDirectory() as tmpdir:
                        tmp_sam = os.path.join(tmpdir, "stream.sam")
                        tmp_lifted_sam = os.path.join(tmpdir, "stream_lifted.sam")

                        self.log(self.i18n.t("log_samtools_stream"))
                        subprocess.run([samtools_bin, "view", "-h", "-o", tmp_sam, inp], check=True)

                        self.log(self.i18n.t("log_lifting_coords"))
                        cli_args = ["CrossMap", "bam"]
                        if self.chrom_map_path.get():
                            cli_args.extend(["-c", self.chrom_map_path.get()])
                        cli_args.extend(["--mapq", str(self.min_mapq.get())])
                        cli_args.extend(["-m", str(self.min_insert_size.get())])
                        cli_args.extend(["-t", str(self.max_insert_size.get())])
                        cli_args.extend(["-s", str(self.stdev_insert_size.get())])
                        if self.append_tags.get():
                            cli_args.append("-a")
                        cli_args.extend([chain, tmp_sam, tmp_lifted_sam])

                        sys.argv = cli_args
                        crossmap_entry()

                        self.log(self.i18n.t("log_indexing_bam"))
                        subprocess.run([samtools_bin, "view", "-b", "-o", out, tmp_lifted_sam], check=True)
                        subprocess.run([samtools_bin, "index", out], check=False)

                    self.log("\n" + self.i18n.t("log_finished"))
                    messagebox.showinfo(self.i18n.t("success_title"), self.i18n.t("success_message", fmt="BAM"))
                    return

            # Стандартная обработка остальных форматов
            cli_args = ["CrossMap", cmd_type]
            if self.chrom_map_path.get():
                cli_args.extend(["-c", self.chrom_map_path.get()])

            if cmd_type == 'bam':
                cli_args.extend(["--mapq", str(self.min_mapq.get())])
                cli_args.extend(["-m", str(self.min_insert_size.get())])
                cli_args.extend(["-t", str(self.max_insert_size.get())])
                cli_args.extend(["-s", str(self.stdev_insert_size.get())])
                if self.append_tags.get():
                    cli_args.append("-a")

            if cmd_type in ['vcf', 'gvcf'] and self.no_comp_as_chrom.get():
                cli_args.append("--no-comp-as-chrom")

            cli_args.extend([chain, inp])
            if cmd_type in ["vcf", "gvcf", "maf"]:
                cli_args.append(fasta)
            cli_args.append(out)

            sys.argv = cli_args

            if crossmap_entry is not None:
                crossmap_entry()
            else:
                self.log(self.i18n.t("err_crossmap_not_found"))

            self.log("\n" + self.i18n.t("log_finished"))
            messagebox.showinfo(self.i18n.t("success_title"), self.i18n.t("success_message", fmt=cmd_type.upper()))

        except SystemExit as e:
            if e.code == 0:
                self.log("\n" + self.i18n.t("log_finished"))
                messagebox.showinfo(self.i18n.t("success_title"), self.i18n.t("success_message", fmt=cmd_type.upper()))
            else:
                self.log(f"\n{self.i18n.t('err_process_failed', code=e.code)}")
                messagebox.showerror(self.i18n.t("err_title"), self.i18n.t("err_process_failed", code=e.code))
        except Exception as e:
            self.log(f"\nERROR: {str(e)}")
            messagebox.showerror(self.i18n.t("err_title"), str(e))
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr
            self.btn_run.config(state='normal')


if __name__ == "__main__":
    root = TkinterDnD.Tk()
    app = LiftOverZApp(root)
    root.mainloop()
