import os
import sys
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


def parse_dnd_path(event_data):
    """Очистка пути к файлу при перетаскивании мышью (Drag & Drop)"""
    path = event_data.strip()
    if path.startswith('{') and path.endswith('}'):
        path = path[1:-1]
    return path


class TextRedirector:
    """Перенаправление логов консоли прямо в текстовое окно GUI"""
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
        self.root.title("LiftOverZ — Universal Genomic Liftover GUI")
        self.root.geometry("760x740")
        self.root.minsize(700, 650)

        # Переменные основных путей
        self.input_path = tk.StringVar()
        self.chain_path = tk.StringVar()
        self.fasta_path = tk.StringVar()
        self.out_path = tk.StringVar()
        self.selected_format = tk.StringVar(value=list(self.FORMAT_MAPPING.keys())[0])

        # Интерактивные переменные опций GUI
        self.chrom_map_path = tk.StringVar()
        self.min_mapq = tk.IntVar(value=0)
        self.min_insert_size = tk.IntVar(value=0)
        self.max_insert_size = tk.IntVar(value=1000)
        self.stdev_insert_size = tk.IntVar(value=300)
        self.append_tags = tk.BooleanVar(value=False)
        self.no_comp_as_chrom = tk.BooleanVar(value=False)

        self._build_ui()

    def _build_ui(self):
        main_frame = ttk.Padding(self.root, padding=12)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # 1. Формат данных
        fmt_frame = ttk.LabelFrame(main_frame, text=" 1. Выбор формата данных ", padding=8)
        fmt_frame.pack(fill=tk.X, pady=(0, 8))

        cb_fmt = ttk.Combobox(
            fmt_frame, 
            textvariable=self.selected_format, 
            values=list(self.FORMAT_MAPPING.keys()), 
            state="readonly"
        )
        cb_fmt.pack(fill=tk.X, expand=True)
        cb_fmt.bind("<<ComboboxSelected>>", self._on_format_changed)

        # 2. Выбор файлов (Drag & Drop)
        files_frame = ttk.LabelFrame(main_frame, text=" 2. Исходные файлы (перетащите файлы мышью) ", padding=8)
        files_frame.pack(fill=tk.X, pady=(0, 8))
        files_frame.columnconfigure(1, weight=1)

        # Входной файл
        ttk.Label(files_frame, text="Входной файл:").grid(row=0, column=0, sticky="w", pady=3)
        ent_inp = ttk.Entry(files_frame, textvariable=self.input_path)
        ent_inp.grid(row=0, column=1, sticky="ew", padx=5, pady=3)
        ttk.Button(files_frame, text="Обзор...", command=self._browse_input).grid(row=0, column=2, pady=3)
        ent_inp.drop_target_register(DND_FILES)
        ent_inp.dnd_bind('<<Drop>>', self._on_drop_input)

        # Chain файл
        ttk.Label(files_frame, text="Chain файл:").grid(row=1, column=0, sticky="w", pady=3)
        ent_chain = ttk.Entry(files_frame, textvariable=self.chain_path)
        ent_chain.grid(row=1, column=1, sticky="ew", padx=5, pady=3)
        ttk.Button(files_frame, text="Обзор...", command=self._browse_chain).grid(row=1, column=2, pady=3)
        ent_chain.drop_target_register(DND_FILES)
        ent_chain.dnd_bind('<<Drop>>', lambda e: self.chain_path.set(parse_dnd_path(e.data)))

        # Выходной файл + Диалог Сохранения
        ttk.Label(files_frame, text="Сохранить результат:").grid(row=2, column=0, sticky="w", pady=3)
        ent_out = ttk.Entry(files_frame, textvariable=self.out_path)
        ent_out.grid(row=2, column=1, sticky="ew", padx=5, pady=3)
        ttk.Button(files_frame, text="Сохранить как...", command=self._save_out_dialog).grid(row=2, column=2, pady=3)

        # 3. Динамическая панель опций (меняется под каждый формат)
        self.opts_frame = ttk.LabelFrame(main_frame, text=" 3. Настройки для выбранного формата ", padding=8)
        self.opts_frame.pack(fill=tk.X, pady=(0, 8))

        # Rebuild dynamic widgets
        self._rebuild_dynamic_options()

        # 4. Запуск
        self.btn_run = ttk.Button(main_frame, text="🚀 Запустить LiftOverZ", command=self._start_processing)
        self.btn_run.pack(fill=tk.X, pady=5, ipady=4)

        # 5. Окно логов
        log_frame = ttk.LabelFrame(main_frame, text=" Процесс выполнения ", padding=5)
        log_frame.pack(fill=tk.BOTH, expand=True, pady=(5, 0))

        self.log_area = tk.Text(log_frame, height=9, state='disabled', wrap='word')
        self.log_area.pack(fill=tk.BOTH, expand=True)

    def _rebuild_dynamic_options(self):
        """Очищает и перестраивает переключатели в зависимости от формата"""
        for child in self.opts_frame.winfo_children():
            child.destroy()

        cmd = self.FORMAT_MAPPING[self.selected_format.get()]
        self.opts_frame.columnconfigure(1, weight=1)

        row_idx = 0

        # Опция Target FASTA (Обязательна для VCF, gVCF, MAF)
        if cmd in ["vcf", "gvcf", "maf"]:
            ttk.Label(self.opts_frame, text="Target FASTA (Геном):").grid(row=row_idx, column=0, sticky="w", pady=3)
            ent_fa = ttk.Entry(self.opts_frame, textvariable=self.fasta_path)
            ent_fa.grid(row=row_idx, column=1, sticky="ew", padx=5, pady=3)
            ttk.Button(self.opts_frame, text="Обзор...", command=self._browse_fasta).grid(row=row_idx, column=2, pady=3)
            ent_fa.drop_target_register(DND_FILES)
            ent_fa.dnd_bind('<<Drop>>', lambda e: self.fasta_path.set(parse_dnd_path(e.data)))
            row_idx += 1

        # Опция файла маппинга хромосом (Для всех форматов)
        ttk.Label(self.opts_frame, text="Chrom Alias Map (необязательно):").grid(row=row_idx, column=0, sticky="w", pady=3)
        ent_cmap = ttk.Entry(self.opts_frame, textvariable=self.chrom_map_path)
        ent_cmap.grid(row=row_idx, column=1, sticky="ew", padx=5, pady=3)
        ttk.Button(self.opts_frame, text="Обзор...", command=self._browse_chrom_map).grid(row=row_idx, column=2, pady=3)
        ent_cmap.drop_target_register(DND_FILES)
        ent_cmap.dnd_bind('<<Drop>>', lambda e: self.chrom_map_path.set(parse_dnd_path(e.data)))
        row_idx += 1

        # Специфичные счетчики для BAM/SAM/CRAM
        if cmd == "bam":
            bam_panel = ttk.Frame(self.opts_frame)
            bam_panel.grid(row=row_idx, column=0, columnspan=3, sticky="ew", pady=5)

            ttk.Label(bam_panel, text="Min MAPQ:").grid(row=0, column=0, sticky="w", padx=(0, 2))
            ttk.Spinbox(bam_panel, from_=0, to=255, textvariable=self.min_mapq, width=5).grid(row=0, column=1, padx=(0, 10))

            ttk.Label(bam_panel, text="Min Insert Size:").grid(row=0, column=2, sticky="w", padx=(0, 2))
            ttk.Spinbox(bam_panel, from_=0, to=100000, textvariable=self.min_insert_size, width=7).grid(row=0, column=3, padx=(0, 10))

            ttk.Label(bam_panel, text="Max Insert Size:").grid(row=0, column=4, sticky="w", padx=(0, 2))
            ttk.Spinbox(bam_panel, from_=0, to=100000, textvariable=self.max_insert_size, width=7).grid(row=0, column=5, padx=(0, 10))

            ttk.Label(bam_panel, text="StDev:").grid(row=1, column=0, sticky="w", padx=(0, 2), pady=5)
            ttk.Spinbox(bam_panel, from_=0, to=10000, textvariable=self.stdev_insert_size, width=5).grid(row=1, column=1, padx=(0, 10), pady=5)

            ttk.Checkbutton(bam_panel, text="Добавлять теги сопряжения (-a / --append-tags)", variable=self.append_tags).grid(row=1, column=2, columnspan=4, sticky="w", pady=5)
            row_idx += 1

        # Специфичная галочка для VCF/gVCF
        if cmd in ["vcf", "gvcf"]:
            ttk.Checkbutton(
                self.opts_frame, 
                text="Отключить автоматическое сопоставление префиксов chr/non-chr (--no-comp-as-chrom)", 
                variable=self.no_comp_as_chrom
            ).grid(row=row_idx, column=0, columnspan=3, sticky="w", pady=4)

    # --- Обработчики GUI ---
    def _on_format_changed(self, event=None):
        self._rebuild_dynamic_options()

    def _on_drop_input(self, event):
        path = parse_dnd_path(event.data)
        self.input_path.set(path)
        self._auto_detect_format(path)

    def _browse_input(self):
        filename = filedialog.askopenfilename(
            title="Выберите входной файл",
            filetypes=[("Genomic files", "*.vcf *.vcf.gz *.gvcf *.gvcf.gz *.bam *.sam *.cram *.bed *.gff *.gtf *.wig *.bw *.bigwig *.maf"), ("All Files", "*.*")]
        )
        if filename:
            self.input_path.set(filename)
            self._auto_detect_format(filename)

    def _browse_chain(self):
        filename = filedialog.askopenfilename(
            title="Выберите Chain файл",
            filetypes=[("Chain Files", "*.chain *.chain.gz"), ("All Files", "*.*")]
        )
        if filename:
            self.chain_path.set(filename)

    def _browse_fasta(self):
        filename = filedialog.askopenfilename(
            title="Выберите Target Genome FASTA",
            filetypes=[("FASTA Files", "*.fa *.fasta *.fa.gz *.fasta.gz"), ("All Files", "*.*")]
        )
        if filename:
            self.fasta_path.set(filename)

    def _browse_chrom_map(self):
        filename = filedialog.askopenfilename(
            title="Выберите файл соответствия хромосом",
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
            title="Сохранить результат как...",
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
            messagebox.showerror("Ошибка", "Укажите существующий входной файл.")
            return
        if not chain or not os.path.exists(chain):
            messagebox.showerror("Ошибка", "Укажите существующий Chain файл.")
            return
        if cmd_type in ["vcf", "gvcf", "maf"] and (not fasta or not os.path.exists(fasta)):
            messagebox.showerror("Ошибка", f"Для формата {cmd_type.upper()} требуется выставить путь к FASTA файлу генома.")
            return
        if not out:
            messagebox.showerror("Ошибка", "Укажите путь для сохранения результата.")
            return

        self.btn_run.config(state='disabled')
        threading.Thread(target=self._run_crossmap, args=(cmd_type, chain, inp, fasta, out), daemon=True).start()

    def _run_crossmap(self, cmd_type, chain, inp, fasta, out):
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        redirector = TextRedirector(self.log_area)

        try:
            sys.stdout = redirector
            sys.stderr = redirector

            self.log(f"=== LiftOverZ: {cmd_type.upper()} ===")
            self.log(f"Входной файл: {inp}")
            self.log(f"Chain файл: {chain}")
            if cmd_type in ["vcf", "gvcf", "maf"]:
                self.log(f"Target FASTA: {fasta}")
            self.log(f"Результат: {out}\n")

            # Формирование аргументов на основе выставленных настроек GUI
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
                self.log("ОШИБКА: Движок CrossMap не найден!")

            self.log("\n=== Завершено успешно! ===")
            messagebox.showinfo("Успех", f"Конвертация {cmd_type.upper()} завершена!")

        except SystemExit as e:
            if e.code == 0:
                self.log("\n=== Завершено успешно! ===")
                messagebox.showinfo("Успех", "Операция успешно завершена!")
            else:
                self.log(f"\nОШИБКА: Процесс завершился с кодом {e.code}")
                messagebox.showerror("Ошибка", f"Код ошибки: {e.code}")
        except Exception as e:
            self.log(f"\nОШИБКА: {str(e)}")
            messagebox.showerror("Ошибка", str(e))
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr
            self.btn_run.config(state='normal')


if __name__ == "__main__":
    root = TkinterDnD.Tk()
    app = LiftOverZApp(root)
    root.mainloop()