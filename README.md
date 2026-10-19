# LiftOverZ - A fast, standalone, cross-platform graphical tool for lifting over genomic coordinates across genome assembly builds.

![Python](https://img.shields.io/badge/Language-Python-blue.svg)
![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey.svg)
![License](https://img.shields.io/badge/License-MIT-orange.svg)

**Languages: [English](#english) | [Русский](#русский)**

---

<a name="english"></a>
## English

**LiftOverZ** is a fast, standalone, cross-platform graphical tool designed for lifting over genomic coordinates across genome assembly builds (e.g., hg19/GRCh37, hg38/GRCh38, T2T-CHM13v2, and mouse builds).

It requires **no installation**, runs as a single portable executable, includes **pre-bundled, gap-corrected UCSC chain files**, supports **Drag & Drop**, handles compressed files (`.gz`) out of the box, preserves chromosome naming styles (UCSC vs. Ensembl/NCBI), supports **optional target FASTA reference validation**, and features hot-swappable multilingual support.

---

## Screenshot

![Screenshot LiftOverZ](Screenshot.png)

---

### Key Features

* **Zero Installation / 100% Portable:** Single executable running out of the box on Windows 11/10, Linux, and macOS.
* **Native Cross-Platform Binaries:**
  * Windows (`LiftOverZ-windows.exe`)
  * Linux x86_64 (`LiftOverZ-linux`)
  * macOS Apple Silicon M1/M2/M3/M4 (`LiftOverZ-macos-arm64`)
  * macOS Intel x86_64 (`LiftOverZ-macos-intel`)
* **Pre-bundled & Gap-Corrected Chain Presets:** Popular human and mouse assembly chains are embedded directly into the binary. Crucially, human mitochondrial chains (`chrM`) are patched to eliminate official UCSC alignment gaps.
* **Chromosome Naming & Header Integrity:** Automatically preserves chromosome prefix style (e.g. Ensembl `1`, `MT` stay `1`, `MT` rather than forcing UCSC `chr1`, `chrM`). In VCF output the `##contig` lines are rebuilt for the **target** assembly (lengths come from the target FASTA, or from the chain file when no FASTA is given) and the stale `##reference` line is dropped, so the header always matches the records.
* **Automatic Coordinate Sorting (bounded memory):** Re-orders output entries by chromosome and genomic coordinate so that results can be indexed with `tabix` or loaded into genome browsers without external sorting tools. Very large inputs are sorted with an external merge sort (temporary files), so whole-genome VCFs do not have to fit into RAM.
* **Optional Target FASTA Reference Validation:** Provide a target reference genome FASTA file (`pyfaidx`) to enable full biological validation for VCF files:
  * **Strand Flips:** On inverted (`-`) alignment blocks `REF`/`ALT` are reverse-complemented (IUPAC-aware), multi-base alleles are anchored at the correct end, indels are re-anchored to the base on their left using the target FASTA, and records are tagged `LIFTOVER_REVERSE_COMPLEMENTED`.
  * **REF/ALT Swapping:** When the reference allele in the new build switches to the variant allele, alleles are swapped (`REF` $\leftrightarrow$ `ALT`); diploid **and haploid** genotypes (`GT`), `PL`/`GL`, `AD`/`ADF`/`ADR` and `FORMAT/AF` are recalculated, INFO `AF`/`MLEAF` become `1-x`, `AC`/`MLEAC` become `AN-x`, and entries are flagged with `LIFTOVER_SWAPPED`. Other custom INFO fields (e.g. frequencies copied from external databases) cannot be recalculated and are left unchanged.
  * **Mismatch Audit:** True reference mismatches are cleanly diverted into a `.rejected.vcf` file with detailed annotations (`LIFTOVER_REJECT_REASON=MismatchedRefAllele`).
* **Intuitive GUI with Drag & Drop & Reset:** Drop input files, chain files, target FASTA, or mapping tables straight into the interface. Includes a one-click Reset button to safely process multiple files consecutively.
* **Broad Genomic Format Support:** Automatically detects and converts 6 major genomic formats (VCF, gVCF, BED, GFF/GTF, Wiggle, MAF).
* **Transparent Gzip / BGZF Compression:** Reads `.gz` / `.bgz` (detected by file signature, not by extension) without pre-decompression. Output files ending in `.gz` or `.bgz` are written as true **BGZF** (block-gzip, built in, no extra dependencies), so they can be indexed directly with `tabix -p vcf out.vcf.gz` or `bcftools index`.
* **Unmapped Records Audit:** Exports unmapped records (including malformed lines and records whose `REF` crosses an alignment gap) into a dedicated `.unmapped` file for downstream verification.
* **Extensible Internationalization (i18n):** Shipped with 7 languages (English, Russian, Latvian, German, Spanish, Lithuanian, Estonian) with automatic OS language detection.

---

### Pre-bundled Chains & Mitochondrial (chrM) Alignment Fixes

LiftOverZ embeds the most common chain files directly into the executable:
* **Human:**
  * hg19 (GRCh37) ⇄ hg38 (GRCh38)
  * hg38 (GRCh38) ⇄ T2T-CHM13v2.0 (hs1)
  * hg19 (GRCh37) ⇄ T2T-CHM13v2.0 (hs1)
* **Mouse:**
  * mm10 (GRCm38) ⇄ mm39 (GRCm39)
  * mm9 (NCBI37) → mm10 (GRCm38)

#### The UCSC Mitochondrial Gap Issue & How LiftOverZ Fixes It
Official UCSC chain files frequently fail to map variants on mitochondrial DNA (`chrM` / `MT`) due to automated whole-genome alignment artifacts:
1. **D-loop / Hypervariable Regions:** High variation in the control region causes alignment algorithms to break chains into fragmented blocks with gaps (e.g., dropping positions 150, 263, 310).
2. **Nuclear Mitochondrial Pseudogenes (NUMTs):** Homology between mtDNA and nuclear pseudogenes causes the aligner to falsely branch or penalize mitochondrial alignments (e.g., dropping positions 11719 and 15326).
3. **The hg19 2-bp InDel (NC_001807 vs. rCRS):** UCSC hg19 contains an artificial 2-bp `NN` insertion near position 3106 (chain offset 3106, 0-based; 16,571 bp total) compared to the standard rCRS sequence used in hg38 and T2T-CHM13v2 (16,569 bp total).

**LiftOverZ bundles custom, gap-free mitochondrial alignments based on Broad Institute / GATK standards:**
* **hg38 / hs1 ⇄ hg19:** A clean, two-block continuous alignment: the first 3106 bases (0-based offsets 0–3105) map 1:1 and everything after them maps with a +2 / −2 bp shift, without any gaps.
* **hg38 ⇄ hs1 (T2T):** A pure 1:1 identity alignment across all 16,569 bp, removing fragmented UCSC chain splits.

---

### Dual Mode: Fast Mathematical Mapping vs. Full FASTA Validation

LiftOverZ provides two flexible operating modes for VCF liftover:

1. **Fast Mode (Target FASTA left empty):**
   * Instantaneous mathematical coordinate translation and automatic chromosome sorting.
   * Does not require downloading multi-gigabyte reference genomes; memory use stays bounded thanks to external sorting.
   * Ideal for targeted sequencing, SNP arrays, mitochondrial analyses, and non-critical research queries.
   * On inverted (`-`) alignment blocks, SNVs and equal-length MNVs are reverse-complemented automatically. Length-changing variants (indels) on inverted blocks cannot be represented correctly without the target sequence, so they are written to `.unmapped` and counted in the log — use Strict Mode for them.
2. **Strict Mode (Target FASTA provided):**
   * Rapidly queries target nucleotide sequence using indexed FASTA (`pyfaidx`).
   * Validates target `REF` sequence against the new reference build.
   * Handles strand inversion (reverse complement, including re-anchoring of indels) and allele flips (REF/ALT swap) for clinical and research rigor.
   * Discards incompatible variants into a valid, annotated `.rejected.vcf` file.

---

### Supported File Formats

| Format | Extensions | Coordinate System | Notes |
| :--- | :--- | :--- | :--- |
| **VCF** | `.vcf`, `.vcf.gz`, `.vcf.bgz` | 1-based | Retains variant metadata and contig prefix style; rebuilds `##contig` for the target build; lifts `INFO/END`; supports FASTA validation |
| **gVCF** | `.gvcf`, `.g.vcf`, `.gvcf.gz`, `.g.vcf.gz` | 1-based | Reference blocks: `POS` **and** `INFO/END` are both lifted |
| **BED / BEDGraph** | `.bed`, `.bedgraph` (+ `.gz`) | 0-based start, 1-based end | Updates interval spans; flips the strand column on inverted blocks; converts BED12 `thickStart`/`thickEnd`/`blockSizes`/`blockStarts` |
| **GFF / GTF / GFF3** | `.gff`, `.gtf`, `.gff3` (+ `.gz`) | 1-based inclusive | Updates feature bounds; flips the strand column on inverted blocks; keeps a trailing `##FASTA` section intact |
| **Wiggle** | `.wig` (+ `.gz`) | 1-based | Supports `variableStep` and `fixedStep`; output is always `variableStep` because lifted positions are no longer evenly spaced; multiple tracks stay separate |
| **MAF** | `.maf` (+ `.gz`) | 1-based | Mutation Annotation Format: updates `Start_Position` / `End_Position`; reverse-complements allele columns on inverted blocks (`NCBI_Build` is not modified) |

---

### How Records Are Handled

* **Chain gaps:** a VCF record whose `REF` span crosses an alignment gap is not lifted (it goes to `.unmapped`), because its allele lengths would no longer be valid. BED / GFF / MAF / Wiggle intervals are lifted when both ends map to the same chromosome and strand.
* **Inverted blocks:** coordinates, strand columns and alleles are flipped as described above. VCF breakend (`BND`) alleles on inverted blocks are sent to `.unmapped`.
* **Malformed and short lines:** unparsable lines are written to `.unmapped`; VCF lines with fewer than 8 columns are padded with `.`; blank lines are ignored. Therefore `Total = Mapped + Unmapped + Rejected`.
* **Line endings:** output always uses `LF`, also on Windows (CRLF breaks several VCF tools).
* **Failures:** if a run fails, partially written output files are removed instead of being left behind. The output file may not be the same as the input file.

---

### Quick Start

1. Go to the **[Releases](../../releases)** tab of this repository.
2. Download the binary for your platform.
3. Launch the application.
4. Drag & drop your genomic data file into **Input File**.
5. Select the desired conversion from the **Built-in Chain Preset** dropdown (or drop your own custom chain file).
6. *(Optional)* Select or drop a target reference FASTA (e.g., `hg38.fa`) into **Target Reference FASTA** for biological REF validation.
7. Confirm the destination in **Save Result To**.
8. Click **🚀 Run LiftOverZ**.
9. Use **🔄 Reset** to clear paths and state when processing consecutive files.

---

### Running the Downloaded Binaries

* **Linux / macOS:** GitHub release assets do not keep the executable bit — run `chmod +x LiftOverZ-linux` (or `LiftOverZ-macos-*`) once.
* **macOS:** the binaries are not code-signed or notarized. On first launch use right-click → *Open*, or run `xattr -d com.apple.quarantine LiftOverZ-macos-*`.
* **Linux:** the binary is built on Ubuntu 24.04 and needs a similarly recent glibc; on older distributions run from source (`python liftover_gui.py`).

---

### Options & Configuration

* **Target Reference FASTA (Optional):** Path to target assembly FASTA file (`.fa`, `.fasta`, `.fa.gz`). Automatically generates a `.fai` index on first run if missing (the folder must be writable). Compressed FASTA must be **bgzip**-compressed — plain gzip is not supported by `pyfaidx`.
* **Auto-match 'chr' prefix:** When checked, reconciles naming variations between NCBI/Ensembl (`1`, `2`, `MT`) and UCSC (`chr1`, `chr2`, `chrM`) while retaining original input prefix format in the output.
* **Save unmapped records:** Generates a `.unmapped` file containing lines that could not be lifted over.
* **Chromosome Mapping File (Optional):** Whitespace-delimited text file mapping custom contig names (lines starting with `#` are ignored). Contigs listed here are written to the output with the target-build name from the chain file (e.g. `NC_000001.10` → `chr1`):
```text
NC_000001.10   chr1
NC_000002.11   chr2
MT             chrM
```

---

### Adding Custom Translations

1. Create a `locales/` directory next to the executable.
2. Copy any existing JSON file (e.g. `en.json`) to your language code (e.g., `locales/fr.json`, `locales/it.json`).
3. Translate the values in any text editor.
4. Restart the app — the new language will be automatically loaded into the language selector.

---

### Local Build & Development

```bash
# Clone the repository
git clone https://github.com/zajakin/LiftOverZ.git
cd LiftOverZ

# Install dependencies
pip install -r requirements.txt

# Run directly
python liftover_gui.py

# Build a standalone binary
pyinstaller --onefile --noconsole --add-data "locales:locales" --add-data "chains:chains" --collect-all tkinterdnd2 --collect-all pyliftover --collect-all pyfaidx --name LiftOverZ liftover_gui.py
```
*(On Windows, replace `:locales` and `:chains` with `;locales` and `;chains` in `--add-data`)*

---

<a name="русский"></a>
## Русский

**LiftOverZ** — это быстрый автономный кроссплатформенный графический инструмент для переноса (liftover) геномных координат между сборками референсного генома (hg19/GRCh37, hg38/GRCh38, T2T-CHM13v2, а также сборками генома мыши).

Программа **не требует установки**, поставляется в виде одного портативного исполняемого файла, содержит **встроенные исправленные chain-файлы UCSC (без пропусков в митохондрии)**, поддерживает **Drag & Drop**, напрямую работает со сжатыми архивами (`.gz`), сохраняет нотацию хромосом (Ensembl против UCSC), автоматически сортирует координаты, поддерживает **опциональную строгую валидацию по целевой FASTA** и включает 7 языков интерфейса с автоопределением языка ОС.

---

## Скриншот

![Screenshot LiftOverZ](Screenshot.png)

---

### Основные возможности

* **Полная портативность:** Один исполняемый файл без внешних зависимостей. Запускается без инсталлятора на Windows 11/10, Linux и macOS.
* **Нативные бинарники под все платформы:**
  * Windows (`LiftOverZ-windows.exe`)
  * Linux x86_64 (`LiftOverZ-linux`)
  * macOS на процессорах Apple Silicon M1/M2/M3/M4 (`LiftOverZ-macos-arm64`)
  * macOS на процессорах Intel (`LiftOverZ-macos-intel`)
* **Встроенные и исправленные цепочки UCSC:** Все популярные цепочки для человека и мыши уже включены в программу. Цепочки человеческой митохондрии (`chrM`) модифицированы для устранения официальных дефектов выравнивания UCSC.
* **Сохранение стиля префиксов и валидность заголовков VCF:** Если исходный файл использует нотацию Ensembl/NCBI (`1`, `2`, `MT`), результат сохраняет её без принудительного добавления `chr`. Строки `##contig` в выходном VCF пересобираются для **целевой** сборки (длины берутся из целевой FASTA, а без неё — из chain-файла), устаревшая строка `##reference` удаляется, поэтому заголовок всегда соответствует записям.
* **Автоматическая сортировка координат (ограниченная память):** Записи в выходном файле упорядочиваются по хромосоме и координатам, поэтому результат можно сразу индексировать с помощью `tabix` или загружать в IGV без внешних утилит сортировки. Очень большие файлы сортируются внешней сортировкой слиянием (временные файлы), поэтому VCF полного генома не обязан помещаться в ОЗУ.
* **Опциональная биологическая валидация по целевой FASTA:** Подключение референсного файла новой сборки (`pyfaidx`) позволяет выполнять полную валидацию VCF:
  * **Инверсия цепей (Strand Flip):** На инвертированных блоках (`-`) аллели `REF`/`ALT` переводятся в обратный комплемент (с учётом кодов IUPAC), многобазовые аллели привязываются к правильному концу, инделы заново якорятся на основание слева по целевой FASTA, а записи получают тег `LIFTOVER_REVERSE_COMPLEMENTED`.
  * **Замена аллелей (REF/ALT Swap):** Если референс в новой сборке сменился на мутантный аллель, программа производит замену `REF` $\leftrightarrow$ `ALT`; пересчитываются диплоидные **и гаплоидные** генотипы `GT`, поля `PL`/`GL`, `AD`/`ADF`/`ADR` и `FORMAT/AF`, в INFO `AF`/`MLEAF` становятся `1-x`, а `AC`/`MLEAC` — `AN-x`; запись маркируется тегом `LIFTOVER_SWAPPED`. Прочие пользовательские поля INFO (например, частоты из внешних баз) пересчитать невозможно, они остаются без изменений.
  * **Отбраковка несовпадений:** Строки, где референс не сошёлся с новой сборкой, выводятся в отдельный валидный файл `.rejected.vcf` с указанием причины `LIFTOVER_REJECT_REASON=MismatchedRefAllele`.
* **Интерфейс с Drag & Drop и кнопкой Reset:** Перетаскивайте входные файлы, сторонние chain-файлы, целевую FASTA и карты хромосом мышью. Кнопка «Сброс» позволяет быстро очистить состояние при пакетной обработке файлов.
* **Поддержка 6 основных геномных форматов:** Автоматическое распознавание VCF, gVCF, BED, GFF/GTF, Wiggle, MAF.
* **Прозрачная работа с Gzip / BGZF:** Чтение `.gz` / `.bgz` (определяется по сигнатуре файла, а не по расширению) без предварительной распаковки. Выходные файлы с расширением `.gz` или `.bgz` записываются как настоящий **BGZF** (блочный gzip, реализован встроенно, без дополнительных зависимостей), поэтому их можно сразу индексировать командами `tabix -p vcf out.vcf.gz` или `bcftools index`.
* **Аудит неперенесённых строк:** Запись пропущенных строк (включая некорректные и те, у которых `REF` пересекает разрыв выравнивания) в файл `.unmapped` для последующего контроля.
* **Мультиязычность (i18n):** Встроенная поддержка 7 языков (английский, русский, латышский, немецкий, испанский, литовский, эстонский) с автоопределением языка системы.

---

### Встроенные цепочки и исправление митохондриальных разрывов (chrM)

В состав LiftOverZ включены готовые цепочки UCSC:
* **Человек:**
  * hg19 (GRCh37) ⇄ hg38 (GRCh38)
  * hg38 (GRCh38) ⇄ T2T-CHM13v2.0 (hs1)
  * hg19 (GRCh37) ⇄ T2T-CHM13v2.0 (hs1)
* **Мышь:**
  * mm10 (GRCm38) ⇄ mm39 (GRCm39)
  * mm9 (NCBI37) → mm10 (GRCm38)

#### Проблема официальных файлов UCSC и её решение в LiftOverZ
В стандартных chain-файлах UCSC многие митохондриальные варианты (`chrM` / `MT`) ошибочно отбрасываются в `.unmapped` из-за артефактов автоматического полногеномного выравнивания:
1. **Гипервариабельные участки D-петли:** Из-за высокой плотности мутаций и поли-C тракта цепочки выравнивания разрываются на куски (выпадают маркерные позиции гаплогрупп: 150, 263, 310).
2. **Ядерные псевдогены митохондрии (NUMTs):** Из-за высокой гомологии участков мтДНК с псевдогенами на хромосомах 1 и 2 алгоритм `lastz` дробит цепочки митохондрии (выпадают позиции 11719 и 15326).
3. **Разница сборок hg19 и hg38:** В UCSC hg19 (NC_001807) рядом с позицией 3106 (смещение в chain-файле 3106, 0-based) присутствуют 2 фиктивных нуклеотида `NN` (длина 16 571 п.н.), тогда как hg38 и T2T hs1 используют стандарт rCRS (16 569 п.н.).

**Во встроенных цепочках LiftOverZ митохондрия заменена на непрерывные блоки по стандартам Broad Institute / GATK:**
* **hg38 / hs1 ⇄ hg19:** Непрерывное выравнивание из двух блоков: первые 3106 оснований (смещения 0–3105, 0-based) совпадают 1-в-1, далее — строгий сдвиг на +2 / −2 п.н. без единого разрыва.
* **hg38 ⇄ hs1 (T2T):** Строгое непрерывное совпадение 1-в-1 по всей длине 16 569 п.н.

---

### Два режима работы: Быстрый перенос vs. Строгая валидация по FASTA

1. **Быстрый режим (поле Target FASTA оставлено пустым):**
   * Мгновенный перенос координат по цепочкам и автоматическая сортировка.
   * Не требует загрузки многогигабайтных референсов и сторонних баз; потребление памяти ограничено благодаря внешней сортировке.
   * Оптимален для таргетных панелей, микрочипов, митохондрий и рутинных исследовательских задач.
   * На инвертированных блоках (`-`) SNV и MNV одинаковой длины автоматически переводятся в обратный комплемент. Варианты с изменением длины (инделы) на таких блоках нельзя корректно записать без целевой последовательности, поэтому они попадают в `.unmapped` и учитываются в логе — для них используйте строгий режим.
2. **Строгий режим (указан файл Target FASTA):**
   * Мгновенное чтение нуклеотидов новой сборки по индексу FASTA (`pyfaidx`).
   * Полная проверка совпадения референсного нуклеотида новой сборки с полем `REF`.
   * Корректная обработка инверсий (обратный комплемент, включая повторное якорение инделов) и смены референса (REF/ALT Swap).
   * Вывод биологически несовместимых вариантов в файл `.rejected.vcf`.

---

### Поддерживаемые форматы данных

| Формат | Расширения | Система координат | Примечание |
| :--- | :--- | :--- | :--- |
| **VCF** | `.vcf`, `.vcf.gz`, `.vcf.bgz` | 1-based | Сохраняет аннотации и стиль префиксов; пересобирает `##contig` для целевой сборки; переносит `INFO/END`; поддерживает проверку по FASTA |
| **gVCF** | `.gvcf`, `.g.vcf`, `.gvcf.gz`, `.g.vcf.gz` | 1-based | Невариантные блоки: переносятся и `POS`, **и** `INFO/END` |
| **BED / BEDGraph** | `.bed`, `.bedgraph` (+ `.gz`) | 0-based старт, 1-based конец | Обновляет границы интервалов; инвертирует столбец strand на инвертированных блоках; пересчитывает `thickStart`/`thickEnd`/`blockSizes`/`blockStarts` у BED12 |
| **GFF / GTF / GFF3** | `.gff`, `.gtf`, `.gff3` (+ `.gz`) | 1-based включительно | Обновляет границы признаков; инвертирует столбец strand на инвертированных блоках; секцию `##FASTA` в конце сохраняет без изменений |
| **Wiggle** | `.wig` (+ `.gz`) | 1-based | Поддерживает `variableStep` и `fixedStep`; на выходе всегда `variableStep`, так как после переноса позиции уже не равноотстоят; несколько дорожек остаются раздельными |
| **MAF** | `.maf` (+ `.gz`) | 1-based | Mutation Annotation Format: обновляет `Start_Position` / `End_Position`; на инвертированных блоках заменяет столбцы аллелей обратным комплементом (`NCBI_Build` не изменяется) |

---

### Как обрабатываются записи

* **Разрывы цепочки:** запись VCF, у которой `REF` пересекает разрыв выравнивания, не переносится (попадает в `.unmapped`), поскольку длины аллелей перестали бы быть корректными. Интервалы BED / GFF / MAF / Wiggle переносятся, если оба конца попадают на одну хромосому и одну цепь.
* **Инвертированные блоки:** координаты, столбцы цепи и аллели переворачиваются, как описано выше. Аллели-брейкенды VCF (`BND`) на инвертированных блоках отправляются в `.unmapped`.
* **Некорректные и короткие строки:** строки, которые не удалось разобрать, пишутся в `.unmapped`; строки VCF короче 8 столбцов дополняются `.`; пустые строки игнорируются. Поэтому всегда `Всего = Перенесено + Пропущено + Отклонено`.
* **Переводы строк:** на выходе всегда `LF`, в том числе в Windows (CRLF ломает ряд VCF-утилит).
* **Сбои:** при ошибке недописанные выходные файлы удаляются, а не остаются на диске. Выходной файл не может совпадать с входным.

---

### Быстрый старт

1. Перейдите во вкладку **[Releases](../../releases)** репозитория.
2. Скачайте файл под свою операционную систему.
3. Запустите программу.
4. Перетащите геномный файл в поле **Входной файл**.
5. Выберите нужное направление в списке **Встроенный пресет Chain** (или укажите свой chain-файл).
6. *(Опционально)* Укажите файл референса целевого генома (например, `hg38.fa`) в поле **Целевой референс FASTA** для биологической валидации REF.
7. Проверьте путь сохранения в поле **Сохранить результат**.
8. Нажмите **🚀 Запустить LiftOverZ**.
9. При необходимости используйте кнопку **🔄 Сброс** для очистки путей и лога перед обработкой следующего файла.

---

### Запуск скачанных бинарников

* **Linux / macOS:** файлы релизов GitHub не сохраняют бит исполнения — один раз выполните `chmod +x LiftOverZ-linux` (или `LiftOverZ-macos-*`).
* **macOS:** бинарники не подписаны и не нотаризованы. При первом запуске используйте правый клик → *Открыть* или выполните `xattr -d com.apple.quarantine LiftOverZ-macos-*`.
* **Linux:** бинарник собран на Ubuntu 24.04 и требует достаточно свежей glibc; на старых дистрибутивах запускайте из исходников (`python liftover_gui.py`).

---

### Настройки и параметры

* **Целевой референс FASTA (Опционально):** Путь к FASTA-файлу целевой сборки (`.fa`, `.fasta`, `.fa.gz`). При первом запуске автоматически создаёт `.fai` индекс, если он отсутствует (папка должна быть доступна для записи). Сжатая FASTA должна быть упакована **bgzip** — обычный gzip `pyfaidx` не поддерживает.
* **Авто-сопоставление префикса 'chr':** Устраняет различия между нотациями NCBI/Ensembl (`1`, `2`, `MT`) и UCSC (`chr1`, `chr2`, `chrM`), сохраняя исходный стиль в итоговом файле.
* **Сохранять неперенесённые записи:** Создаёт файл `.unmapped` со строками, для которых в цепочке переноса не нашлось пары.
* **Файл соответствия хромосом (Опционально):** Текстовый файл с соответствиями нестандартных имен (разделитель — пробел или табуляция, строки с `#` игнорируются). Перечисленные здесь контиги записываются в результат под именем целевой сборки из chain-файла (например, `NC_000001.10` → `chr1`):
```text
NC_000001.10   chr1
NC_000002.11   chr2
MT             chrM
```

---

### Добавление новых языков

1. Создайте папку `locales/` рядом с бинарником.
2. Скопируйте любой JSON файл (например, `ru.json` или `en.json`) с новым кодом языка (например, `locales/fr.json`).
3. Переведите строки в любом текстовом редакторе.
4. Перезапустите LiftOverZ — новый язык появится в выпадающем меню вверху.

---

### Сборка из исходного кода

```bash
# Клонирование репозитория
git clone https://github.com/zajakin/LiftOverZ.git
cd LiftOverZ

# Установка зависимостей
pip install -r requirements.txt

# Запуск скрипта
python liftover_gui.py

# Сборка автономного бинарника
pyinstaller --onefile --noconsole --add-data "locales:locales" --add-data "chains:chains" --collect-all tkinterdnd2 --collect-all pyliftover --collect-all pyfaidx --name LiftOverZ liftover_gui.py
```
*(В Windows замените разделители двоеточия `:locales` и `:chains` на точку с запятой `;locales` и `;chains`)*
