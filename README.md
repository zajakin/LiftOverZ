# LiftOverZ

[English](#english) | [Русский](#русский)

---

<a name="english"></a>
## English

**LiftOverZ** is a fast, standalone, cross-platform graphical tool designed for lifting over genomic coordinates across genome assembly builds (e.g., hg19/GRCh37, hg38/GRCh38, T2T-CHM13v2, and mouse builds).

It requires **no installation**, runs as a single portable executable, includes **pre-bundled popular UCSC chain files**, supports **Drag & Drop**, handles compressed files (`.gz`) out of the box, and features hot-swappable language support.

---

### Key Features

* **Zero Installation / 100% Portable:** Single executable running out of the box on Windows 11/10, Linux, and macOS.
* **Native Cross-Platform Binaries:**
  * Windows (`LiftOverZ-windows.exe`)
  * Linux x86_64 (`LiftOverZ-linux`)
  * macOS Apple Silicon M1/M2/M3/M4 (`LiftOverZ-macos-arm64`)
  * macOS Intel x86_64 (`LiftOverZ-macos-intel`)
* **Pre-bundled Chain Presets:** Ready-to-use human and mouse assembly conversions embedded directly into the app (no manual downloads required).
* **Intuitive GUI with Drag & Drop:** Drop input genomic files, custom `.chain` files, or mapping tables straight into the interface.
* **Broad Genomic Format Support:** Automatically detects and converts 6 major genomic formats.
* **Transparent Gzip Compression:** Direct read and write support for `.gz` archives without pre-decompression.
* **Chromosome Prefix Tolerance:** Automatically resolves nomenclature discrepancies (e.g., `chr1` vs `1`).
* **Custom Chromosome Mapping:** Supply custom alias tables for non-standard contig names.
* **Unmapped Records Audit:** Exports unmapped records into a dedicated `.unmapped` file for verification.
* **Extensible Internationalization (i18n):** Shipped with English and Russian; expandable by dropping JSON files into `locales/`.

---

### Pre-bundled Chain Presets

LiftOverZ comes packaged with the most widely used UCSC chain files:
* **Human:**
  * hg19 (GRCh37) ⇄ hg38 (GRCh38)
  * hg38 (GRCh38) ⇄ T2T-CHM13v2.0 (hs1)
  * hg19 (GRCh37) ⇄ T2T-CHM13v2.0 (hs1)
* **Mouse:**
  * mm10 (GRCm38) ⇄ mm39 (GRCm39)
  * mm9 (NCBI37) → mm10 (GRCm38)

*You can also pick "Custom Chain File" to load any other `.chain` or `.chain.gz` file.*

---

### Supported File Formats

| Format | Extensions | Coordinate System | Notes |
| :--- | :--- | :--- | :--- |
| **VCF** | `.vcf`, `.vcf.gz` | 1-based | Retains headers and variant metadata |
| **gVCF** | `.gvcf`, `.gvcf.gz` | 1-based | Accurately handles non-variant blocks |
| **BED / BEDGraph** | `.bed`, `.bedgraph` | 0-based start, 1-based end | Updates interval spans cleanly |
| **GFF / GTF / GFF3** | `.gff`, `.gtf`, `.gff3` | 1-based inclusive | Updates feature start and end bounds |
| **Wiggle** | `.wig` | 1-based | Supports `variableStep` and `fixedStep` |
| **MAF** | `.maf` | 1-based | Updates `Start_Position` and `End_Position` |

---

### Quick Start

1. Go to the **[Releases](../../releases)** tab of this repository.
2. Download the binary for your platform:
   * **Windows 11 / 10:** `LiftOverZ-windows.exe`
   * **Linux:** `LiftOverZ-linux` (run `chmod +x LiftOverZ-linux` before launch)
   * **macOS (M1/M2/M3/M4):** `LiftOverZ-macos-arm64`
   * **macOS (Intel):** `LiftOverZ-macos-intel`
3. Launch the application.
4. Drag & drop your genomic data into **Input File**.
5. Select a conversion from the **Built-in Chain Preset** dropdown (or drop your own chain file).
6. Confirm or change the output location in **Save Result To**.
7. Click **🚀 Run LiftOverZ**.

---

### Options & Configuration

* **Auto-match 'chr' prefix:** When checked, reconciles naming variations between NCBI/Ensembl (`1`, `2`) and UCSC (`chr1`, `chr2`).
* **Save unmapped records:** Generates a `.unmapped` file containing records that could not be lifted over.
* **Chromosome Mapping File (Optional):** Whitespace-delimited text file mapping custom contig names:
```text
NC_000001.10   chr1
NC_000002.11   chr2
MT             chrM
```

---

### Adding Custom Translations

1. Create a `locales/` directory next to the executable.
2. Copy `locales/en.json` to your language code (e.g., `locales/de.json`, `locales/es.json`).
3. Translate the values in any text editor.
4. Restart the app — the new language will be available in the language picker.

---

### Local Build & Development

```bash
# Clone the repository
git clone https://github.com/your-username/LiftOverZ.git
cd LiftOverZ

# Install dependencies
pip install -r requirements.txt

# Run directly
python liftover_gui.py

# Build a standalone binary
pyinstaller --onefile --noconsole \
  --add-data "locales:locales" \
  --add-data "chains:chains" \
  --collect-all tkinterdnd2 \
  --collect-all pyliftover \
  --name LiftOverZ liftover_gui.py
```
*(On Windows, replace `:locales` and `:chains` with `;locales` and `;chains` in `--add-data`)*

---

<a name="русский"></a>
## Русский

**LiftOverZ** — это быстрый автономный графический инструмент для переноса (liftover) геномных координат между сборками референсного генома (hg19/GRCh37, hg38/GRCh38, T2T-CHM13v2, а также сборками генома мыши).

Программа **не требует установки**, поставляется в виде одного портативного исполняемого файла, содержит **встроенные популярные chain-файлы UCSC**, поддерживает **Drag & Drop**, напрямую работает с архивами `.gz` и позволяет переключать язык на лету.

---

### Основные возможности

* **Полная портативность:** Один исполняемый файл без внешних зависимостей. Запускается без инсталлятора на Windows 11/10, Linux и macOS.
* **Нативные бинарники под все платформы:**
  * Windows (`LiftOverZ-windows.exe`)
  * Linux x86_64 (`LiftOverZ-linux`)
  * macOS на процессорах Apple Silicon M1/M2/M3/M4 (`LiftOverZ-macos-arm64`)
  * macOS на процессорах Intel (`LiftOverZ-macos-intel`)
* **Встроенные пресеты Chain:** Самые востребованные цепочки для человека и мыши уже включены в программу (не нужно искать и скачивать вручную).
* **Графический интерфейс с Drag & Drop:** Перетаскивайте исходные геномные файлы, сторонние chain-файлы и таблицы соответствия хромосом прямо в окно программы.
* **Поддержка 6 основных форматов данных:** Автоматическое определение формата по расширению.
* **Прозрачная работа с Gzip:** Чтение и запись архивов `.gz` на лету без предварительной распаковки.
* **Устойчивость к префиксу 'chr':** Автоматическое сопоставление хромосом `chr1` и `1`.
* **Пользовательский маппинг хромосом:** Поддержка таблиц соответствия нестандартных идентификаторов контигов.
* **Аудит неперенесённых строк:** Запись пропущенных строк в файл `.unmapped` для последующего контроля.
* **Расширяемая локализация (i18n):** Встроенные русский и английский интерфейсы с возможностью добавления новых языков через JSON.

---

### Встроенные цепочки переноса

В дистрибутив LiftOverZ включены следующие официальные цепочки UCSC:
* **Человек:**
  * hg19 (GRCh37) ⇄ hg38 (GRCh38)
  * hg38 (GRCh38) ⇄ T2T-CHM13v2.0 (hs1)
  * hg19 (GRCh37) ⇄ T2T-CHM13v2.0 (hs1)
* **Мышь:**
  * mm10 (GRCm38) ⇄ mm39 (GRCm39)
  * mm9 (NCBI37) → mm10 (GRCm38)

*Для других сборок выберите вариант «Пользовательский файл» и укажите любой сторонний `.chain` или `.chain.gz` файл.*

---

### Поддерживаемые форматы данных

| Формат | Расширения | Система координат | Примечание |
| :--- | :--- | :--- | :--- |
| **VCF** | `.vcf`, `.vcf.gz` | 1-based | Сохраняет метаданные заголовка и аннотации |
| **gVCF** | `.gvcf`, `.gvcf.gz` | 1-based | Корректно переносит невариантные интервалы |
| **BED / BEDGraph** | `.bed`, `.bedgraph` | 0-based старт, 1-based конец | Сохраняет длину и границы интервалов |
| **GFF / GTF / GFF3** | `.gff`, `.gtf`, `.gff3` | 1-based включительно | Обновляет начальные и конечные координаты |
| **Wiggle** | `.wig` | 1-based | Поддерживает дорожки `variableStep` и `fixedStep` |
| **MAF** | `.maf` | 1-based | Обновляет колонки `Start_Position` и `End_Position` |

---

### Быстрый старт

1. Перейдите во вкладку **[Releases](../../releases)** репозитория.
2. Скачайте файл под свою операционную систему:
   * **Windows 11 / 10:** `LiftOverZ-windows.exe`
   * **Linux:** `LiftOverZ-linux` (выполните `chmod +x LiftOverZ-linux` перед первым запуском)
   * **macOS (M1/M2/M3/M4):** `LiftOverZ-macos-arm64`
   * **macOS (Intel):** `LiftOverZ-macos-intel`
3. Запустите программу.
4. Перетащите геномный файл в поле **Входной файл**.
5. Выберите нужное направление в выпадающем списке **Встроенный пресет Chain** (или укажите свой chain-файл).
6. Проверьте путь сохранения в поле **Сохранить результат**.
7. Нажмите **🚀 Запустить LiftOverZ**.

---

### Настройки и параметры

* **Авто-сопоставление префикса 'chr':** Устраняет различия между нотациями NCBI/Ensembl (`1`, `2`) и UCSC (`chr1`, `chr2`).
* **Сохранять неперенесённые записи:** Создаёт файл `.unmapped` со строками, для которых в цепочке переноса не нашлось пары.
* **Файл соответствия хромосом (Опционально):** Текстовый файл с соответствиями нестандартных имен (разделитель — пробел или табуляция):
```text
NC_000001.10   chr1
NC_000002.11   chr2
MT             chrM
```

---

### Добавление новых языков

1. Создайте папку `locales/` рядом с бинарником.
2. Скопируйте файл `locales/ru.json` или `locales/en.json` с новым кодом языка (например, `locales/de.json`, `locales/fr.json`).
3. Переведите строки в любом текстовом редакторе.
4. Перезапустите LiftOverZ — новый язык появится в меню выбора вверху.

---

### Сборка из исходного кода

```bash
# Клонирование репозитория
git clone https://github.com/your-username/LiftOverZ.git
cd LiftOverZ

# Установка зависимостей
pip install -r requirements.txt

# Запуск скрипта
python liftover_gui.py

# Сборка автономного бинарника
pyinstaller --onefile --noconsole \
  --add-data "locales:locales" \
  --add-data "chains:chains" \
  --collect-all tkinterdnd2 \
  --collect-all pyliftover \
  --name LiftOverZ liftover_gui.py
```
*(В Windows замените разделители двоеточия `:locales` и `:chains` на точку с запятой `;locales` и `;chains`)*
