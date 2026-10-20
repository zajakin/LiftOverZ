import os
import gzip
import sys
import shutil
import urllib.request

# Ссылки для автоматической загрузки исходных chain-файлов с сервера UCSC
CHAIN_URLS = {
    "hg38ToHg19.over.chain.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
    "hg19ToHg38.over.chain.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
    "hs1ToHg19.over.chain.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg19.over.chain.gz",
    "hg19ToHs1.over.chain.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHs1.over.chain.gz",
    "hg38ToHs1.over.chain.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz",
    "hs1ToHg38.over.chain.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg38.over.chain.gz",
    "mm9ToMm10.over.chain.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz",
    "mm10ToMm39.over.chain.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz",
    "mm39ToMm10.over.chain.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/liftOver/mm39ToMm10.over.chain.gz",
}

PATCHES = {
    # 1. hg38 -> hg19: source tName=chrM (16,569 bp) maps to the legacy
    #    target qName=chrM (16,571 bp). The custom shift is name-specific.
    "hg38ToHg19.over.chain.gz": (
        "chain 100000000 chrM 16569 + 0 16569 chrM 16571 + 0 16571 9000001\n"
        "3106 0 2\n"
        "13463\n"
    ),

    # 2. hg19 -> hg38: patch only the legacy source tName=chrM (16,571 bp).
    "hg19ToHg38.over.chain.gz": (
        "chain 100000000 chrM 16571 + 0 16571 chrM 16569 + 0 16569 9000002\n"
        "3106 2 0\n"
        "13463\n"
    ),

    # --- МЫШЬ (NC_005089.1, длина 16299 п.н., тождественное совпадение 1-в-1) ---
    "mm9ToMm10.over.chain.gz": (
        "chain 100000000 chrM 16299 + 0 16299 chrM 16299 + 0 16299 9000010\n"
        "16299\n"
    ),
    "mm10ToMm39.over.chain.gz": (
        "chain 100000000 chrM 16299 + 0 16299 chrM 16299 + 0 16299 9000011\n"
        "16299\n"
    ),
    "mm39ToMm10.over.chain.gz": (
        "chain 100000000 chrM 16299 + 0 16299 chrM 16299 + 0 16299 9000012\n"
        "16299\n"
    ),
}

# Дополнительные тождественные привязки (например, актуальный chrMT hg19)
EXTRA_PATCHES = {
    "hg19ToHg38.over.chain.gz": [
        "chain 100000001 chrMT 16569 + 0 16569 chrM 16569 + 0 16569 9000020\n"
        "16569\n"
    ],
}


def _verify_gzip(path):
    found_chain = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("chain "):
                continue
            parts = line.split()
            if len(parts) < 13:
                raise ValueError("повреждённый заголовок chain-файла UCSC")
            try:
                int(parts[1])
                int(parts[3]); int(parts[5]); int(parts[6])
                int(parts[8]); int(parts[10]); int(parts[11]); int(parts[12])
            except ValueError as exc:
                raise ValueError("некорректные числовые поля заголовка chain-файла UCSC") from exc
            found_chain = True
            break
        for _ in f:
            pass
    if not found_chain:
        raise ValueError("файл не похож на chain-файл UCSC")


def download_chain(url, target_path):
    print(f"[СКАЧИВАНИЕ] Загрузка {os.path.basename(target_path)}...")
    tmp_path = target_path + ".download"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"})
    try:
        with urllib.request.urlopen(req, timeout=60) as response, open(tmp_path, "wb") as out_file:
            shutil.copyfileobj(response, out_file)
        _verify_gzip(tmp_path)
        os.replace(tmp_path, target_path)
        print(f"[СКАЧИВАНИЕ] Файл {os.path.basename(target_path)} успешно получен.")
        return True
    except Exception as e:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        print(f"[ОШИБКА] Не удалось скачать {os.path.basename(target_path)} с {url}: {e}")
        return False


def _validate_patch_block(new_block):
    lines = [line for line in new_block.strip().splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError("пустой или неполный блок chain для патча")
    header = lines[0].split()
    if len(header) < 13 or header[0] != "chain":
        raise ValueError("некорректный заголовок chain в PATCHES")
    try:
        t_start, t_end = int(header[5]), int(header[6])
        q_start, q_end = int(header[10]), int(header[11])
    except ValueError as exc:
        raise ValueError("некорректные координаты в заголовке PATCHES") from exc
    t_covered = q_covered = 0
    terminal_seen = False
    for index, line in enumerate(lines[1:], start=1):
        cols = line.split()
        if len(cols) == 1:
            if index != len(lines) - 1 or terminal_seen:
                raise ValueError("терминальный блок chain должен быть последним")
            size = int(cols[0])
            if size <= 0:
                raise ValueError("размер блока chain должен быть положительным")
            t_covered += size
            q_covered += size
            terminal_seen = True
        elif len(cols) == 3:
            if terminal_seen:
                raise ValueError("обнаружены данные после терминального блока chain")
            size, dt, dq = map(int, cols)
            if size <= 0 or dt < 0 or dq < 0:
                raise ValueError("некорректные размеры блока/разрывов chain")
            t_covered += size + dt
            q_covered += size + dq
        else:
            raise ValueError("строка блока chain должна содержать 1 или 3 числа")
    if not terminal_seen:
        raise ValueError("в chain отсутствует терминальный блок")
    if t_covered != t_end - t_start or q_covered != q_end - q_start:
        raise ValueError(
            f"длины блоков не совпадают с заголовком: t={t_covered}/{t_end - t_start}, "
            f"q={q_covered}/{q_end - q_start}"
        )
    return header


def patch_chain_file(filepath, new_block, *, require_match=True):
    new_hdr = _validate_patch_block(new_block)
    expected_source_name = new_hdr[2].lower()
    expected_source_size = int(new_hdr[3])

    with gzip.open(filepath, "rt", encoding="utf-8", errors="strict") as f_in:
        lines = f_in.readlines()

    clean_lines = []
    skipping_chain = False
    removed_count = 0
    saw_matching_source = False

    for line in lines:
        if line.startswith("chain "):
            parts = line.split()
            skipping_chain = False
            if len(parts) < 13:
                raise ValueError("обнаружен некорректный заголовок chain")
            try:
                source_size = int(parts[3])
            except ValueError as exc:
                raise ValueError("некорректный tSize в chain-заголовке") from exc

            if parts[2].lower() == expected_source_name and source_size == expected_source_size:
                saw_matching_source = True
                skipping_chain = True
                removed_count += 1
                continue
            clean_lines.append(line)
        elif skipping_chain:
            continue
        else:
            clean_lines.append(line)

    if require_match and not saw_matching_source:
        raise ValueError(
            f"не найдены исходные цепочки для контига {new_hdr[2]} длиной "
            f"{expected_source_size} п.н."
        )

    final_content = "".join(clean_lines).rstrip() + "\n\n" + new_block.strip() + "\n\n"
    tmp_path = filepath + ".tmp"
    try:
        with gzip.open(tmp_path, "wt", encoding="utf-8", newline="\n") as f_out:
            f_out.write(final_content)
        _verify_gzip(tmp_path)
        os.replace(tmp_path, filepath)
    except BaseException:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise
    return removed_count


def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    chains_dir = os.path.join(base_dir, "chains")
    os.makedirs(chains_dir, exist_ok=True)

    print(f"Обработка файлов в папке '{chains_dir}/':\n" + "-" * 50)
    ok_count, failed = 0, []

    for filename, url in CHAIN_URLS.items():
        filepath = os.path.join(chains_dir, filename)
        if not os.path.isfile(filepath):
            if not download_chain(url, filepath):
                print(f"[ОШИБКА] {filename} — не удалось получить файл.")
                failed.append(filename)
                continue

        try:
            _verify_gzip(filepath)
            patches_applied = []
            
            legacy_block = PATCHES.get(filename)
            if legacy_block is not None:
                patches_applied.append((legacy_block, True))
            
            extra_blocks = EXTRA_PATCHES.get(filename)
            if extra_blocks is not None:
                # Если передан как строка, преобразуем в список во избежание итерации по символам
                if isinstance(extra_blocks, str):
                    extra_blocks = [extra_blocks]
                patches_applied.extend((block, False) for block in extra_blocks)

            removed_total = 0
            for block, require_match in patches_applied:
                removed_total += patch_chain_file(filepath, block, require_match=require_match)

            size_kb = os.path.getsize(filepath) / 1024
            if patches_applied:
                print(
                    f"[УСПЕХ]  {filename}: удалено {removed_total} старых цепочек по точному "
                    f"имени исходного контига, применено патч-блоков: {len(patches_applied)} "
                    f"({size_kb:.1f} КБ)"
                )
            else:
                print(
                    f"[БЕЗ ПАТЧА] {filename}: оставлена официальная UCSC-цепочка; "
                    "для hs1 не подменяем митохондрию синтетическим 1:1-выравниванием"
                )
            ok_count += 1
        except Exception as e:
            print(f"[ОШИБКА] Не удалось обработать {filename}: {e}")
            failed.append(filename)

    print("-" * 50)
    if failed:
        print(f"Готово с ошибками: обработано {ok_count} из {len(CHAIN_URLS)}. Проблемные файлы: {', '.join(failed)}")
        return 1
    patched_files = sum(1 for name in CHAIN_URLS if name in PATCHES or name in EXTRA_PATCHES)
    print(f"Готово. Проверено цепочек: {ok_count}; файлов с патчами: {patched_files}.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
