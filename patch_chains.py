import os
import gzip
import shutil

PATCHES = {
    # 1. hg38 -> hg19: 16569 п.н. -> 16571 п.н. (вставка 2 п.н. NN в hg19 после 3106-й позиции; в chain-файле смещение 0-based = 3106)
    "hg38ToHg19.over.chain.gz": (
        "chain 100000000 chrM 16569 + 0 16569 chrM 16571 + 0 16571 9000001\n"
        "3106 0 2\n"
        "13463\n"
    ),

    # 2. hg19 -> hg38: 16571 п.н. -> 16569 п.н. (пропуск 2 п.н. NN в hg19 после 3106-й позиции)
    "hg19ToHg38.over.chain.gz": (
        "chain 100000000 chrM 16571 + 0 16571 chrM 16569 + 0 16569 9000002\n"
        "3106 2 0\n"
        "13463\n"
    ),

    # 3. hs1 (T2T) -> hg19: 16569 п.н. -> 16571 п.н.
    "hs1ToHg19.over.chain.gz": (
        "chain 100000000 chrM 16569 + 0 16569 chrM 16571 + 0 16571 9000003\n"
        "3106 0 2\n"
        "13463\n"
    ),

    # 4. hg19 -> hs1 (T2T): 16571 п.н. -> 16569 п.н.
    "hg19ToHs1.over.chain.gz": (
        "chain 100000000 chrM 16571 + 0 16571 chrM 16569 + 0 16569 9000004\n"
        "3106 2 0\n"
        "13463\n"
    ),

    # 5. hg38 <-> hs1: обе сборки rCRS (16569 п.н.), чистое выравнивание без разрывов
    "hg38ToHs1.over.chain.gz": (
        "chain 100000000 chrM 16569 + 0 16569 chrM 16569 + 0 16569 9000005\n"
        "16569\n"
    ),
    "hs1ToHg38.over.chain.gz": (
        "chain 100000000 chrM 16569 + 0 16569 chrM 16569 + 0 16569 9000006\n"
        "16569\n"
    ),
}


def patch_chain_file(filepath, new_block):
    """Заменяет все блоки chrM в chain-файле на new_block.
    Запись идёт во временный файл и атомарно подменяет оригинал: при сбое исходный файл не страдает."""
    with gzip.open(filepath, "rt", encoding="utf-8", errors="replace") as f_in:
        lines = f_in.readlines()

    clean_lines = []
    skipping_chrM = False
    removed_count = 0

    for line in lines:
        if line.startswith("chain "):
            parts = line.split()
            if len(parts) >= 8 and parts[2] == "chrM":
                skipping_chrM = True
                removed_count += 1
                continue
            else:
                skipping_chrM = False
                clean_lines.append(line)
        elif skipping_chrM:
            if not line.strip():
                skipping_chrM = False
            continue
        else:
            clean_lines.append(line)

    final_content = "".join(clean_lines).rstrip() + "\n\n" + new_block.strip() + "\n\n"

    tmp_path = filepath + ".tmp"
    try:
        with gzip.open(tmp_path, "wt", encoding="utf-8") as f_out:
            f_out.write(final_content)
        # проверяем, что новый файл читается, и только потом подменяем оригинал
        with gzip.open(tmp_path, "rt", encoding="utf-8") as f_chk:
            for _ in f_chk:
                pass
        os.replace(tmp_path, filepath)
    except BaseException:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise
    return removed_count


def main():
    chains_dir = "chains"
    if not os.path.isdir(chains_dir):
        print(f"Ошибка: папка '{chains_dir}' не найдена в текущей директории.")
        return

    print(f"Обработка файлов в папке '{chains_dir}/':\n" + "-" * 50)

    for filename, block in PATCHES.items():
        filepath = os.path.join(chains_dir, filename)
        if not os.path.exists(filepath):
            print(f"[ПРОПУСК] {filename} — файл не найден.")
            continue

        removed = patch_chain_file(filepath, block)
        size_kb = os.path.getsize(filepath) / 1024
        print(f"[УСПЕХ]  {filename}: удалено {removed} старых блоков chrM, внедрена единая непрерывная цепь ({size_kb:.1f} КБ)")

    print("-" * 50)
    print("Готово. Цепочки успешно исправлены.")


if __name__ == "__main__":
    main()