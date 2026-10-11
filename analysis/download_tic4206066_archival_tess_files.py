"""Download the archival TESS files of TIC 4206066 that Rabtsevich used, and check each sha256.

The file list is `01_inputs/INPUT_PRODUCTS.csv` of the Rabtsevich data package.
This script downloads the light curves, the target pixel files and the 21x21 TESScut cutouts.
It skips the MAST listings, the PRF models, the documentation snapshots and the 31x31 cutouts.

Run from the repository root:
    .venv/bin/python analysis/download_tic4206066_archival_tess_files.py
"""

import csv
import hashlib
import sys
import urllib.request
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
INPUT_PRODUCTS_CSV = (
    REPOSITORY_ROOT
    / "reference/rabtsevich_data_zenodo_22999133/unpacked/01_inputs/INPUT_PRODUCTS.csv"
)
LIGHT_CURVE_DIRECTORY = REPOSITORY_ROOT / "data/tic4206066_mast"
TESSCUT_DIRECTORY = REPOSITORY_ROOT / "data/tic4206066_tesscut"


def sha256_of_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as file:
        for block in iter(lambda: file.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def selected_input_products():
    for row in csv.DictReader(open(INPUT_PRODUCTS_CSV)):
        is_fits = row["file_name"].endswith(".fits")
        if row["category"].startswith("TESS light curve") and is_fits:
            yield row, LIGHT_CURVE_DIRECTORY
        elif row["category"] == "TESScut FFI cutout" and "_21x21_" in row["file_name"]:
            yield row, TESSCUT_DIRECTORY


def main():
    mismatched_files = []
    for row, directory in selected_input_products():
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / row["file_name"]
        if not path.exists():
            print(f"download {row['file_name']} ({int(row['bytes']) / 1e6:.1f} MB)", flush=True)
            urllib.request.urlretrieve(row["url"], path)
        if sha256_of_file(path) == row["sha256"]:
            print(f"OK       {row['file_name']}")
        else:
            print(f"MISMATCH {row['file_name']}")
            mismatched_files.append(row["file_name"])
    if mismatched_files:
        sys.exit(f"{len(mismatched_files)} files do not match the sha256 in INPUT_PRODUCTS.csv")


if __name__ == "__main__":
    main()
