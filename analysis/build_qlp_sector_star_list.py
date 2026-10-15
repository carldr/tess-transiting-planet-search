"""Build the list of stars with a QLP light curve in one TESS sector, with their TIC stellar parameters.

The script reads the TIC IDs from the MAST bulk download script of the sector, and downloads that script
when it is not yet in data/qlp_bulk_download_scripts/.
It queries the TESS Input Catalog (TIC) at MAST in batches of 2000 IDs.

Run from the repository root, for example for Sector 98:
    .venv/bin/python analysis/build_qlp_sector_star_list.py 98

Output: data/qlp_sector_<nnnn>/star_list.csv, one row per star, with the columns
tic_id, ra, dec, tmag, teff, radius, mass, logg, contamination_ratio, light_curve_url.
"""

import argparse
import csv
import re
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from astroquery.mast import Catalogs

from pipeline_status import PipelineStatus

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
BULK_DOWNLOAD_SCRIPT_DIRECTORY = REPOSITORY_ROOT / "data/qlp_bulk_download_scripts"
BULK_DOWNLOAD_SCRIPT_URL = "https://archive.stsci.edu/hlsps/qlp/download_scripts/{name}"
TIC_QUERY_BATCH_SIZE = 2000
TIC_QUERY_THREADS = 4
TIC_QUERY_ATTEMPTS = 3

# TIC column name -> star list column name.
TIC_COLUMNS = {
    "ID": "tic_id",
    "ra": "ra",
    "dec": "dec",
    "Tmag": "tmag",
    "Teff": "teff",
    "rad": "radius",
    "mass": "mass",
    "logg": "logg",
    "contratio": "contamination_ratio",
}


def sector_directory(sector):
    return REPOSITORY_ROOT / f"data/qlp_sector_{sector:04d}"


def bulk_download_script_path(sector):
    name = f"hlsp_qlp_tess_ffi_s{sector:04d}_tess_v01_llc-fits.sh"
    path = BULK_DOWNLOAD_SCRIPT_DIRECTORY / name
    if not path.exists():
        BULK_DOWNLOAD_SCRIPT_DIRECTORY.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(BULK_DOWNLOAD_SCRIPT_URL.format(name=name), path)
    return path


def light_curve_url_by_tic_id(sector):
    """Return {TIC ID: light curve URL} from the bulk download script. The script is only read, never run."""
    pattern = re.compile(rf"'(https://\S+s{sector:04d}-(\d{{16}})_tess_v01_llc\.fits)'")
    urls = {}
    for line in open(bulk_download_script_path(sector)):
        match = pattern.search(line)
        if match:
            urls[int(match.group(2))] = match.group(1)
    return urls


def query_tic_batch(tic_ids):
    for attempt in range(1, TIC_QUERY_ATTEMPTS + 1):
        try:
            table = Catalogs.query_criteria(catalog="Tic", ID=tic_ids)
            return [{name: row[column] for column, name in TIC_COLUMNS.items()} for row in table]
        except Exception as error:
            if attempt == TIC_QUERY_ATTEMPTS:
                raise
            print(f"TIC query failed ({error}); attempt {attempt + 1} of {TIC_QUERY_ATTEMPTS}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sector", type=int)
    sector = parser.parse_args().sector

    urls = light_curve_url_by_tic_id(sector)
    tic_ids = sorted(urls)
    print(f"Sector {sector}: {len(tic_ids)} QLP light curves", flush=True)
    batches = [tic_ids[i : i + TIC_QUERY_BATCH_SIZE] for i in range(0, len(tic_ids), TIC_QUERY_BATCH_SIZE)]

    rows = []
    status = PipelineStatus(sector, 1, "star_list", "build_qlp_sector_star_list.py", len(tic_ids), "stars")
    with ThreadPoolExecutor(TIC_QUERY_THREADS) as executor, status:
        for batch_number, batch_rows in enumerate(executor.map(query_tic_batch, batches), start=1):
            rows.extend(batch_rows)
            status.set_items_done(min(batch_number * TIC_QUERY_BATCH_SIZE, len(tic_ids)))
            if batch_number % 50 == 0:
                print(f"  {batch_number} of {len(batches)} TIC batches", flush=True)

    for row in rows:
        row["tic_id"] = int(row["tic_id"])
        row["light_curve_url"] = urls[row["tic_id"]]
    rows.sort(key=lambda row: row["tic_id"])
    missing = len(tic_ids) - len(rows)
    if missing:
        print(f"WARNING: {missing} TIC IDs have no TIC row", flush=True)

    output_path = sector_directory(sector) / "star_list.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(TIC_COLUMNS.values()) + ["light_curve_url"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} stars to {output_path.relative_to(REPOSITORY_ROOT)}")


if __name__ == "__main__":
    main()
