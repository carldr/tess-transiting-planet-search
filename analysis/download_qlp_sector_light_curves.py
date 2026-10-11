"""Download the QLP light curves of one TESS sector and keep them as compact arrays.

The script reads data/qlp_sector_<nnnn>/star_list.csv (from build_qlp_sector_star_list.py) and keeps the stars at
or brighter than --max-tmag. For each star it downloads the QLP FITS file into memory and saves only
TIME, SAP_FLUX and QUALITY. A QLP FITS file of Sector 98 is about 1.5 MB; the saved arrays are about 0.25 MB.

The script skips stars that it already saved, so it continues an interrupted run.
It writes the stars that failed to download into failed_downloads.csv; run the script again to retry them.

Run from the repository root, for example for Sector 98:
    .venv/bin/python analysis/download_qlp_sector_light_curves.py 98 --max-tmag 11.5

Output: data/qlp_sector_<nnnn>/light_curves/<last three digits of TIC ID>/<TIC ID>.npz,
with the arrays time (BTJD), sap_flux and quality.
"""

import argparse
import csv
import io
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from astropy.io import fits

from build_qlp_sector_star_list import sector_directory

DOWNLOAD_THREADS = 16
DOWNLOAD_ATTEMPTS = 3
DOWNLOAD_TIMEOUT_SECONDS = 60


def light_curve_path(sector, tic_id):
    return sector_directory(sector) / "light_curves" / f"{tic_id % 1000:03d}" / f"{tic_id}.npz"


def download_and_save(sector, tic_id, url):
    path = light_curve_path(sector, tic_id)
    if path.exists():
        return tic_id, "skipped", 0
    for attempt in range(1, DOWNLOAD_ATTEMPTS + 1):
        try:
            with urllib.request.urlopen(url, timeout=DOWNLOAD_TIMEOUT_SECONDS) as response:
                content = response.read()
            break
        except Exception as error:
            if attempt == DOWNLOAD_ATTEMPTS:
                return tic_id, f"failed: {error}", 0
            time.sleep(2 * attempt)
    data = fits.open(io.BytesIO(content))[1].data
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(".partial.npz")
    np.savez_compressed(
        temporary_path,
        time=np.asarray(data["TIME"], dtype=np.float64),
        sap_flux=np.asarray(data["SAP_FLUX"], dtype=np.float32),
        quality=np.asarray(data["QUALITY"], dtype=np.int32),
    )
    temporary_path.rename(path)
    return tic_id, "saved", len(content)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sector", type=int)
    parser.add_argument("--max-tmag", type=float, default=11.5)
    arguments = parser.parse_args()
    sector = arguments.sector

    stars = [
        row
        for row in csv.DictReader(open(sector_directory(sector) / "star_list.csv"))
        if row["tmag"] and float(row["tmag"]) <= arguments.max_tmag
    ]
    print(f"Sector {sector}: {len(stars)} stars with Tmag <= {arguments.max_tmag}", flush=True)

    failures = []
    downloaded_bytes = 0
    start = time.time()
    with ThreadPoolExecutor(DOWNLOAD_THREADS) as executor:
        futures = [
            executor.submit(download_and_save, sector, int(row["tic_id"]), row["light_curve_url"]) for row in stars
        ]
        for count, future in enumerate(as_completed(futures), start=1):
            tic_id, status, size = future.result()
            downloaded_bytes += size
            if status.startswith("failed"):
                failures.append({"tic_id": tic_id, "error": status})
            if count % 1000 == 0:
                elapsed = time.time() - start
                print(
                    f"  {count} of {len(stars)} stars, {downloaded_bytes / 1e9:.1f} GB downloaded, "
                    f"{len(failures)} failed, {elapsed / 60:.0f} min",
                    flush=True,
                )

    failures_path = sector_directory(sector) / "failed_downloads.csv"
    with open(failures_path, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=["tic_id", "error"])
        writer.writeheader()
        writer.writerows(failures)
    print(f"Done: {len(stars) - len(failures)} stars saved or already present, {len(failures)} failed")


if __name__ == "__main__":
    main()
