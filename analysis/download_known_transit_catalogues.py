"""Download the catalogues of known transit signals, to tell new candidates from known ones.

The script downloads:
- the TESS Objects of Interest (TOI) list from ExoFOP;
- the Community TESS Objects of Interest (CTOI) list from ExoFOP;
- every SPOC threshold crossing event (TCE) table from MAST, one per sector or range of sectors;
- the TESS eclipsing binary (EB) catalogue of Prsa et al. (2022) from VizieR, with the EB periods and epochs.

Run it again to get the current lists. It downloads every file again, about 300 MB.

Run from the repository root:
    .venv/bin/python analysis/download_known_transit_catalogues.py

Output: data/known_transit_catalogues/toi.csv, ctoi.csv, tess_eb.csv and spoc_tce/<MAST file name>.csv,
and download_date.txt with the UTC date of the download.
"""

import datetime
import re
import urllib.request
from pathlib import Path

from astroquery.vizier import Vizier

from pipeline_status import PipelineStatus

# Every VizieR query of the pipeline uses this server. The IDIA mirror holds the same catalogues as CDS Strasbourg.
# From 2026-10-28 CDS was unreachable, and the CfA mirror answered every query with 504, but this mirror answered.
VIZIER_SERVER = "vizier.idia.ac.za"
TESS_EB_VIZIER_CATALOGUE = "J/ApJS/258/16/tess-ebs"
REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
CATALOGUE_DIRECTORY = REPOSITORY_ROOT / "data/known_transit_catalogues"
TOI_URL = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
CTOI_URL = "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv"
TCE_INDEX_URL = "https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_tce.html"
MAST_ROOT_URL = "https://archive.stsci.edu"


def download(url, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, path)
    print(f"{path.relative_to(REPOSITORY_ROOT)} ({path.stat().st_size / 1e6:.1f} MB)", flush=True)


def download_tess_eb_catalogue(path):
    """Write the TIC ID, position, period (d) and epoch (BTJD) of each EB in the TESS EB catalogue."""
    table = Vizier(row_limit=-1, vizier_server=VIZIER_SERVER).get_catalogs(TESS_EB_VIZIER_CATALOGUE)[0]
    table = table["TIC", "RAJ2000", "DEJ2000", "Per", "BJD0", "Morph"]
    table.rename_columns(
        ["TIC", "RAJ2000", "DEJ2000", "Per", "BJD0", "Morph"],
        ["tic_id", "ra", "dec", "period_days", "epoch_btjd", "morphology"],
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    table.write(path, format="csv", overwrite=True)
    print(f"{path.relative_to(REPOSITORY_ROOT)} ({len(table)} EBs)", flush=True)


def main():
    index = urllib.request.urlopen(TCE_INDEX_URL).read().decode()
    tce_paths = sorted(set(re.findall(r'href="(/missions/tess/catalogs/tce/[^"]*tcestats[^"]*\.csv)"', index)))
    with PipelineStatus(None, 5, "known_transit_catalogue_download", "download_known_transit_catalogues.py",
                        3 + len(tce_paths), item_name="catalogue files") as status:
        download(TOI_URL, CATALOGUE_DIRECTORY / "toi.csv")
        download(CTOI_URL, CATALOGUE_DIRECTORY / "ctoi.csv")
        download_tess_eb_catalogue(CATALOGUE_DIRECTORY / "tess_eb.csv")
        status.set_items_done(3)
        for count, tce_path in enumerate(tce_paths, start=4):
            download(MAST_ROOT_URL + tce_path, CATALOGUE_DIRECTORY / "spoc_tce" / tce_path.split("/")[-1])
            status.set_items_done(count)
    today = datetime.datetime.now(datetime.timezone.utc).date().isoformat()
    (CATALOGUE_DIRECTORY / "download_date.txt").write_text(today + "\n")
    print(f"Downloaded the TOI and CTOI lists and {len(tce_paths)} SPOC TCE tables on {today}")


if __name__ == "__main__":
    main()
