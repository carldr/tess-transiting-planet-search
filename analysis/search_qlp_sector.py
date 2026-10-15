"""Search every downloaded QLP light curve of one TESS sector for periodic transits with box least squares (BLS).

For each star the script:
1. removes cadences with a QUALITY bit other than QUALITY_BITS_KEPT, and normalises the SAP flux;
2. removes trends with a wotan biweight filter with a 0.5 d window, and removes points more than 5 sigma high;
3. runs BLS over periods of 0.5 to 14 d and durations of 1 to 4.8 h;
4. records the strongest signal, masks its transits, and searches again, up to SIGNALS_PER_STAR times.
   It stops early when a signal has a BLS SNR below MINIMUM_SNR_TO_CONTINUE.

Run from the repository root, for example for Sector 98:
    .venv/bin/python analysis/search_qlp_sector.py 98

The script works in chunks of 1000 stars and skips chunks that are already done, so it continues an interrupted run.
--restart deletes the chunks of an earlier run and searches all stars again.

Output: results/qlp_sector_<nnnn>_bls_search/chunk_<nnnn>.csv, one row per signal, with the columns
tic_id, signal_number, period_days, epoch_btjd, duration_hours, depth_ppm, snr, sde, transit_count, point_count,
noise_ppm. A row with signal_number 0 records a star that the search could not process, in the column error.
"""

import argparse
import csv
import os
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from astropy.timeseries import BoxLeastSquares
from wotan import flatten

from build_qlp_sector_star_list import REPOSITORY_ROOT, sector_directory
from pipeline_status import PipelineStatus

# QUALITY bits that do not remove a cadence. Every other non-zero QUALITY bit removes it.
# QLP does not document bit 30 (2^29). In Sector 98 it flags 31% of the cadences, and their noise is 1.1 times that
# of unflagged cadences. Keeping them doubles the confirmed planets that pass vetting, from 13 to 25 of 43.
QUALITY_BITS_KEPT = 2**29
BIWEIGHT_WINDOW_LENGTH_DAYS = 0.5
MINIMUM_PERIOD_DAYS = 0.5
MAXIMUM_PERIOD_DAYS = 14.0
DURATIONS_DAYS = np.array([0.04, 0.06, 0.08, 0.11, 0.15, 0.2])
FREQUENCY_FACTOR = 2.0
SIGNALS_PER_STAR = 3
MINIMUM_SNR_TO_CONTINUE = 7.0
MINIMUM_POINT_COUNT = 1000
CHUNK_SIZE = 1000

OUTPUT_COLUMNS = [
    "tic_id", "signal_number", "period_days", "epoch_btjd", "duration_hours", "depth_ppm", "snr", "sde",
    "transit_count", "point_count", "noise_ppm", "error",
]


def output_directory(sector):
    return REPOSITORY_ROOT / f"results/qlp_sector_{sector:04d}_bls_search"


def robust_standard_deviation(values):
    return 1.4826 * np.median(np.abs(values - np.median(values)))


def select_usable_cadences(time, flux, quality):
    """Return time and normalised flux, without the cadences that QUALITY_BITS_KEPT does not allow."""
    flux = np.asarray(flux, dtype=np.float64)
    keep = (quality & ~QUALITY_BITS_KEPT == 0) & np.isfinite(time) & np.isfinite(flux) & (flux > 0)
    return time[keep], flux[keep] / np.median(flux[keep])


def load_light_curve(light_curve_path):
    """Return time and normalised SAP flux of the usable cadences of one downloaded QLP light curve."""
    arrays = np.load(light_curve_path)
    return select_usable_cadences(arrays["time"], arrays["sap_flux"], arrays["quality"])


def detrend(time, flux):
    """Return time and detrended flux, with points more than 5 sigma above the trend removed."""
    detrended_flux = flatten(
        time, flux, method="biweight", window_length=BIWEIGHT_WINDOW_LENGTH_DAYS, edge_cutoff=0.0
    )
    keep = np.isfinite(detrended_flux)
    time, detrended_flux = time[keep], detrended_flux[keep]
    keep = detrended_flux < 1 + 5 * robust_standard_deviation(detrended_flux)
    return time[keep], detrended_flux[keep]


def search_light_curve(time, flux):
    """Return a list of signals, each a dict of the output columns except tic_id."""
    # Some QLP files have a finite flux at only a few cadences, and the detrending fails on so few points.
    if len(time) < MINIMUM_POINT_COUNT:
        return []
    time, flux = detrend(time, flux)
    signals = []
    for signal_number in range(1, SIGNALS_PER_STAR + 1):
        if len(time) < MINIMUM_POINT_COUNT:
            break
        noise = robust_standard_deviation(flux)
        model = BoxLeastSquares(time, flux, np.full_like(flux, noise))
        result = model.autopower(
            DURATIONS_DAYS,
            minimum_period=MINIMUM_PERIOD_DAYS,
            maximum_period=MAXIMUM_PERIOD_DAYS,
            frequency_factor=FREQUENCY_FACTOR,
            method="fast",
            objective="likelihood",
        )
        power = result.power
        best = int(np.argmax(power))
        period, epoch, duration = result.period[best], result.transit_time[best], result.duration[best]
        transit_numbers = np.round((time - epoch) / period)
        in_transit = np.abs(time - epoch - transit_numbers * period) < 0.5 * duration
        signals.append(
            {
                "signal_number": signal_number,
                "period_days": round(float(period), 6),
                "epoch_btjd": round(float(epoch), 5),
                "duration_hours": round(float(duration * 24), 3),
                "depth_ppm": round(float(result.depth[best] * 1e6), 1),
                "snr": round(float(result.depth_snr[best]), 2),
                "sde": round(float((power[best] - np.median(power)) / robust_standard_deviation(power)), 2),
                "transit_count": int(len(np.unique(transit_numbers[in_transit]))),
                "point_count": int(len(time)),
                "noise_ppm": round(float(noise * 1e6), 1),
                "error": "",
            }
        )
        if result.depth_snr[best] < MINIMUM_SNR_TO_CONTINUE:
            break
        keep = np.abs(time - epoch - transit_numbers * period) > duration
        time, flux = time[keep], flux[keep]
    return signals


def search_star(light_curve_path):
    tic_id = int(Path(light_curve_path).stem)
    try:
        time, flux = load_light_curve(light_curve_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            signals = search_light_curve(time, flux)
        if not signals:
            raise ValueError(f"fewer than {MINIMUM_POINT_COUNT} usable points")
    except Exception as error:
        return [{"tic_id": tic_id, "signal_number": 0, "error": f"{type(error).__name__}: {error}"}]
    return [{"tic_id": tic_id, **signal} for signal in signals]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sector", type=int)
    parser.add_argument("--workers", type=int, default=max(1, os.cpu_count() - 2))
    parser.add_argument("--restart", action="store_true", help="delete the chunks of an earlier run and search all stars")
    arguments = parser.parse_args()
    sector = arguments.sector
    if arguments.restart:
        for chunk_path in output_directory(sector).glob("chunk_*.csv"):
            chunk_path.unlink()

    light_curve_paths = sorted(
        (sector_directory(sector) / "light_curves").glob("*/*.npz"), key=lambda path: int(path.stem)
    )
    light_curve_paths = [path for path in light_curve_paths if not path.name.endswith(".partial.npz")]
    chunks = [light_curve_paths[i : i + CHUNK_SIZE] for i in range(0, len(light_curve_paths), CHUNK_SIZE)]
    output_directory(sector).mkdir(parents=True, exist_ok=True)
    print(f"Sector {sector}: {len(light_curve_paths)} light curves in {len(chunks)} chunks", flush=True)

    status = PipelineStatus(sector, 3, "bls_search", "search_qlp_sector.py", len(light_curve_paths), "stars")
    with Pool(arguments.workers) as pool, status:
        for chunk_number, chunk in enumerate(chunks):
            status.set_items_done(chunk_number * CHUNK_SIZE)
            chunk_path = output_directory(sector) / f"chunk_{chunk_number:04d}.csv"
            if chunk_path.exists():
                continue
            rows = [row for star_rows in pool.map(search_star, chunk, chunksize=4) for row in star_rows]
            temporary_path = chunk_path.with_suffix(".partial.csv")
            with open(temporary_path, "w", newline="") as file:
                writer = csv.DictWriter(file, fieldnames=OUTPUT_COLUMNS)
                writer.writeheader()
                writer.writerows(rows)
            temporary_path.rename(chunk_path)
            print(f"  chunk {chunk_number + 1} of {len(chunks)} done", flush=True)


if __name__ == "__main__":
    main()
