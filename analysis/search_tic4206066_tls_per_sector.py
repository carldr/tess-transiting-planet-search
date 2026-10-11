"""Search each archival sector of TIC 4206066 for transits without using Rabtsevich's period.

The script runs Transit Least Squares (TLS) on Sectors 6, 32 and 98 separately, over periods of 0.5 to 14 d.
For each sector it finds the strongest signal, masks its transits, and searches again, up to three times.
The test passes when the 3.18 d period appears in a sector without any prior knowledge of it.

Run from the repository root:
    .venv/bin/python analysis/search_tic4206066_tls_per_sector.py

Outputs, in results/tls_per_sector/:
- tls_signals_by_sector.json: period, epoch, depth, duration, SDE and transit count of each signal found;
- tls_periodogram_sector_<n>.png: the first periodogram of each sector.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from transitleastsquares import transit_mask, transitleastsquares

from tic4206066_light_curves import REPOSITORY_ROOT, detrend_sector_light_curve

OUTPUT_DIRECTORY = REPOSITORY_ROOT / "results/tls_per_sector"
SECTORS = (6, 32, 98)
SIGNALS_PER_SECTOR = 3
RABTSEVICH_PERIOD_SIGNAL_A_DAYS = 3.182785
RABTSEVICH_PERIOD_SIGNAL_B_DAYS = 11.13274

# Stellar values from Rabtsevich section 1 (Mann et al. relations); TLS uses them to set its duration grid.
STELLAR_RADIUS_SOLAR = 0.623
STELLAR_MASS_SOLAR = 0.606


def search_sector(sector):
    time, flux, flux_error = detrend_sector_light_curve(sector)
    signals = []
    for search_number in range(1, SIGNALS_PER_SECTOR + 1):
        model = transitleastsquares(time, flux, flux_error)
        result = model.power(
            period_min=0.5,
            period_max=14.0,
            R_star=STELLAR_RADIUS_SOLAR,
            R_star_min=0.5,
            R_star_max=0.75,
            M_star=STELLAR_MASS_SOLAR,
            M_star_min=0.5,
            M_star_max=0.75,
            oversampling_factor=5,
            duration_grid_step=1.05,
            use_threads=8,
            show_progress_bar=False,
        )
        signals.append(
            {
                "search_number": search_number,
                "period_days": float(result.period),
                "epoch_btjd": float(result.T0),
                "depth_ppm": float((1 - result.depth) * 1e6),
                "duration_hours": float(result.duration * 24),
                "signal_detection_efficiency": float(result.SDE),
                "transit_count_with_data": int(result.distinct_transit_count),
                "ratio_to_rabtsevich_period_signal_a": float(result.period / RABTSEVICH_PERIOD_SIGNAL_A_DAYS),
                "ratio_to_rabtsevich_period_signal_b": float(result.period / RABTSEVICH_PERIOD_SIGNAL_B_DAYS),
            }
        )
        if search_number == 1:
            plot_periodogram(sector, result)
        in_transit = transit_mask(time, result.period, 2 * result.duration, result.T0)
        time, flux, flux_error = time[~in_transit], flux[~in_transit], flux_error[~in_transit]
    return signals


def plot_periodogram(sector, result):
    figure, axes = plt.subplots(figsize=(10, 3.5))
    axes.plot(result.periods, result.power, color="black", linewidth=0.6)
    for period, label in (
        (RABTSEVICH_PERIOD_SIGNAL_A_DAYS, "Rabtsevich A 3.18 d"),
        (RABTSEVICH_PERIOD_SIGNAL_B_DAYS, "Rabtsevich B 11.13 d"),
    ):
        axes.axvline(period, color="tab:red", alpha=0.4, linestyle="--", label=label)
    axes.set_xscale("log")
    axes.set_xlabel("period (d)")
    axes.set_ylabel("SDE")
    axes.set_title(f"TIC 4206066, Sector {sector}: TLS periodogram (first search)")
    axes.legend(loc="upper left", fontsize=8)
    figure.tight_layout()
    figure.savefig(OUTPUT_DIRECTORY / f"tls_periodogram_sector_{sector}.png", dpi=130)
    plt.close(figure)


def main():
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    signals_by_sector = {}
    for sector in SECTORS:
        signals_by_sector[sector] = search_sector(sector)
        for signal in signals_by_sector[sector]:
            print(
                f"Sector {sector:3d} search {signal['search_number']}: "
                f"P = {signal['period_days']:.5f} d, depth = {signal['depth_ppm']:.0f} ppm, "
                f"duration = {signal['duration_hours']:.2f} h, SDE = {signal['signal_detection_efficiency']:.1f}, "
                f"transits = {signal['transit_count_with_data']}",
                flush=True,
            )
    with open(OUTPUT_DIRECTORY / "tls_signals_by_sector.json", "w") as file:
        json.dump(signals_by_sector, file, indent=1)


if __name__ == "__main__":
    main()
