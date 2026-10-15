"""Run the ellipsoidal variation fit of the all-sector verification on known transiting planets of hot stars.

The known planets show how often the fit gives a false detection, and the known massive companions show whether it
recovers a real variation. The script:
1. downloads from the NASA Exoplanet Archive the transiting planets with a mass, a host hotter than
   MINIMUM_STELLAR_TEMPERATURE_KELVIN, a TESS magnitude below MAXIMUM_TESS_MAGNITUDE and a period below
   MAXIMUM_PERIOD_DAYS;
2. computes the ellipsoidal amplitude that each planet gives, from its mass, a/R* and inclination;
3. loads up to MAXIMUM_SECTORS_PER_PLANET sectors of each host, with the loader of the all-sector verification, and
   fits the ellipsoidal amplitude of each sector with ellipsoidal_variation.py;
4. prints the spread of the significance of the planets with a predicted amplitude below
   LOW_PREDICTED_AMPLITUDE_PPM, and the measured and predicted amplitudes of the others.

Run from the repository root:
    .venv/bin/python analysis/calibrate_ellipsoidal_variation_test_on_known_planets.py

Output, in results/ellipsoidal_variation_test_calibration/:
- known_hot_star_planets.csv: the Exoplanet Archive query result;
- ellipsoidal_variation_of_known_planets.csv: one row per planet and sector.
"""

import csv
import io
import urllib.parse
import urllib.request
import warnings

import numpy as np

from build_qlp_sector_star_list import REPOSITORY_ROOT
from check_qlp_sector_new_candidates import load_other_sector_light_curve, other_sector_products
from ellipsoidal_variation import ELLIPSOIDAL_COEFFICIENT, ellipsoidal_amplitude

MINIMUM_STELLAR_TEMPERATURE_KELVIN = 6500
MAXIMUM_TESS_MAGNITUDE = 11
MAXIMUM_PERIOD_DAYS = 12
MAXIMUM_SECTORS_PER_PLANET = 4
LOW_PREDICTED_AMPLITUDE_PPM = 15
JUPITER_MASS_IN_SOLAR_MASSES = 9.546e-4
EXOPLANET_ARCHIVE_TAP_URL = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
EXOPLANET_ARCHIVE_QUERY = (
    "select pl_name, tic_id, sy_tmag, st_teff, st_mass, pl_orbper, pl_tranmid, pl_trandur, pl_bmassj, pl_ratdor, "
    "pl_orbincl, pl_orbeccen from pscomppars where tran_flag = 1 and pl_bmassj is not null "
    f"and st_teff > {MINIMUM_STELLAR_TEMPERATURE_KELVIN} and sy_tmag < {MAXIMUM_TESS_MAGNITUDE} "
    f"and pl_orbper < {MAXIMUM_PERIOD_DAYS}"
)
OUTPUT_DIRECTORY = REPOSITORY_ROOT / "results/ellipsoidal_variation_test_calibration"
OUTPUT_COLUMNS = [
    "planet_name", "tic_id", "sector", "author", "period_days", "orbital_eccentricity", "predicted_amplitude_ppm",
    "measured_amplitude_ppm", "measured_amplitude_error_ppm", "significance",
]


def known_hot_star_planets():
    url = f"{EXOPLANET_ARCHIVE_TAP_URL}?{urllib.parse.urlencode({'query': EXOPLANET_ARCHIVE_QUERY, 'format': 'csv'})}"
    with urllib.request.urlopen(url, timeout=120) as response:
        text = response.read().decode()
    (OUTPUT_DIRECTORY / "known_hot_star_planets.csv").write_text(text)
    return list(csv.DictReader(io.StringIO(text)))


def predicted_amplitude(planet):
    inclination = float(planet["pl_orbincl"] or 90)
    mass_ratio = float(planet["pl_bmassj"]) * JUPITER_MASS_IN_SOLAR_MASSES / float(planet["st_mass"])
    return (ELLIPSOIDAL_COEFFICIENT * mass_ratio / float(planet["pl_ratdor"]) ** 3
            * np.sin(np.radians(inclination)) ** 2)


def measure_planet(planet, writer, rng):
    tic_id = int(planet["tic_id"].split()[-1])
    period, epoch = float(planet["pl_orbper"]), float(planet["pl_tranmid"]) - 2457000
    duration = (float(planet["pl_trandur"]) / 24 if planet["pl_trandur"]
                else period / (np.pi * float(planet["pl_ratdor"])))
    predicted = predicted_amplitude(planet)
    # Sector 0 does not exist, so other_sector_products returns every sector.
    for sector, author, product in other_sector_products(tic_id, 0)[:MAXIMUM_SECTORS_PER_PLANET]:
        time, flux = load_other_sector_light_curve(product, author)
        result = ellipsoidal_amplitude(time, flux, period, epoch, duration, rng)
        if result is None:
            continue
        amplitude, error = result
        row = {
            "planet_name": planet["pl_name"], "tic_id": tic_id, "sector": sector, "author": author,
            "period_days": period, "orbital_eccentricity": planet["pl_orbeccen"],
            "predicted_amplitude_ppm": f"{predicted * 1e6:.1f}", "measured_amplitude_ppm": f"{amplitude * 1e6:.1f}",
            "measured_amplitude_error_ppm": f"{error * 1e6:.1f}", "significance": f"{amplitude / error:.2f}",
        }
        writer.writerow(row)
        print(f"{planet['pl_name']} S{sector} {author}: predicted {row['predicted_amplitude_ppm']} ppm, measured "
              f"{row['measured_amplitude_ppm']} +/- {row['measured_amplitude_error_ppm']} ppm", flush=True)


def print_summary(rows):
    low = [row for row in rows if float(row["predicted_amplitude_ppm"]) < LOW_PREDICTED_AMPLITUDE_PPM]
    significances = np.array([float(row["significance"]) for row in low])
    print(f"\n{len(low)} sectors of {len({row['planet_name'] for row in low})} planets with a predicted amplitude below "
          f"{LOW_PREDICTED_AMPLITUDE_PPM} ppm: significance mean {significances.mean():.2f}, "
          f"standard deviation {significances.std():.2f}")
    for limit in (2, 3, 4):
        planets = sorted({row["planet_name"] for row, value in zip(low, significances) if value >= limit})
        print(f"  {np.sum(significances >= limit)} sectors at {limit} sigma or more: {', '.join(planets) or 'none'}")
    print(f"Planets with a predicted amplitude of {LOW_PREDICTED_AMPLITUDE_PPM} ppm or more:")
    for row in rows:
        if float(row["predicted_amplitude_ppm"]) >= LOW_PREDICTED_AMPLITUDE_PPM:
            print(f"  {row['planet_name']} S{row['sector']}: predicted {row['predicted_amplitude_ppm']}, measured "
                  f"{row['measured_amplitude_ppm']} +/- {row['measured_amplitude_error_ppm']} ppm, "
                  f"eccentricity {row['orbital_eccentricity'] or 'unknown'}")


def main():
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    planets = known_hot_star_planets()
    print(f"{len(planets)} known transiting planets of hot stars", flush=True)
    rng = np.random.default_rng(1)
    output_path = OUTPUT_DIRECTORY / "ellipsoidal_variation_of_known_planets.csv"
    with open(output_path, "w", newline="") as file:
        writer = csv.DictWriter(file, OUTPUT_COLUMNS)
        writer.writeheader()
        for planet in planets:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    measure_planet(planet, writer, rng)
            except Exception as error:
                print(f"{planet['pl_name']}: {type(error).__name__}: {error}", flush=True)
            file.flush()
    print_summary(list(csv.DictReader(open(output_path))))
    print(f"Wrote {output_path.relative_to(REPOSITORY_ROOT)}")


if __name__ == "__main__":
    main()
