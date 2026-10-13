# TIC 4200338, period 0.96726 d

| date (UTC) | sector | step | script | code version | verdict | result |
|---|---|---|---|---|---|---|
| 2026-10-12 21:32 | 98 | search, vetting and known-signal cross-match | cross_match_qlp_sector_signals_with_known_catalogues.py | 49142f8 | new candidate | BLS SNR 15.3, MES 12.4, depth 173 ppm, 43 transits; passes the LEO-Vetter false alarm and false positive tests; no TOI, CTOI, TCE or TESS EB match in the catalogues of 2026-10-09 |
| 2026-10-13 09:46 | 98 | independent data checks | check_qlp_sector_new_candidates.py | 9674c6b | refuted | other sectors: S32 TESS-SPOC: SNR 26.4, expected SNR 11.8; joint search SNR 25.7 against highest null SNR 6.4; Gaia RUWE 1.04; companion radius 3.94494 Earth radii; refuted by: nearby star is the likely source; cautions: neighbour able to cause dip |
