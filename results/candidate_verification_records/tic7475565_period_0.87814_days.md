# TIC 7475565, period 0.87814 d

| date (UTC) | sector | step | script | code version | verdict | result |
|---|---|---|---|---|---|---|
| 2026-10-28 01:51 | 102 | search, vetting and known-signal cross-match | cross_match_qlp_sector_signals_with_known_catalogues.py | 378b532 | new candidate | BLS SNR 17.6, MES 10.4, depth 2292 ppm, 24 transits; passes the LEO-Vetter false alarm and false positive tests; no TOI, CTOI, TCE or TESS EB match in the catalogues of 2026-10-28 |
| 2026-10-30 11:13 | 102 | independent data checks | check_qlp_sector_new_candidates.py | 492b76f | refuted | other sectors: S12 QLP: SNR 10.1, expected SNR 5.2; S39 QLP: SNR 16.7, expected SNR 9.7; S66 QLP: SNR 15.5, expected SNR 8.9; S93 QLP: SNR 13.8, expected SNR 8.4; joint search SNR 30.8 against highest null SNR 7.2; Gaia RUWE 4.19; companion radius 6.56929 Earth radii; refuted by: nearby star is the likely source; VSX binary; cautions: RUWE above 1.4; neighbour able to cause dip |
