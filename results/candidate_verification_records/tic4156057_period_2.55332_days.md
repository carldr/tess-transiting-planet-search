# TIC 4156057, period 2.55332 d

| date (UTC) | sector | step | script | code version | verdict | result |
|---|---|---|---|---|---|---|
| 2026-10-12 21:32 | 98 | search, vetting and known-signal cross-match | cross_match_qlp_sector_signals_with_known_catalogues.py | 49142f8 | new candidate | BLS SNR 40.6, MES 21.2, depth 2283 ppm, 16 transits; passes the LEO-Vetter false alarm and false positive tests; no TOI, CTOI, TCE or TESS EB match in the catalogues of 2026-10-09 |
| 2026-10-13 09:12 | 98 | independent data checks | check_qlp_sector_new_candidates.py | 9674c6b | refuted | other sectors: S32 QLP: SNR 68.6, expected SNR 5.3; joint search SNR 129.3 against highest null SNR 33.3; Gaia RUWE 1.3; companion radius 21.3051 Earth radii; refuted by: nearby star is the likely source; Gaia eclipsing binary neighbour; VSX binary; cautions: TIC stellar parameters missing; neighbour able to cause dip |
