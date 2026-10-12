# ICEtool: QGIS vs. Grasshopper — implementation comparison

Side-by-side comparison of the upstream QGIS plugin ([Art-Ev/ICEtool](https://github.com/Art-Ev/ICEtool), analyzed in [icetool-qgis.md](icetool-qgis.md)) and this repo's Rhino/Grasshopper port. Both implement the same physical model — `Q_R = Q_H + Q_L + Q_C + δQ_S`, solved per point as `A + B·T + C·T⁴ = 0`.

**Status (2026-09-01): numerical parity achieved — including against real plugin output.**

- *Identical inputs*: `tests/compare_qgis_example.py --json` (local dev-only, not tracked in this repo) compares `qgis_compat` mode with the actual QGIS result shipped upstream (`Scripts/Example/Step_4/ComputedPoints.csv`, INSA-Toulouse, 68 949 points, computed for July 1): fully-sunlit points match to **0.005 °C mean** (the CSV's own rounding). Remaining spread for partially-shaded points (median 0.75 °C) is exactly the per-hour shadow information the example does not ship (Step_3 rasters are absent), not physics error.
- *Full GH pipeline*: the same site was rebuilt as `data/qgis-example.3dm` (gitignored, 14 MB, generated — rebuilt locally from the upstream GeoPackage via a dev-only script; note the ground breps must face **up**, or Ladybug's offset grid sees no sun; fixed 2026-09-01) and solved end-to-end in `ee_icetool.gh` with the shipped Toulouse EPW: 47 414 mesh points → 5 568 unique groups, 9 s. Agreement with the QGIS output: **hottest 24 h curve within 0.009 °C**, min/mean/max distribution percentiles within ~0.3–1.9 °C, most-shaded curve within 0.44 °C. Residuals are Ladybug-vs-UMEP shading differences (binary canopy occlusion vs raster transmissivity, mesh vs 1 m grid), not energy-balance differences.

 The port was rewritten as the [icetool_gh](../icetool_gh) package running in-canvas (Rhino 8 GH Python 3, ScriptEditor v3 file-synced components — no Flask/Hops). Its `qgis_compat` mode reproduces upstream Step 4 to ≤0.05 °C/hour, verified by `tests/test_parity.py` (local dev-only) against a vendored copy of the upstream physics (`tests/reference/qgis_step4_ref.py`, drift-guarded by `test_reference_guard.py`).

## Historical divergences of the old port (`app.py`, since deleted) — all fixed

| Aspect | ICEtool (QGIS) | old `app.py` (deleted) | now (`icetool_gh`) |
|---|---|---|---|
| Stefan-Boltzmann `σ` | `5.67e-08` | `5.57e-08` (typo) | **fixed**: `5.67e-08` |
| Evapotranspiration `ETO` | full hourly Penman-Monteith | raw precipitation mm fed into the `ETo` slot | **fixed**: Penman-Monteith ported (`physics.evapo_penman_monteith`) |
| Sky temperature `Tsky` | Fuentes (1987) per hour | one annual scalar `mean(Tair)+1` | **fixed**: hourly Fuentes (`physics.tsky_fuentes`) |
| Diffuse radiation albedo | `a0 = −(0.8·Gh·sh·(1−alb) + 0.2·Gh·(1−alb))` | `(1−alb)` missing on the 0.2 diffuse term | **fixed**: both terms carry `(1−alb)` |
| Convergence | iterate the 24 h cycle until `\|T[23]−T0\| < 0.5 °C`, 2–25 iterations, shadow-aware seeds | compared two adjacent *hours* (not iterations) at a loose 1 % relative threshold — effectively always exited after 3 passes | **fixed**: upstream loop reproduced exactly (`physics.solve_periodic_day`) |
| Ground temperature `Tint` | Cableizer annual sinusoid at 20 cm, per hour, in Kelvin | single user constant (default 35 °C) | **fixed**: Cableizer model ported (`physics.tint_cableizer`) |
| `FixedTemp[degC]` override (water…) | skips the solve | absent | **fixed**: supported via the v2026.07 material DB |
| Radiation units | Wh/m² per hour from EPW | comment claimed a kWh→Wh conversion the code didn't perform | **fixed**: solver takes Wh/m²; GH glue converts Ladybug kWh explicitly (`rad_units`) |
| Material database | v2026.07 schema (`category`, `FixedTemp[degC]`, provenance columns, revised values) | pre-2026 CSV fetched live from GitHub, parsed positionally | **fixed**: v2026.07 vendored locally, parsed by column name, best `database_version` per ID |

Point deduplication note: QGIS dedups by `(material, shadow-pattern)` inside Step 4; the GH canvas does the equivalent in the `'reduced'` group *before* the solver and remaps after (`'remap from reduced'` / the solver's `remap` input) — same optimization, different layer.

## Intentional differences that remain (by design)

All weather now comes from the EPW file itself (`epw_file` input on both components); the "gh"-mode refinements below are derived from EPW columns the QGIS plugin reads but does not use:

| Aspect                    | ICEtool (QGIS) & compat component                                                                                           | ICEtool-gh solver component ("gh" mode)                                                                                                                                                                                                                          |
| ------------------------- | --------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Convection `hc`           | fixed 5 W/m²K                                                                                                               | hourly wind-derived (McAdams `5.7 + 3.8·u₂`, EPW 10 m wind reduced to 2 m via FAO-56); the Python API also accepts an explicit `hc`                                                                                                                              |
| Sky temperature `Tsky`    | Fuentes (1987) empirical correlation from air temperature — upstream's own comment flags it "to replace in future versions" | `(IR/σ)^0.25` from the EPW's measured **Horizontal Infrared Radiation Intensity**; falls back to Fuentes clamped at 0 °C (the correlation is NaN for sub-zero air)                                                                                               |
| Evapotranspiration wind   | fixed 1 km/h assumption inside Penman-Monteith                                                                              | real hourly EPW wind (same 2 m reduction)                                                                                                                                                                                                                        |
| Time coverage             | one representative day (24 h)                                                                                               | that same day by default (Ladybug analysis period), or the full 8760 h year when fed year-long sun trees (day-1 spin-up, then a single forward pass)                                                                                                             |
| Radiation                 | EPW global horizontal × shadow fraction                                                                                     | same; the Python API additionally accepts per-point `incident_rad_Wh_m2`                                                                                                                                                                                         |
| Ground temperature `Tint` | Cableizer annual sinusoid, always                                                                                           | same by default; the Python API accepts `tint_override_C` (e.g. measured ground temperature — successor of the old app.py `Tn` input, now an override instead of the only model)                                                                                 |
| Precipitation             | not used by Step 4                                                                                                          | not used either. The old app.py fed raw precipitation mm into the `ETo` slot — that was a bug, not a feature; latent heat is now the actual Penman-Monteith `ET0`. A soil-moisture/rain modulation of `kc` would be a genuine future improvement, in both tools. |

The **QGIS-compat component** (`gh/icetool_qgis_compat.py` in [gh](../gh)) disables all of these and matches QGIS exactly — use it for validation runs.

## Architecture (current)

| | ICEtool (QGIS) | ICEtool-gh |
|---|---|---|
| Runtime | QGIS Processing plugin (in-process) | Rhino 8 GH Python 3 script components (file-synced from [gh](../gh)), importing [icetool_gh](../icetool_gh) — no server |
| Workflow | 5 manual Processing steps with intermediate files | one `.gh` definition, one recompute |
| Shadow/solar | UMEP ShadowGenerator rasters, one day | Ladybug `LB Direct Sun Hours` per point |
| Output | `ComputedPoints.csv` (min/mean/max per point) | trees/lists in-graph (`temp_c`, `min_c/mean_c/max_c`) |

The previous architecture (Flask + ghhops server called through two nested Hops components, `app.py`) has since been deleted outright, not archived — see the divergences table above for what it got wrong.

### Canvas note: sun-tree hour gap-fill

`ee_icetool.gh`'s "shadows" group reconciles two Ladybug outputs before feeding the solver's `sun` input: `LB Analysis Period.hoys` (all 24 Hour-Of-Year values for the analysed day) and `LB SunPath.hoys` (only the hours the sun is actually above the horizon — e.g. 15 of 24 on the Toulouse example, HOY 4350–4364 = 06:00–20:00 on July 1). `LB Direct Sun Hours` is only ever given the daylight sun vectors, so its `results` tree has one **branch per hour actually tested**, not one per point — it never "drops" points, it was simply never asked about the missing hours (00:00–05:00, 21:00–23:00 in that example). The group fills those missing hour-branches with an all-zero sun-fraction row (correct — there is no sun at night) and merges them back in by their HOY value, restoring a complete 24-hour tree in the order the solver expects. (An earlier draft of this note wrongly described this as points being dropped for having zero sun all day; verified against the live `hoys` data 2026-09-01 — it is hours, not points, that are missing, and the fill targets nighttime, not unlit geometry.)

See [icetool-qgis.md](icetool-qgis.md) for the upstream source analysis.
