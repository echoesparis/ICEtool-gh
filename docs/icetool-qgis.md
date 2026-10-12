# ICEtool for QGIS — source analysis

Notes on the upstream plugin, [Art-Ev/ICEtool](https://github.com/Art-Ev/ICEtool) at commit [`46f31db`](https://github.com/Art-Ev/ICEtool/tree/46f31db) / v2026.7.2, 2026-08-31 (vendored locally at commit time, not tracked in this repo). This is the reference implementation the Grasshopper port (`ee_icetool.gh` + [icetool_gh](../icetool_gh)) is adapting — see [How it works in README.md](../README.md#how-it-works) for the Rhino/GH side, and [qgis-vs-gh.md](qgis-vs-gh.md) for a full side-by-side comparison of the two implementations.

## What it is

A QGIS **Processing provider** plugin (GPL-3.0, author Arthur Evrard, Egis/INSA Toulouse). Registers 9 algorithms ([`ICEtool_provider.py`](https://github.com/Art-Ev/ICEtool/blob/46f31db/ICEtool_provider.py)) grouped into steps, run in order from the QGIS Processing Toolbox:

| Group | Algorithm | File | Purpose |
|---|---|---|---|
| Help | Generate an example project | [SeeEXAMPLE.py](https://github.com/Art-Ev/ICEtool/blob/46f31db/Scripts/SeeEXAMPLE.py) | Unpacks the bundled INSA-Toulouse example project |
| Help | Load an EPW file | [Display_EPW.py](https://github.com/Art-Ev/ICEtool/blob/46f31db/Scripts/Display_EPW.py) | Weather-file preview/import helper |
| Help | (fr/en help) | [Display_HELP_fr.py](https://github.com/Art-Ev/ICEtool/blob/46f31db/Scripts/Display_HELP_fr.py) / [_en.py](https://github.com/Art-Ev/ICEtool/blob/46f31db/Scripts/Display_HELP_en.py) | Opens the user manual PDFs in `Scripts/Docs/` |
| — | Create Project | [Create_Project.py](https://github.com/Art-Ev/ICEtool/blob/46f31db/Scripts/Create_Project.py) | Scaffolds a new project folder (Step_1..4 subfolders) by copying [`Scripts/ProjectModel/`](https://github.com/Art-Ev/ICEtool/tree/46f31db/Scripts/ProjectModel) (a template `.qgz` + geopackage), reprojected to the chosen CRS |
| Step_2 | Create rasters (tree points) | [Step2_TreePoints.py](https://github.com/Art-Ev/ICEtool/blob/46f31db/Scripts/Step2_TreePoints.py) | Rasterizes buildings (height field) + point-based trees (height/radius fields) into a combined canopy/obstruction raster for shadow-casting |
| Step_2 | Create rasters (tree poly) | [Step2_TreePoly.py](https://github.com/Art-Ev/ICEtool/blob/46f31db/Scripts/Step2_TreePoly.py) | Same, but trees given as polygons instead of points |
| Step_3 | Compute shadows [UMEP] | [Step3_UMEPshadows.py](https://github.com/Art-Ev/ICEtool/blob/46f31db/Scripts/Step3_UMEPshadows.py) | Thin wrapper that launches the bundled [`ShadowGenerator`](https://github.com/Art-Ev/ICEtool/tree/46f31db/Scripts/ShadowGenerator) dialog (itself adapted from [UMEP](https://github.com/UMEP-dev/UMEP)'s daily shading tool) to produce one shadow-fraction `.tif` per hour into `Step_3/` |
| Step_4 | Compute ground temperature (weather data as .epw) | [Step4_ComputeTemperatureEPW.py](https://github.com/Art-Ev/ICEtool/blob/46f31db/Scripts/Step4_ComputeTemperatureEPW.py) | **The core calculation** — registered and the one actually exposed in the Processing Toolbox |
| Step_4 | Compute ground temperature (weather data as csv) | [Step4_ComputeTemperatureCSV.py](https://github.com/Art-Ev/ICEtool/blob/46f31db/Scripts/Step4_ComputeTemperatureCSV.py) | CSV-weather variant of the same algorithm — **present in the repo but not imported/registered in `ICEtool_provider.py`**, so currently unreachable from the QGIS UI |

Supporting code: [`Scripts/Utilities/misc.py`](https://github.com/Art-Ev/ICEtool/blob/46f31db/Scripts/Utilities/misc.py) and [`shadowingfunctions.py`](https://github.com/Art-Ev/ICEtool/blob/46f31db/Scripts/Utilities/shadowingfunctions.py) (raster/shadow helpers used by `ShadowGenerator`), [`Scripts/Fuseaux_horaires.csv`](https://github.com/Art-Ev/ICEtool/blob/46f31db/Scripts/Fuseaux_horaires.csv) (timezone → longitude lookup used for solar-time correction), [`Scripts/point_style.qml`](https://github.com/Art-Ev/ICEtool/blob/46f31db/Scripts/point_style.qml) (QGIS styling applied to the result layer).

## Step 4 in detail — `Step4_ComputeTemperatureEPW.py`

This is the algorithm the port's [icetool_gh](../icetool_gh) solver is a transposition of. Inputs: a ground-description polygon layer (fields `Material, alb, em, Cv, lambd, ep, kc, FixedTemp[degC]`), a buildings polygon layer, an `.epw` weather file, target `day`/`month`, a timezone enum, and two advanced params (`altitude`, `spatial accuracy`).

Pipeline (`processAlgorithm`, 6 feedback steps):

1. **Point grid** — `native:creategrid` over the ground layer's extent at `2/spatialaccuracy` m spacing.
2. **Material join** — spatial join against the ground-description polygons to tag each grid point with its material + physical properties.
3. **Exclude buildings** — `native:extractbylocation` keeps only points disjoint from the buildings layer; a field-calculator step stamps each point with local (`x,y`) and reprojected WGS84 (`Long,Lat`) coordinates.
4. **Shadow sampling** — reads every hourly `.tif` from `Step_3/` (produced by the UMEP shadow step, named so the hour is parsed out of the filename), and samples each grid point's shadow fraction per hour via `gdal.InvGeoTransform` + array indexing (no interpolation) into a 24-value `Shadow1` tuple per point.
5. **Problem simplification** — points are grouped (deduplicated) by a `key = Material + shadow-pattern` composite, so the expensive iterative solve below only runs once per unique (material, shadow-history) combination rather than once per grid point — a perf optimization [icetool_gh](../icetool_gh) replicates too (see [qgis-vs-gh.md](qgis-vs-gh.md)).
6. **Weather** — parses the `.epw` (skips the header block until the first numeric row), computes sky temperature via the empirical Fuentes (1987) correlation, and extracts the single `day`/`month`'s 24 hourly values for air temp, global horizontal radiation, and relative humidity.
7. **Deep-soil temperature (`Tint`)** — instead of a fixed constant, a sinusoidal annual-cycle model estimates ground temperature at 20 cm depth per hour from the year's mean temperature (`Tyear`) and amplitude (`deltaT`) at that hour, the material's thermal diffusivity `Dh = lambd/Cv * 86400`, and day-of-year `t` (the "Cableizer" approximation referenced in the main README).
8. **Evapotranspiration (`ETO`)** — a full **Penman-Monteith** implementation (`Evapo()`), day/night branches, using solar geometry (hour angle, declination, extraterrestrial radiation), wind fixed at 1 km/h, humidity, and albedo. This feeds the energy balance as a genuine latent-heat term — see below for how this differs from the GH port.
9. **Thermal equilibrium solve** — same physical model as [icetool_gh](../icetool_gh): `A + B·T + C·T⁴ = 0`, split as
   - `B = hc + lambd/ep + Cv·ep/3600` (constant per point)
   - `C = em · σ`
   - `A = -(0.8·Gh·Shadow·(1-alb) + 0.2·Gh·(1-alb)) - em·σ·Tsky⁴ - hc·Tair - (lambd/ep)·Tint - (Cv·ep/3600)·T0 + ETO·kc`

   solved hour-by-hour with `scipy.optimize.root`, chained (each hour's solution seeds the next as `T0`), for **one representative day only** (not a full year), iterated at least 2 and at most 25 times around the 24-hour loop until `|T[23h] - T0| < 0.5°C`. A per-material `FixedTemp[degC]` column can bypass the solve entirely and hard-set all 24 hours to a fixed value (used for e.g. water bodies).
10. **Output** — merges results back onto every original grid point (not just the deduplicated set), writes `Step_4/ComputedPoints.csv` (`id,x,y,Temp_DegC,min_DegC,mean_DegC,max_DegC`), and loads it into QGIS as a styled point layer via `point_style.qml`.

For how this compares to the Grasshopper port's calculation, see [qgis-vs-gh.md](qgis-vs-gh.md).

## Version / provenance

- Upstream repo: [Art-Ev/ICEtool](https://github.com/Art-Ev/ICEtool.git), HEAD [`46f31dbf06f279246b476ea63876bf73bc3859b9`](https://github.com/Art-Ev/ICEtool/commit/46f31dbf06f279246b476ea63876bf73bc3859b9) ("Update README.md", 2026-08-31), matching `metadata.txt` version `2026.7.0` / `CITATION.cff` version `2026.7.2`.
- Companion asset repo (material databases, diagrams): [Art-Ev/ICEtool_sources](https://github.com/Art-Ev/ICEtool_sources).
- Runtime deps for the calculation code itself (excluding QGIS/GDAL, which come from the QGIS Python environment): `pandas`, `scipy`, `numpy` ([requirements.txt](https://github.com/Art-Ev/ICEtool/blob/46f31db/requirements.txt)).
