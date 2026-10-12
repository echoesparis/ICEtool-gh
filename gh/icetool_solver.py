#! python3
# r: numpy, scipy, pandas
# icetool_gh component glue v3 — full-improvements mode ("gh")
"""GH ScriptComponent V3 (file-synced). Grasshopper only supplies geometry,
material tags and the Ladybug shadow study; ALL weather comes from the EPW file.

Inputs (set on the component):

  run       bool  item   solve toggle
  mat_ids   str   list   material tag per point, FULL resolution (one per mesh point/face)
  sun       num   tree   sun-exposure fraction per point per hour (1 = sunlit);
                         hour-major (Ladybug native) or point-major, 24 columns for the
                         analysed day or 8760 for a full-year run. Dedup of identical
                         (tag, sun-pattern) points happens inside the solver.
  epw_file  str   item   path to the .epw weather file (same file the shadow study uses)
  month,day num   item   the analysed day when sun has 24 columns (defaults 7/21)
  db_path   str   item   optional materials CSV override (empty -> bundled v2026.07 db)

Refinements over the QGIS-compat component (all derived from the EPW itself):
  - hc: hourly wind-derived (McAdams, EPW wind at 10 m -> 2 m) instead of constant 5
  - Tsky: from EPW measured horizontal infrared radiation instead of the Fuentes correlation
  - ET0: Penman-Monteith with real hourly wind instead of a fixed 1 km/h
  - Tsky Fuentes fallback clamped at 0 degC so winter/full-year runs survive

Outputs (index 0 of the component MUST stay the reserved 'out' stdout param):
  temp_c  tree  one branch per point, one item per hour, surface temperature degC
  min_c, mean_c, max_c  lists (N)
  report  str
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)  # repo root
for _p in (_ROOT, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
# fresh import so edits to the package are picked up without restarting Rhino
for _k in [k for k in list(sys.modules)
           if k.split(".")[0] in ("icetool_gh", "_glue_common")]:
    del sys.modules[_k]

from _glue_common import run_component

temp_c, min_c, mean_c, max_c, report = run_component(
    "gh", run, mat_ids, sun, epw_file, month, day, db_path)
