#! python3
# r: numpy, scipy, pandas
# icetool_gh component glue v3 — QGIS-compat validation mode (exact upstream Step 4 physics)
"""GH ScriptComponent V3 (file-synced). Grasshopper only supplies geometry,
material tags and the Ladybug shadow study; ALL weather comes from the EPW file.
Reproduces ICEtool QGIS Step 4 exactly: hc = 5 W/m2K, one 24 h representative
day, Fuentes sky temperature, Cableizer ground temperature, Penman-Monteith
ET0 with fixed 1 km/h wind, |dT| < 0.5 degC convergence over 2-25 daily cycles.

Inputs (set on the component):

  run       bool  item   solve toggle
  mat_ids   str   list   material tag per point, FULL resolution (one per mesh point/face)
  sun       num   tree   sun-exposure fraction per point per hour of the day (24 h);
                         hour-major (Ladybug native) or point-major. Dedup of identical
                         (tag, sun-pattern) points happens inside the solver.
  epw_file  str   item   path to the .epw weather file (same file the shadow study uses)
  month,day num   item   representative day to solve (as in QGIS Step 4; defaults 7/21)
  db_path   str   item   optional materials CSV override (empty -> bundled v2026.07 db)

Outputs (index 0 of the component MUST stay the reserved 'out' stdout param):
  temp_c  tree  one branch per point, 24 items, surface temperature degC
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
    "qgis_compat", run, mat_ids, sun, epw_file, month, day, db_path)
