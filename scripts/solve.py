# /// script
# requires-python = ">=3.9"
# dependencies = [
#   "numpy>=1.22,<3",
#   "scipy>=1.8,<2",
#   "pandas>=1.4,<3",
# ]
# ///
"""Headless ICEtool-gh solve: surface temperatures without Rhino/Grasshopper.

Same solver and EPW handling as the two GH components; the sun-exposure
fractions the Ladybug shadow study supplies on the canvas come here from
--sun (one constant) or --points (one row per point).
"""

import argparse
import csv
import json
import os
import sys
from typing import NoReturn

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

EXIT_SOLVE, EXIT_USAGE, EXIT_INPUT, EXIT_NOT_CONVERGED = 1, 2, 3, 4

EPILOG = """\
points file (--points): CSV without header, one row per point:
  <material ID>,<sun h1>,...,<sun h24>      (24 values: one day, needs --month/--day)
  <material ID>,<sun h1>,...,<sun h8760>    (8760 values: full year, --mode gh only)
sun values are exposure fractions, 1 = sunlit, 0 = shaded.

stdout: one JSON object (run metadata, overall min/mean/max, per-point
min/mean/max limited by --limit/--offset). Hourly curves go to --output only
(CSV: point,material,h1..hN in degC). Diagnostics go to stderr.

exit codes:
  0  solved
  1  solver error
  2  invalid arguments
  3  input file missing or unreadable (EPW, points file, materials CSV)
  4  solved, but some points did not converge (see "non_converged")

examples:
  uv run scripts/solve.py --list-materials
  uv run scripts/solve.py --epw data/ToulouseWeatherData.epw --materials CC3,VG1 --month 7 --day 1
  uv run scripts/solve.py --epw data/ToulouseWeatherData.epw --materials CC3 --sun 0.25 \\
      --month 7 --day 1 --mode qgis_compat --output curves.csv
  uv run scripts/solve.py --epw data/ToulouseWeatherData.epw --points points.csv --month 7 --day 1
"""


def fail(code, msg) -> NoReturn:
    print("Error: " + msg, file=sys.stderr)
    sys.exit(code)


def build_parser():
    p = argparse.ArgumentParser(
        prog="scripts/solve.py", description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--epw", metavar="FILE", help="EPW weather file (all weather and site location are read from it)")
    src = p.add_mutually_exclusive_group()
    src.add_argument("--materials", metavar="IDS", help="comma-separated material IDs, one point each (see --list-materials)")
    src.add_argument("--points", metavar="FILE", help="CSV of points with per-hour sun fractions (format below)")
    p.add_argument("--sun", type=float, default=1.0, metavar="FRACTION",
                   help="constant sun fraction 0..1 for --materials points, over one 24 h day (default: 1.0)")
    p.add_argument("--month", type=int, choices=range(1, 13), metavar="1-12", help="analysed month (required for a 24 h solve)")
    p.add_argument("--day", type=int, choices=range(1, 32), metavar="1-31", help="analysed day (required for a 24 h solve)")
    p.add_argument("--mode", choices=["gh", "qgis_compat"], default="gh",
                   help="gh: EPW-derived refinements; qgis_compat: exact upstream Step 4 (default: gh)")
    p.add_argument("--db", metavar="FILE", help="materials CSV override (default: bundled v2026.07 database)")
    p.add_argument("--output", metavar="FILE", help="write hourly curves as CSV to FILE ('-' for stdout, replaces the JSON)")
    p.add_argument("--limit", type=int, default=50, metavar="N", help="per-point rows in the JSON summary (default: 50)")
    p.add_argument("--offset", type=int, default=0, metavar="N", help="first per-point row in the JSON summary (default: 0)")
    p.add_argument("--list-materials", action="store_true", help="print the material IDs of the database as JSON and exit")
    return p


def read_points(path):
    ids, rows = [], []
    with open(path, newline="") as f:
        for n, row in enumerate(csv.reader(f), 1):
            if not row or not row[0].strip():
                continue
            try:
                vals = [float(v) for v in row[1:]]
            except ValueError:
                fail(EXIT_INPUT, "%s line %d: sun values must be numbers (no header row expected)." % (path, n))
            ids.append(row[0].strip())
            rows.append(vals)
    if not ids:
        fail(EXIT_INPUT, "%s contains no points." % path)
    widths = sorted(set(len(r) for r in rows))
    if len(widths) != 1 or widths[0] not in (24, 8760):
        fail(EXIT_INPUT, "%s: every row needs the same number of sun values, 24 or 8760. Found: %s."
             % (path, ", ".join(map(str, widths))))
    return ids, rows


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    for label, path in (("--epw", args.epw), ("--points", args.points), ("--db", args.db)):
        if path and not os.path.isfile(path):
            fail(EXIT_INPUT, "%s file not found: %s" % (label, path))

    import numpy as np
    from icetool_gh import __version__, load_materials, read_epw, solve_surface_temperature
    from icetool_gh.epw import read_epw_location
    from icetool_gh.solver import Weather

    try:
        materials = load_materials(args.db)
    except Exception as e:
        fail(EXIT_INPUT, "cannot read materials database %s: %s" % (args.db or "(bundled)", e))

    if args.list_materials:
        json.dump([{"id": m.id, "name": m.name, "category": m.category,
                    "fixed_temp_degC": m.fixed_temp_degC or None}
                   for m in materials.values()], sys.stdout, indent=1)
        print()
        return 0

    if not args.epw:
        parser.error("--epw is required. Example: --epw data/ToulouseWeatherData.epw")
    if not args.materials and not args.points:
        parser.error("one of --materials or --points is required. Example: --materials CC3,VG1")
    if args.limit < 0 or args.offset < 0:
        parser.error("--limit and --offset must be >= 0.")

    if args.points:
        ids, rows = read_points(args.points)
        sun = np.array(rows, dtype=float)
    else:
        if not 0.0 <= args.sun <= 1.0:
            parser.error("--sun must be between 0 and 1. Received: %s" % args.sun)
        ids = [s.strip() for s in args.materials.split(",") if s.strip()]
        if not ids:
            parser.error("--materials is empty. Example: --materials CC3,VG1")
        sun = np.full((len(ids), 24), args.sun)
    if sun.min() < 0 or sun.max() > 1:
        fail(EXIT_INPUT, "sun fractions must be between 0 and 1 (found %g..%g)." % (sun.min(), sun.max()))

    hours = sun.shape[1]
    if hours == 24 and (args.month is None or args.day is None):
        parser.error("--month and --day are required for a 24 h solve. Example: --month 7 --day 1")
    if hours != 24 and args.mode == "qgis_compat":
        parser.error("--mode qgis_compat solves one 24 h day; the points file has %d sun values per row." % hours)
    unknown = sorted(set(ids) - set(materials))
    if unknown:
        parser.error("unknown material ID(s): %s. Valid IDs: %s"
                     % (", ".join(unknown), ", ".join(sorted(materials))))

    try:
        lat, lon, tz_lon, alt = read_epw_location(args.epw)
        weather = Weather.from_epw(read_epw(args.epw))
    except Exception as e:
        fail(EXIT_INPUT, "cannot read EPW file %s: %s" % (args.epw, e))

    print("solving %d points x %d h, mode=%s ..." % (len(ids), hours, args.mode), file=sys.stderr)
    try:
        res = solve_surface_temperature(
            ids, sun, materials, weather,
            latitude=lat, longitude=lon, tz_meridian_longitude=tz_lon, altitude_m=alt,
            month=args.month or 7, day=args.day or 21, mode=args.mode)
    except Exception as e:
        fail(EXIT_SOLVE, "solve failed: %s" % e)

    tc = res.temp_C
    if not np.all(np.isfinite(tc)):
        fail(EXIT_SOLVE, "solve produced non-finite temperatures (qgis_compat mode fails on days "
                         "with air temperature below 0 degC; use --mode gh).")

    if args.output:
        out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
        w = csv.writer(out)
        w.writerow(["point", "material"] + ["h%d" % (h + 1) for h in range(hours)])
        for n, row in enumerate(tc):
            w.writerow([n, ids[n]] + ["%.3f" % v for v in row])
        if out is not sys.stdout:
            out.close()

    rows = slice(args.offset, args.offset + args.limit)
    summary = {
        "icetool_gh": __version__, "mode": args.mode,
        "epw": os.path.basename(args.epw),
        "site": {"latitude": lat, "longitude": lon, "altitude_m": alt},
        "day": None if hours != 24 else "%02d-%02d" % (args.month, args.day),
        "hours": hours, "points": len(ids), "groups": res.n_groups,
        "non_converged": [{"point": int(n), "iterations": int(c)} for n, c in res.non_converged],
        "overall_C": {"min": round(float(tc.min()), 3), "mean": round(float(tc.mean()), 3),
                      "max": round(float(tc.max()), 3)},
        "output": None if not args.output else args.output,
        "points_offset": args.offset, "points_returned": len(ids[rows]),
        "per_point": [{"point": n, "material": ids[n],
                       "min_C": round(float(res.min_C[n]), 3),
                       "mean_C": round(float(res.mean_C[n]), 3),
                       "max_C": round(float(res.max_C[n]), 3)}
                      for n in range(len(ids))[rows]],
    }
    json.dump(summary, sys.stderr if args.output == "-" else sys.stdout, indent=1)
    print(file=sys.stderr if args.output == "-" else sys.stdout)
    return EXIT_NOT_CONVERGED if res.non_converged else 0


if __name__ == "__main__":
    sys.exit(main())
