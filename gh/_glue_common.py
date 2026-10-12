"""Shared helpers for the two ICEtool GH script components (imported by both)."""

import os
import time

import numpy as np

try:
    from Grasshopper import DataTree
    from Grasshopper.Kernel.Data import GH_Path
except ImportError:  # allows py_compile / harness import outside Rhino
    DataTree = None
    GH_Path = None


def unwrap(x):
    """Unwrap an IGH_Goo (script inputs without a type hint arrive wrapped)."""
    return x.Value if hasattr(x, "Value") else x


def num(x):
    """Coerce a (possibly goo-wrapped) item to float."""
    x = unwrap(x)
    if isinstance(x, (int, float)):
        return float(x)
    return float(str(x))


def _branch_floats(branch):
    """Fast branch -> list[float], falling back to per-item goo coercion."""
    try:
        return np.fromiter(branch, dtype=float, count=branch.Count)
    except (TypeError, ValueError):
        return np.array([num(v) for v in branch], dtype=float)


def tree_to_array(tree):
    """Grasshopper DataTree -> 2D numpy array (one row per branch); lists pass through."""
    if tree is None:
        return None
    if hasattr(tree, "Paths") and hasattr(tree, "Branch"):
        rows = [_branch_floats(tree.Branch(p)) for p in tree.Paths]
        if not rows:
            return None
        return np.asarray(rows, dtype=float)
    return np.asarray([num(v) for v in tree], dtype=float)


def orient(arr, n_points):
    """Ensure a (points, hours) 2D array. GH trees from Ladybug are often
    hour-major (one branch per hour, one item per point) — transpose those.
    A square tree (n_points equals the hour count) cannot be told apart and is
    rejected rather than guessed."""
    if arr is None:
        return None
    arr = np.atleast_2d(arr)
    if arr.shape[0] == arr.shape[1] == n_points and n_points > 1:
        raise ValueError(
            "sun tree is %dx%d: with as many points as hours its orientation "
            "(hour-major or point-major) is ambiguous; change the number of "
            "analysis points (e.g. the grid size)" % arr.shape)
    if arr.shape[0] != n_points and arr.shape[1] == n_points:
        arr = arr.T
    if arr.shape[0] != n_points:
        raise ValueError(
            "sun tree is %dx%d but mat_ids has %d points: one of its two "
            "dimensions must equal the point count" % (arr.shape + (n_points,)))
    return arr


def array_to_tree(arr):
    t = DataTree[float]()
    for i, row in enumerate(np.atleast_2d(arr)):
        t.AddRange(row.tolist(), GH_Path(i))
    return t


def run_component(mode, run, mat_ids, sun, epw_file, month, day, db_path):
    """Shared body of the two ICEtool script components.

    Returns (temp_c_tree, min_c, mean_c, max_c, report). The mode string is the
    ONLY difference between the components — everything else lives here once.
    """
    if not run:
        return None, None, None, None, "run = False"
    from icetool_gh import load_materials, solve_surface_temperature, read_epw
    from icetool_gh.epw import read_epw_location
    from icetool_gh.solver import Weather
    import icetool_gh

    t0 = time.time()
    ids = [str(x) for x in mat_ids]
    epw_path = str(unwrap(epw_file))
    if not os.path.isfile(epw_path):
        raise ValueError("epw_file not found: %r" % epw_path)
    lat, lon, tz_lon, alt = read_epw_location(epw_path)
    weather = Weather.from_epw(read_epw(epw_path))

    month_used = int(month) if month is not None else 7
    day_used = int(day) if day is not None else 21
    res = solve_surface_temperature(
        ids, orient(tree_to_array(sun), len(ids)),
        load_materials(db_path if db_path else None), weather,
        latitude=lat, longitude=lon, tz_meridian_longitude=tz_lon, altitude_m=alt,
        month=month_used, day=day_used, mode=mode)

    tc = res.temp_C
    label = "gh (EPW wind hc, IR Tsky, wind ET0)" if mode == "gh" else "qgis_compat"
    if tc.shape[1] != 24:
        when = "full year"
    else:
        when = "day %d/%d%s" % (day_used, month_used,
                                " (default, month/day not set)" if month is None or day is None else "")
    report = "icetool_gh v%s | mode=%s | %d points -> %d groups x %d h | %s | %s | %.1f s | non-converged: %s" % (
        icetool_gh.__version__, label, len(ids), res.n_groups, tc.shape[1],
        when, os.path.basename(epw_path), time.time() - t0,
        res.non_converged or "none")
    return (array_to_tree(tc), res.min_C.tolist(), res.mean_C.tolist(),
            res.max_C.tolist(), report)
