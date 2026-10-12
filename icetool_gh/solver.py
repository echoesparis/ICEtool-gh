"""Surface-temperature solver: composes physics.py over points and hours.

Two modes:
- "qgis_compat": reproduces upstream ICEtool Step 4 exactly — hc=5, one
  representative day (24 h), periodic-steady-state iteration (2-25 passes,
  |T[23]-T0| < 0.5 degC). Used for validation against QGIS.
- "gh": full-year 8760 h. Per-day Tint/ET0/solar geometry as upstream, hc may
  be wind-derived per hour (scalar, (H,), or (N,H)). Convergence: the first
  day is iterated to periodic steady state (spin-up), then a single forward
  pass runs the whole year seeded from it — the yearly "cycle" replaces
  upstream's 24 h cycle.

Deduplication is built in (dedup=True): identical (material, sun-pattern, hc,
radiation) points are solved once and broadcast back, mirroring the QGIS
plugin's internal "simplified" groupby — feed the solver full-resolution
points; no external reduce/remap is needed.
"""

from dataclasses import dataclass, field

import numpy as np

from . import physics

KELVIN = 273.15


@dataclass
class Weather:
    """Full-year hourly weather (8760 rows), EPW conventions."""
    month: np.ndarray        # 1-12
    day: np.ndarray          # 1-31
    hour: np.ndarray         # 1-24 (EPW convention)
    tair_C: np.ndarray       # dry bulb, degC
    gh_Wh_m2: np.ndarray     # global horizontal radiation, Wh/m2 per hour
    rh_pct: np.ndarray       # relative humidity, %
    wind_mps: np.ndarray = None      # wind speed at 10 m, m/s (gh-mode refinements)
    ir_horiz_Wh_m2: np.ndarray = None  # horizontal infrared radiation intensity, Wh/m2

    @classmethod
    def from_epw(cls, df):
        return cls(
            month=df["month"].to_numpy(), day=df["day"].to_numpy(),
            hour=df["hour"].to_numpy(),
            tair_C=df["Dry Bulb Temperature [DegC]"].to_numpy(dtype=float),
            gh_Wh_m2=df["Global Horizontal Radiation [Wh/m2]"].to_numpy(dtype=float),
            rh_pct=df["Relative Humidity"].to_numpy(dtype=float),
            wind_mps=df["Wind Speed"].to_numpy(dtype=float),
            ir_horiz_Wh_m2=df["Horizontal Infrared Radiation Intensity"].to_numpy(dtype=float),
        )

    def day_slice(self, month, day):
        return (self.month == month) & (self.day == day)


@dataclass
class SolveResult:
    temp_K: np.ndarray               # (N, H)
    min_C: np.ndarray                # (N,)
    mean_C: np.ndarray               # (N,)
    max_C: np.ndarray                # (N,)
    non_converged: list = field(default_factory=list)  # (point_index, iterations)
    n_groups: int = 0                # unique (material, sun, hc, rad) groups solved

    @property
    def temp_C(self):
        return self.temp_K - KELVIN

    @classmethod
    def from_temp_K(cls, temp_K, non_converged, n_groups):
        temp_C = temp_K - KELVIN
        return cls(
            temp_K=temp_K,
            min_C=temp_C.min(axis=1),
            mean_C=np.round(temp_C.mean(axis=1), 2),
            max_C=temp_C.max(axis=1),
            non_converged=non_converged,
            n_groups=n_groups,
        )


def _hc_array(hc, mode, weather, sel, N, H):
    """Resolve the (N, H) convection-coefficient view for one solve window.

    compat -> upstream's constant 5; gh with no explicit hc -> hourly McAdams
    from the EPW wind (falling back to 5); explicit hc -> broadcast as given.
    """
    if mode == "qgis_compat":
        hc_in = physics.HC_QGIS
    elif hc is None:
        if weather.wind_mps is not None:
            wind = weather.wind_mps if sel is None else weather.wind_mps[sel]
            hc_in = physics.hc_mcadams(wind)
        else:
            hc_in = physics.HC_QGIS
    else:
        hc_in = hc
    return np.broadcast_to(np.asarray(hc_in, dtype=float), (N, H))


def solve_surface_temperature(
    point_material_ids,
    sun_fraction,                 # (N, H): 1 = fully sunlit (upstream "Shadow" semantics)
    materials,                    # dict[str, Material]
    weather,                      # Weather (always the full 8760 series)
    *,
    latitude,
    longitude,
    tz_meridian_longitude,        # longz: timezone meridian, deg (Fuseaux)
    altitude_m=100.0,
    day=21, month=7,              # representative day (compat mode)
    mode="gh",
    hc=None,                      # gh mode: scalar | (H,) | (N,H); ignored in compat (hc=5)
    incident_rad_Wh_m2=None,      # (N, H) per-point radiation incl. shading; replaces Gh*sun path
    tint_override_C=None,         # scalar or 24-length degC: replaces the Cableizer ground-temp
                                  # model (e.g. measured ground temperature); None = Cableizer
    threshold_degC=0.5,
    max_iter=25,
    dedup=True,
):
    if mode not in ("gh", "qgis_compat"):
        raise ValueError("mode must be 'gh' or 'qgis_compat'")
    ids = list(point_material_ids)
    N = len(ids)
    unknown = sorted({str(i) for i in ids if i not in materials})
    if unknown:
        raise ValueError("unknown material id(s) %s; valid IDs: %s"
                         % (", ".join(unknown), ", ".join(sorted(materials))))
    sun = np.asarray(sun_fraction, dtype=float).reshape(N, -1)

    if dedup and N > 1:
        # Group identical (material, sun-pattern, hc-row, rad-row) points, solve
        # one representative per group and broadcast — the same optimization the
        # QGIS plugin does internally (Step4 :232-238) and that the GH canvas
        # used to do with its 'reduce and remap' groups.
        rad_full = None if incident_rad_Wh_m2 is None else np.asarray(incident_rad_Wh_m2, dtype=float).reshape(N, -1)
        hc_rows = np.asarray(hc, dtype=float) if (hc is not None and np.ndim(hc) == 2) else None
        groups = {}
        inverse = np.empty(N, dtype=int)
        reps = []
        for n in range(N):
            key = (ids[n], sun[n].round(6).tobytes(),
                   None if hc_rows is None else hc_rows[n].round(6).tobytes(),
                   None if rad_full is None else rad_full[n].round(6).tobytes())
            g = groups.get(key)
            if g is None:
                g = groups[key] = len(reps)
                reps.append(n)
            inverse[n] = g
        rep = np.array(reps)
        sub = solve_surface_temperature(
            [ids[n] for n in reps], sun[rep], materials, weather,
            latitude=latitude, longitude=longitude,
            tz_meridian_longitude=tz_meridian_longitude, altitude_m=altitude_m,
            day=day, month=month, mode=mode,
            hc=hc_rows[rep] if hc_rows is not None else hc,
            incident_rad_Wh_m2=None if rad_full is None else rad_full[rep],
            tint_override_C=tint_override_C,
            threshold_degC=threshold_degC, max_iter=max_iter, dedup=False)
        return SolveResult.from_temp_K(sub.temp_K[inverse], sub.non_converged, len(reps))

    mats = [materials[i] for i in ids]
    l = physics.longitude_correction(tz_meridian_longitude, longitude)

    tair_K_year = weather.tair_C + KELVIN
    # Sky temperature: compat = upstream's exact Fuentes correlation (NaN-prone
    # below 0 degC). gh mode prefers the EPW's measured horizontal infrared
    # radiation (Tsky = (IR/sigma)^0.25) and falls back to clamped Fuentes.
    ir = weather.ir_horiz_Wh_m2
    if mode == "gh" and ir is not None and np.any(np.asarray(ir) > 0):
        tsky_year = physics.tsky_from_ir(ir)
    else:
        tsky_year = physics.tsky_fuentes(weather.tair_C, clamp_freezing=(mode == "gh"))
    tyear, deltaT = physics.annual_tyear_deltaT(tair_K_year, weather.hour)
    # gh-mode ET0 uses real hourly wind (compat keeps upstream's fixed 1 km/h)
    wind_year = weather.wind_mps if (mode == "gh" and weather.wind_mps is not None) else None

    # per-point radiation incl. shading: a0 = -(rad*(1-alb)) via sun=1
    rad = None if incident_rad_Wh_m2 is None else np.asarray(incident_rad_Wh_m2, dtype=float).reshape(N, -1)
    H = sun.shape[1]

    tint_fixed = None
    if tint_override_C is not None:
        tint_fixed = np.broadcast_to(np.asarray(tint_override_C, dtype=float) + KELVIN, (24,))
    # Tint/ET0 depend only on (material properties, day) / (albedo, sun row, day):
    # memoize across the group loop instead of recomputing per group
    tint_cache = {}
    et0_cache = {}

    def tint_for(m, t):
        if tint_fixed is not None:
            return tint_fixed
        key = (m.Dh, t)
        if key not in tint_cache:
            tint_cache[key] = physics.tint_cableizer(tyear, deltaT, m.Dh, t)
        return tint_cache[key]

    def et0_for(m, tair, gh_n, ha, sh, t, wind):
        key = (m.alb, t, sh.tobytes(), gh_n.tobytes())
        if key not in et0_cache:
            et0_cache[key] = physics.evapo_penman_monteith(
                tair, gh_n, ha, latitude, m.alb, sh, t, altitude_m, l,
                wind10_mps=wind)
        return et0_cache[key]

    if H == 24:
        # one representative day. Compat = upstream exactly (hc=5); gh mode may
        # bring wind-derived hc and per-point radiation to the same day solve.
        sel = weather.day_slice(month, day)
        t = physics.day_of_year(month, day)
        tair = tair_K_year[sel]
        gh = weather.gh_Wh_m2[sel]
        tsky = tsky_year[sel]
        ha = weather.rh_pct[sel]
        if tair.size != 24:
            raise ValueError("weather has %d rows for day %d/%d" % (tair.size, day, month))
        hc_arr = _hc_array(hc, mode, weather, sel, N, 24)
        wind_day = None if wind_year is None else wind_year[sel]
        temp = np.empty((N, 24))
        non_conv = []
        for n, m in enumerate(mats):
            if m.fixed_temp_degC != 0:
                temp[n, :] = m.fixed_temp_degC + KELVIN
                continue
            sh = sun[n] if rad is None else np.ones(24)
            gh_n = gh if rad is None else rad[n]
            tint = tint_for(m, t)
            et0 = et0_for(m, tair, gh_n, ha, sh, t, wind_day)
            C = physics.coeff_C(m.em)
            temp[n], count, ok = physics.solve_periodic_day(
                sh, C, gh_n, m.alb, m.em, tsky, tair, m.lambd, m.ep, m.Cv,
                et0, tint, m.kc, hc_arr[n],
                threshold=threshold_degC, max_iter=max_iter)
            if not ok:
                non_conv.append((n, count))
        # upstream reports degC rounded to 2; we keep full precision in K
    elif mode == "gh":
        if H != weather.tair_C.size:
            raise ValueError("gh mode expects sun_fraction with 24 or %d hours, got %d" % (weather.tair_C.size, H))
        hc_arr = _hc_array(hc, mode, weather, None, N, H)
        # per-day boundaries in the EPW series
        day_starts = np.flatnonzero(weather.hour == 1)
        temp = np.empty((N, H))
        non_conv = []
        for n, m in enumerate(mats):
            if m.fixed_temp_degC != 0:
                temp[n, :] = m.fixed_temp_degC + KELVIN
                continue
            sh = sun[n] if rad is None else np.ones(H)
            gh_n = weather.gh_Wh_m2 if rad is None else rad[n]
            C = physics.coeff_C(m.em)
            # per-day derived series
            tint = np.empty(H)
            et0 = np.empty(H)
            for s in day_starts:
                e = s + 24
                t = physics.day_of_year(int(weather.month[s]), int(weather.day[s]))
                tint[s:e] = tint_for(m, t)
                et0[s:e] = et0_for(
                    m, tair_K_year[s:e], gh_n[s:e], weather.rh_pct[s:e], sh[s:e], t,
                    None if wind_year is None else wind_year[s:e])
            # spin-up: iterate day 1 to periodic steady state
            s0, e0 = day_starts[0], day_starts[0] + 24
            day1, count, ok = physics.solve_periodic_day(
                sh[s0:e0], C, gh_n[s0:e0], m.alb, m.em, tsky_year[s0:e0],
                tair_K_year[s0:e0], m.lambd, m.ep, m.Cv, et0[s0:e0], tint[s0:e0],
                m.kc, hc_arr[n, s0:e0],
                threshold=threshold_degC, max_iter=max_iter)
            if not ok:
                non_conv.append((n, count))
            temp[n] = physics.solve_forward(
                sh, C, gh_n, m.alb, m.em, tsky_year, tair_K_year,
                m.lambd, m.ep, m.Cv, et0, tint, m.kc, hc_arr[n], T_start=day1[-1])
    else:
        raise ValueError("compat mode requires 24-hour sun_fraction, got %d columns" % H)

    return SolveResult.from_temp_K(temp, non_conv, N)
