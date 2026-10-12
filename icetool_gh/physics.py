"""Energy-balance physics, ported from ICEtool (QGIS) Step 4.

Reference: archive/sources/ICEtool/Scripts/Step4_ComputeTemperatureEPW.py @ 46f31db.
Each function mirrors the upstream nested function of the same purpose, with
every closed-over variable promoted to an explicit parameter. Math is kept
line-identical to upstream so results match; do not "clean up" formulas here.

All functions operate on one 24-hour day, like upstream (which solves a single
representative day). The solver composes them over 365 days for full-year runs.

Units: temperatures in Kelvin unless suffixed _C, radiation in Wh/m2 per hour.
"""

import numpy as np
from scipy import optimize

SIGMA = 5.67e-08      # Stefan-Boltzmann, W.m-2.K-4 (upstream :401)
SIGMA_H = 2.043e-10   # Stefan-Boltzmann per hour (upstream :324)
Z_DEPTH = 0.2         # ground temperature burial depth, m (upstream :300)
W_ANNUAL = 2 * np.pi / 365
HC_QGIS = 5.0         # upstream constant wind coefficient, W.m-2.K-1 (:400)

DAYS_PER_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def day_of_year(month, day):
    """Upstream's `t` (:274-279): days of full months elapsed + day."""
    return sum(DAYS_PER_MONTH[: month - 1]) + day


def tsky_fuentes(tair_C, clamp_freezing=False):
    """Hourly sky temperature [K], Fuentes (1987) correlation (upstream :251).

    The correlation's T^1.5 term is NaN for tair_C < 0 — upstream has the same
    behavior (it is only run on summer days). clamp_freezing=True evaluates the
    T^1.5 term with tair clamped at 0 degC so full-year runs survive winter;
    a documented deviation used by the solver's "gh" mode only.
    """
    tair_C = np.asarray(tair_C, dtype=float)
    t15 = np.clip(tair_C, 0, None) if clamp_freezing else tair_C
    return np.round(0.037536 * t15**1.5 + 0.32 * tair_C + 273.15, 2)


def wind_10m_to_2m(wind10_mps):
    """EPW 10 m wind speed reduced to 2 m via the FAO-56 log profile."""
    return np.asarray(wind10_mps, dtype=float) * (4.87 / np.log(67.8 * 10 - 5.42))


def tsky_from_ir(ir_horiz_Wh_m2):
    """Sky temperature [K] from EPW Horizontal Infrared Radiation Intensity:
    Tsky = (IR / sigma)^0.25. Measured/modelled sky emission — more precise
    than the Fuentes correlation ("gh"-mode refinement; upstream's own comment
    at :250 flags Fuentes as "to replace in future versions")."""
    ir = np.asarray(ir_horiz_Wh_m2, dtype=float)
    return (np.clip(ir, 1e-6, None) / SIGMA) ** 0.25


def hc_mcadams(wind10_mps):
    """Convective coefficient [W/m2K] from wind speed (McAdams: 5.7 + 3.8 v),
    with EPW 10 m wind reduced to 2 m. "gh"-mode refinement of upstream's
    constant hc = 5."""
    return 5.7 + 3.8 * wind_10m_to_2m(wind10_mps)


def annual_tyear_deltaT(tair_K_8760, hours_1_24):
    """Per-hour-of-day annual mean temperature and max deviation (upstream :281-295).

    tair_K_8760: full-year dry-bulb series in Kelvin (upstream averages AFTER
    the +273.15 conversion, :254). hours_1_24: EPW 'hour' column (1..24).
    Returns (Tyear (24,), deltaT (24,)) in Kelvin.
    """
    tair = np.asarray(tair_K_8760, dtype=float)
    hrs = np.asarray(hours_1_24)
    tyear = np.array([tair[hrs == i].mean() for i in range(1, 25)])
    tmax, tmin = tair.max(), tair.min()
    deltaT = np.maximum(tmax - tyear, tyear - tmin)
    return tyear, deltaT


def tint_cableizer(tyear_K, deltaT, Dh, t):
    """Ground temperature at 20 cm, per hour of day [K] (upstream Tint(), :303-309).

    Dh = (lambd / Cv) * 86400 (upstream :298); t = day of year.
    """
    Zo = np.sqrt(2 * Dh / W_ANNUAL)
    tyear_K = np.asarray(tyear_K, dtype=float)
    deltaT = np.asarray(deltaT, dtype=float)
    return tyear_K - deltaT * np.e ** (-Z_DEPTH / Zo) * np.cos(W_ANNUAL * t - Z_DEPTH / Zo)


def evapo_penman_monteith(tair_K, gh, ha, lat, alb, shadow, t, alt, l, wind10_mps=None):
    """Hourly reference evapotranspiration term for one day (upstream Evapo(), :331-382).

    tair_K, gh, ha, shadow: 24-length arrays (air temp K, global horizontal
    Wh/m2, relative humidity %, sun-exposure fraction 0-1).
    lat: latitude deg; t: day of year; alt: altitude m;
    l: longitude correction = -(longz - longm) with upstream's sign fix (:316-321).
    wind10_mps: optional 24-length EPW wind speed (10 m); None reproduces
    upstream's fixed 1 km/h assumption, an array is the "gh"-mode refinement
    (real hourly wind, reduced to 2 m via the same FAO-56 profile).
    Returns 24-length array (already scaled by 2260000/3600, clamped at 0).
    """
    if wind10_mps is None:
        Vvent_h = np.full(24, 0.27 * (4.87 / np.log(67.8 * 2 - 5.42)))  # wind = 1 km/h (:325)
    else:
        Vvent_h = wind_10m_to_2m(wind10_mps)
    Pa = 101.3 * ((293 - 0.0065 * alt) / 293) ** 5.26  # (:326)
    gamma = 0.000665 * Pa                              # (:327)
    dr = 1 + 0.033 * np.cos((2 * np.pi / 365) * t)     # (:328)
    d = 0.409 * np.sin((2 * np.pi / 365) * t - 1.39)   # (:329)

    ETO = []
    th = [i + 0.5 for i in range(24)]
    phi = (np.pi / 180) * lat
    b = 2 * np.pi * (t - 81) / 364
    Sc = 0.1645 * np.sin(2 * b) - 0.1255 * np.cos(b) - 0.025 * np.sin(b)
    for i in range(24):
        Vvent = Vvent_h[i]
        Tmean = tair_K[i] - 273.15
        Rs = gh[i] * shadow[i] * 0.0036
        Rns = (1 - alb) * Rs
        delta = 4098 * (0.6108 * np.e ** (17.27 * Tmean / (Tmean + 237.3))) / (Tmean + 237.3) ** 2
        DTjour = delta / (delta + gamma * (1 + 0.24 * Vvent))
        DTnuit = delta / (delta + gamma * (1 + 0.96 * Vvent))
        PTjour = gamma / (delta + gamma * (1 + 0.24 * Vvent))
        PTnuit = gamma / (delta + gamma * (1 + 0.96 * Vvent))
        TT = (37 / (Tmean + 273.15)) * Vvent

        es = 0.6108 * np.e ** ((17.27 * Tmean) / (Tmean + 237.3))
        ea = es * (ha[i] / 100)

        w = (np.pi / 12) * ((th[i] + 0.06667 * l + Sc) - 12)
        ws = np.arccos(-np.tan(phi) * np.tan(d))
        w1 = w - np.pi / 24
        w2 = w + np.pi / 24
        if -ws < w < ws:  # day
            Ra = (12 * 60 / np.pi) * 0.0820 * dr * (
                (w2 - w1) * np.sin(phi) * np.sin(d)
                + np.cos(phi) * np.cos(d) * (np.sin(w2) - np.sin(w1))
            )
            Rso = (0.75 + (2 * 10 ** -5) * alt) * Ra
            Rnl = SIGMA_H * (Tmean + 273.15) ** 4 * (0.34 - 0.14 * np.sqrt(ea)) * (1.35 * (Rs / Rso) - 0.35)
            Rn = Rns - Rnl
            G = 0.1 * Rn
            Rng = 0.408 * Rn - G
            ETrad = DTjour * Rng
            ETwind = PTjour * TT * (es - ea)
        else:  # night
            Rnl = SIGMA_H * (Tmean + 273.15) ** 4 * (0.34 - 0.14 * np.sqrt(ea)) * (1.35 * 0.8 - 0.35)
            Rn = Rns - Rnl
            G = 0.5 * Rn
            Rng = 0.408 * Rn - G
            ETrad = DTnuit * Rng
            ETwind = PTnuit * TT * (es - ea)

        ET0 = (ETwind + ETrad) * 2260000 / 3600
        ETO.append(max(ET0, 0.0))
    return np.array(ETO)


def longitude_correction(longz, longm):
    """Upstream :316-321: l = -(longz - longm), with sign fix for longz>170 west."""
    if longz > 170 and longm < 0:
        longz = longz * -1
    return (longz - longm) * -1


def coeff_B(lambd, ep, Cv, hc):
    """Upstream :404."""
    return hc + (lambd / ep) + (Cv * ep / 3600)


def coeff_C(em):
    """Upstream :405."""
    return em * SIGMA


def compute_A(shadow, gh, alb, em, tsky, tair, lambd, ep, Cv, et0, kc, tint, T0, hc):
    """Upstream compute_A (:407-415). Note (1-alb) on BOTH direct and diffuse."""
    a0 = -((0.8 * gh * shadow * (1 - alb)) + (0.2 * gh * (1 - alb)))
    a1 = (em * SIGMA) * -(tsky) ** 4
    a2 = -hc * tair
    a3 = -(lambd / ep) * tint
    a4 = -(Cv * ep / 3600) * T0
    a5 = et0 * kc
    return a0 + a1 + a2 + a3 + a4 + a5


def thermal_equation(x, A, B, C):
    """Upstream :387-388: A + B*T + C*T^4."""
    return A + (B * x) + (C * x ** 4)


def solve_periodic_day(shadow, C, gh, alb, em, tsky, tair, lambd, ep, Cv, et0,
                       tint, kc, hc, threshold=0.5, max_iter=25, T0=28 + 273.15):
    """Upstream compute_temp solve branch (:430-453): iterate a 24h cycle to
    periodic steady state. All series are 24-length. hc may be a scalar
    (upstream: 5) or a 24-length array (wind-derived; B recomputed per hour).
    Returns (temp_K (24,), iterations, converged).
    """
    H = len(tair)
    hc_arr = np.broadcast_to(np.asarray(hc, dtype=float), (H,))
    temp = [0.0] * H
    count = 0
    equilibrium = False
    while not equilibrium:
        for h in range(H):
            prev = T0 if h == 0 else temp[h - 1]
            hch = hc_arr[h]
            B = coeff_B(lambd, ep, Cv, hch)
            A = compute_A(shadow[h], gh[h], alb, em, tsky[h], tair[h],
                          lambd, ep, Cv, et0[h], kc, tint[h], prev, hch)
            if h == 0:
                guess = T0 - 0.5
            elif shadow[h] > 0.4:
                guess = temp[h - 1] + 1.0
            else:
                guess = temp[h - 1] - 0.5
            temp[h] = optimize.root(thermal_equation, guess, (A, B, C)).x[0]
        count += 1
        if count >= 2:
            if abs(temp[H - 1] - T0) < threshold:
                equilibrium = True
            elif count >= max_iter:
                return np.array(temp), count, False
        T0 = temp[H - 1]
    return np.array(temp), count, True


def solve_forward(shadow, C, gh, alb, em, tsky, tair, lambd, ep, Cv, et0,
                  tint, kc, hc, T_start):
    """Single forward pass over H hours (full-year mode), seeded with T_start [K].

    hc may be a scalar or an H-length array (wind-derived). Returns temp_K (H,).
    """
    H = len(tair)
    hc_arr = np.broadcast_to(np.asarray(hc, dtype=float), (H,))
    temp = np.empty(H)
    prev = T_start
    for h in range(H):
        hch = hc_arr[h]
        Bh = coeff_B(lambd, ep, Cv, hch)
        A = compute_A(shadow[h], gh[h], alb, em, tsky[h], tair[h],
                      lambd, ep, Cv, et0[h], kc, tint[h], prev, hch)
        guess = prev + 1.0 if shadow[h] > 0.4 else prev - 0.5
        temp[h] = optimize.root(thermal_equation, guess, (A, Bh, C)).x[0]
        prev = temp[h]
    return temp
