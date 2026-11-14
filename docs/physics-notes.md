---
type: report
description: "Literature notes on the convective heat transfer coefficient hc: linear wind law (Palyvos 2008), Nusselt/Grashof correlation (Morille 2012), and what the solver uses."
applies_to:
  - icetool_gh/physics.py
---

# Convective heat transfer coefficient — notes

What the solver uses is in [[qgis-vs-gh|qgis-vs-gh.md]] (fixed 5 W/m²K in `qgis_compat`, hourly McAdams in `gh`) and `icetool_gh/physics.py` → `hc_mcadams`. This page keeps the literature behind that choice.

## Convective flux

`H = hc · (T_air − T_surface)`, with `hc` in W·m⁻²·K⁻¹. Weak wind: `hc` ≈ 5.

## Linear law in wind speed

`hc = a + b · V_air` (V in m/s). Palyvos (2008) lists about forty `(a, b)` pairs depending on surface texture, velocity, windward or leeward face; several are proposed for a horizontal surface with `V_air` < 5 m/s.

- J.A. Palyvos (2008), A survey of wind convection coefficient correlations for building envelope energy systems' modeling, Applied Thermal Engineering 28(8-9), 801–808, doi:10.1016/j.applthermaleng.2007.12.005.

## Dimensionless correlation

Morille (2012), A new urban soil model for SOLENE-microclimat: review, sensitivity analysis and validation on a car park (https://hal.archives-ouvertes.fr/hal-01629430v2/file/post_print.pdf), gives for flat surfaces coefficients depending on the flow regime:

`Nu = e · (a·Re^b + c·Gr^d)^f`

Natural (free) convection, air driven by buoyancy, laminar regime: `Nu = 0.49 · Gr^(1/4)`, with `Gr` the Grashof number.

`Nu = hc · Lc / λ_fluid`, with `Lc` the characteristic length (m) and `λ_fluid` the thermal conductivity of the fluid (W·m⁻¹·K⁻¹). With `Lc` = 1 m: `hc = Nu · λ_fluid`.

The coefficient table of the paper (regime by regime) is in the paper.
