# Material_database.csv provenance

- Source: `Scripts/Example/Step_1/Material_database.csv` of [Art-Ev/ICEtool](https://github.com/Art-Ev/ICEtool)
- Upstream: https://github.com/Art-Ev/ICEtool, commit `46f31dbf06f279246b476ea63876bf73bc3859b9` (2026-08-31), release v2026.07 schema
- Vendored: 2026-09-01, byte-identical copy
- License: upstream ICEtool is GPL-3.0; this repo's license must remain compatible.
- Schema: `ID;name;french_name;color;category;rugosity_age;FixedTemp[degC];alb;alb_source;em;em_source;Cv;Cv_source;lambd;lambd_source;ep;kc;CR;CR_source;CBS_base;CBS_source;Artificial;database_version;comments`
  (semicolon-separated, cp1252, multi-line quoted citation cells; duplicate IDs resolved by highest `database_version` in `materials.load_materials`)
