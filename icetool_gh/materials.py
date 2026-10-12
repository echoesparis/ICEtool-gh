"""v2026.07 material database parser (upstream Scripts/Example/Step_1/Material_database.csv)."""

import os
from dataclasses import dataclass

import pandas as pd

DEFAULT_DB = os.path.join(os.path.dirname(__file__), "data", "Material_database.csv")

NUMERIC_COLS = ["FixedTemp[degC]", "alb", "em", "Cv", "lambd", "ep", "kc"]


@dataclass
class Material:
    id: str
    alb: float
    em: float
    Cv: float
    lambd: float
    ep: float
    kc: float
    fixed_temp_degC: float = 0.0  # 0 = solve normally (upstream :426 `if not FixedTemp==0`)
    name: str = ""
    category: str = ""

    @property
    def Dh(self):
        """Thermal diffusivity of soil, upstream :298."""
        return (self.lambd / self.Cv) * 86400


def load_materials(path=None):
    """Parse the semicolon-separated v2026.07 database by column name.

    Where an ID appears with several database_version entries, the highest
    version wins. Returns dict[str, Material].
    """
    df = pd.read_csv(path or DEFAULT_DB, sep=";", encoding="cp1252")
    for col in NUMERIC_COLS:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if "database_version" in df.columns:
        df["database_version"] = pd.to_numeric(df["database_version"], errors="coerce").fillna(1)
        df = df.sort_values("database_version").drop_duplicates("ID", keep="last")
    out = {}
    for _, row in df.iterrows():
        out[str(row["ID"])] = Material(
            id=str(row["ID"]),
            alb=row["alb"], em=row["em"], Cv=row["Cv"], lambd=row["lambd"],
            ep=row["ep"], kc=row["kc"],
            fixed_temp_degC=0.0 if pd.isna(row["FixedTemp[degC]"]) else float(row["FixedTemp[degC]"]),
            name=str(row.get("name", "")),
            category=str(row.get("category", "")),
        )
    return out
