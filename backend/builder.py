import os
from pathlib import Path

import pandas as pd

from .parsing import parse_exp_file

EXCLUDE_DIRS = {"venv", "outputs", "optinome_outputs", "model_outputs", "__pycache__"}

# You can put specific priority filenames here if you want them scanned first
PRIORITY: set[str] = set()
# Or limit number of files for debugging, e.g. PROCESS_LIMIT = 5
PROCESS_LIMIT: int | None = None


def iter_excels(root: Path, priority_names: set[str], limit: int | None) -> list[Path]:
    """
    Find Excel files under root, skipping venv / output dirs.
    """
    files: list[Path] = []
    for dirpath, dirnames, fs in os.walk(root):
        # Skip noisy or non-data dirs
        dirnames[:] = [d for d in dirnames if d not in EXCLUDE_DIRS]
        for f in fs:
            if f.lower().endswith((".xls", ".xlsx")):
                files.append(Path(dirpath) / f)

    # Priority ordering by filename if needed
    files = sorted(
        files,
        key=lambda p: (0 if p.name in priority_names else 1, p.name.lower()),
    )
    if limit is not None:
        files = files[:limit]
    return files


def _load_compilation_if_exists(project_dir: Path) -> pd.DataFrame:
    """
    Try to load a pre-compiled master table if it exists:
    - online_data_compilation.xlsx
    - data_compilation.xlsx

    This is the preferred source for the dashboard.
    """
    candidates = [
        project_dir / "online_data_compilation.xlsx",
        project_dir / "data_compilation.xlsx",
    ]

    for p in candidates:
        if p.exists():
            df = pd.read_excel(p)
            # Normalize column names (strip whitespace)
            df.columns = [str(c).strip() for c in df.columns]
            return df

    return pd.DataFrame()


def _build_from_raw_excels(project_dir: Path) -> pd.DataFrame:
    """
    Fallback: try to build from raw EXP Excel files using parse_exp_file.
    This is best-effort and may return empty if the formats are too custom.
    """
    files = iter_excels(project_dir, PRIORITY, PROCESS_LIMIT)
    out: list[pd.DataFrame] = []

    for p in files:
        df = parse_exp_file(p)
        if df is None or df.empty:
            continue
        out.append(df)

    if not out:
        return pd.DataFrame()

    master = pd.concat(out, ignore_index=True)
        # ---------- normalise vessel IDs ----------
    # We want logical V1..V4 per exp, not raw volumes like 620.0, 650.292...
    # Use Vreactor as the underlying identifier if present.
    if "vessel" not in master.columns or master["vessel"].isna().all():
        master["vessel"] = np.nan

    if "Vreactor" in master.columns:
        # Work exp-by-exp
        for exp_id, idx in master.groupby("exp_id").groups.items():
            sub = master.loc[idx]

            # unique reactor volumes for that EXP
            vals = (
                sub["Vreactor"]
                .dropna()
                .unique()
            )
            vals = np.sort(vals)

            # map first 4 distinct reactors to V1..V4
            for i, v in enumerate(vals[:4]):
                mask = (master.index.isin(idx)) & (master["Vreactor"] == v)
                master.loc[mask, "vessel"] = f"V{i+1}"

    # If still missing (no Vreactor etc.), fall back to a single V1 per exp
    if master["vessel"].isna().all():
        for exp_id, idx in master.groupby("exp_id").groups.items():
            master.loc[idx, "vessel"] = "V1"


    sort_cols = [c for c in ["exp_id", "vessel", "time_hours"] if c in master.columns]
    if sort_cols:
        master = master.sort_values(sort_cols).reset_index(drop=True)

    return master


def build_master(project_dir: Path) -> pd.DataFrame:
    """
    Preferred pipeline for the Optinome dashboard:

    1. Try to load the compiled dataset (online_data_compilation.xlsx or data_compilation.xlsx).
    2. If not found, fall back to attempting a raw parse of EXP Excel files.
    """
    # 1) Prefer the compiled dataset
    compiled = _load_compilation_if_exists(project_dir)
    if not compiled.empty:
        return compiled

    # 2) Fallback to raw parsing
    return _build_from_raw_excels(project_dir)
