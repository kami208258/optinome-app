from __future__ import annotations

from pathlib import Path
import re
import numpy as np
import pandas as pd

# keep your existing clean_column_name() and parse_exp_file() below


def parse_data_sheet(df_raw: pd.DataFrame, sheet_name: str, exp_id: str) -> pd.DataFrame:
    """
    Parse a single EXP data sheet into tidy format.

    Special-cases Liven/DASGIP naming:
    - 'InoculationTime [h]'  (-> inoculationtime_h)
    - 'Duration'
    - 'Timestamp'
    and produces a numeric 'time_hours' column.
    """
    # ---------- find header row ----------
    header_idx = None
    for i, row in df_raw.iterrows():
        if not row.isna().all():
            header_idx = i
            break
    if header_idx is None:
        return pd.DataFrame()

    df = df_raw.iloc[header_idx:].copy()

    # use first non-NaN row as column names, clean them
    df.columns = df.iloc[0].astype(str).apply(clean_column_name)
    df = df.iloc[1:].reset_index(drop=True)

    cols = list(df.columns)

    # ---------- build time_hours ----------
    time_hours = None

    # 1) strongest signal: inoculation time already in hours
    #   "InoculationTime [h]" -> "inoculationtime_h" after cleaning
    for cand in ["inoculationtime_h", "inoculation_time_h"]:
        if cand in cols:
            time_hours = pd.to_numeric(df[cand], errors="coerce")
            break

    # 2) Duration – may be "4 days", "12:00:00", or numeric
    if time_hours is None and "duration" in cols:
        dur_raw = df["duration"].astype(str).str.strip()
        if dur_raw.str.contains("day", case=False, na=False).any():
            td = pd.to_timedelta(dur_raw, errors="coerce")
            time_hours = td.dt.total_seconds() / 3600.0
        else:
            # if looks like H:MM:SS or H:MM
            if dur_raw.str.contains(":", na=False).any():
                td = pd.to_timedelta(dur_raw, errors="coerce")
                time_hours = td.dt.total_seconds() / 3600.0
            else:
                time_hours = pd.to_numeric(dur_raw, errors="coerce")

    # 3) Existing time column
    if time_hours is None:
        for cand in ["time_hours", "time_hr", "time_h", "time"]:
            if cand in cols:
                time_hours = pd.to_numeric(df[cand], errors="coerce")
                break

    # 4) Fallback: Timestamp
    if time_hours is None and "timestamp" in cols:
        ts = pd.to_datetime(df["timestamp"], errors="coerce")
        t0 = ts.dropna().min()
        time_hours = (ts - t0).dt.total_seconds() / 3600.0

    # If we still couldn't make time, give up on this sheet
    if time_hours is None:
        return pd.DataFrame()

    df["time_hours"] = time_hours
    df = df.dropna(subset=["time_hours"])

    # ---------- metadata ----------
    df["exp_id"] = exp_id
    # Do NOT set 'vessel' here – builder.py will normalise it.

    return df
def parse_exp_file(path: Path) -> pd.DataFrame:
    """
    Parse a single EXP workbook (e.g. 24-0048-EXP1.xlsx) into a tidy DataFrame.

    - Figures out exp_id from the filename (e.g. "24-0048-EXP1.xlsx" → "EXP1")
    - Loops over all sheets
    - Uses parse_data_sheet(...) to turn each sheet into a tidy chunk
    - Concatenates everything.

    NOTE: We do NOT set 'vessel' here – builder.py can normalise that later
    if needed.
    """
    name = path.name
    m = re.search(r"(EXP\\d+)", name)
    exp_id = m.group(1) if m else name

    try:
        xls = pd.ExcelFile(path)
    except Exception:
        # If we can't open the workbook, skip it
        return pd.DataFrame()

    dfs: list[pd.DataFrame] = []

    for sheet in xls.sheet_names:
        # Read raw sheet with no header (parse_data_sheet will find it)
        df_raw = xls.parse(sheet, header=None)
        sub = parse_data_sheet(df_raw, sheet, exp_id)
        if sub is not None and not sub.empty:
            dfs.append(sub)

    if not dfs:
        return pd.DataFrame()

    df_all = pd.concat(dfs, ignore_index=True)

    # Sort nicely if time_hours exists
    if "time_hours" in df_all.columns:
        df_all = df_all.sort_values("time_hours").reset_index(drop=True)

    return df_all

