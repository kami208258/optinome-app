# optinome_app.py
# Clean version – DASGIP data1–data4 → V1–V4 + HPLC analysis

from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st
import plotly.graph_objects as go
import textwrap

# ---------------------------------------------------------------------
# Paths / constants
# ---------------------------------------------------------------------

PROJECT_DIR = Path(__file__).resolve().parent
ASSETS_DIR = PROJECT_DIR / "assets"

# ---------------------------------------------------------------------
# HPLC helpers
# ---------------------------------------------------------------------


@st.cache_data(show_spinner="Loading HPLC data...")
def load_hplc_for_exp(exp_id: str) -> pd.DataFrame:
    """
    Load HPLC-like data for an experiment.

    First tries to find a separate EXP workbook with an HPLC sheet.
    If not found, derives data from Vglycerol/Vmeoh columns in compilation files.

    Columns returned:
      - unit        (A/B/C/D)
      - sample      (CM, T0..T10)
      - eft_hours   (h)
      - glycerol_gL
      - methanol_gL
      - exp_id
    """
    # Pull the numeric part, e.g. "EXP-1" or "EXP1" -> "1"
    m = re.search(r"(\d+)", str(exp_id))
    if not m:
        return pd.DataFrame()
    exp_num = m.group(1)

    # First, try the original approach: look for EXP workbook with HPLC sheet
    candidates: list[Path] = []
    for p in PROJECT_DIR.rglob(f"*EXP{exp_num}*.xlsx"):
        name = p.name.lower()
        if any(skip in name for skip in ["raw", "sds", "summary", "anthrone"]):
            continue
        candidates.append(p)

    if candidates:
        wb_path = sorted(candidates)[0]
        try:
            xls = pd.ExcelFile(wb_path)
            if "HPLC" in xls.sheet_names:
                df_raw = xls.parse("HPLC", header=None)
                df = tidy_hplc_sheet(df_raw)
                df["exp_id"] = exp_id
                return df
        except Exception:
            pass

    # Fallback: derive HPLC data from compilation files
    return derive_hplc_from_compilation(exp_id, exp_num)


def tidy_hplc_sheet(df_raw: pd.DataFrame) -> pd.DataFrame:
    """
    Parse the 'HPLC' sheet into tidy rows.

    Assumes structure like:
      row: ..., 'Unit', 'Sample', 'EFT', 'DF',
           'area', 'conc(g/L)', 'area', 'conc(g/L)', 'comment'
      then blocks for units A, B, C, D with CM, T0..T10 samples.
    """
    header_idx = None
    for i in range(len(df_raw)):
        val = str(df_raw.iloc[i, 1]).strip()
        if val == "Sample":
            header_idx = i
            break

    if header_idx is None:
        return pd.DataFrame()

    # Keep the first 9 columns (unit .. comment)
    df = df_raw.iloc[header_idx + 1 :, :9].copy()
    df.columns = [
        "unit",        # A/B/C/D (first row, then NaN below)
        "sample",      # CM, T0, T1...
        "eft",         # EFT (h)
        "df",          # dilution factor
        "area_gly",
        "glycerol_gL",
        "area_meoh",
        "methanol_gL",
        "notes",
    ]

    # Forward-fill A/B/C/D down each block
    df["unit"] = df["unit"].ffill()
    
    # Drop completely empty rows
    df = df[~df["sample"].isna()].copy()
    
    # Numeric conversions
    df["eft_hours"] = pd.to_numeric(df["eft"], errors="coerce")
    df["glycerol_gL"] = pd.to_numeric(df["glycerol_gL"], errors="coerce")
    df["methanol_gL"] = pd.to_numeric(df["methanol_gL"], errors="coerce")

    # Keep only rows with a time value (CM, T0..T10)
    df = df.dropna(subset=["eft_hours"])

    return df


def derive_hplc_from_compilation(exp_id: str, exp_num: str) -> pd.DataFrame:
    """
    Derive HPLC-like data from the compilation files.
    
    Looks for:
    1. Sheets matching EXP{exp_num}_X
    2. Sheets with 'exp_id' column matching EXP{exp_num}_X
    
    Extracts Vglycerol/Vmeoh columns as concentration proxies.
    Returns unit as A, B, C, D (unmapped).
    """
    compilation_files = [
        PROJECT_DIR / "online_data_compilation.xlsx",
        PROJECT_DIR / "data_compilation.xlsx",
    ]
    
    runs: list[pd.DataFrame] = []
    
    # No mapping needed if we want A, B, C, D
    
    for cpath in compilation_files:
        if not cpath.exists():
            continue
            
        try:
            xf = pd.ExcelFile(cpath)
        except Exception:
            continue
        
        for sh in xf.sheet_names:
            # Check strategy: sheet name vs content
            target_unit = None
            df = None
            
            # 1. Check sheet name
            m = re.match(rf"EXP{exp_num}_([A-D])", sh, re.IGNORECASE)
            if m:
                target_unit = m.group(1).upper()
                try:
                    df = xf.parse(sh, header=0)
                except Exception:
                    continue
            else:
                # 2. Check content (e.g. data_compilation.xlsx Sheet1)
                try:
                    # quick check
                    df_check = xf.parse(sh, header=0, nrows=5)
                    if "exp_id" in df_check.columns:
                        # Read full sheet
                        df_full = xf.parse(sh, header=0)
                        # Filter for this experiment number
                        # We look for exp_id containing "EXP{exp_num}"
                        mask = df_full["exp_id"].astype(str).str.contains(rf"EXP{exp_num}", case=False, na=False)
                        if not mask.any():
                            continue
                        
                        # We might have multiple vessels (EXP1_B, EXP1_C) in this sheet
                        # We need to process each vessel separately
                        sub_groups = df_full[mask].groupby("exp_id")
                        for eid, group in sub_groups:
                            m_id = re.match(r"(EXP\d+)_([A-D])", str(eid), re.IGNORECASE)
                            if m_id:
                                unit = m_id.group(2).upper()
                                # Process this group as if it were a sheet
                                processed = _process_hplc_group(group, exp_id, unit)
                                if not processed.empty:
                                    runs.append(processed)
                        continue # Done with this sheet
                except Exception:
                    pass
            
            # If we got a dataframe from Strategy 1 (sheet name match)
            if df is not None and target_unit:
                processed = _process_hplc_group(df, exp_id, target_unit)
                if not processed.empty:
                    runs.append(processed)

    if not runs:
        return pd.DataFrame()
    
    result = pd.concat(runs, ignore_index=True)
    return result[["unit", "sample", "eft_hours", "glycerol_gL", "methanol_gL", "exp_id"]]


def _process_hplc_group(df: pd.DataFrame, exp_id: str, unit: str) -> pd.DataFrame:
    """Helper to extract sampling points from a dataframe."""
    # Check for required columns
    # Note: Column names might vary (Vglycerol vs V_glycerol etc)
    # columns in data_compilation: 'cGly [g/L]', 'cMeOH [g/L]' ? 
    # Let's check the columns I saw earlier:
    # 'cGly [g/L]', 'cMeOH [g/L]' exist in data_compilation.xlsx!
    # 'Vglycerol', 'Vmeoh' exist in online_data_compilation.xlsx!
    
    cols = list(df.columns)
    
    gly_col = None
    meoh_col = None
    
    # Try to find best match
    for c in cols:
        cl = c.lower()
        if "vglycerol" in cl: gly_col = c
        elif "cgly" in cl: gly_col = c
        
        if "vmeoh" in cl: meoh_col = c
        elif "cmeoh" in cl: meoh_col = c
    
    # If we strictly need both? Or partial?
    # Let's be lenient
    if not gly_col and not meoh_col:
        return pd.DataFrame()

    # Get time column
    time_col = None
    if "EFT" in cols:
        time_col = "EFT"
    elif "Duration" in cols:
        time_col = "Duration"
    else:
        for c in cols:
            if "InoculationTime" in str(c):
                time_col = c
                break
    
    if not time_col:
        # Fallback to Time/Timestamp
        if "Time" in cols: time_col = "Time"
        elif "Timestamp" in cols: time_col = "Timestamp"
        else: return pd.DataFrame()
    
    # Convert time to hours
    time_hours = get_time_hours_from_series(df[time_col])
    
    # Extract values
    glycerol_vals = pd.to_numeric(df[gly_col], errors="coerce") if gly_col else pd.Series(0, index=df.index)
    methanol_vals = pd.to_numeric(df[meoh_col], errors="coerce") if meoh_col else pd.Series(0, index=df.index)
    
    # Create temp df
    temp_df = pd.DataFrame({
        "eft_hours": time_hours,
        "glycerol_gL": glycerol_vals,
        "methanol_gL": methanol_vals,
    }).dropna(subset=["eft_hours"])
    
    if temp_df.empty:
        return pd.DataFrame()
    
    # If there are very few points (manual sampling data?), just return them
    if len(temp_df) < 20:
        sampled_rows = []
        # Sort by time
        temp_df = temp_df.sort_values("eft_hours")
        for i in range(len(temp_df)):
            row = temp_df.iloc[i].copy()
            # Try to infer sample name from 'Sample' column if it exists
            sample_name = f"S{i}"
            if "Sample" in df.columns:
                 # Re-align index to fetch sample name
                 # This is tricky because temp_df re-indexed?
                 # Actually temp_df has original index.
                 orig_idx = temp_df.index[i]
                 val = str(df.loc[orig_idx, "Sample"])
                 if val and val.lower() != "nan":
                     sample_name = val
            
            row["sample"] = sample_name
            sampled_rows.append(row)
        
        sampled_df = pd.DataFrame(sampled_rows)
        sampled_df["unit"] = unit
        sampled_df["exp_id"] = exp_id
        return sampled_df

    # Otherwise (online data with many points), sample at intervals
    temp_df = temp_df.sort_values("eft_hours")
    max_time = temp_df["eft_hours"].max()
    
    sample_times = [0]
    t = 10
    while t < max_time:
        sample_times.append(t)
        t += 10
    sample_times.append(max_time)
    
    sampled_rows = []
    for i, target_t in enumerate(sample_times):
        idx = (temp_df["eft_hours"] - target_t).abs().idxmin()
        row = temp_df.loc[idx].copy()
        if i == 0:
            row["sample"] = "CM"
        else:
            row["sample"] = f"T{i-1}"
        sampled_rows.append(row)
    
    sampled_df = pd.DataFrame(sampled_rows)
    sampled_df["unit"] = unit
    sampled_df["exp_id"] = exp_id
    return sampled_df


# ---------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------


def iter_excels(project_dir: Path) -> list[Path]:
    """Return all Excel files under project_dir (excluding temp ~ files)."""
    files: list[Path] = []
    for p in project_dir.rglob("*.xls*"):
        if p.name.startswith("~$"):
            continue
        files.append(p)
    files.sort()
    return files


# ---------------------------------------------------------------------
# Helpers for time & columns
# ---------------------------------------------------------------------


def get_time_hours_from_series(s: pd.Series) -> pd.Series:
    """
    Convert a DASGIP time-like column to hours.

    Handles:
    - strings like '0 days 03:00:00'
    - '03:00:00' or '3:00'
    - numeric (already hours)
    - timestamps (absolute time → hours since first sample)
    """
    s = s.copy()
    s = s.replace("", np.nan)

    # If numeric → assume hours
    if pd.api.types.is_numeric_dtype(s):
        return s.astype(float)

    s_str = s.astype(str).str.strip()

    # 1) 'X days HH:MM:SS' → Timedelta
    if s_str.str.contains("day", case=False, na=False).any():
        td = pd.to_timedelta(s_str, errors="coerce")
        return td.dt.total_seconds() / 3600.0

    # 2) Try direct numeric
    num = pd.to_numeric(s_str, errors="coerce")
    if num.notna().sum() >= 2:
        return num.astype(float)

    # 3) Try as time / timestamp using to_timedelta
    td = pd.to_timedelta(s_str, errors="coerce")
    if td.notna().sum() >= 2:
        td0 = td.dropna().iloc[0]
        return (td - td0).dt.total_seconds() / 3600.0

    # 4) Last resort: treat as datetime
    dt = pd.to_datetime(s_str, errors="coerce")
    if dt.notna().sum() >= 2:
        dt0 = dt.dropna().iloc[0]
        return (dt - dt0) / np.timedelta64(1, "h")

    # 5) Fallback: monotonically increasing index
    idx = np.arange(len(s), dtype=float)
    out = pd.Series(idx, index=s.index, dtype=float)
    return out


# Regex for control parameters in DASGIP headers
RE_DO_PV = re.compile(r"DO\d*\.PV", re.IGNORECASE)
RE_FAIR_PV = re.compile(r"FAir\d*\.PV", re.IGNORECASE)
RE_FA_PV = re.compile(r"FA\d*\.PV", re.IGNORECASE)
RE_FB_PV = re.compile(r"FB\d*\.PV", re.IGNORECASE)
RE_FC_PV = re.compile(r"FC\d*\.PV", re.IGNORECASE)
RE_RPM = re.compile(r"N1\.PV", re.IGNORECASE)
RE_PH = re.compile(r"pH1\.PV", re.IGNORECASE)
RE_TEMP = re.compile(r"T1\.PV", re.IGNORECASE)


def _detect_col(cols, pattern: re.Pattern[str]) -> str | None:
    for c in cols:
        if pattern.search(str(c)):
            return c
    return None


# ---------------------------------------------------------------------
# Parsing DASGIP data sheets
# ---------------------------------------------------------------------


def parse_data_sheet(df_raw: pd.DataFrame, sheet_name: str, exp_id: str, headers_in_row0: bool = False) -> pd.DataFrame:
    """
    Parse one DASGIP 'dataX' sheet → tidy dataframe with core signals.

    - Picks a proper time column (Duration preferred).
    - Uses get_time_hours_from_series to build `time_hours`.
    - Extracts DO, air, FA1, FB1, FC1, pH, Temp, RPM.
    - Adds `exp_id` and `vessel`.
    
    Args:
        headers_in_row0: If True, headers are in row 0 (compilation files).
                         If False, headers are in row 1 (original EXP files).
    """
    if df_raw is None or df_raw.empty:
        return pd.DataFrame()

    df = df_raw.copy()
    
    if headers_in_row0:
        # Compilation files: headers in row 0, data starts row 1
        if len(df) < 2:
            return pd.DataFrame()
        df.columns = df.iloc[0].astype(str)
        df = df.iloc[1:].reset_index(drop=True)
    else:
        # Original DASGIP export: row 0 metadata, row 1 headers, row 2+ data
        if len(df) < 3:
            return pd.DataFrame()
        df.columns = df.iloc[1].astype(str)
        df = df.iloc[2:].reset_index(drop=True)

    cols = list(df.columns)

    # ----- choose time column -----
    t_col = None
    if "Duration" in cols:
        t_col = "Duration"
    else:
        for c in cols:
            if "InoculationTime" in str(c):
                t_col = c
                break
        if t_col is None:
            if "Timestamp" in cols:
                t_col = "Timestamp"
            else:
                t_col = cols[0]  # fallback

    time_hours = get_time_hours_from_series(df[t_col])

    # ----- detect control signals -----
    do_col = _detect_col(cols, RE_DO_PV)
    air_col = _detect_col(cols, RE_FAIR_PV)
    fa_col = _detect_col(cols, RE_FA_PV)
    fb_col = _detect_col(cols, RE_FB_PV)
    fc_col = _detect_col(cols, RE_FC_PV)
    rpm_col = _detect_col(cols, RE_RPM)
    ph_col = _detect_col(cols, RE_PH)
    temp_col = _detect_col(cols, RE_TEMP)

    out = pd.DataFrame({"time_hours": time_hours})
    out["do_data"] = pd.to_numeric(df.get(do_col), errors="coerce") if do_col else np.nan
    out["air_slph"] = pd.to_numeric(df.get(air_col), errors="coerce") if air_col else np.nan
    out["feed_fa_ml_h"] = pd.to_numeric(df.get(fa_col), errors="coerce") if fa_col else np.nan
    out["feed_fb_ml_h"] = pd.to_numeric(df.get(fb_col), errors="coerce") if fb_col else np.nan
    out["feed_fc_ml_h"] = pd.to_numeric(df.get(fc_col), errors="coerce") if fc_col else np.nan
    out["rpm"] = pd.to_numeric(df.get(rpm_col), errors="coerce") if rpm_col else np.nan
    out["ph"] = pd.to_numeric(df.get(ph_col), errors="coerce") if ph_col else np.nan
    out["temp_c"] = pd.to_numeric(df.get(temp_col), errors="coerce") if temp_col else np.nan

    out = out.dropna(subset=["time_hours"])
    if out.empty:
        return out

    # ----- metadata: exp_id, vessel -----
    sheet_key = sheet_name.strip().lower()
    # Handle "data1", "data2" style
    m = re.search(r"data\s*(\d+)", sheet_key)
    if m:
        vessel = f"V{m.group(1)}"
    # Handle pure numeric sheet names "1", "2", etc.
    elif sheet_key.isdigit():
        vessel = f"V{sheet_key}"
    else:
        vessel = sheet_name  # fallback

    out["exp_id"] = exp_id
    out["vessel"] = vessel

    return out


def parse_exp_file(path: Path) -> pd.DataFrame:
    """
    Parse a single EXP Excel workbook.

    Only sheets whose name starts with 'data' are used (data1–data4 → V1–V4).
    """
    name = path.name
    if "EXP" not in name.upper():
        return pd.DataFrame()

    # EXP id, e.g. 24-0048-EXP1 → EXP1
    m = re.search(r"EXP(\d+)", name.upper())
    exp_id = f"EXP{m.group(1)}" if m else name

    try:
        xf = pd.ExcelFile(path)
    except Exception:
        return pd.DataFrame()

    sheets = [s for s in xf.sheet_names if s.lower().startswith("data")]
    runs: list[pd.DataFrame] = []

    for sh in sheets:
        df_raw = xf.parse(sh, header=None)
        sub = parse_data_sheet(df_raw, sh, exp_id, headers_in_row0=False)
        if sub is None or sub.empty:
            continue
        sub["source_file"] = name
        runs.append(sub)

    if not runs:
        return pd.DataFrame()

    return pd.concat(runs, ignore_index=True)


def parse_compilation_file(path: Path) -> pd.DataFrame:
    """
    Parse a compilation Excel workbook (data_compilation.xlsx or online_data_compilation.xlsx).

    Handle two formats:
    1. Sheets named "EXP1_B", "EXP3_A" etc.
    2. Sheets with an "exp_id" column (e.g. "EXP1_C") that contains multiple experiments.
    
    Maps A->1, B->2, C->3, D->4 for vessel IDs.
    """
    name = path.name

    try:
        xf = pd.ExcelFile(path)
    except Exception:
        return pd.DataFrame()

    runs: list[pd.DataFrame] = []
    
    # Helper to map letter to number
    def _map_vessel(letter: str) -> str:
        letter = letter.upper()
        mapping = {"A": "1", "B": "2", "C": "3", "D": "4"}
        return f"V{mapping.get(letter, letter)}"

    for sh in xf.sheet_names:
        # --- Strategy A: Sheet name contains info (EXP1_B) ---
        m = re.match(r"(EXP\d+)_([A-D])", sh, re.IGNORECASE)
        if m:
            exp_id = m.group(1).upper()
            vessel = _map_vessel(m.group(2))
            
            df_raw = xf.parse(sh, header=None)
            sub = parse_data_sheet_simple(df_raw, exp_id, vessel)
            if sub is not None and not sub.empty:
                sub["source_file"] = name
                runs.append(sub)
            continue
            
        # --- Strategy B: Check content for "exp_id" column ---
        # Try reading with header=0 to check for exp_id column
        try:
            df_check = xf.parse(sh, header=0, nrows=5)
            if "exp_id" in df_check.columns:
                # Read full sheet
                df = xf.parse(sh, header=0)
                
                # Split by exp_id
                groups = df.groupby("exp_id")
                for eid, group in groups:
                    eid_str = str(eid).strip() # e.g. "EXP1_C"
                    
                    # Parse EXP and Vessel from the ID
                    # EXP1_C -> exp_id="EXP1", vessel="VC" -> "V3"
                    m_id = re.match(r"(EXP\d+)_([A-D])", eid_str, re.IGNORECASE)
                    if m_id:
                        real_exp_id = m_id.group(1).upper()
                        real_vessel = _map_vessel(m_id.group(2))
                    else:
                        real_exp_id = eid_str
                        real_vessel = "V1" # fallback
                    
                    sub = _parse_dataframe_group(group, real_exp_id, real_vessel)
                    if sub is not None and not sub.empty:
                        sub["source_file"] = name
                        runs.append(sub)
                continue
                
        except Exception:
            pass

        # --- Strategy C: Numeric sheet names or "Sheet1" without exp_id col ---
        if sh.isdigit():
            exp_id = "Compilation"
            vessel = f"V{sh}"
        elif sh.lower() == "sheet1":
            exp_id = "Compilation"
            vessel = "V1"
        else:
            continue

        df_raw = xf.parse(sh, header=None)
        sub = parse_data_sheet_simple(df_raw, exp_id, vessel)
        if sub is not None and not sub.empty:
            sub["source_file"] = name
            runs.append(sub)

    if not runs:
        return pd.DataFrame()

    return pd.concat(runs, ignore_index=True)


def _parse_dataframe_group(df: pd.DataFrame, exp_id: str, vessel: str) -> pd.DataFrame:
    """
    Parse a pre-filtered dataframe (from data_compilation splitting).
    Expected to have columns like 'EFT', 'Time', 'pH', etc.
    """
    cols = list(df.columns)
    
    # ----- choose time column -----
    t_col = None
    if "EFT" in cols:
        t_col = "EFT" # compilation files often use EFT
    elif "Duration" in cols:
        t_col = "Duration"
    else:
         for c in cols:
            if "InoculationTime" in str(c):
                t_col = c
                break
    
    if t_col is None:
        if "Time" in cols:
             t_col = "Time"
        elif "Timestamp" in cols:
            t_col = "Timestamp"
        else:
            return pd.DataFrame() # No time column found

    time_hours = get_time_hours_from_series(df[t_col])
    
    # ----- detect control signals -----
    # Mapping for data_compilation.xlsx specific columns if needed
    
    out = pd.DataFrame({"time_hours": time_hours})
    
    # Helper to find column
    def get_col(candidates):
        for c in candidates:
            if c in cols:
                return df[c]
            # Try regex/fuzzy matching?
            for col_name in cols:
                if c.lower() == str(col_name).lower():
                    return df[col_name]
        return np.nan

    out["do_data"] = pd.to_numeric(get_col(["DO", "DO1.PV", "DO1.Out"]), errors="coerce")
    out["air_slph"] = pd.to_numeric(get_col(["Air", "FAir1.PV", "FAir"]), errors="coerce")
    out["feed_fa_ml_h"] = pd.to_numeric(get_col(["FA1.PV", "FA1", "Feed A"]), errors="coerce")
    out["feed_fb_ml_h"] = pd.to_numeric(get_col(["FB1.PV", "FB1", "Feed B"]), errors="coerce")
    out["feed_fc_ml_h"] = pd.to_numeric(get_col(["FC1.PV", "FC1", "Feed C"]), errors="coerce")
    out["rpm"] = pd.to_numeric(get_col(["N1.PV", "RPM", "Stirrer"]), errors="coerce")
    out["ph"] = pd.to_numeric(get_col(["pH1.PV", "pH"]), errors="coerce")
    out["temp_c"] = pd.to_numeric(get_col(["T1.PV", "T", "Temp"]), errors="coerce")

    out = out.dropna(subset=["time_hours"])
    out["exp_id"] = exp_id
    out["vessel"] = vessel
    
    return out


def parse_data_sheet_simple(df_raw: pd.DataFrame, exp_id: str, vessel: str) -> pd.DataFrame:
    """
    Parse a data sheet with headers in row 0.
    
    Returns a tidy dataframe with time_hours and control signals.
    """
    if df_raw is None or df_raw.empty:
        return pd.DataFrame()

    df = df_raw.copy()
    
    # Headers in row 0, data starts row 1
    if len(df) < 2:
        return pd.DataFrame()
    df.columns = df.iloc[0].astype(str)
    df = df.iloc[1:].reset_index(drop=True)

    cols = list(df.columns)

    # ----- choose time column -----
    t_col = None
    if "Duration" in cols:
        t_col = "Duration"
    else:
        for c in cols:
            if "InoculationTime" in str(c):
                t_col = c
                break
        if t_col is None:
            if "Timestamp" in cols:
                t_col = "Timestamp"
            else:
                t_col = cols[0]  # fallback

    time_hours = get_time_hours_from_series(df[t_col])

    # ----- detect control signals -----
    do_col = _detect_col(cols, RE_DO_PV)
    air_col = _detect_col(cols, RE_FAIR_PV)
    fa_col = _detect_col(cols, RE_FA_PV)
    fb_col = _detect_col(cols, RE_FB_PV)
    fc_col = _detect_col(cols, RE_FC_PV)
    rpm_col = _detect_col(cols, RE_RPM)
    ph_col = _detect_col(cols, RE_PH)
    temp_col = _detect_col(cols, RE_TEMP)

    out = pd.DataFrame({"time_hours": time_hours})
    out["do_data"] = pd.to_numeric(df.get(do_col), errors="coerce") if do_col else np.nan
    out["air_slph"] = pd.to_numeric(df.get(air_col), errors="coerce") if air_col else np.nan
    out["feed_fa_ml_h"] = pd.to_numeric(df.get(fa_col), errors="coerce") if fa_col else np.nan
    out["feed_fb_ml_h"] = pd.to_numeric(df.get(fb_col), errors="coerce") if fb_col else np.nan
    out["feed_fc_ml_h"] = pd.to_numeric(df.get(fc_col), errors="coerce") if fc_col else np.nan
    out["rpm"] = pd.to_numeric(df.get(rpm_col), errors="coerce") if rpm_col else np.nan
    out["ph"] = pd.to_numeric(df.get(ph_col), errors="coerce") if ph_col else np.nan
    out["temp_c"] = pd.to_numeric(df.get(temp_col), errors="coerce") if temp_col else np.nan

    out = out.dropna(subset=["time_hours"])
    if out.empty:
        return out

    out["exp_id"] = exp_id
    out["vessel"] = vessel

    return out


def build_master(project_dir: Path) -> pd.DataFrame:
    """Scan project dir for EXP workbooks or compilation files and build the master table."""
    files = iter_excels(project_dir)
    out: list[pd.DataFrame] = []

    # First, try to find compilation files
    compilation_files = [
        "data_compilation.xlsx",
        "online_data_compilation.xlsx",
    ]
    
    for cname in compilation_files:
        cpath = project_dir / cname
        if cpath.exists():
            df = parse_compilation_file(cpath)
            if not df.empty:
                out.append(df)

    # Then, look for EXP files (original behavior)
    for p in files:
        if "EXP" not in p.name.upper():
            continue
        df = parse_exp_file(p)
        if not df.empty:
            out.append(df)

    if not out:
        return pd.DataFrame()

    master = pd.concat(out, ignore_index=True)
    master = master.sort_values(["exp_id", "vessel", "time_hours"]).reset_index(drop=True)
    return master


# ---------------------------------------------------------------------
# Plotting helpers (matplotlib → Streamlit)
# ---------------------------------------------------------------------


def _style(ax, x_label: str, y_label: str):
    """Apply common styling to axes."""
    ax.set_xlabel(x_label, fontsize=9, fontweight="bold", color="#333")
    ax.set_ylabel(y_label, fontsize=9, fontweight="bold", color="#333")
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.3, color="#ccc")
    ax.tick_params(axis="both", which="major", labelsize=8)
    # Remove top and right spines for a cleaner look
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_ph_temp(df: pd.DataFrame, title_suffix: str):
    fig, ax1 = plt.subplots(figsize=(6, 4))  # Smaller default size for grid
    ax2 = ax1.twinx()

    for v, sub in df.groupby("vessel"):
        sub = sub.dropna(subset=["time_hours"])
        if sub["ph"].notna().any():
            ax1.plot(sub["time_hours"], sub["ph"], label=f"{v} pH", linewidth=1.5)
        if sub["temp_c"].notna().any():
            ax2.plot(
                sub["time_hours"],
                sub["temp_c"],
                linestyle="--",
                label=f"{v} Temp",
                linewidth=1.5,
                alpha=0.7,
            )

    _style(ax1, "Time (h)", "pH")
    ax2.set_ylabel("Temp (°C)", fontsize=9, fontweight="bold", color="#333")
    ax2.tick_params(axis="both", labelsize=8)
    ax2.spines["top"].set_visible(False)

    ax1.set_title(f"pH & Temp", fontsize=10, fontweight="bold", pad=10)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(
        lines1 + lines2,
        labels1 + labels2,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=4,
        fontsize=8,
        frameon=False,
    )
    plt.tight_layout()
    return fig


def plot_do_air_rpm(df: pd.DataFrame, title_suffix: str):
    fig, ax1 = plt.subplots(figsize=(6, 4))
    ax2 = ax1.twinx()

    for v, sub in df.groupby("vessel"):
        sub = sub.dropna(subset=["time_hours"])
        if sub["do_data"].notna().any():
            ax1.plot(sub["time_hours"], sub["do_data"], label=f"{v} DO", linewidth=1.5)
        if sub["air_slph"].notna().any():
            ax1.plot(
                sub["time_hours"],
                sub["air_slph"],
                linestyle="--",
                label=f"{v} Air",
                linewidth=1.5,
                alpha=0.7,
            )
        if sub["rpm"].notna().any():
            ax2.plot(
                sub["time_hours"],
                sub["rpm"],
                linestyle=":",
                label=f"{v} RPM",
                linewidth=1.5,
                alpha=0.6,
            )

    _style(ax1, "Time (h)", "DO (%) / Air (sL/h)")
    ax2.set_ylabel("RPM", fontsize=9, fontweight="bold", color="#333")
    ax2.tick_params(labelsize=8)
    ax2.spines["top"].set_visible(False)
    
    ax1.set_title(f"DO / Air / RPM", fontsize=10, fontweight="bold", pad=10)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(
        lines1 + lines2,
        labels1 + labels2,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=4,
        fontsize=8,
        frameon=False,
    )
    plt.tight_layout()
    return fig


def plot_feeds(df: pd.DataFrame, title_suffix: str):
    fig, ax = plt.subplots(figsize=(6, 4))
    for v, sub in df.groupby("vessel"):
        if sub["feed_fa_ml_h"].notna().any():
            ax.plot(
                sub["time_hours"], sub["feed_fa_ml_h"], label=f"{v} FA1", linewidth=1.5
            )

    _style(ax, "Time (h)", "FA1 feed (mL/h)")
    ax.set_title(f"Glycerol/MeOH feed (FA1)", fontsize=10, fontweight="bold", pad=10)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=4,
        fontsize=8,
        frameon=False,
    )
    plt.tight_layout()
    return fig


def plot_base(df: pd.DataFrame, title_suffix: str):
    fig, ax = plt.subplots(figsize=(6, 4))
    for v, sub in df.groupby("vessel"):
        if sub["feed_fb_ml_h"].notna().any():
            ax.plot(
                sub["time_hours"], sub["feed_fb_ml_h"], label=f"{v} FB1", linewidth=1.5
            )

    _style(ax, "Time (h)", "Base FB1 (mL/h)")
    ax.set_title(f"Base feed", fontsize=10, fontweight="bold", pad=10)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=4,
        fontsize=8,
        frameon=False,
    )
    plt.tight_layout()
    return fig


def plot_antifoam(df: pd.DataFrame, title_suffix: str):
    fig, ax = plt.subplots(figsize=(6, 4))
    for v, sub in df.groupby("vessel"):
        if sub["feed_fc_ml_h"].notna().any():
            ax.plot(
                sub["time_hours"], sub["feed_fc_ml_h"], label=f"{v} FC1", linewidth=1.5
            )

    _style(ax, "Time (h)", "Antifoam FC1 (mL/h)")
    ax.set_title(f"Antifoam", fontsize=10, fontweight="bold", pad=10)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=4,
        fontsize=8,
        frameon=False,
    )
    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------
# Styling & Theme
# ---------------------------------------------------------------------

def get_custom_css():
    return """
    <style>
        /* --- Root Variables --- */
        :root {
            --primary-color: #0E76A8;
            --background-color: #F8F9FA;
            --card-bg: #FFFFFF;
            --text-color: #2C3E50;
            --text-secondary: #6C757D;
            --border-radius: 12px;
            --shadow-sm: 0 2px 4px rgba(0,0,0,0.05);
            --shadow-md: 0 4px 12px rgba(0,0,0,0.08);
            --font-main: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
        }

        /* --- Global Reset / Defaults --- */
        body {
            font-family: var(--font-main);
            color: var(--text-color);
            background-color: var(--background-color);
        }
        
        /* Hide standard Streamlit chrome */
        #MainMenu {visibility: hidden;}
        footer {visibility: hidden;}
        header[data-testid="stHeader"] {
            display: none;
        }
        
        .stApp {
            background-color: var(--background-color);
        }

        /* --- Layout & Spacing --- */
        .block-container {
            padding-top: 1rem; /* Reduced since header is hidden */
            padding-bottom: 3rem;
            max-width: 1600px;
        }

        /* --- Sidebar Styling (Always visible, fixed) --- */
        section[data-testid="stSidebar"] {
            background-color: #FFFFFF;
            border-right: 1px solid #E9ECEF;
            width: 280px !important;
            min-width: 280px !important;
            transform: none !important;
            transition: none !important;
        }
        section[data-testid="stSidebar"] > div:first-child {
            position: sticky;
            top: 0;
            height: 100vh;
            overflow-y: auto;
            padding-top: 1rem;
            scrollbar-width: none; /* Firefox */
            -ms-overflow-style: none; /* IE/Edge */
        }
        section[data-testid="stSidebar"] > div:first-child::-webkit-scrollbar {
            display: none; /* Chrome/Safari */
        }
        /* Hide ALL sidebar collapse/toggle controls */
        button[data-testid="stSidebarCollapseButton"],
        button[kind="headerNoPadding"],
        [data-testid="collapsedControl"],
        [data-testid="stSidebarNavCollapseButton"],
        .stSidebar button[kind="header"],
        [data-testid="stSidebarCollapsedControl"] {
            display: none !important;
            visibility: hidden !important;
            pointer-events: none !important;
        }
        /* Prevent sidebar from collapsing */
        section[data-testid="stSidebar"][aria-expanded="false"] {
            width: 280px !important;
            min-width: 280px !important;
            transform: none !important;
            margin-left: 0 !important;
        }
        [data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2 {
            font-size: 1.5rem;
            color: var(--primary-color);
        }
        
        /* --- Filter Bar (Top) --- */
        .filter-bar {
            background: var(--card-bg);
            padding: 1rem 1.5rem;
            border-radius: var(--border-radius);
            box-shadow: var(--shadow-sm);
            margin-bottom: 2rem;
            border: 1px solid #E9ECEF;
            display: flex;
            align-items: center;
            gap: 20px;
        }

        /* --- Cards --- */
        .stCard {
            background-color: var(--card-bg);
            padding: 1.5rem;
            border-radius: var(--border-radius);
            box-shadow: var(--shadow-sm);
            border: 1px solid #F1F3F5;
            margin-bottom: 1.5rem;
            transition: transform 0.2s, box-shadow 0.2s;
        }
        .stCard:hover {
            box-shadow: var(--shadow-md);
        }
        
        /* Fix for overlaying plots/content */
        .main .block-container {
            z-index: 1;
        }

        /* --- Metrics / KPI --- */
        [data-testid="stMetric"] {
            background-color: #FFFFFF;
            padding: 10px 15px;
            border-radius: 8px;
            border: 1px solid #E9ECEF;
            box-shadow: 0 1px 3px rgba(0,0,0,0.04);
            text-align: center;
            min-height: 110px; /* Ensure equal height */
            display: flex;
            flex-direction: column;
            justify-content: center;
        }
        [data-testid="stMetricLabel"] {
            color: var(--text-secondary);
            font-size: 0.85rem;
            font-weight: 500;
        }
        [data-testid="stMetricValue"] {
            color: var(--primary-color);
            font-weight: 700;
            font-size: 1.5rem;
        }

        /* --- Tabs --- */
        .stTabs [data-baseweb="tab-list"] {
            gap: 8px;
            background-color: transparent;
            padding-bottom: 5px;
        }
        .stTabs [data-baseweb="tab"] {
            height: 40px;
            border-radius: 20px;
            padding: 0 20px;
            background-color: #E9ECEF;
            border: none;
            color: var(--text-secondary);
            font-weight: 600;
        }
        .stTabs [data-baseweb="tab"][aria-selected="true"] {
            background-color: var(--primary-color);
            color: white;
        }

        /* --- Headings --- */
        h1, h2, h3 {
            font-family: var(--font-main);
            font-weight: 700;
            color: var(--text-color);
        }
        h3 {
            font-size: 1.1rem;
            margin-top: 0;
            margin-bottom: 1rem;
        }

        /* --- Tables --- */
        [data-testid="stDataFrame"] {
            border: 1px solid #E9ECEF;
            border-radius: 8px;
            overflow: hidden;
        }

        /* --- Accordion (Expander) --- */
        .streamlit-expanderHeader {
            background-color: #fff;
            border-radius: 8px;
            font-weight: 600;
        }
    </style>
    """

# ---------------------------------------------------------------------
# Suggestions / DOE
# ---------------------------------------------------------------------


def render_suggestions(exp_id: str):
    st.markdown("### Process Optimization Suggestions")
    
    with st.expander("A. Core set-points (Induction Phase)", expanded=True):
        st.info("🎯 **Target**: Balance growth & protein stability.")
        c1, c2 = st.columns(2)
        with c1:
            st.markdown("""
            - **DO set-point**: `25–35%` (Avoid <15%)
            - **Temperature**: `25–26°C` (Promotes folding)
            """)
        with c2:
            st.markdown("""
            - **pH**: `5.2–5.8`
            - **Agitation**: Use Airflow as main lever; minimize RPM shear.
            """)

    with st.expander("B. Carbon & Induction Strategy"):
        st.markdown("""
        1. **Glycerol Batch**: Target OD600 ~50–100.
        2. **Derepression**: 2–4h low glycerol feed.
        3. **Adaptation**: Slowly ramp FA1. Residual `0.5–1.0 g/L`.
        4. **Induction**: Residual MeOH `0.5–2.0 g/L`.
        """)
        st.caption("Future: Consider Sorbitol co-feed (20–40% carbon).")

    with st.expander("C. Control Heuristics"):
        st.success("💡 **Tip**: Base feed (FB1) is a good proxy for metabolic activity.")
        st.markdown("""
        - Keep DO curve smooth.
        - Make small FA1 changes frequently.
        - Once collagen assay is online, use titer/qP/STY to tune profiles.
        """)


def generate_default_doe() -> pd.DataFrame:
    df = pd.DataFrame(
        [
            (25, 22, 0.5, 0.0, 5.5),
            (35, 28, 2.0, 0.5, 5.5),
            (25, 22, 2.0, 0.5, 5.5),
            (25, 28, 0.5, 0.5, 5.5),
        ],
        columns=[
            "DO (%)",
            "Temp (°C)",
            "MeOH Res (g/L)",
            "Sorbitol Frac",
            "pH",
        ],
        index=["Run 1", "Run 2", "Run 3", "Run 4"],
    )
    # Add strategy explanation column
    # df["Strategy"] = ["Base", "High Stress", "High MeOH", "Mixed"] 
    return df


def render_doe():
    st.markdown("### DOE – Next runs (Protein Production)")
    
    df = generate_default_doe()
    
    # Toolbar
    c1, c2 = st.columns([1, 4])
    with c1:
        csv = df.to_csv().encode()
        st.download_button(
            "⬇️ Download CSV", 
            csv, 
            "optinome_doe.csv", 
            "text/csv", 
            type="primary",
            use_container_width=True
        )
    
    st.dataframe(
        df.style.background_gradient(cmap="Blues", subset=["DO (%)", "MeOH Res (g/L)"]),
        use_container_width=True,
        height=200
    )


# ---------------------------------------------------------------------
# Main app
# ---------------------------------------------------------------------


@st.cache_data(show_spinner="Building master table from EXP files…")
def load_master_cached() -> pd.DataFrame:
    return build_master(PROJECT_DIR)


def main():
    st.set_page_config(
        page_title="Optinome Dashboard",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )

    # Inject CSS
    st.markdown(get_custom_css(), unsafe_allow_html=True)

    # --- Data Loading ---
    master = load_master_cached()
    if master.empty:
        st.error("❌ No EXP data found under this project directory.")
        st.stop()

    # --- Native Sidebar ---
    with st.sidebar:
        if (ASSETS_DIR / "fermenter.png").exists():
            st.image(str(ASSETS_DIR / "fermenter.png"), width=180)
        st.markdown("---")
        st.markdown("#### Navigation")
        st.info("Use the tabs on the right to switch views.")
        st.markdown("")
        st.markdown("#### Quick Links")
        st.markdown("📄 [Documentation](#)")
        st.markdown("🐛 [Report Issue](#)")
        st.markdown("---")
        if st.button("🔄 Refresh Data", use_container_width=True):
            st.cache_data.clear()
            st.rerun()

    # --- Header & Filter Bar ---
    c_title, c_filters = st.columns([2, 4], gap="large")

    with c_title:
        st.markdown("# Optinome")
        st.caption("Bioreactor Optimization Dashboard")

    with c_filters:
        exps = sorted(master["exp_id"].unique())
        f1, f2 = st.columns([1, 2])
        with f1:
            exp_id = st.selectbox("Experiment", exps, index=0, label_visibility="collapsed")

        df_exp = master[master["exp_id"] == exp_id].copy()
        vessels = sorted(df_exp["vessel"].unique())

        with f2:
            vessels_selected = st.multiselect(
                "Select Vessels",
                vessels,
                default=vessels,
                label_visibility="collapsed",
                placeholder="Select vessels..."
            )

    st.markdown("---")

    # Filtered Data
    if not vessels_selected:
        st.warning("Please select at least one vessel.")
        st.stop()
        
    df_view = df_exp[df_exp["vessel"].isin(vessels_selected)].copy()

    # --- Main Navigation Tabs ---
    tab_overview, tab_suggestions, tab_doe, tab_analysis = st.tabs(
        ["Overview & Trends", "Suggestions", "Doe Design", "HPLC Analysis"]
    )

    # -------------------- Overview tab --------------------
    with tab_overview:
        
        # KPI Cards
        kpi_cols = st.columns(4)
        with kpi_cols[0]:
            duration = df_view["time_hours"].max() if df_view["time_hours"].notna().any() else 0
            st.metric("Run Duration", f"{duration:.1f} h", delta=None)
        with kpi_cols[1]:
            st.metric("Vessels Active", f"{len(vessels_selected)}", delta=f"{len(vessels)} total")
        
        st.markdown("### Time Course Trends")
        
        # Chart Grid
        # Row 1
        r1c1, r1c2 = st.columns(2)
        with r1c1:
            fig = plot_ph_temp(df_view, exp_id)
            st.pyplot(fig)
            
        with r1c2:
            fig = plot_do_air_rpm(df_view, exp_id)
            st.pyplot(fig)

        # Row 2
        r2c1, r2c2 = st.columns(2)
        with r2c1:
            fig = plot_feeds(df_view, exp_id)
            st.pyplot(fig)
            
        with r2c2:
            fig = plot_base(df_view, exp_id)
            st.pyplot(fig)
            
        # Row 3 (Antifoam + Raw Data)
        r3c1, r3c2 = st.columns(2)
        with r3c1:
            fig = plot_antifoam(df_view, exp_id)
            st.pyplot(fig)
        
        with r3c2:
            pass # Placeholder or leave empty for balance

        # Raw Data Expander
        st.markdown("### Raw Data Export")
        with st.expander("View Filtered Data Table", expanded=False):
            st.dataframe(
                df_view.sort_values(["vessel", "time_hours"]),
                use_container_width=True,
                height=300,
            )
            csv = df_view.to_csv(index=False).encode("utf-8")
            st.download_button(
                "⬇️ Download Filtered CSV",
                csv,
                f"{exp_id}_filtered.csv",
                "text/csv",
            )

    # -------------------- Suggestions tab --------------------
    with tab_suggestions:
        render_suggestions(exp_id)

    # -------------------- DOE tab --------------------
    with tab_doe:
        render_doe()

    # -------------------- Analysis (HPLC) tab --------------------
    with tab_analysis:
        st.markdown("### HPLC Analysis")
        
        hplc_all = load_hplc_for_exp(exp_id)

        if hplc_all.empty:
            st.info("No HPLC data found for this experiment.")
        else:
            # Layout: Mapping Table (Left) + Chart (Right)
            col_map, col_chart = st.columns([1, 2])
            
            # Map logic
            vessels_order = sorted(df_exp["vessel"].unique())
            units = sorted(hplc_all["unit"].dropna().unique())
            unit_for_vessel = {v: units[i] for i, v in enumerate(vessels_order) if i < len(units)}
            
            with col_map:
                st.markdown("#### Vessel Mapping")
                if not unit_for_vessel:
                    st.warning("Mapping failed.")
                else:
                    map_df = pd.DataFrame([
                         {"Vessel": v, "HPLC Unit": unit_for_vessel.get(v, "—")}
                         for v in vessels_order
                    ]).set_index("Vessel")
                    st.dataframe(map_df, use_container_width=True)

            with col_chart:
                 # Plot plotly
                fig = go.Figure()
                for v in vessels_selected:
                    unit = unit_for_vessel.get(v)
                    if not unit: continue
                    sub = hplc_all[hplc_all["unit"] == unit].sort_values("eft_hours")
                    if sub.empty: continue
                    
                    fig.add_trace(go.Scatter(x=sub["eft_hours"], y=sub["glycerol_gL"], mode="lines+markers", name=f"{v} Gly"))
                    fig.add_trace(go.Scatter(x=sub["eft_hours"], y=sub["methanol_gL"], mode="lines+markers", name=f"{v} MeOH", yaxis="y2"))
                
                fig.update_layout(
                    template="plotly_white",
                    height=450,
                    margin=dict(l=20, r=20, t=30, b=20),
                    xaxis_title="EFT (h)",
                    yaxis_title="Glycerol (g/L)",
                    yaxis2=dict(title="Methanol (g/L)", overlaying="y", side="right"),
                    legend=dict(orientation="h", y=1.1)
                )
                st.plotly_chart(fig, use_container_width=True)


if __name__ == "__main__":
    main()
