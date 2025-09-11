# -*- coding: utf-8 -*-
"""
Utilities for collating, analyzing, and plotting docking and training results.
This module provides helpers to load/merge data, find key columns, convert SDF to CSV,
and plot training metrics.
"""

from __future__ import annotations

# Standard library imports
from pathlib import Path
from datetime import datetime
from typing import Optional, List, Iterable, Tuple, Dict, Any
import os
import re
import logging
import subprocess

# Third-party imports
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import yaml

# ---------------------------------------------------------------------------
# Configuration and constants
# ---------------------------------------------------------------------------

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

# Centralized defaults for commonly auto-detected columns
DEFAULT_SMILES_COL_CANDIDATES: List[str] = ["SMILES", "smiles", "Smiles", "SMILES_staged"]
DEFAULT_SCORE_COL_CANDIDATES: List[str] = [
    "m_score__min__energy",
    "maize (raw)",
    "maize",
    "Score",
    "score",
]


# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------

def _find_first_col(df: pd.DataFrame, candidates: Optional[List[str]] = None) -> Optional[str]:
    """
    Return the first column present in df from the provided candidates list.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe to inspect.
    candidates : Optional[List[str]]
        Ordered list of candidate column names to look for. If None, a union of
        defaults for SMILES and score columns is used.

    Returns
    -------
    Optional[str]
        Name of the first matching column, or None if not found.
    """
    if df is None or df.empty:
        return None
    _cands = (
        candidates
        if candidates is not None
        else list(dict.fromkeys(DEFAULT_SMILES_COL_CANDIDATES + DEFAULT_SCORE_COL_CANDIDATES))
    )
    for c in _cands:
        if c in df.columns:
            return c
    return None


def _find_smiles_col(df: pd.DataFrame) -> str:
    """Detect the SMILES column name."""
    col = _find_first_col(df, DEFAULT_SMILES_COL_CANDIDATES)
    if not col:
        raise KeyError(f"SMILES column not found in dataframe. Tried: {DEFAULT_SMILES_COL_CANDIDATES}")
    return col


def _find_score_col(df: pd.DataFrame) -> str:
    """Detect a score/energy column name."""
    col = _find_first_col(df, DEFAULT_SCORE_COL_CANDIDATES)
    if not col:
        raise KeyError(f"Score column not found in dataframe. Tried: {DEFAULT_SCORE_COL_CANDIDATES}")
    return col


def _lineplot(y_col: str, title: str, fname: str, data: pd.DataFrame) -> None:
    """
    Simple wrapper around seaborn lineplot, saving to file.
    Expects a column 'step' for x-axis if present; otherwise uses data index.
    """
    if y_col not in data.columns:
        logger.warning("Column '%s' not found, skipping plot '%s'", y_col, title)
        return
    plt.figure(figsize=(8, 4.5))
    x = "step" if "step" in data.columns else None
    sns.lineplot(data=data, x=x, y=y_col, marker="o")
    plt.title(title)
    plt.tight_layout()
    out = Path(fname)
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out.as_posix(), dpi=200)
    plt.close()
    logger.info("Saved plot: %s", out)


# ---------------------------------------------------------------------------
# Data loading and merging
# ---------------------------------------------------------------------------

def load_data(path: str | os.PathLike) -> pd.DataFrame:
    """PEP 8-compliant data loader. Reads CSV/TSV into a DataFrame with dtype inference."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() in {".tsv", ".tab"}:
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path)
    logger.info("Loaded %d rows from %s", len(df), path.name)
    return df


def merge_staged_with_data_by_smiles(staged_df: pd.DataFrame, data_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge staged learning results with external data using the SMILES string.

    The function tries to auto-detect the SMILES column in both frames, normalizes
    whitespace, and performs an inner join.
    """
    staged_smiles = _find_smiles_col(staged_df)
    data_smiles = _find_smiles_col(data_df)

    lhs = staged_df.copy()
    rhs = data_df.copy()

    lhs[staged_smiles] = lhs[staged_smiles].astype(str).str.strip()
    rhs[data_smiles] = rhs[data_smiles].astype(str).str.strip()

    merged = lhs.merge(rhs, left_on=staged_smiles, right_on=data_smiles, how="inner", suffixes=("", "_ext"))
    logger.info("Merged frames: %d (staged) x %d (data) -> %d (inner)", len(lhs), len(rhs), len(merged))
    return merged


# Backward-compatibility aliases (retain until all call-sites migrate)
def loadData(path: str | os.PathLike) -> pd.DataFrame:  # noqa: N802 - legacy alias
    return load_data(path)


def merge_staged_and_data_matching_smiles(staged_df: pd.DataFrame, data_df: pd.DataFrame) -> pd.DataFrame:  # noqa: N802
    return merge_staged_with_data_by_smiles(staged_df, data_df)


# ---------------------------------------------------------------------------
# File system helpers
# ---------------------------------------------------------------------------

def get_files_after_datetime(directory: str | os.PathLike, target_datetime: datetime) -> list[Path]:
    """
    Return files in 'directory' whose modification time is after 'target_datetime'.
    """
    base = Path(directory)
    if not base.exists():
        raise FileNotFoundError(base)
    out: list[Path] = []
    cutoff = target_datetime.timestamp()
    for p in base.iterdir():
        if p.is_file():
            try:
                if p.stat().st_mtime > cutoff:
                    out.append(p)
            except OSError:
                continue
    out.sort()
    logger.info("Found %d files newer than %s in %s", len(out), target_datetime.isoformat(), base)
    return out


# ---------------------------------------------------------------------------
# Data analysis helpers
# ---------------------------------------------------------------------------

def pick_lowest_energy(df: pd.DataFrame) -> pd.Series:
    """
    Pick the row with the lowest energy/score from df.

    Uses the first detectable score column from DEFAULT_SCORE_COL_CANDIDATES.
    """
    score_col = _find_score_col(df)
    row = df.loc[df[score_col].astype(float).idxmin()]
    logger.info("Lowest energy row has %s = %s", score_col, row[score_col])
    return row


def run_collate(staged_path: str | os.PathLike, data_path: str | os.PathLike, out_csv: str | os.PathLike) -> Path:
    """
    Load staged and external data, merge on SMILES, and write the result to CSV.
    Returns the output path.
    """
    staged_df = load_data(staged_path)
    data_df = load_data(data_path)
    merged = merge_staged_with_data_by_smiles(staged_df, data_df)
    out = Path(out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out, index=False)
    logger.info("Wrote merged CSV (%d rows): %s", len(merged), out)
    return out


def plot_training_metrics(training_csv: str | os.PathLike, out_dir: str | os.PathLike = "plots") -> None:
    """
    Generate a few standard training metric plots from a staged learning CSV.
    Expects columns like: step, Score, maize, 'maize (raw)', QED, etc.
    """
    df = load_data(training_csv)
    out_dir = Path(out_dir)

    # Choose some common metrics to visualize if present
    candidates = ["Score", "maize", "maize (raw)", "QED", "QED (raw)"]
    for metric in candidates:
        if metric in df.columns:
            _lineplot(metric, f"{metric} over steps", out_dir / f"{metric}_over_steps.png", df)

    # If both 'maize' and 'Score' exist, plot them together for comparison
    if {"maize", "Score"}.issubset(df.columns):
        plt.figure(figsize=(8, 4.5))
        x = "step" if "step" in df.columns else None
        sns.lineplot(data=df, x=x, y="maize", label="maize", marker="o")
        sns.lineplot(data=df, x=x, y="Score", label="Score", marker="o")
        plt.title("maize and Score over steps")
        plt.legend()
        plt.tight_layout()
        out = out_dir / "maize_and_Score_over_steps.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out.as_posix(), dpi=200)
        plt.close()
        logger.info("Saved plot: %s", out)


# ---------------------------------------------------------------------------
# SDF conversion (pure Python minimal parser)
# ---------------------------------------------------------------------------

def _parse_sdf_blocks(lines: Iterable[str]) -> Iterable[List[str]]:
    """Yield SDF molecule blocks (list of lines) delimited by '$$$$'."""
    buf: List[str] = []
    for line in lines:
        if line.strip() == "$$$$":
            if buf:
                yield buf
                buf = []
        else:
            buf.append(line.rstrip("\n"))
    if buf:
        # Some SDFs may not end with '$$$$'
        yield buf


def _parse_sdf_properties(block: List[str]) -> Dict[str, Any]:
    """
    Extract property fields from an SDF molecule block. Coordinates are ignored.
    Returns a dict of property_name -> value (string).
    """
    props: Dict[str, Any] = {}
    # Properties start at lines with '> <Name>'
    i = 0
    n = len(block)
    while i < n:
        line = block[i]
        if line.startswith(">"):
            m = re.match(r">+\s*<([^>]+)>", line)
            if m:
                key = m.group(1).strip()
                i += 1
                values: List[str] = []
                # Read until empty line or next property header
                while i < n and not block[i].startswith(">") and block[i].strip() != "":
                    values.append(block[i])
                    i += 1
                props[key] = "\n".join(values).strip()
                continue
        i += 1
    return props


def sdf_to_csv(sdf_file: str | os.PathLike) -> pd.DataFrame:
    """
    Convert an SDF file to a pandas DataFrame by extracting property blocks.
    Writes a CSV next to the input with the same stem.

    Notes
    -----
    - This parser does not extract atom/bond coordinates; it focuses on SDF data fields.
    - If a property appears multiple times across molecules, columns will be unioned.
    """
    sdf_path = Path(sdf_file)
    if not sdf_path.exists():
        raise FileNotFoundError(sdf_path)

    records: List[Dict[str, Any]] = []
    with sdf_path.open("r", encoding="utf-8", errors="ignore") as f:
        for block in _parse_sdf_blocks(f):
            props = _parse_sdf_properties(block)
            records.append(props)

    if not records:
        logger.warning("No molecule records parsed from SDF: %s", sdf_path)

    df = pd.DataFrame.from_records(records)
    out_csv = sdf_path.with_suffix(".csv")
    df.to_csv(out_csv, index=False)
    logger.info("Parsed %d molecules from %s -> %s", len(df), sdf_path.name, out_csv.name)
    return df


# ---------------------------------------------------------------------------
# Text parsing and YAML writing
# ---------------------------------------------------------------------------

def _parse_ic50_from_text(text: str) -> Optional[Tuple[float, str]]:
    """
    Parse an IC50-like expression from text and return (value_in_nM, original_unit).

    Supports units: nM, uM/µM, mM. If multiple matches exist, returns the first.
    """
    if not text:
        return None
    pattern = r"(?i)\bIC50\s*[:=]?\s*([0-9]*\.?[0-9]+)\s*(nM|µM|uM|mM)\b"
    m = re.search(pattern, text)
    if not m:
        return None
    val = float(m.group(1))
    unit = m.group(2)
    unit_norm = unit.replace("µ", "u")
    # Convert to nM
    factor = {"nM": 1.0, "uM": 1e3, "mM": 1e6}[unit_norm]
    return val * factor, unit


def _write_boltz_yaml(yaml_path: Path, prot_id: str, prot_seq: str, lig_id: str, smiles: str, msa: str) -> None:
    """
    Write a minimal YAML configuration for an external 'boltz' predictor.
    """
    payload = {
        "version": 1,
        "sequences": [
            {"protein": {"id": prot_id, "sequence": prot_seq, "msa": msa}},
            {"ligand": {"id": lig_id, "smiles": smiles}},
        ],
        "properties": [{"affinity": {"binder": lig_id}}],
    }
    yaml_path.parent.mkdir(parents=True, exist_ok=True)
    with yaml_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(payload, f, sort_keys=False)
    logger.info("Wrote boltz YAML: %s", yaml_path)


def predict_ic50_with_boltz_conda(
        boltz_exe: str,
        prot_id: str,
        prot_seq: str,
        lig_id: str,
        smiles: str,
        msa: str,
        work_dir: str | os.PathLike = ".",
) -> Dict[str, Any]:
    """
    Prepare a boltz prediction job and return the command information.

    This function writes the YAML needed by 'boltz' and returns the command to run.
    Execution is not performed here to keep IO side-effects explicit for callers.
    """
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    yaml_path = work / f"{prot_id}_{lig_id}.yaml"
    _write_boltz_yaml(yaml_path, prot_id, prot_seq, lig_id, smiles, msa)

    cmd = [boltz_exe, "predict", yaml_path.as_posix()]
    logger.info("Prepared boltz command: %s", " ".join(cmd))
    return {"yaml": yaml_path, "cmd": cmd}


# ---------------------------------------------------------------------------
# Example / convenience functions
# ---------------------------------------------------------------------------

def test() -> pd.DataFrame:
    """
    Convenience test: load 'test.csv' from the current directory if present.
    Returns the DataFrame (or empty DataFrame if not found).
    """
    p = Path("test.csv")
    if p.exists():
        return load_data(p)
    logger.warning("test.csv not found in current directory")
    return pd.DataFrame()


def quick_prediction(training_csv: str | os.PathLike) -> Dict[str, Any]:
    """
    Quick summary over a staged learning CSV: returns simple statistics.
    """
    df = load_data(training_csv)
    out: Dict[str, Any] = {"rows": len(df), "columns": list(df.columns)}
    for col in ["Score", "maize", "maize (raw)", "QED", "QED (raw)"]:
        if col in df.columns:
            s = pd.to_numeric(df[col], errors="coerce").dropna()
            if not s.empty:
                out[f"{col}_mean"] = float(s.mean())
                out[f"{col}_min"] = float(s.min())
                out[f"{col}_max"] = float(s.max())
    score_col = None
    try:
        score_col = _find_score_col(df)
    except KeyError:
        pass
    if score_col:
        row = pick_lowest_energy(df)
        out["lowest_energy_score_col"] = score_col
        out["lowest_energy_row_index"] = int(row.name) if hasattr(row, "name") else None
        out["lowest_energy_value"] = float(row[score_col])
    return out




