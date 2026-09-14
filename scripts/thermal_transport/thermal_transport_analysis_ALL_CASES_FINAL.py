#!/usr/bin/env python3
"""Unified NEMD thermal-transport analysis for h-AlN project datasets.

Supported input layouts
=======================

1) Ordinary seed replicates or one-seed dopant realizations (wide)

   Length (nm) | Seed 1 | Seed 2 | ... | Seed 5
   Length (nm) | r1     | r2     | ... | r5

2) Hierarchical dopant realizations x seeds (wide)

   Length (nm) | R1_S1 | R1_S2 | R1_S3 | ... | R5_S3

   One seed per realization is also accepted:
   Length (nm) | R1_S1 | R2_S1 | ... | R5_S1

3) Hierarchical long format

   Length (nm) | Realization | Seed | Kappa

The first column may instead be Width, Delta, Strain, Concentration, Porosity,
Temperature, or Direction. In those cases the script calculates point-wise mean
conductivity and 95% Student-t confidence intervals but NEVER performs an
infinite-length extrapolation.

Automatic analysis
==================

* Genuine transport-length data with >=3 lengths:
  - length-wise mean kappa(L), SD, SE, two-sided 95% Student-t CI;
  - reciprocal fit 1/kappa(L) = a + b/L;
  - kappa_inf = 1/a;
  - lambda_eff = b/a;
  - 10,000-iteration bootstrap by default;
  - hierarchical bootstrap for realization x seed data;
  - residual, shortest-length, and leave-one-length-out diagnostics.

* One-length data:
  - finite-length kappa(L), SD, SE, and 95% Student-t CI;
  - no kappa_inf or lambda_eff.

* Width series:
  - point-wise statistics;
  - consecutive-width percent changes and CI-overlap convergence checks.

* Delta series:
  - point-wise statistics;
  - comparison against Delta=0.05 by default;
  - percent-change and CI-overlap checks;
  - reminds the user that temperature-profile and cumulative-energy linearity
    must still be checked from the underlying NEMD outputs.

* Direction series:
  - point-wise statistics;
  - two-direction percent difference and CI-overlap isotropy check.

* Hierarchical realization x seed data with >=2 seeds:
  - seed-noise audit at every supplied x value;
  - compares the preselected first-seed result with the all-seed realization result;
  - compares pooled within-realization seed SD with between-realization SD.

The extracted lambda_eff is an effective characteristic length of the selected
finite-size suppression model, not a mode-resolved phonon MFP spectrum.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.drawing.image import Image as XLImage
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from scipy.stats import t

DEFAULT_ITERATIONS = 10_000
DEFAULT_RANDOM_SEED = 20260728
DEFAULT_DPI = 1200


@dataclass(frozen=True)
class FitResult:
    intercept: float
    slope: float
    r2: float
    adjusted_r2: float
    rmse: float
    intercept_se: float
    slope_se: float
    fitted_y: np.ndarray
    residuals: np.ndarray
    standardized_residuals: np.ndarray

    @property
    def kappa_inf(self) -> float:
        with np.errstate(divide="ignore", invalid="ignore"):
            return float(np.divide(1.0, self.intercept))

    @property
    def lambda_eff(self) -> float:
        with np.errstate(divide="ignore", invalid="ignore"):
            return float(np.divide(self.slope, self.intercept))


@dataclass(frozen=True)
class ParsedData:
    mode: str  # simple or hierarchical
    layout: str  # wide or long
    selected_sheet: str
    x_header: str
    x_kind: str
    x_values: list[Any]
    x_numeric: np.ndarray | None
    simple_names: list[str] | None
    simple_values: list[np.ndarray] | None
    realization_ids: list[str] | None
    seed_ids: list[str] | None
    hierarchical_values: list[np.ndarray] | None
    raw_dataframe: pd.DataFrame


@dataclass(frozen=True)
class Analysis:
    case_name: str
    analysis_type: str  # length_scaling, finite_single, series
    x_kind: str
    data_mode: str
    x_header: str
    x_values: list[Any]
    x_numeric: np.ndarray | None
    unit_names: list[list[str]]
    units: list[np.ndarray]
    realization_means: list[np.ndarray] | None
    n_units: np.ndarray
    mean: np.ndarray
    sd: np.ndarray
    se: np.ndarray
    t_critical: np.ndarray
    ci_half: np.ndarray
    ci_lower: np.ndarray
    ci_upper: np.ndarray
    point_bootstrap_summary: pd.DataFrame
    point_bootstrap_samples: list[np.ndarray]
    fit: FitResult | None
    reciprocal_x: np.ndarray | None
    reciprocal_y: np.ndarray | None
    reciprocal_neg_error: np.ndarray | None
    reciprocal_pos_error: np.ndarray | None
    bootstrap_intercept: np.ndarray | None
    bootstrap_slope: np.ndarray | None
    bootstrap_kappa_inf: np.ndarray | None
    bootstrap_lambda_eff: np.ndarray | None
    bootstrap_physical: np.ndarray | None
    fit_bootstrap_summary: pd.DataFrame
    shortest_fit: FitResult | None
    leave_one_out: pd.DataFrame
    seed_audit: pd.DataFrame
    special_checks: pd.DataFrame
    diagnostics: pd.DataFrame
    thickness_nm: float | None
    exact_concentration_pct: float | None


def normalize_header(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    # Normalize common scientific symbols used in worksheet headers.
    replacements = {"Δ": "Delta", "δ": "delta", "κ": "kappa", "λ": "lambda", "∞": "inf", "%": " percent "}
    for old, new in replacements.items():
        text = text.replace(old, new)
    return re.sub(r"\s+", " ", text)


def compact_header(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", normalize_header(value).lower())


def natural_key(value: object) -> list[Any]:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", str(value))]


def detect_x_kind(header: str) -> str:
    compact = compact_header(header)
    if compact in {"l", "lnm", "length", "lengthnm", "transportlength", "transportlengthnm"} or "length" in compact:
        return "length"
    if "width" in compact:
        return "width"
    if compact in {"delta", "deltat", "temperaturebias", "thermostatbias"} or "delta" in compact:
        return "delta"
    if "strain" in compact:
        return "strain"
    if "concentration" in compact or "doping" in compact or "dopant" in compact:
        return "concentration"
    if "porosity" in compact or compact in {"pore", "porepercent", "porepct"}:
        return "porosity"
    if "temperature" in compact or compact in {"t", "tk"}:
        return "temperature"
    if "direction" in compact or compact in {"dir", "orientation"}:
        return "direction"
    return "screening"




def is_metadata_column(header: str) -> bool:
    compact = compact_header(header)
    metadata_tokens = {
        "nominalconcentration", "exactconcentration", "realizedconcentration",
        "realiseddopantconcentration", "dopantcount", "hostsites", "alsites",
        "temperaturek", "caseid", "casename", "notes", "comment", "units",
    }
    return compact in metadata_tokens or compact.startswith("exactrealized")


def find_column(columns: Sequence[str], candidates: set[str]) -> str | None:
    for column in columns:
        if compact_header(column) in candidates:
            return column
    return None


def parse_realization_seed_header(header: str) -> tuple[int, int] | None:
    text = normalize_header(header).lower()
    patterns = [
        r"(?:realization|realisation|real|r)\s*[_\- ]*(\d+)\s*[_\- ]*(?:seed|s)\s*[_\- ]*(\d+)",
        r"(?:seed|s)\s*[_\- ]*(\d+)\s*[_\- ]*(?:realization|realisation|real|r)\s*[_\- ]*(\d+)",
    ]
    match = re.fullmatch(patterns[0], text)
    if match:
        return int(match.group(1)), int(match.group(2))
    match = re.fullmatch(patterns[1], text)
    if match:
        return int(match.group(2)), int(match.group(1))
    compact = compact_header(header)
    for pattern in [r"r(\d+)s(\d+)", r"realization(\d+)seed(\d+)", r"realisation(\d+)seed(\d+)"]:
        match = re.fullmatch(pattern, compact)
        if match:
            return int(match.group(1)), int(match.group(2))
    return None


def load_sheet(path: Path, sheet: str | None, header_row: int) -> tuple[pd.DataFrame, str]:
    if not path.exists():
        raise FileNotFoundError(f"Input file does not exist: {path}")
    suffix = path.suffix.lower()
    if suffix == ".csv":
        frame = pd.read_csv(path, header=header_row - 1)
        selected_sheet = "CSV"
    elif suffix in {".xlsx", ".xlsm"}:
        excel = pd.ExcelFile(path)
        selected_sheet = sheet or excel.sheet_names[0]
        if selected_sheet not in excel.sheet_names:
            raise ValueError(f"Sheet '{selected_sheet}' not found. Available sheets: {excel.sheet_names}")
        frame = pd.read_excel(path, sheet_name=selected_sheet, header=header_row - 1)
    else:
        raise ValueError("Supported input formats are .xlsx, .xlsm, and .csv")

    frame.columns = [normalize_header(c) for c in frame.columns]
    frame = frame.loc[:, [c for c in frame.columns if c and not c.lower().startswith("unnamed")]]
    frame = frame.dropna(how="all").reset_index(drop=True)
    if frame.empty:
        raise ValueError("The selected worksheet contains no data.")
    return frame, selected_sheet


def choose_x_column(columns: list[str], requested: str | None) -> str:
    if requested:
        for column in columns:
            if column.lower() == requested.strip().lower():
                return column
        raise ValueError(f"X column '{requested}' was not found.")
    non_data = {
        "realization", "realisation", "realizationid", "realisationid", "real", "r",
        "seed", "seedid", "s", "kappa", "thermalconductivity", "conductivity",
        "kappawmk", "thermalconductivitywmk",
    }
    candidates = [c for c in columns if compact_header(c) not in non_data]
    if not candidates:
        raise ValueError("Could not determine the independent-variable column.")
    known = [c for c in candidates if detect_x_kind(c) != "screening"]
    return known[0] if known else candidates[0]


def convert_x_values(series: pd.Series, kind: str) -> tuple[list[Any], np.ndarray | None, pd.Series]:
    if kind == "direction" or kind == "screening":
        numeric = pd.to_numeric(series, errors="coerce")
        if numeric.notna().all():
            order = np.argsort(numeric.to_numpy(dtype=float))
            values = numeric.iloc[order].to_numpy(dtype=float).tolist()
            return values, np.asarray(values, dtype=float), pd.Series(order, index=series.index)
        values = series.astype(str).str.strip().tolist()
        return values, None, pd.Series(np.arange(len(series)), index=series.index)
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.isna().any():
        bad = series[numeric.isna()].tolist()
        raise ValueError(f"The '{kind}' x column must be numeric. Invalid values: {bad[:5]}")
    values = numeric.to_numpy(dtype=float)
    order = np.argsort(values)
    return values[order].tolist(), values[order], pd.Series(order, index=series.index)


def parse_long_hierarchical(
    frame: pd.DataFrame, selected_sheet: str, x_column: str | None
) -> ParsedData | None:
    columns = list(frame.columns)
    realization_col = find_column(columns, {
        "realization", "realisation", "realizationid", "realisationid", "real", "r"
    })
    seed_col = find_column(columns, {"seed", "seedid", "s"})
    kappa_col = find_column(columns, {
        "kappa", "thermalconductivity", "conductivity", "kappawmk", "thermalconductivitywmk"
    })
    if not all([realization_col, seed_col, kappa_col]):
        return None
    x_col = choose_x_column(columns, x_column)
    kind = detect_x_kind(x_col)

    clean = pd.DataFrame({
        "X": frame[x_col],
        "Realization": frame[realization_col].astype(str).str.strip(),
        "Seed": frame[seed_col].astype(str).str.strip(),
        "Kappa": pd.to_numeric(frame[kappa_col], errors="coerce"),
    }).dropna(subset=["X", "Kappa"])
    if clean.empty:
        raise ValueError("No valid hierarchical long-format rows were found.")
    if (clean["Kappa"] <= 0).any():
        raise ValueError("All conductivity values must be strictly positive.")
    if clean.duplicated(["X", "Realization", "Seed"]).any():
        raise ValueError("Duplicate x-realization-seed entries were found.")

    if kind in {"length", "width", "delta", "strain", "concentration", "porosity", "temperature"}:
        clean["X_numeric"] = pd.to_numeric(clean["X"], errors="coerce")
        if clean["X_numeric"].isna().any():
            raise ValueError(f"The '{kind}' x values must be numeric.")
        x_values = sorted(clean["X_numeric"].unique().astype(float).tolist())
        x_numeric = np.asarray(x_values, dtype=float)
        clean["X_key"] = clean["X_numeric"]
    else:
        x_values = list(dict.fromkeys(clean["X"].astype(str).str.strip().tolist()))
        x_numeric = None
        clean["X_key"] = clean["X"].astype(str).str.strip()

    realization_ids = sorted(clean["Realization"].unique().tolist(), key=natural_key)
    seed_ids = sorted(clean["Seed"].unique().tolist(), key=natural_key)
    matrices: list[np.ndarray] = []
    for x in x_values:
        subset = clean.loc[clean["X_key"] == x]
        pivot = subset.pivot(index="Realization", columns="Seed", values="Kappa")
        pivot = pivot.reindex(index=realization_ids, columns=seed_ids)
        if pivot.isna().any().any():
            raise ValueError(
                f"X={x} has missing realization-seed values. A complete rectangular block is required."
            )
        matrices.append(pivot.to_numpy(dtype=float))

    raw = clean[["X", "Realization", "Seed", "Kappa"]].rename(columns={"X": x_col})
    return ParsedData(
        mode="hierarchical", layout="long", selected_sheet=selected_sheet,
        x_header=x_col, x_kind=kind, x_values=x_values, x_numeric=x_numeric,
        simple_names=None, simple_values=None,
        realization_ids=[str(v) for v in realization_ids],
        seed_ids=[str(v) for v in seed_ids], hierarchical_values=matrices,
        raw_dataframe=raw,
    )


def parse_wide(
    frame: pd.DataFrame, selected_sheet: str, x_column: str | None
) -> ParsedData:
    columns = list(frame.columns)
    x_col = choose_x_column(columns, x_column)
    kind = detect_x_kind(x_col)

    work = frame.loc[frame[x_col].notna()].copy().reset_index(drop=True)
    if work.empty:
        raise ValueError("No valid x values were found.")

    # Sort numeric x types, preserve categorical order otherwise.
    if kind in {"length", "width", "delta", "strain", "concentration", "porosity", "temperature"}:
        work[x_col] = pd.to_numeric(work[x_col], errors="coerce")
        if work[x_col].isna().any():
            raise ValueError(f"The '{kind}' x column must be numeric.")
        if kind == "length" and (work[x_col] <= 0).any():
            raise ValueError("All transport lengths must be strictly positive.")
        work = work.sort_values(x_col).reset_index(drop=True)
        x_values = work[x_col].to_numpy(dtype=float).tolist()
        x_numeric = np.asarray(x_values, dtype=float)
    else:
        x_values = work[x_col].astype(str).str.strip().tolist()
        numeric_try = pd.to_numeric(work[x_col], errors="coerce")
        x_numeric = numeric_try.to_numpy(dtype=float) if numeric_try.notna().all() else None

    if len(set(map(str, x_values))) != len(x_values):
        raise ValueError("Duplicate x rows were found. Combine replicates into one row per condition.")

    parsed_columns: dict[str, tuple[int, int]] = {}
    for column in columns:
        if column == x_col:
            continue
        parsed = parse_realization_seed_header(column)
        if parsed is not None:
            parsed_columns[column] = parsed

    non_x_columns = [c for c in columns if c != x_col and not is_metadata_column(c)]
    # Use hierarchical mode when every non-x data column follows R#_S# notation.
    if parsed_columns and len(parsed_columns) == len(non_x_columns):
        realization_numbers = sorted({pair[0] for pair in parsed_columns.values()})
        seed_numbers = sorted({pair[1] for pair in parsed_columns.values()})
        if len(realization_numbers) < 2:
            raise ValueError("At least two dopant realizations are required.")
        # One seed per realization is valid and is treated as realization-level sampling.
        expected = {(r, s) for r in realization_numbers for s in seed_numbers}
        observed = set(parsed_columns.values())
        missing = sorted(expected - observed)
        if missing:
            raise ValueError(f"Missing realization-seed columns: {missing}")
        pair_to_column = {pair: column for column, pair in parsed_columns.items()}
        matrices: list[np.ndarray] = []
        for _, row in work.iterrows():
            matrix = np.empty((len(realization_numbers), len(seed_numbers)), dtype=float)
            for i, realization in enumerate(realization_numbers):
                for j, seed in enumerate(seed_numbers):
                    value = pd.to_numeric(row[pair_to_column[(realization, seed)]], errors="coerce")
                    if pd.isna(value) or float(value) <= 0:
                        raise ValueError(
                            f"Missing or nonpositive conductivity at x={row[x_col]}, "
                            f"realization {realization}, seed {seed}."
                        )
                    matrix[i, j] = float(value)
            matrices.append(matrix)
        raw_columns = [x_col] + [pair_to_column[p] for p in sorted(expected)]
        return ParsedData(
            mode="hierarchical", layout="wide", selected_sheet=selected_sheet,
            x_header=x_col, x_kind=kind, x_values=x_values, x_numeric=x_numeric,
            simple_names=None, simple_values=None,
            realization_ids=[f"R{r}" for r in realization_numbers],
            seed_ids=[f"S{s}" for s in seed_numbers], hierarchical_values=matrices,
            raw_dataframe=work.copy(),
        )

    # Ordinary seeds or one-seed realization columns such as r1...r5.
    replicate_columns: list[str] = []
    for column in non_x_columns:
        numeric = pd.to_numeric(work[column], errors="coerce")
        if numeric.notna().any():
            replicate_columns.append(column)
    if len(replicate_columns) < 2:
        raise ValueError("At least two independent seed/realization columns are required.")

    values: list[np.ndarray] = []
    for _, row in work.iterrows():
        row_values = pd.to_numeric(row[replicate_columns], errors="coerce").dropna().to_numpy(dtype=float)
        if len(row_values) < 2:
            raise ValueError(f"Condition {row[x_col]} contains fewer than two valid values.")
        if np.any(row_values <= 0):
            raise ValueError("All conductivity values must be strictly positive.")
        values.append(row_values)

    return ParsedData(
        mode="simple", layout="wide", selected_sheet=selected_sheet,
        x_header=x_col, x_kind=kind, x_values=x_values, x_numeric=x_numeric,
        simple_names=replicate_columns, simple_values=values,
        realization_ids=None, seed_ids=None, hierarchical_values=None,
        raw_dataframe=work.copy(),
    )


def parse_input(path: Path, sheet: str | None, header_row: int, x_column: str | None) -> ParsedData:
    frame, selected_sheet = load_sheet(path, sheet, header_row)
    long_data = parse_long_hierarchical(frame, selected_sheet, x_column)
    return long_data if long_data is not None else parse_wide(frame, selected_sheet, x_column)


def ordinary_fit(x: np.ndarray, y: np.ndarray) -> FitResult:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3:
        raise ValueError("At least three lengths are mathematically required for extrapolation.")
    design = np.column_stack([np.ones_like(x), x])
    beta, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
    if rank < 2:
        raise ValueError("The reciprocal fit is rank-deficient.")
    intercept, slope = map(float, beta)
    fitted = design @ beta
    residuals = y - fitted
    sse = float(np.sum(residuals**2))
    sst = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - sse / sst if sst > 0 else float("nan")
    n = len(x)
    dof = n - 2
    adjusted = 1.0 - (1.0 - r2) * (n - 1) / dof if dof > 0 and np.isfinite(r2) else float("nan")
    rmse = math.sqrt(sse / dof) if dof > 0 else float("nan")
    if dof > 0:
        covariance = (sse / dof) * np.linalg.inv(design.T @ design)
        intercept_se = math.sqrt(max(float(covariance[0, 0]), 0.0))
        slope_se = math.sqrt(max(float(covariance[1, 1]), 0.0))
        standardized = residuals / rmse if rmse > 0 else np.zeros_like(residuals)
    else:
        intercept_se = slope_se = float("nan")
        standardized = np.full_like(residuals, np.nan)
    return FitResult(
        intercept=intercept, slope=slope, r2=r2, adjusted_r2=adjusted,
        rmse=rmse, intercept_se=intercept_se, slope_se=slope_se,
        fitted_y=fitted, residuals=residuals, standardized_residuals=standardized,
    )


def vectorized_fit(x: np.ndarray, y_matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    y_matrix = np.asarray(y_matrix, dtype=float)
    xc = x - np.mean(x)
    denominator = float(np.sum(xc**2))
    if denominator <= 0:
        raise ValueError("Reciprocal transport lengths are identical.")
    ym = np.mean(y_matrix, axis=1)
    slopes = ((y_matrix - ym[:, None]) @ xc) / denominator
    intercepts = ym - slopes * np.mean(x)
    return intercepts, slopes


def t_statistics(units: list[np.ndarray]) -> tuple[np.ndarray, ...]:
    n = np.array([len(v) for v in units], dtype=int)
    if np.any(n < 2):
        raise ValueError("At least two independent statistical units are required per condition.")
    mean = np.array([np.mean(v) for v in units], dtype=float)
    sd = np.array([np.std(v, ddof=1) for v in units], dtype=float)
    se = sd / np.sqrt(n)
    critical = t.ppf(0.975, n - 1)
    half = critical * se
    return n, mean, sd, se, critical, half, mean - half, mean + half


def bootstrap_point_simple(values: np.ndarray, iterations: int, rng: np.random.Generator) -> np.ndarray:
    n = len(values)
    indices = rng.integers(0, n, size=(iterations, n))
    return np.mean(values[indices], axis=1)


def bootstrap_point_hierarchical(matrix: np.ndarray, iterations: int, rng: np.random.Generator) -> np.ndarray:
    n_real, n_seed = matrix.shape
    output = np.empty(iterations, dtype=float)
    for b in range(iterations):
        sampled_real = rng.integers(0, n_real, size=n_real)
        realization_means = np.empty(n_real, dtype=float)
        for i, real_idx in enumerate(sampled_real):
            sampled_seed = rng.integers(0, n_seed, size=n_seed)
            realization_means[i] = float(np.mean(matrix[real_idx, sampled_seed]))
        output[b] = float(np.mean(realization_means))
    return output


def bootstrap_scaling_simple(
    values_by_length: list[np.ndarray], x: np.ndarray, iterations: int,
    rng: np.random.Generator, paired_units: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    means = np.empty((iterations, len(values_by_length)), dtype=float)
    if paired_units:
        sizes = {len(v) for v in values_by_length}
        if len(sizes) != 1:
            raise ValueError("Paired-unit bootstrap requires the same number of units at every length.")
        n = next(iter(sizes))
        indices = rng.integers(0, n, size=(iterations, n))
        for j, values in enumerate(values_by_length):
            means[:, j] = np.mean(values[indices], axis=1)
    else:
        for j, values in enumerate(values_by_length):
            n = len(values)
            indices = rng.integers(0, n, size=(iterations, n))
            means[:, j] = np.mean(values[indices], axis=1)
    intercept, slope = vectorized_fit(x, 1.0 / means)
    with np.errstate(divide="ignore", invalid="ignore"):
        kappa_inf = 1.0 / intercept
        lambda_eff = slope / intercept
    physical = (
        np.isfinite(intercept) & np.isfinite(slope) & np.isfinite(kappa_inf)
        & np.isfinite(lambda_eff) & (intercept > 0) & (slope > 0)
        & (kappa_inf > 0) & (lambda_eff > 0)
    )
    return intercept, slope, kappa_inf, lambda_eff, physical


def bootstrap_scaling_hierarchical(
    matrices: list[np.ndarray], x: np.ndarray, iterations: int,
    rng: np.random.Generator, paired_realizations: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    shapes = {matrix.shape for matrix in matrices}
    if len(shapes) != 1:
        raise ValueError("Every length must have the same realization x seed dimensions.")
    n_real, n_seed = next(iter(shapes))
    means = np.empty((iterations, len(matrices)), dtype=float)
    paired_indices = rng.integers(0, n_real, size=(iterations, n_real)) if paired_realizations else None
    for b in range(iterations):
        for j, matrix in enumerate(matrices):
            sampled_real = paired_indices[b] if paired_indices is not None else rng.integers(0, n_real, size=n_real)
            realization_means = np.empty(n_real, dtype=float)
            for i, real_idx in enumerate(sampled_real):
                sampled_seed = rng.integers(0, n_seed, size=n_seed)
                realization_means[i] = np.mean(matrix[real_idx, sampled_seed])
            means[b, j] = np.mean(realization_means)
    intercept, slope = vectorized_fit(x, 1.0 / means)
    with np.errstate(divide="ignore", invalid="ignore"):
        kappa_inf = 1.0 / intercept
        lambda_eff = slope / intercept
    physical = (
        np.isfinite(intercept) & np.isfinite(slope) & np.isfinite(kappa_inf)
        & np.isfinite(lambda_eff) & (intercept > 0) & (slope > 0)
        & (kappa_inf > 0) & (lambda_eff > 0)
    )
    return intercept, slope, kappa_inf, lambda_eff, physical


def percentile_row(values: np.ndarray, quantity: str, units: str) -> dict[str, object]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return {
            "Quantity": quantity, "Median": np.nan, "Lower_95": np.nan,
            "Upper_95": np.nan, "Lower_Error": np.nan, "Upper_Error": np.nan,
            "Mean": np.nan, "SD": np.nan, "Units": units, "Finite_Samples": 0,
        }
    lower, median, upper = np.percentile(finite, [2.5, 50, 97.5])
    return {
        "Quantity": quantity, "Median": float(median), "Lower_95": float(lower),
        "Upper_95": float(upper), "Lower_Error": float(median - lower),
        "Upper_Error": float(upper - median), "Mean": float(np.mean(finite)),
        "SD": float(np.std(finite, ddof=1)) if finite.size > 1 else 0.0,
        "Units": units, "Finite_Samples": int(finite.size),
    }


def confidence_intervals_overlap(low1: float, high1: float, low2: float, high2: float) -> bool:
    return max(low1, low2) <= min(high1, high2)


def build_seed_audit(
    parsed: ParsedData, mean: np.ndarray, ci_lower: np.ndarray, ci_upper: np.ndarray,
    preselected_seed: int,
) -> pd.DataFrame:
    if parsed.mode != "hierarchical" or parsed.hierarchical_values is None:
        return pd.DataFrame()
    if parsed.hierarchical_values[0].shape[1] < 2:
        return pd.DataFrame()
    n_seed_available = parsed.hierarchical_values[0].shape[1]
    if preselected_seed < 1 or preselected_seed > n_seed_available:
        raise ValueError(
            f"--preselected-seed must be between 1 and {n_seed_available} for this dataset."
        )
    seed_index = preselected_seed - 1
    rows: list[dict[str, object]] = []
    for j, matrix in enumerate(parsed.hierarchical_values):
        n_real, n_seed = matrix.shape
        realization_means = np.mean(matrix, axis=1)
        selected_seed_values = matrix[:, seed_index]
        selected_mean = float(np.mean(selected_seed_values))
        selected_sd = float(np.std(selected_seed_values, ddof=1))
        selected_half = float(t.ppf(0.975, n_real - 1) * selected_sd / math.sqrt(n_real))
        within_variances = np.var(matrix, axis=1, ddof=1)
        pooled_within_sd = float(math.sqrt(np.mean(within_variances)))
        between_sd = float(np.std(realization_means, ddof=1))
        inside = float(ci_lower[j]) <= selected_mean <= float(ci_upper[j])
        variance_pass = pooled_within_sd < between_sd
        rows.append({
            "X": parsed.x_values[j],
            "Realizations": n_real,
            "Seeds_per_realization": n_seed,
            "All_seed_mean": float(mean[j]),
            "All_seed_CI_lower": float(ci_lower[j]),
            "All_seed_CI_upper": float(ci_upper[j]),
            "Preselected_seed": preselected_seed,
            "Preselected_seed_mean": selected_mean,
            "Preselected_seed_CI_half": selected_half,
            "Preselected_seed_mean_inside_all_seed_CI": inside,
            "Pooled_within_realization_seed_SD": pooled_within_sd,
            "Between_realization_SD": between_sd,
            "Seed_variation_smaller_than_realization_variation": variance_pass,
            "Audit_status": "PASS" if inside and variance_pass else "FAIL",
        })
    return pd.DataFrame(rows)


def build_special_checks(
    parsed: ParsedData, mean: np.ndarray, ci_lower: np.ndarray, ci_upper: np.ndarray,
    delta_reference: float,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    kind = parsed.x_kind
    if kind == "width" and parsed.x_numeric is not None and len(mean) >= 2:
        for i in range(len(mean) - 1):
            change = 100.0 * abs(mean[i + 1] - mean[i]) / abs(mean[i])
            overlap = confidence_intervals_overlap(ci_lower[i], ci_upper[i], ci_lower[i + 1], ci_upper[i + 1])
            rows.append({
                "Check": "Width convergence",
                "Comparison": f"{parsed.x_values[i]:g} to {parsed.x_values[i+1]:g}",
                "Percent_change": change,
                "CI_overlap": overlap,
                "Status": "PASS" if change < 5.0 and overlap else "FAIL",
                "Interpretation": "Pass requires <5% mean change and overlapping 95% CIs.",
            })
    elif kind == "delta" and parsed.x_numeric is not None:
        index = int(np.argmin(np.abs(parsed.x_numeric - delta_reference)))
        if not np.isclose(parsed.x_numeric[index], delta_reference, atol=1e-12):
            rows.append({
                "Check": "Delta reference", "Comparison": delta_reference,
                "Percent_change": np.nan, "CI_overlap": np.nan, "Status": "REVIEW",
                "Interpretation": f"Requested reference Delta={delta_reference:g} was not present.",
            })
        else:
            for i, x in enumerate(parsed.x_numeric):
                if i == index:
                    continue
                change = 100.0 * abs(mean[i] - mean[index]) / abs(mean[index])
                overlap = confidence_intervals_overlap(ci_lower[i], ci_upper[i], ci_lower[index], ci_upper[index])
                rows.append({
                    "Check": "Delta sensitivity",
                    "Comparison": f"Delta {x:g} vs {delta_reference:g}",
                    "Percent_change": change,
                    "CI_overlap": overlap,
                    "Status": "PASS" if change < 5.0 and overlap else "FAIL",
                    "Interpretation": (
                        "Statistical criterion only. Interior-profile linearity and cumulative-energy "
                        "linearity require manual checks from NEMD outputs."
                    ),
                })
    elif kind == "direction" and len(mean) == 2:
        difference = 100.0 * abs(mean[1] - mean[0]) / ((abs(mean[0]) + abs(mean[1])) / 2.0)
        overlap = confidence_intervals_overlap(ci_lower[0], ci_upper[0], ci_lower[1], ci_upper[1])
        rows.append({
            "Check": "Directional isotropy",
            "Comparison": f"{parsed.x_values[0]} vs {parsed.x_values[1]}",
            "Percent_change": difference,
            "CI_overlap": overlap,
            "Status": "PASS" if difference < 5.0 and overlap else "FAIL",
            "Interpretation": "Pass supports effective isotropy within statistical uncertainty.",
        })
    return pd.DataFrame(rows)


def determine_analysis_type(parsed: ParsedData, requested: str) -> tuple[str, str]:
    kind = parsed.x_kind
    if requested != "auto":
        if requested == "length":
            kind = "length"
        elif requested in {"screening", "width", "delta", "direction"}:
            kind = requested
    n = len(parsed.x_values)
    if kind == "length":
        if parsed.x_numeric is None:
            raise ValueError("Length scaling requires numeric lengths.")
        if n == 1:
            return "finite_single", kind
        if n == 2:
            raise ValueError("Two lengths cannot support a defensible extrapolation. Use one or >=3 lengths.")
        return "length_scaling", kind
    return ("finite_single" if n == 1 else "series"), kind


def build_diagnostics(
    analysis_type: str, fit: FitResult | None, shortest_fit: FitResult | None,
    fit_bootstrap_summary: pd.DataFrame, physical: np.ndarray | None,
    leave_one_out: pd.DataFrame, n_points: int, ci_lower: np.ndarray,
    special_checks: pd.DataFrame, seed_audit: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    def add(check: str, value: object, status: str, interpretation: str) -> None:
        rows.append({"Check": check, "Value": value, "Status": status, "Interpretation": interpretation})

    if np.any(ci_lower <= 0):
        add("Point confidence intervals", "At least one lower CI <= 0", "CAUTION",
            "Reciprocal error bars are undefined or nonphysical when a kappa CI crosses zero.")
    else:
        add("Point confidence intervals", "All lower limits > 0", "PASS",
            "All reciprocal confidence-limit transformations are positive.")

    if analysis_type != "length_scaling":
        add("Infinite-length extrapolation", "Not calculated", "NOT APPLICABLE",
            "kappa_inf and lambda_eff are calculated only from genuine multi-length data.")
        if not special_checks.empty:
            fails = int((special_checks["Status"] == "FAIL").sum())
            add("Special comparison checks", f"{fails} failed", "PASS" if fails == 0 else "REVIEW",
                "See the Special_Checks worksheet.")
        if not seed_audit.empty:
            fails = int((seed_audit["Audit_status"] == "FAIL").sum())
            add("Seed-noise audit", f"{fails} failed", "PASS" if fails == 0 else "REVIEW",
                "See the Seed_Noise_Audit worksheet.")
        return pd.DataFrame(rows)

    assert fit is not None
    add("Positive intercept", fit.intercept, "PASS" if fit.intercept > 0 else "FAIL",
        "A positive intercept is required for positive kappa_inf.")
    add("Positive slope", fit.slope, "PASS" if fit.slope > 0 else "FAIL",
        "A positive slope is expected for the adopted suppression model.")
    r2_status = "PASS" if fit.r2 >= 0.95 else ("REVIEW" if fit.r2 >= 0.90 else "CAUTION")
    add("Reciprocal-fit R^2", fit.r2, r2_status,
        "R^2 is descriptive; residual and sensitivity checks remain necessary.")
    max_std_resid = float(np.nanmax(np.abs(fit.standardized_residuals)))
    add("Maximum absolute standardized residual", max_std_resid,
        "PASS" if max_std_resid <= 2.0 else "REVIEW",
        "A value above about 2 suggests an influential or poorly fitted length.")
    if physical is not None:
        rate = 100.0 * float(np.mean(physical))
        status = "PASS" if rate >= 99.0 else ("REVIEW" if rate >= 95.0 else "FAIL")
        add("Physical bootstrap fits", f"{rate:.3f}%", status,
            "Nonphysical fits are counted and retained as an instability diagnostic.")
    if shortest_fit is not None and not fit_bootstrap_summary.empty:
        krow = fit_bootstrap_summary.loc[fit_bootstrap_summary["Quantity"] == "kappa_inf"].iloc[0]
        mrow = fit_bootstrap_summary.loc[fit_bootstrap_summary["Quantity"] == "lambda_eff"].iloc[0]
        kin = float(krow["Lower_95"]) <= shortest_fit.kappa_inf <= float(krow["Upper_95"])
        minside = float(mrow["Lower_95"]) <= shortest_fit.lambda_eff <= float(mrow["Upper_95"])
        add("Shortest-length exclusion: kappa_inf", shortest_fit.kappa_inf,
            "PASS" if kin else "REVIEW", "Checks sensitivity to the shortest bar.")
        add("Shortest-length exclusion: lambda_eff", shortest_fit.lambda_eff,
            "PASS" if minside else "REVIEW", "Checks MFP sensitivity to the shortest bar.")
    if not leave_one_out.empty:
        max_k = float(leave_one_out["Kappa_inf_change_pct"].abs().max())
        max_l = float(leave_one_out["Lambda_eff_change_pct"].abs().max())
        influential = float(leave_one_out.loc[
            leave_one_out["Kappa_inf_change_pct"].abs().idxmax(), "Excluded_Length"
        ])
        add("Maximum leave-one-out kappa_inf change", f"{max_k:.3f}% (exclude {influential:g})",
            "PASS" if max_k <= 5 else ("REVIEW" if max_k <= 10 else "CAUTION"),
            "Large changes indicate that one length controls the intercept.")
        add("Maximum leave-one-out lambda_eff change", f"{max_l:.3f}%",
            "PASS" if max_l <= 10 else ("REVIEW" if max_l <= 20 else "CAUTION"),
            "lambda_eff is usually more sensitive than kappa_inf.")
    add("Number of simulated lengths", n_points, "PASS" if n_points >= 5 else "REVIEW",
        "Five or more lengths generally provide a more defensible extrapolation.")
    if not seed_audit.empty:
        fails = int((seed_audit["Audit_status"] == "FAIL").sum())
        add("Seed-noise audit", f"{fails} failed", "PASS" if fails == 0 else "REVIEW",
            "See the Seed_Noise_Audit worksheet.")
    return pd.DataFrame(rows)


def analyze_data(
    parsed: ParsedData, case_name: str, iterations: int, random_seed: int,
    paired_realizations: bool, paired_units: bool, requested_kind: str,
    delta_reference: float, thickness_nm: float | None,
    dopant_count: int | None, host_sites: int | None, preselected_seed: int,
) -> Analysis:
    if iterations < 100:
        raise ValueError("Use at least 100 bootstrap iterations; 10,000 is recommended.")
    if thickness_nm is not None and thickness_nm <= 0:
        raise ValueError("Thickness must be positive.")
    if (dopant_count is None) ^ (host_sites is None):
        raise ValueError("Provide both --dopant-count and --host-sites, or neither.")
    exact_concentration = None
    if dopant_count is not None and host_sites is not None:
        if dopant_count < 0 or host_sites <= 0 or dopant_count > host_sites:
            raise ValueError("Invalid dopant and host-site counts.")
        exact_concentration = 100.0 * dopant_count / host_sites

    analysis_type, effective_kind = determine_analysis_type(parsed, requested_kind)
    if effective_kind != parsed.x_kind:
        parsed = ParsedData(**{**parsed.__dict__, "x_kind": effective_kind})
    rng = np.random.default_rng(random_seed)

    realization_means: list[np.ndarray] | None = None
    units: list[np.ndarray] = []
    unit_names: list[list[str]] = []
    if parsed.mode == "simple":
        assert parsed.simple_values is not None
        units = [np.asarray(v, dtype=float) for v in parsed.simple_values]
        names = parsed.simple_names or [f"Unit {i+1}" for i in range(len(units[0]))]
        unit_names = [[str(v) for v in names[:len(row)]] for row in units]
    else:
        assert parsed.hierarchical_values is not None
        assert parsed.realization_ids is not None
        realization_means = [np.mean(matrix, axis=1) for matrix in parsed.hierarchical_values]
        units = [np.asarray(v, dtype=float) for v in realization_means]
        unit_names = [[str(v) for v in parsed.realization_ids] for _ in units]

    n, mean, sd, se, tcrit, half, lower, upper = t_statistics(units)

    point_samples: list[np.ndarray] = []
    point_rows: list[dict[str, object]] = []
    for i in range(len(parsed.x_values)):
        if parsed.mode == "simple":
            sample = bootstrap_point_simple(units[i], iterations, rng)
        else:
            assert parsed.hierarchical_values is not None
            sample = bootstrap_point_hierarchical(parsed.hierarchical_values[i], iterations, rng)
        point_samples.append(sample)
        row = percentile_row(sample, "kappa_point", "W m^-1 K^-1")
        row[parsed.x_header] = parsed.x_values[i]
        point_rows.append(row)
    point_bootstrap_summary = pd.DataFrame(point_rows)

    seed_audit = build_seed_audit(parsed, mean, lower, upper, preselected_seed)
    special_checks = build_special_checks(parsed, mean, lower, upper, delta_reference)

    fit = None
    reciprocal_x = reciprocal_y = reciprocal_neg = reciprocal_pos = None
    bi = bs = bk = bl = bp = None
    fit_bootstrap_summary = pd.DataFrame()
    shortest_fit = None
    leave_one_out = pd.DataFrame()

    if analysis_type == "length_scaling":
        assert parsed.x_numeric is not None
        reciprocal_x = 1.0 / parsed.x_numeric
        reciprocal_y = 1.0 / mean
        reciprocal_neg = np.full_like(reciprocal_y, np.nan)
        reciprocal_pos = np.full_like(reciprocal_y, np.nan)
        positive = lower > 0
        reciprocal_low = np.full_like(reciprocal_y, np.nan)
        reciprocal_high = np.full_like(reciprocal_y, np.nan)
        reciprocal_low[positive] = 1.0 / upper[positive]
        reciprocal_high[positive] = 1.0 / lower[positive]
        reciprocal_neg[positive] = reciprocal_y[positive] - reciprocal_low[positive]
        reciprocal_pos[positive] = reciprocal_high[positive] - reciprocal_y[positive]
        fit = ordinary_fit(reciprocal_x, reciprocal_y)
        if len(parsed.x_values) >= 4:
            shortest_fit = ordinary_fit(reciprocal_x[1:], reciprocal_y[1:])
            loo_rows: list[dict[str, float]] = []
            for i, excluded in enumerate(parsed.x_numeric):
                mask = np.ones(len(parsed.x_numeric), dtype=bool)
                mask[i] = False
                subfit = ordinary_fit(reciprocal_x[mask], reciprocal_y[mask])
                loo_rows.append({
                    "Excluded_Length": float(excluded), "Intercept_a": subfit.intercept,
                    "Slope_b": subfit.slope, "R2": subfit.r2,
                    "Kappa_inf": subfit.kappa_inf, "Lambda_eff_nm": subfit.lambda_eff,
                    "Kappa_inf_change_pct": 100.0 * (subfit.kappa_inf - fit.kappa_inf) / fit.kappa_inf,
                    "Lambda_eff_change_pct": 100.0 * (subfit.lambda_eff - fit.lambda_eff) / fit.lambda_eff,
                })
            leave_one_out = pd.DataFrame(loo_rows)
        if parsed.mode == "simple":
            bi, bs, bk, bl, bp = bootstrap_scaling_simple(
                units, reciprocal_x, iterations, rng, paired_units=paired_units
            )
        else:
            assert parsed.hierarchical_values is not None
            bi, bs, bk, bl, bp = bootstrap_scaling_hierarchical(
                parsed.hierarchical_values, reciprocal_x, iterations, rng,
                paired_realizations=paired_realizations,
            )
        fit_bootstrap_summary = pd.DataFrame([
            percentile_row(bk, "kappa_inf", "W m^-1 K^-1"),
            percentile_row(bl, "lambda_eff", "nm"),
            percentile_row(bi, "intercept_a", "m K W^-1"),
            percentile_row(bs, "slope_b", "nm m K W^-1"),
        ])

    diagnostics = build_diagnostics(
        analysis_type, fit, shortest_fit, fit_bootstrap_summary, bp,
        leave_one_out, len(parsed.x_values), lower, special_checks, seed_audit,
    )

    return Analysis(
        case_name=case_name, analysis_type=analysis_type, x_kind=effective_kind,
        data_mode=parsed.mode, x_header=parsed.x_header, x_values=parsed.x_values,
        x_numeric=parsed.x_numeric, unit_names=unit_names, units=units,
        realization_means=realization_means, n_units=n, mean=mean, sd=sd, se=se,
        t_critical=tcrit, ci_half=half, ci_lower=lower, ci_upper=upper,
        point_bootstrap_summary=point_bootstrap_summary,
        point_bootstrap_samples=point_samples, fit=fit,
        reciprocal_x=reciprocal_x, reciprocal_y=reciprocal_y,
        reciprocal_neg_error=reciprocal_neg, reciprocal_pos_error=reciprocal_pos,
        bootstrap_intercept=bi, bootstrap_slope=bs, bootstrap_kappa_inf=bk,
        bootstrap_lambda_eff=bl, bootstrap_physical=bp,
        fit_bootstrap_summary=fit_bootstrap_summary, shortest_fit=shortest_fit,
        leave_one_out=leave_one_out, seed_audit=seed_audit,
        special_checks=special_checks, diagnostics=diagnostics,
        thickness_nm=thickness_nm, exact_concentration_pct=exact_concentration,
    )


def safe_x_for_plot(analysis: Analysis) -> tuple[np.ndarray, list[str] | None]:
    if analysis.x_numeric is not None:
        return analysis.x_numeric, None
    return np.arange(len(analysis.x_values), dtype=float), [str(v) for v in analysis.x_values]


def create_figures(analysis: Analysis, figures_dir: Path, dpi: int) -> dict[str, Path]:
    figures_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    x, labels = safe_x_for_plot(analysis)

    fig, ax = plt.subplots(figsize=(6.4, 4.5))
    ax.errorbar(x, analysis.mean, yerr=analysis.ci_half, fmt="o-" if len(x) > 1 else "o", capsize=4)
    if labels is not None:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_xlabel(analysis.x_header)
    ax.set_ylabel(r"Thermal conductivity, $\kappa$ (W m$^{-1}$ K$^{-1}$)")
    ax.set_title(analysis.case_name)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    point_path = figures_dir / "kappa_with_95CI.png"
    fig.savefig(point_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    paths["kappa_with_95CI"] = point_path

    if analysis.analysis_type == "length_scaling":
        assert analysis.fit is not None
        assert analysis.x_numeric is not None
        assert analysis.reciprocal_x is not None and analysis.reciprocal_y is not None
        fig, ax = plt.subplots(figsize=(6.4, 4.5))
        ax.errorbar(analysis.x_numeric, analysis.mean, yerr=analysis.ci_half, fmt="o", capsize=4,
                    label="Mean and 95% CI")
        dense_l = np.linspace(np.min(analysis.x_numeric), np.max(analysis.x_numeric), 500)
        predicted = 1.0 / (analysis.fit.intercept + analysis.fit.slope / dense_l)
        ax.plot(dense_l, predicted, label="Finite-size fit")
        ax.set_xlabel("Transport length, L (nm)")
        ax.set_ylabel(r"Thermal conductivity, $\kappa(L)$ (W m$^{-1}$ K$^{-1}$)")
        ax.set_title(analysis.case_name)
        ax.legend()
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        path = figures_dir / "kappa_vs_length_fit.png"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        paths["kappa_vs_length_fit"] = path

        fig, ax = plt.subplots(figsize=(6.4, 4.5))
        valid_errors = np.isfinite(analysis.reciprocal_neg_error).all() and np.isfinite(analysis.reciprocal_pos_error).all()
        yerr = None
        if valid_errors:
            yerr = np.vstack([analysis.reciprocal_neg_error, analysis.reciprocal_pos_error])
        ax.errorbar(analysis.reciprocal_x, analysis.reciprocal_y, yerr=yerr, fmt="o", capsize=4,
                    label="Transformed means" + (" and 95% CI" if valid_errors else ""))
        dense_x = np.linspace(0, np.max(analysis.reciprocal_x) * 1.05, 500)
        ax.plot(dense_x, analysis.fit.intercept + analysis.fit.slope * dense_x, label="Linear fit")
        equation = f"1/kappa = {analysis.fit.intercept:.6g} + {analysis.fit.slope:.6g}(1/L)\nR² = {analysis.fit.r2:.4f}"
        ax.text(0.04, 0.96, equation, transform=ax.transAxes, va="top")
        ax.set_xlabel(r"$1/L$ (nm$^{-1}$)")
        ax.set_ylabel(r"$1/\kappa$ (m K W$^{-1}$)")
        ax.legend()
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        path = figures_dir / "reciprocal_fit.png"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        paths["reciprocal_fit"] = path

        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        ax.axhline(0, linewidth=1)
        ax.scatter(analysis.reciprocal_x, analysis.fit.residuals)
        ax.set_xlabel(r"$1/L$ (nm$^{-1}$)")
        ax.set_ylabel("Residual in 1/kappa")
        ax.set_title("Reciprocal-fit residuals")
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        path = figures_dir / "residuals.png"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        paths["residuals"] = path

        for values, filename, xlabel, title in [
            (analysis.bootstrap_kappa_inf, "bootstrap_kappa_inf.png",
             r"$\kappa_\infty$ (W m$^{-1}$ K$^{-1}$)", "Bootstrap infinite-length conductivity"),
            (analysis.bootstrap_lambda_eff, "bootstrap_lambda_eff.png",
             r"$\lambda_{eff}$ (nm)", "Bootstrap effective characteristic MFP"),
        ]:
            assert values is not None
            finite = values[np.isfinite(values)]
            fig, ax = plt.subplots(figsize=(6.4, 4.4))
            ax.hist(finite, bins=50)
            lower, median, upper = np.percentile(finite, [2.5, 50, 97.5])
            ax.axvline(median, linestyle="--", label="Median")
            ax.axvline(lower, linestyle=":", label="95% limits")
            ax.axvline(upper, linestyle=":")
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Count")
            ax.set_title(title)
            ax.legend()
            fig.tight_layout()
            path = figures_dir / filename
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
            plt.close(fig)
            paths[filename.removesuffix(".png")] = path

        if not analysis.leave_one_out.empty:
            fig, ax = plt.subplots(figsize=(6.4, 4.4))
            ax.plot(analysis.leave_one_out["Excluded_Length"], analysis.leave_one_out["Kappa_inf"], marker="o")
            ax.axhline(analysis.fit.kappa_inf, linestyle="--", label="All-length fit")
            ax.set_xlabel("Excluded length (nm)")
            ax.set_ylabel(r"Refitted $\kappa_\infty$ (W m$^{-1}$ K$^{-1}$)")
            ax.set_title("Leave-one-length-out sensitivity")
            ax.legend()
            ax.grid(True, alpha=0.25)
            fig.tight_layout()
            path = figures_dir / "leave_one_out_kappa_inf.png"
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
            plt.close(fig)
            paths["leave_one_out_kappa_inf"] = path
    return paths


HEADER_FILL = PatternFill("solid", fgColor="1F4E78")
HEADER_FONT = Font(color="FFFFFF", bold=True)
PASS_FILL = PatternFill("solid", fgColor="D9EAD3")
CAUTION_FILL = PatternFill("solid", fgColor="FCE4D6")
ERROR_FILL = PatternFill("solid", fgColor="F4CCCC")
THIN_GRAY = Side(style="thin", color="B7B7B7")


def style_header(ws, row: int = 1) -> None:
    for cell in ws[row]:
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = Border(bottom=THIN_GRAY)


def auto_width(ws, maximum: int = 55) -> None:
    for column_cells in ws.columns:
        letter = column_cells[0].column_letter
        width = max((len("" if c.value is None else str(c.value)) for c in column_cells), default=10)
        ws.column_dimensions[letter].width = min(max(width + 2, 10), maximum)


def write_dataframe(ws, frame: pd.DataFrame, start_row: int = 1, start_col: int = 1) -> None:
    for j, column in enumerate(frame.columns, start_col):
        ws.cell(start_row, j, column)
    for i, row in enumerate(frame.itertuples(index=False, name=None), start_row + 1):
        for j, value in enumerate(row, start_col):
            if isinstance(value, np.generic):
                value = value.item()
            ws.cell(i, j, None if pd.isna(value) else value)


def reporting_lines(analysis: Analysis) -> list[tuple[str, str]]:
    lines: list[tuple[str, str]] = []
    if analysis.analysis_type == "length_scaling":
        krow = analysis.fit_bootstrap_summary.loc[
            analysis.fit_bootstrap_summary["Quantity"] == "kappa_inf"
        ].iloc[0]
        mrow = analysis.fit_bootstrap_summary.loc[
            analysis.fit_bootstrap_summary["Quantity"] == "lambda_eff"
        ].iloc[0]
        lines.append((
            "kappa_inf",
            f"{float(krow['Median']):.6g} (+{float(krow['Upper_Error']):.6g}/-"
            f"{float(krow['Lower_Error']):.6g}) W m^-1 K^-1, bootstrap 95% CI",
        ))
        lines.append((
            "lambda_eff",
            f"{float(mrow['Median']):.6g} (+{float(mrow['Upper_Error']):.6g}/-"
            f"{float(mrow['Lower_Error']):.6g}) nm, bootstrap 95% CI",
        ))
        if analysis.bootstrap_physical is not None:
            physical_rate = 100.0 * float(np.mean(analysis.bootstrap_physical))
            if physical_rate < 95.0 or analysis.fit is None or analysis.fit.intercept <= 0 or analysis.fit.slope <= 0:
                lines.append((
                    "Reporting warning",
                    f"Extrapolation is unstable ({physical_rate:.3f}% physical bootstrap fits). "
                    "Do not report kappa_inf or lambda_eff without resolving the diagnostics.",
                ))
        if analysis.thickness_nm is not None:
            t_m = analysis.thickness_nm * 1e-9
            median_g = float(krow["Median"]) * t_m
            lower_g = float(krow["Lower_95"]) * t_m
            upper_g = float(krow["Upper_95"]) * t_m
            lines.append((
                "G_sheet_inf",
                f"{median_g:.6g} (+{upper_g-median_g:.6g}/-{median_g-lower_g:.6g}) "
                "W K^-1 m^-1, bootstrap 95% CI",
            ))
    else:
        for i, x in enumerate(analysis.x_values):
            lines.append((
                f"kappa at {analysis.x_header}={x}",
                f"{analysis.mean[i]:.6g} +/- {analysis.ci_half[i]:.6g} W m^-1 K^-1, "
                f"two-sided 95% Student-t CI (n={analysis.n_units[i]})",
            ))
        lines.append(("kappa_inf", "Not calculated: no genuine multi-length extrapolation was performed."))
        lines.append(("lambda_eff", "Not calculated: effective MFP requires genuine multi-length data."))
    if analysis.exact_concentration_pct is not None:
        lines.append(("Exact realized concentration", f"{analysis.exact_concentration_pct:.6g}%"))
    return lines


def create_workbook(
    analysis: Analysis, parsed: ParsedData, source_path: Path, output_path: Path,
    figures: dict[str, Path], iterations: int, random_seed: int,
    paired_realizations: bool, paired_units: bool, include_bootstrap_samples: bool,
    dpi: int,
) -> None:
    wb = Workbook()
    wb.remove(wb.active)
    wb.calculation.fullCalcOnLoad = True
    wb.calculation.forceFullCalc = True

    ws_info = wb.create_sheet("README")
    info_rows = [
        ["Field", "Value"],
        ["Source file", str(source_path)],
        ["Source sheet", parsed.selected_sheet],
        ["Case", analysis.case_name],
        ["Detected input mode", f"{parsed.mode} ({parsed.layout})"],
        ["Detected x kind", analysis.x_kind],
        ["Analysis type", analysis.analysis_type],
        ["Bootstrap iterations", iterations],
        ["Bootstrap random seed", random_seed],
        ["Figure DPI", dpi],
        ["Paired realizations", paired_realizations],
        ["Paired simple units", paired_units],
        ["MFP definition", "lambda_eff=b/a from 1/kappa=a+b/L; not a spectral MFP distribution"],
    ]
    if analysis.thickness_nm is not None:
        info_rows.append(["Effective thickness (nm)", analysis.thickness_nm])
    if analysis.exact_concentration_pct is not None:
        info_rows.append(["Exact realized dopant concentration (%)", analysis.exact_concentration_pct])
    for row in info_rows:
        ws_info.append(row)
    style_header(ws_info)

    ws_raw = wb.create_sheet("Raw_Data")
    write_dataframe(ws_raw, parsed.raw_dataframe)
    style_header(ws_raw)

    if analysis.realization_means is not None:
        ws_rm = wb.create_sheet("Realization_Means")
        headers = [analysis.x_header] + [f"Realization_mean_{i+1}" for i in range(len(analysis.realization_means[0]))]
        ws_rm.append(headers)
        for x, values in zip(analysis.x_values, analysis.realization_means):
            ws_rm.append([x] + [float(v) for v in values])
        style_header(ws_rm)

    ws_stats = wb.create_sheet("Point_Statistics")
    headers = [
        analysis.x_header, "N_independent_units", "Mean_kappa", "Sample_SD", "SE",
        "t_critical_95", "CI_half_95", "CI_lower_95", "CI_upper_95",
    ]
    if analysis.thickness_nm is not None:
        headers += ["Mean_G_sheet_W_K-1_m-1", "G_sheet_CI_half"]
    ws_stats.append(headers)
    for i, x in enumerate(analysis.x_values):
        row = [
            x, int(analysis.n_units[i]), float(analysis.mean[i]), float(analysis.sd[i]),
            float(analysis.se[i]), float(analysis.t_critical[i]), float(analysis.ci_half[i]),
            float(analysis.ci_lower[i]), float(analysis.ci_upper[i]),
        ]
        if analysis.thickness_nm is not None:
            factor = analysis.thickness_nm * 1e-9
            row += [float(analysis.mean[i] * factor), float(analysis.ci_half[i] * factor)]
        ws_stats.append(row)
    style_header(ws_stats)

    ws_point_bs = wb.create_sheet("Point_Bootstrap_Summary")
    write_dataframe(ws_point_bs, analysis.point_bootstrap_summary)
    style_header(ws_point_bs)

    ws_fit = wb.create_sheet("Fit_Summary")
    if analysis.fit is None:
        ws_fit.append(["Quantity", "Value", "Interpretation"])
        ws_fit.append(["kappa_inf", "Not calculated", "Only genuine multi-length data support extrapolation."])
        ws_fit.append(["lambda_eff", "Not calculated", "Only genuine multi-length data support MFP extraction."])
    else:
        f = analysis.fit
        rows = [
            ["Quantity", "Value", "Units/notes"],
            ["Fit equation", "1/kappa = a + b/L", "Unweighted ordinary least squares"],
            ["Intercept_a", f.intercept, "m K W^-1"],
            ["Slope_b", f.slope, "nm m K W^-1"],
            ["R2", f.r2, ""],
            ["Adjusted_R2", f.adjusted_r2, ""],
            ["RMSE", f.rmse, "reciprocal-kappa units"],
            ["Direct_kappa_inf", f.kappa_inf, "W m^-1 K^-1"],
            ["Direct_lambda_eff", f.lambda_eff, "nm"],
        ]
        if analysis.shortest_fit is not None:
            rows += [
                ["Without_shortest_kappa_inf", analysis.shortest_fit.kappa_inf, "Sensitivity check"],
                ["Without_shortest_lambda_eff", analysis.shortest_fit.lambda_eff, "Sensitivity check"],
            ]
        for row in rows:
            ws_fit.append(row)
    style_header(ws_fit)

    ws_fit_bs = wb.create_sheet("Fit_Bootstrap_Summary")
    if analysis.fit_bootstrap_summary.empty:
        ws_fit_bs.append(["Result", "Not applicable"])
    else:
        write_dataframe(ws_fit_bs, analysis.fit_bootstrap_summary)
    style_header(ws_fit_bs)

    ws_report = wb.create_sheet("Recommended_Reporting")
    ws_report.append(["Quantity", "Recommended reporting"])
    for quantity, line in reporting_lines(analysis):
        ws_report.append([quantity, line])
    style_header(ws_report)

    ws_sens = wb.create_sheet("Sensitivity")
    if analysis.leave_one_out.empty:
        ws_sens.append(["Check", "Not applicable"])
    else:
        write_dataframe(ws_sens, analysis.leave_one_out)
    style_header(ws_sens)

    ws_audit = wb.create_sheet("Seed_Noise_Audit")
    if analysis.seed_audit.empty:
        ws_audit.append(["Result", "Not applicable: hierarchical data with at least two seeds per realization were not supplied."])
    else:
        write_dataframe(ws_audit, analysis.seed_audit)
    style_header(ws_audit)

    ws_special = wb.create_sheet("Special_Checks")
    if analysis.special_checks.empty:
        ws_special.append(["Result", "No width, Delta, or two-direction automatic check was applicable."])
    else:
        write_dataframe(ws_special, analysis.special_checks)
    style_header(ws_special)

    ws_diag = wb.create_sheet("Diagnostics")
    write_dataframe(ws_diag, analysis.diagnostics)
    style_header(ws_diag)
    if "Status" in analysis.diagnostics.columns:
        status_col = list(analysis.diagnostics.columns).index("Status") + 1
        for row in range(2, ws_diag.max_row + 1):
            cell = ws_diag.cell(row, status_col)
            if cell.value == "PASS":
                cell.fill = PASS_FILL
            elif cell.value in {"REVIEW", "CAUTION", "NOT APPLICABLE"}:
                cell.fill = CAUTION_FILL
            elif cell.value == "FAIL":
                cell.fill = ERROR_FILL

    ws_origin = wb.create_sheet("OriginPro_Data")
    if analysis.analysis_type == "length_scaling":
        ws_origin.append([
            "Length_nm", "Mean_kappa", "CI_half_95", "CI_lower_95", "CI_upper_95",
            "1_over_L", "1_over_kappa", "Reciprocal_negative_error",
            "Reciprocal_positive_error", "Fitted_1_over_kappa", "Residual", "Predicted_kappa",
        ])
        assert analysis.fit is not None
        for i, x in enumerate(analysis.x_values):
            ws_origin.append([
                x, float(analysis.mean[i]), float(analysis.ci_half[i]),
                float(analysis.ci_lower[i]), float(analysis.ci_upper[i]),
                float(analysis.reciprocal_x[i]), float(analysis.reciprocal_y[i]),
                None if not np.isfinite(analysis.reciprocal_neg_error[i]) else float(analysis.reciprocal_neg_error[i]),
                None if not np.isfinite(analysis.reciprocal_pos_error[i]) else float(analysis.reciprocal_pos_error[i]),
                float(analysis.fit.fitted_y[i]), float(analysis.fit.residuals[i]),
                float(1.0 / analysis.fit.fitted_y[i]),
            ])
    else:
        ws_origin.append([analysis.x_header, "Mean_kappa", "CI_half_95", "CI_lower_95", "CI_upper_95"])
        for i, x in enumerate(analysis.x_values):
            ws_origin.append([x, float(analysis.mean[i]), float(analysis.ci_half[i]),
                              float(analysis.ci_lower[i]), float(analysis.ci_upper[i])])
    style_header(ws_origin)

    if include_bootstrap_samples:
        ws_samples = wb.create_sheet("Bootstrap_Samples")
        if analysis.analysis_type == "length_scaling":
            ws_samples.append(["Iteration", "Intercept_a", "Slope_b", "Kappa_inf", "Lambda_eff_nm", "Physical_fit"])
            assert analysis.bootstrap_intercept is not None
            for i in range(len(analysis.bootstrap_intercept)):
                ws_samples.append([
                    i + 1, float(analysis.bootstrap_intercept[i]), float(analysis.bootstrap_slope[i]),
                    float(analysis.bootstrap_kappa_inf[i]), float(analysis.bootstrap_lambda_eff[i]),
                    bool(analysis.bootstrap_physical[i]),
                ])
        else:
            headers = ["Iteration"] + [f"Bootstrap_mean_{i+1}" for i in range(len(analysis.x_values))]
            ws_samples.append(headers)
            for b in range(iterations):
                ws_samples.append([b + 1] + [float(samples[b]) for samples in analysis.point_bootstrap_samples])
        style_header(ws_samples)
        ws_samples.freeze_panes = "A2"

    ws_fig = wb.create_sheet("Figures")
    ws_fig.sheet_view.showGridLines = False
    anchor_rows = [1, 28, 55, 82, 109, 136, 163]
    for anchor_row, (name, path) in zip(anchor_rows, figures.items()):
        ws_fig.cell(anchor_row, 1, name).font = Font(bold=True)
        image = XLImage(str(path))
        image.width = 640
        image.height = 450
        ws_fig.add_image(image, f"A{anchor_row + 1}")

    for ws in wb.worksheets:
        ws.sheet_view.showGridLines = False
        for row in ws.iter_rows():
            for cell in row:
                cell.alignment = Alignment(vertical="center", wrap_text=True)
        auto_width(ws)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(output_path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Unified thermal-conductivity statistics, extrapolation, MFP, and robustness analysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input", type=Path, help="Input .xlsx, .xlsm, or .csv file")
    parser.add_argument("--sheet", help="Worksheet name; defaults to first sheet")
    parser.add_argument("--header-row", type=int, default=1, help="1-based header row")
    parser.add_argument("--x-column", help="Exact independent-variable column header")
    parser.add_argument("--case-name", help="Case label; defaults to input filename stem")
    parser.add_argument(
        "--analysis-kind", choices=["auto", "length", "screening", "width", "delta", "direction"],
        default="auto", help="Override automatic interpretation of the x column",
    )
    parser.add_argument("--iterations", type=int, default=DEFAULT_ITERATIONS, help="Bootstrap iterations")
    parser.add_argument("--random-seed", type=int, default=DEFAULT_RANDOM_SEED, help="Bootstrap RNG seed")
    parser.add_argument("--paired-realizations", action="store_true",
                        help="Use matched realization IDs jointly across lengths in hierarchical bootstrap")
    parser.add_argument("--paired-units", action="store_true",
                        help="Use matched simple unit IDs jointly across lengths")
    parser.add_argument(
        "--preselected-seed", type=int, default=1,
        help="1-based seed used for the reduced one-seed dopant audit; default is S1",
    )
    parser.add_argument("--delta-reference", type=float, default=0.05, help="Reference Delta for sensitivity checks")
    parser.add_argument("--thickness-nm", type=float, help="Effective thickness for G_sheet=kappa*t")
    parser.add_argument("--dopant-count", type=int, help="Number of substituted host atoms")
    parser.add_argument("--host-sites", type=int, help="Number of eligible host sites before substitution")
    parser.add_argument("--output", type=Path, help="Output workbook")
    parser.add_argument("--figures-dir", type=Path, help="Output figure directory")
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI, help="PNG resolution")
    parser.add_argument("--omit-bootstrap-samples", action="store_true",
                        help="Do not include all bootstrap samples in the workbook")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.dpi <= 0:
            raise ValueError("Figure DPI must be positive.")
        if args.preselected_seed < 1:
            raise ValueError("--preselected-seed must be a positive 1-based index.")
        parsed = parse_input(args.input, args.sheet, args.header_row, args.x_column)
        case_name = args.case_name or args.input.stem
        output = args.output or args.input.with_name(f"{args.input.stem}_analyzed.xlsx")
        figures_dir = args.figures_dir or args.input.with_name(f"{args.input.stem}_analysis_figures")
        analysis = analyze_data(
            parsed, case_name, args.iterations, args.random_seed,
            args.paired_realizations, args.paired_units, args.analysis_kind,
            args.delta_reference, args.thickness_nm, args.dopant_count, args.host_sites,
            args.preselected_seed,
        )
        figures = create_figures(analysis, figures_dir, args.dpi)
        create_workbook(
            analysis, parsed, args.input, output, figures, args.iterations,
            args.random_seed, args.paired_realizations, args.paired_units,
            not args.omit_bootstrap_samples, args.dpi,
        )
        print("Analysis completed successfully.")
        print(f"Detected structure: {parsed.mode} ({parsed.layout})")
        print(f"Detected x kind: {analysis.x_kind}")
        print(f"Analysis type: {analysis.analysis_type}")
        print(f"Default/requested figure resolution: {args.dpi} dpi")
        print("Recommended reporting:")
        for quantity, line in reporting_lines(analysis):
            print(f"  {quantity}: {line}")
        print(f"Workbook: {output.resolve()}")
        print(f"Figures: {figures_dir.resolve()}")
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
