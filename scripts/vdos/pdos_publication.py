#!/usr/bin/env python3
"""
Publication-grade MD-derived VDOS/PDOS from a LAMMPS velocity dump.

Designed for monolayers and small equilibrium VDOS cells.

Input LAMMPS dump must contain at least:
    id type vx vy vz
Coordinates x y z are optional and ignored.

Recommended LAMMPS dump command:
    dump VEL all custom N pdos.lammpstrj id type vx vy vz
    dump_modify VEL sort id

What this script does
---------------------
1. Reads a LAMMPS dump with velocity columns.
2. Verifies constant sampled timestep spacing.
3. Removes center-of-mass velocity frame-by-frame using all atoms in the dump
   before storing sampled/selected velocities.
4. Computes a block-averaged velocity power spectrum using rFFT.
5. Supports optional mass weighting by applying sqrt(m_i) to velocities before FFT.
6. Supports atom-type projected spectra and allows repeated labels, e.g.
      --species_map "1:S,2:Mo,3:S"
   which automatically combines both S atom types into one S curve.
7. If atom sampling is used, corrects each sampled atom's spectral contribution by
   N_type_total / N_type_sampled to preserve expected species weights.
8. Normalizes the total integrated spectrum to 3 * primitive_atoms.

Important wording for papers
----------------------------
This is an MD-derived velocity power spectrum / VDOS / PDOS. It is not a
Phonopy calculation and should not be described as a harmonic eigenmode PDOS.

Output files
------------
    <prefix>_THz_3N_raw.dat
    <prefix>_THz_3N_smoothed.dat
    <prefix>_metadata.txt
    <prefix>_total_pub.png / .pdf
    <prefix>_components_pub.png / .pdf
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import OrderedDict, defaultdict
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

import matplotlib as mpl
mpl.use("Agg")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt


def _trapz(y: np.ndarray, x: np.ndarray) -> float:
    """Compatibility wrapper for NumPy trapezoidal integration."""
    fn = getattr(np, "trapezoid", None)
    if fn is None:
        fn = np.trapz
    return float(fn(y, x))


def fail(message: str) -> None:
    raise RuntimeError(message)


def parse_atoms_header(line: str) -> Dict[str, int]:
    parts = line.strip().split()
    if len(parts) < 3 or parts[0] != "ITEM:" or parts[1] != "ATOMS":
        raise ValueError(f"Malformed ITEM: ATOMS header: {line!r}")
    columns = parts[2:]
    required = ["id", "type", "vx", "vy", "vz"]
    missing = [c for c in required if c not in columns]
    if missing:
        raise ValueError(
            f"Dump is missing required columns {missing}. Need at least: id type vx vy vz. "
            f"Found columns: {columns}"
        )
    return {name: i for i, name in enumerate(columns)}


def read_first_snapshot(dump_path: str) -> Tuple[int, np.ndarray, np.ndarray, int]:
    """Return natoms, atom IDs, atom types, and first timestep from first frame."""
    with open(dump_path, "r") as f:
        line = f.readline()
        while line and not line.startswith("ITEM: TIMESTEP"):
            line = f.readline()
        if not line:
            fail("No ITEM: TIMESTEP frame found in dump file.")

        ts_line = f.readline()
        if not ts_line:
            fail("Unexpected EOF after ITEM: TIMESTEP.")
        timestep = int(ts_line.strip())

        if not f.readline().startswith("ITEM: NUMBER OF ATOMS"):
            fail("Malformed dump: expected ITEM: NUMBER OF ATOMS.")
        natoms = int(f.readline().strip())

        if not f.readline().startswith("ITEM: BOX BOUNDS"):
            fail("Malformed dump: expected ITEM: BOX BOUNDS.")
        f.readline(); f.readline(); f.readline()

        atoms_header = f.readline()
        cidx = parse_atoms_header(atoms_header)

        ids = np.empty(natoms, dtype=np.int64)
        types = np.empty(natoms, dtype=np.int32)
        for i in range(natoms):
            parts = f.readline().split()
            ids[i] = int(parts[cidx["id"]])
            types[i] = int(float(parts[cidx["type"]]))

    if len(np.unique(ids)) != natoms:
        fail("Atom IDs are not unique in first frame.")

    return natoms, ids, types, timestep


def read_frame_after_timestep(fh, natoms: int) -> Tuple[Optional[Dict[str, int]], Optional[List[str]]]:
    """
    fh is positioned immediately after the timestep value line.
    Return column index and atom lines, or (None, None) at EOF/malformed end.
    """
    line = fh.readline()
    if not line:
        return None, None
    if not line.startswith("ITEM: NUMBER OF ATOMS"):
        return None, None

    n_atoms_frame = int(fh.readline().strip())
    if n_atoms_frame != natoms:
        fail(f"Number of atoms changed between frames: first={natoms}, current={n_atoms_frame}.")

    if not fh.readline().startswith("ITEM: BOX BOUNDS"):
        return None, None
    fh.readline(); fh.readline(); fh.readline()

    atoms_header = fh.readline()
    if not atoms_header.startswith("ITEM: ATOMS"):
        return None, None
    cidx = parse_atoms_header(atoms_header)
    atom_lines = [fh.readline().rstrip("\n") for _ in range(natoms)]
    return cidx, atom_lines


def parse_species_map(text: str) -> "OrderedDict[int, str]":
    """
    Parse type-to-label map. Duplicate labels are allowed and are combined.
    Example: "1:S,2:Mo,3:S" -> {1:"S", 2:"Mo", 3:"S"}.
    """
    out: "OrderedDict[int, str]" = OrderedDict()
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        if ":" not in token:
            raise ValueError(f"Bad species_map token {token!r}. Expected format type:Label.")
        k, v = token.split(":", 1)
        atom_type = int(k.strip())
        label = v.strip()
        if not label:
            raise ValueError(f"Empty label in species_map token {token!r}.")
        out[atom_type] = label
    if not out:
        raise ValueError("species_map is empty.")
    return out


def parse_mass_map(text: Optional[str]) -> Dict[str, float]:
    if text is None:
        return {}
    out: Dict[str, float] = {}
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        if "=" not in token:
            raise ValueError(f"Bad mass_map token {token!r}. Expected Label=mass.")
        k, v = token.split("=", 1)
        label = k.strip()
        mass = float(v.strip())
        if mass <= 0:
            raise ValueError(f"Mass must be positive for {label}; got {mass}.")
        out[label] = mass
    return out


def unique_preserve_order(items: Iterable[str]) -> List[str]:
    seen = OrderedDict()
    for item in items:
        if item not in seen:
            seen[item] = None
    return list(seen.keys())


def gaussian_smooth(y: np.ndarray, sigma_bins: float) -> np.ndarray:
    """Gaussian smoothing with reflected edge padding; sigma is in frequency bins."""
    if sigma_bins <= 0:
        return np.array(y, copy=True)
    radius = int(max(2, math.ceil(4.0 * sigma_bins)))
    radius = min(radius, max(2, len(y) - 2)) if len(y) > 4 else 1
    x = np.arange(-radius, radius + 1, dtype=np.float64)
    kernel = np.exp(-0.5 * (x / sigma_bins) ** 2)
    kernel /= kernel.sum()
    if len(y) <= 2 * radius + 1:
        pad = radius
        ypad = np.pad(y, pad, mode="edge")
    else:
        left = y[1:radius + 1][::-1]
        right = y[-radius - 1:-1][::-1]
        ypad = np.concatenate([left, y, right])
    return np.convolve(ypad, kernel, mode="same")[radius:-radius]


def validate_args(args: argparse.Namespace) -> None:
    if args.timestep_fs <= 0:
        raise ValueError("--timestep_fs must be positive.")
    if args.dump_interval_steps <= 0:
        raise ValueError("--dump_interval_steps must be a positive integer.")
    if args.stride <= 0:
        raise ValueError("--stride must be a positive integer.")
    if args.block_len < 32:
        raise ValueError("--block_len is too small. Use at least 32; for publication use 4096 or larger.")
    if args.max_sample_per_type < 0:
        raise ValueError("--max_sample_per_type must be >= 0. Use 0 for all atoms.")
    if args.primitive_atoms <= 0:
        raise ValueError("--primitive_atoms must be positive.")
    if args.smooth_sigma_thz < 0:
        raise ValueError("--smooth_sigma_thz must be >= 0.")
    if args.xmax <= args.xmin:
        raise ValueError("--xmax must be larger than --xmin.")
    if args.min_partial_block_fraction <= 0 or args.min_partial_block_fraction > 1:
        raise ValueError("--min_partial_block_fraction must be in (0, 1].")


def build_selection(
    ids: np.ndarray,
    types: np.ndarray,
    present_types: List[int],
    max_sample_per_type: int,
    seed: int,
) -> Tuple[List[int], Dict[int, int], Dict[int, int], Dict[int, float]]:
    """Select atom IDs and return total/selected counts and per-type sampling correction factors."""
    rng = np.random.default_rng(seed)
    ids_list = [int(x) for x in ids.tolist()]
    types_list = [int(x) for x in types.tolist()]

    total_count_by_type: Dict[int, int] = {}
    selected_count_by_type: Dict[int, int] = {}
    selected_ids: List[int] = []

    for tval in present_types:
        pool = [aid for aid, typ in zip(ids_list, types_list) if typ == tval]
        total_count_by_type[tval] = len(pool)
        if max_sample_per_type > 0 and len(pool) > max_sample_per_type:
            pool = rng.choice(pool, size=max_sample_per_type, replace=False).astype(int).tolist()
        selected_count_by_type[tval] = len(pool)
        selected_ids.extend(pool)

    selected_ids = sorted(selected_ids)

    sample_factor_by_type: Dict[int, float] = {}
    for tval in present_types:
        n_total = total_count_by_type[tval]
        n_sel = selected_count_by_type[tval]
        if n_sel <= 0:
            fail(f"No atoms selected for atom type {tval}.")
        sample_factor_by_type[tval] = float(n_total) / float(n_sel)

    return selected_ids, total_count_by_type, selected_count_by_type, sample_factor_by_type


def process_block(
    data_block: np.ndarray,
    window: np.ndarray,
    sqrt_mass_dof: Optional[np.ndarray],
    sample_factor_dof: np.ndarray,
    masks: Dict[str, np.ndarray],
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """Return total and label-resolved periodogram contribution for one block."""
    data = np.array(data_block, dtype=np.float64, copy=True)

    # Detrend each DOF, apply window, then remove any tiny DC introduced by windowing.
    data -= data.mean(axis=0, keepdims=True)
    data *= window[:, None]
    data -= data.mean(axis=0, keepdims=True)

    if sqrt_mass_dof is not None:
        data *= sqrt_mass_dof[None, :]

    X = np.fft.rfft(data, axis=0)
    X[0, :] = 0.0
    P = np.abs(X) ** 2

    # Correct for type-wise atom sampling, if used.
    P *= sample_factor_dof[None, :]

    P_total = P.sum(axis=1)
    P_parts = {
        label: (P[:, mask].sum(axis=1) if mask.any() else np.zeros_like(P_total))
        for label, mask in masks.items()
    }
    return P_total, P_parts


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Publication-grade MD-derived VDOS/PDOS from a LAMMPS velocity dump."
    )
    parser.add_argument("--dump", required=True, help="LAMMPS dump file with id type vx vy vz columns.")
    parser.add_argument("--timestep_fs", type=float, required=True, help="MD timestep in fs, e.g. 0.5.")
    parser.add_argument("--dump_interval_steps", type=int, default=1, help="LAMMPS dump interval in MD steps.")
    parser.add_argument("--stride", type=int, default=1, help="Use every Nth dumped frame. Default 1.")
    parser.add_argument("--block_len", type=int, default=16384, help="Block length in sampled frames. Default 16384.")
    parser.add_argument(
        "--max_sample_per_type",
        type=int,
        default=0,
        help="Max atoms sampled per atom type. Use 0 for all atoms; recommended for final publication.",
    )
    parser.add_argument(
        "--sample_seed",
        type=int,
        default=12345,
        help="Random seed used only if --max_sample_per_type > 0 and sampling is required.",
    )
    parser.add_argument(
        "--species_map",
        required=True,
        help='Type-to-label map, e.g. "1:Al,2:B,3:Ga,4:N" or "1:S,2:Mo,3:S".',
    )
    parser.add_argument(
        "--mass_map",
        default=None,
        help='Optional label-to-mass map in amu, e.g. "Al=26.9815385,N=14.0067".',
    )
    parser.add_argument(
        "--primitive_atoms",
        type=int,
        default=2,
        help="Atoms in the primitive/reference cell for total-area normalization to 3N.",
    )
    parser.add_argument("--smooth_sigma_thz", type=float, default=0.25, help="Gaussian smoothing sigma in THz.")
    parser.add_argument("--xmin", type=float, default=0.0, help="Minimum x-axis frequency in THz.")
    parser.add_argument("--xmax", type=float, default=45.0, help="Maximum x-axis frequency in THz.")
    parser.add_argument("--prefix", default="PDOS", help="Output file prefix.")
    parser.add_argument(
        "--use_partial_block",
        action="store_true",
        help="Also process final incomplete block if it is sufficiently full. Disabled by default.",
    )
    parser.add_argument(
        "--min_partial_block_fraction",
        type=float,
        default=0.50,
        help="Minimum filled fraction for processing an incomplete final block when --use_partial_block is set.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=1200,
        help="PNG output resolution. Default 1200 dpi.",
    )

    args = parser.parse_args()
    validate_args(args)

    species_map = parse_species_map(args.species_map)
    mass_map = parse_mass_map(args.mass_map)
    use_mass = args.mass_map is not None

    print(f"Reading first frame: {args.dump}")
    natoms, ids, types, first_timestep = read_first_snapshot(args.dump)
    present_types = sorted(int(t) for t in np.unique(types))

    missing_type_labels = [t for t in present_types if t not in species_map]
    if missing_type_labels:
        print(
            "[warning] Some atom types are not in --species_map and will be labeled by type number: "
            + ", ".join(str(t) for t in missing_type_labels),
            file=sys.stderr,
        )
        for t in missing_type_labels:
            species_map[t] = str(t)

    label_by_type: Dict[int, str] = {int(t): species_map[int(t)] for t in present_types}
    labels = unique_preserve_order(label_by_type[t] for t in present_types)

    if use_mass:
        missing_masses = sorted({label_by_type[t] for t in present_types if label_by_type[t] not in mass_map})
        if missing_masses:
            raise ValueError(f"--mass_map is missing masses for labels: {missing_masses}")

    selected_ids, total_count_by_type, selected_count_by_type, sample_factor_by_type = build_selection(
        ids=ids,
        types=types,
        present_types=present_types,
        max_sample_per_type=args.max_sample_per_type,
        seed=args.sample_seed,
    )
    n_selected = len(selected_ids)
    id_to_selected_pos = {aid: i for i, aid in enumerate(selected_ids)}
    id_to_type = {int(aid): int(t) for aid, t in zip(ids.tolist(), types.tolist())}

    if args.max_sample_per_type > 0:
        print(
            "[warning] Atom sampling is enabled. Sampling factors will correct species weights, "
            "but using all atoms is preferred for final figures.",
            file=sys.stderr,
        )

    # Build label masks and mass/sample-factor arrays for selected atoms.
    masks: Dict[str, np.ndarray] = OrderedDict((label, np.zeros(3 * n_selected, dtype=bool)) for label in labels)
    selected_type = np.empty(n_selected, dtype=np.int32)
    selected_label: List[str] = []
    sample_factor_atom = np.empty(n_selected, dtype=np.float64)
    mass_atom: Optional[np.ndarray] = np.empty(n_selected, dtype=np.float64) if use_mass else None

    for pos, aid in enumerate(selected_ids):
        atom_type = id_to_type[aid]
        label = label_by_type[atom_type]
        selected_type[pos] = atom_type
        selected_label.append(label)
        masks[label][3 * pos:3 * pos + 3] = True
        sample_factor_atom[pos] = sample_factor_by_type[atom_type]
        if use_mass and mass_atom is not None:
            mass_atom[pos] = mass_map[label]

    sample_factor_dof = np.repeat(sample_factor_atom, 3)
    sqrt_mass_dof = np.repeat(np.sqrt(mass_atom), 3) if use_mass and mass_atom is not None else None

    # Full-system type/mass lookup for COM removal.
    if use_mass:
        mass_by_type = {t: mass_map[label_by_type[t]] for t in present_types}
        total_mass_all = sum(total_count_by_type[t] * mass_by_type[t] for t in present_types)
    else:
        mass_by_type = {}
        total_mass_all = float(natoms)

    window = np.hanning(args.block_len)
    block = np.zeros((args.block_len, 3 * n_selected), dtype=np.float64)

    acc_total: Optional[np.ndarray] = None
    acc_parts: Dict[str, np.ndarray] = OrderedDict()
    processed_blocks = 0
    filled = 0
    dump_frame_index = -1
    used_timesteps: List[int] = []
    all_timesteps: List[int] = []

    print(f"Atoms in dump: {natoms}")
    print(f"Selected atoms: {n_selected}")
    print("Present atom types:")
    for t in present_types:
        print(
            f"  type {t:>3d} -> {label_by_type[t]:>12s}: "
            f"total={total_count_by_type[t]}, selected={selected_count_by_type[t]}, "
            f"sampling_factor={sample_factor_by_type[t]:.6g}"
        )

    with open(args.dump, "r") as fh:
        line = fh.readline()
        while line and not line.startswith("ITEM: TIMESTEP"):
            line = fh.readline()

        while line:
            if not line.startswith("ITEM: TIMESTEP"):
                line = fh.readline()
                continue

            ts_line = fh.readline()
            if not ts_line:
                break
            timestep = int(ts_line.strip())
            all_timesteps.append(timestep)

            cidx, atom_lines = read_frame_after_timestep(fh, natoms)
            if cidx is None or atom_lines is None:
                break

            dump_frame_index += 1
            if dump_frame_index % args.stride != 0:
                line = fh.readline()
                continue

            row = np.zeros(3 * n_selected, dtype=np.float64)
            if use_mass:
                sum_mv = np.zeros(3, dtype=np.float64)
            else:
                sum_v = np.zeros(3, dtype=np.float64)

            for atom_line in atom_lines:
                parts = atom_line.split()
                aid = int(parts[cidx["id"]])
                atom_type = int(float(parts[cidx["type"]]))
                v = np.array(
                    [float(parts[cidx["vx"]]), float(parts[cidx["vy"]]), float(parts[cidx["vz"]])],
                    dtype=np.float64,
                )

                if use_mass:
                    mi = mass_by_type[atom_type]
                    sum_mv += mi * v
                else:
                    sum_v += v

                pos = id_to_selected_pos.get(aid)
                if pos is not None:
                    row[3 * pos:3 * pos + 3] = v

            if use_mass:
                vcom = sum_mv / total_mass_all
            else:
                vcom = sum_v / float(natoms)

            # Remove full-system COM velocity from selected velocities.
            vrow = row.reshape(n_selected, 3)
            vrow -= vcom[None, :]
            row = vrow.reshape(3 * n_selected)

            block[filled, :] = row
            filled += 1
            used_timesteps.append(timestep)

            if filled == args.block_len:
                p_total, p_parts = process_block(
                    data_block=block,
                    window=window,
                    sqrt_mass_dof=sqrt_mass_dof,
                    sample_factor_dof=sample_factor_dof,
                    masks=masks,
                )
                if acc_total is None:
                    acc_total = p_total.copy()
                    acc_parts = OrderedDict((label, p_parts[label].copy()) for label in labels)
                else:
                    acc_total += p_total
                    for label in labels:
                        acc_parts[label] += p_parts[label]
                processed_blocks += 1
                filled = 0

            line = fh.readline()

    if len(used_timesteps) < 2:
        fail("Fewer than two sampled frames were read. Check dump, stride, and input settings.")

    all_dt = np.diff(np.asarray(all_timesteps, dtype=np.int64))
    if len(all_dt) > 0:
        if not np.all(all_dt == all_dt[0]):
            fail("The raw dump timestep spacing is not constant. FFT/PDOS requires uniformly sampled data.")
        if int(all_dt[0]) != int(args.dump_interval_steps):
            fail(
                f"Provided --dump_interval_steps={args.dump_interval_steps}, but the dump file shows "
                f"an interval of {int(all_dt[0])} MD steps. Fix the command-line argument."
            )

    used_dt = np.diff(np.asarray(used_timesteps, dtype=np.int64))
    if not np.all(used_dt == used_dt[0]):
        fail("The sampled frame spacing is not constant after applying stride. FFT/PDOS requires uniform sampling.")

    expected_used_dt = int(args.dump_interval_steps) * int(args.stride)
    if int(used_dt[0]) != expected_used_dt:
        fail(
            f"Expected sampled interval {expected_used_dt} MD steps from dump_interval_steps*stride, "
            f"but found {int(used_dt[0])} MD steps."
        )

    if filled > 0:
        frac = filled / float(args.block_len)
        if args.use_partial_block and frac >= args.min_partial_block_fraction:
            # Use only the filled part and a matching-length Hann window. This avoids zero-padding bias.
            partial_window = np.hanning(filled)
            p_total, p_parts = process_block(
                data_block=block[:filled, :],
                window=partial_window,
                sqrt_mass_dof=sqrt_mass_dof,
                sample_factor_dof=sample_factor_dof,
                masks=masks,
            )
            # Partial block has different frequency grid if filled != block_len; do not mix it.
            print(
                "[warning] Final partial block was not included because different block length changes frequency grid. "
                "Use a trajectory length that contains an integer number of full blocks for final figures.",
                file=sys.stderr,
            )
        else:
            print(
                f"[info] Skipped final incomplete block: {filled}/{args.block_len} frames. "
                "This avoids zero-padding bias.",
                file=sys.stderr,
            )

    if processed_blocks == 0 or acc_total is None:
        fail(
            "No complete blocks were processed. Reduce --block_len, reduce --stride, or use a longer trajectory. "
            "For final figures, a longer trajectory is preferred."
        )

    # Convert to averaged spectra.
    s_total = acc_total / float(processed_blocks)
    s_parts = OrderedDict((label, acc_parts[label] / float(processed_blocks)) for label in labels)

    dt_eff_fs = float(used_dt[0]) * args.timestep_fs
    dt_eff_ps = dt_eff_fs / 1000.0
    freqs_thz = np.fft.rfftfreq(args.block_len, d=dt_eff_ps)
    nyquist_thz = 1.0 / (2.0 * dt_eff_ps)
    df_thz = freqs_thz[1] - freqs_thz[0] if len(freqs_thz) > 1 else float("nan")
    block_duration_ps = args.block_len * dt_eff_ps

    if args.xmax > nyquist_thz:
        raise ValueError(
            f"--xmax={args.xmax} THz exceeds Nyquist frequency {nyquist_thz:.6g} THz. "
            "Decrease xmax or dump velocities more frequently."
        )

    # 3N normalization.
    target_area = 3.0 * float(args.primitive_atoms)
    area_raw = _trapz(s_total, freqs_thz)
    if area_raw <= 0:
        fail("Raw total spectral area is non-positive. Check velocity dump and settings.")
    scale = target_area / area_raw

    g_parts = OrderedDict((label, s_parts[label] * scale) for label in labels)
    g_total = np.zeros_like(s_total)
    for label in labels:
        g_total += g_parts[label]
    g_total[0] = 0.0
    for label in labels:
        g_parts[label][0] = 0.0

    add_err_raw = float(np.max(np.abs(g_total - sum(g_parts[label] for label in labels))))

    # Smooth partials, force origin to zero, rebuild total, then renormalize once.
    sigma_bins = args.smooth_sigma_thz / df_thz if df_thz > 0 and args.smooth_sigma_thz > 0 else 0.0
    g_parts_s = OrderedDict((label, gaussian_smooth(g_parts[label], sigma_bins)) for label in labels)
    for label in labels:
        g_parts_s[label][0] = 0.0
        g_parts_s[label] = np.maximum(g_parts_s[label], 0.0)

    g_total_s = np.zeros_like(g_total)
    for label in labels:
        g_total_s += g_parts_s[label]

    area_s = _trapz(g_total_s, freqs_thz)
    if area_s > 0:
        scale_s = target_area / area_s
        g_total_s *= scale_s
        for label in labels:
            g_parts_s[label] *= scale_s

    add_err_s = float(np.max(np.abs(g_total_s - sum(g_parts_s[label] for label in labels))))

    # Save raw and smoothed data.
    raw_cols = [freqs_thz, g_total] + [g_parts[label] for label in labels]
    raw_header = (
        "freq_THz  G_total  " + "  ".join(f"G_{label}" for label in labels) + "\n"
        "MD-derived velocity power spectrum. Common normalization: integral(G_total) = 3*primitive_atoms.\n"
        f"mass_weighted={use_mass}; primitive_atoms={args.primitive_atoms}; "
        f"block_duration_ps={block_duration_ps:.8f}; df_THz={df_thz:.8f}; "
        "total_equals_sum_of_partials=True"
    )
    np.savetxt(f"{args.prefix}_THz_3N_raw.dat", np.column_stack(raw_cols), header=raw_header, fmt="%.10e")

    sm_cols = [freqs_thz, g_total_s] + [g_parts_s[label] for label in labels]
    sm_header = (
        "freq_THz  G_total_s  " + "  ".join(f"G_{label}_s" for label in labels) + "\n"
        "Smoothed MD-derived velocity power spectrum. Common normalization: integral(G_total_s) = 3*primitive_atoms.\n"
        f"Gaussian_sigma_THz={args.smooth_sigma_thz}; mass_weighted={use_mass}; "
        f"primitive_atoms={args.primitive_atoms}; total_equals_sum_of_partials=True"
    )
    np.savetxt(f"{args.prefix}_THz_3N_smoothed.dat", np.column_stack(sm_cols), header=sm_header, fmt="%.10e")

    # Metadata file.
    with open(f"{args.prefix}_metadata.txt", "w") as meta:
        meta.write("MD-derived VDOS/PDOS metadata\n")
        meta.write("================================\n")
        meta.write(f"dump_file: {args.dump}\n")
        meta.write(f"natoms_dump: {natoms}\n")
        meta.write(f"selected_atoms: {n_selected}\n")
        meta.write(f"mass_weighted: {use_mass}\n")
        meta.write(f"primitive_atoms: {args.primitive_atoms}\n")
        meta.write(f"timestep_fs: {args.timestep_fs}\n")
        meta.write(f"dump_interval_steps: {args.dump_interval_steps}\n")
        meta.write(f"stride: {args.stride}\n")
        meta.write(f"effective_dt_fs: {dt_eff_fs:.12g}\n")
        meta.write(f"nyquist_THz: {nyquist_thz:.12g}\n")
        meta.write(f"block_len_frames: {args.block_len}\n")
        meta.write(f"block_duration_ps: {block_duration_ps:.12g}\n")
        meta.write(f"frequency_spacing_THz: {df_thz:.12g}\n")
        meta.write(f"processed_blocks: {processed_blocks}\n")
        meta.write(f"used_frames: {len(used_timesteps)}\n")
        meta.write(f"skipped_final_partial_frames: {filled}\n")
        meta.write(f"smooth_sigma_THz: {args.smooth_sigma_thz}\n")
        meta.write(f"raw_additivity_error: {add_err_raw:.12e}\n")
        meta.write(f"smoothed_additivity_error: {add_err_s:.12e}\n")
        meta.write("\nAtom type summary:\n")
        for t in present_types:
            meta.write(
                f"type {t}: label={label_by_type[t]}, total={total_count_by_type[t]}, "
                f"selected={selected_count_by_type[t]}, sampling_factor={sample_factor_by_type[t]:.12g}\n"
            )

    # Plots.
    def pub_axes():
        ax = plt.gca()
        ax.spines["top"].set_visible(True)
        ax.spines["right"].set_visible(True)
        ax.tick_params(direction="in", which="both")
        ax.minorticks_on()
        return ax

    plt.figure(figsize=(7.0, 5.0))
    plt.plot(freqs_thz, g_total_s, linewidth=1.6, label="Total")
    pub_axes()
    plt.xlim(args.xmin, args.xmax)
    plt.ylim(bottom=0)
    plt.xlabel("Frequency (THz)")
    plt.ylabel("VDOS / PDOS (states cell$^{-1}$ THz$^{-1}$)")
    plt.tight_layout()
    plt.savefig(f"{args.prefix}_total_pub.png", dpi=args.dpi)
    plt.savefig(f"{args.prefix}_total_pub.pdf")
    plt.close()

    plt.figure(figsize=(7.0, 5.0))
    plt.plot(freqs_thz, g_total_s, label="Total", linewidth=1.6)
    for label in labels:
        plt.plot(freqs_thz, g_parts_s[label], label=label, linewidth=1.2)
    pub_axes()
    plt.xlim(args.xmin, args.xmax)
    plt.ylim(bottom=0)
    plt.xlabel("Frequency (THz)")
    plt.ylabel("VDOS / PDOS (states cell$^{-1}$ THz$^{-1}$)")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(f"{args.prefix}_components_pub.png", dpi=args.dpi)
    plt.savefig(f"{args.prefix}_components_pub.pdf")
    plt.close()

    print("[ok] PDOS completed")
    print(f"  output prefix: {args.prefix}")
    print(f"  processed full blocks: {processed_blocks}")
    print(f"  effective sampling interval: {dt_eff_fs:.6g} fs")
    print(f"  block duration: {block_duration_ps:.6g} ps")
    print(f"  frequency spacing: {df_thz:.6g} THz")
    print(f"  Nyquist frequency: {nyquist_thz:.6g} THz")
    print(f"  mass weighted: {use_mass}")
    print(f"  max raw additivity error: {add_err_raw:.3e}")
    print(f"  max smoothed additivity error: {add_err_s:.3e}")


if __name__ == "__main__":
    main()
