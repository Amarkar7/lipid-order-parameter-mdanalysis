#!/usr/bin/env python3

"""
calc_order_mdanalysis_strict_stride.py

MDAnalysis-based lipid tail order parameter (SCD) calculator
with strict frame-skipping behavior comparable to GROMACS gmx order.

Features
--------
- Strict stride frame selection:
      0, stride, 2*stride, ...

- Explicit trajectory frame jumps for reproducible processing

- Computes averaged lipid tail order parameters (SCD)

- Supports separate sn1 and sn2 lipid tail index groups

- Uses bonded hydrogen detection with distance-based fallback

Inputs
------
- Topology file (.tpr, .gro, .pdb, etc.)
- Trajectory file (.xtc, .trr, etc.)
- sn1 index file
- sn2 index file

Outputs
-------
- lipids_sn1_mda.xvg
- lipids_sn2_mda.xvg

Author
------
Amar Kar
Bioinformatics Centre
Indian Institute of Technology Kharagpur
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element
from collections import OrderedDict
import argparse
import sys

# Default output filenames
DEFAULT_OUT_SN1 = "lipids_sn1_mda.xvg"
DEFAULT_OUT_SN2 = "lipids_sn2_mda.xvg"

# Maximum allowed C-H bond distance (Å)
FALLBACK_DIST = 1.6


def parse_ndx_keep_order(filename):
    """
    Parse GROMACS .ndx file while preserving group order.

    Returns
    -------
    OrderedDict
        group_name -> list of atom indices
    """

    groups = OrderedDict()
    current = None

    with open(filename, 'r') as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            if line.startswith('[') and line.endswith(']'):

                current = line.strip('[]').strip()

                key = current
                kcount = 1

                while key in groups:
                    kcount += 1
                    key = f"{current}__{kcount}"

                groups[key] = []

            else:

                if current is None:
                    continue

                for token in line.split():

                    try:
                        idx = int(token)
                        groups[next(reversed(groups))].append(idx)

                    except ValueError:
                        pass

    return groups


def find_hydrogens_for_carbon(atom):
    """
    Identify hydrogen atoms bonded to a carbon atom.

    Topology bonds are preferred.
    If unavailable, fallback residue-based hydrogen search is used.
    """

    hydrogens = []

    try:

        for bond in atom.bonds:

            a0, a1 = bond.atoms
            other = a1 if a0 == atom else a0

            elem = getattr(other, 'element', None)

            if elem is None:
                elem = guess_atom_element(other.name)

            if elem is not None and elem.upper() == 'H':
                hydrogens.append(other)

            elif other.name.startswith('H') or other.name.startswith('h'):
                hydrogens.append(other)

    except Exception:
        pass

    if not hydrogens:

        for other in atom.residue.atoms:

            if other is atom:
                continue

            elem = getattr(other, 'element', None)

            if elem is None:
                elem = guess_atom_element(other.name)

            if (
                (elem and elem.upper() == 'H')
                or other.name.startswith('H')
                or other.name.startswith('h')
            ):
                hydrogens.append(other)

    return hydrogens


def compute_order_for_groups(
    universe,
    groups_ordered_dict,
    output_file,
    frame_indices,
    verbose=False
):
    """
    Compute averaged lipid order parameters for all index groups.
    """

    group_items = list(groups_ordered_dict.items())

    n_groups = len(group_items)

    if verbose:
        print(f"[INFO] {n_groups} groups found in '{output_file}'")

    groups_atoms = []

    for group_name, atom_indices_1based in group_items:

        atom_inds0 = [i - 1 for i in atom_indices_1based]

        if len(atom_inds0) == 0:
            groups_atoms.append([])
            continue

        carbons = universe.atoms[atom_inds0]

        per_carbon = []

        for carbon in carbons:

            h_candidates = find_hydrogens_for_carbon(carbon)

            per_carbon.append({
                'C': carbon,
                'H_cand': h_candidates
            })

        groups_atoms.append(per_carbon)

    if verbose:

        for i, (gname, atomlist) in enumerate(group_items, start=1):

            print(
                f"  Group {i:3d}: "
                f"'{gname}'  "
                f"atoms={len(atomlist)}"
            )

    S_sum = np.zeros(n_groups, dtype=float)
    S_count = np.zeros(n_groups, dtype=int)

    zvec = np.array([0.0, 0.0, 1.0])

    total_frames = len(frame_indices)

    print(f"[INFO] Processing {total_frames} frames")

    for frame_counter, frame_idx in enumerate(frame_indices, start=1):

        try:
            universe.trajectory[frame_idx]

        except Exception as e:

            print(f"[ERROR] Cannot access frame {frame_idx}: {e}")
            sys.exit(1)

        if verbose and (
            frame_counter % max(1, total_frames // 10) == 0
        ):
            print(
                f"  processed {frame_counter}/{total_frames} "
                f"(frame {frame_idx})"
            )

        for gi, per_carbon in enumerate(groups_atoms):

            if not per_carbon:
                continue

            for entry in per_carbon:

                carbon = entry['C']
                cpos = carbon.position

                h_candidates = entry['H_cand']

                used = False

                for hydrogen in h_candidates:

                    try:
                        hpos = hydrogen.position

                    except Exception:
                        continue

                    dist = np.linalg.norm(hpos - cpos)

                    if dist > FALLBACK_DIST:
                        continue

                    vec = hpos - cpos

                    norm = np.linalg.norm(vec)

                    if norm == 0:
                        continue

                    cos_theta = np.dot(vec, zvec) / norm

                    S = 0.5 * (
                        3.0 * (cos_theta ** 2) - 1.0
                    )

                    S_sum[gi] += S
                    S_count[gi] += 1

                    used = True

                # Last-resort global hydrogen search

                if not used:

                    for other in universe.atoms:

                        if other is carbon:
                            continue

                        elem = getattr(other, 'element', None)

                        if elem is None:
                            elem = guess_atom_element(other.name)

                        if not (
                            (elem and elem.upper() == 'H')
                            or other.name.startswith('H')
                            or other.name.startswith('h')
                        ):
                            continue

                        hpos = other.position

                        dist = np.linalg.norm(hpos - cpos)

                        if dist <= FALLBACK_DIST:

                            vec = hpos - cpos

                            norm = np.linalg.norm(vec)

                            if norm == 0:
                                continue

                            cos_theta = np.dot(vec, zvec) / norm

                            S = 0.5 * (
                                3.0 * (cos_theta ** 2) - 1.0
                            )

                            S_sum[gi] += S
                            S_count[gi] += 1

    S_avg = np.full(n_groups, np.nan)

    for i in range(n_groups):

        if S_count[i] > 0:
            S_avg[i] = S_sum[i] / float(S_count[i])

    with open(output_file, 'w') as fh:

        fh.write(
            "# MDAnalysis lipid order parameters\n"
        )

        fh.write(
            "# Columns: index    S_CD\n"
        )

        for i, sval in enumerate(S_avg, start=1):

            if np.isnan(sval):

                fh.write(f"{i:12d}    nan\n")

            else:

                fh.write(
                    f"{i:12d}    {sval: .6f}\n"
                )

    print(f"[INFO] Wrote output: {output_file}")


def build_frame_indices_strict(
    n_total_frames,
    stride
):
    """
    Generate strict frame sequence:
        0, stride, 2*stride ...
    """

    s = int(stride)

    if s <= 0:
        raise ValueError("stride must be positive")

    return list(range(0, n_total_frames, s))


def main():

    parser = argparse.ArgumentParser(
        description=(
            "MDAnalysis strict-stride "
            "lipid order parameter calculator"
        )
    )

    parser.add_argument(
        "--top",
        required=True,
        help="Topology file"
    )

    parser.add_argument(
        "--traj",
        required=True,
        help="Trajectory file"
    )

    parser.add_argument(
        "--ndx1",
        required=True,
        help="sn1 index file"
    )

    parser.add_argument(
        "--ndx2",
        required=True,
        help="sn2 index file"
    )

    parser.add_argument(
        "--out1",
        default=DEFAULT_OUT_SN1,
        help="sn1 output filename"
    )

    parser.add_argument(
        "--out2",
        default=DEFAULT_OUT_SN2,
        help="sn2 output filename"
    )

    parser.add_argument(
        "--stride",
        type=int,
        default=None,
        help=(
            "Process frames as "
            "0, stride, 2*stride ..."
        )
    )

    parser.add_argument(
        "--verbose",
        action="store_true"
    )

    args = parser.parse_args()

    print("[INFO] Loading trajectory")

    universe = mda.Universe(
        args.top,
        args.traj
    )

    total_frames = universe.trajectory.n_frames

    print(
        f"[INFO] Total frames: {total_frames}"
    )

    if args.stride is not None:

        frame_indices = build_frame_indices_strict(
            total_frames,
            args.stride
        )

        print(
            f"[INFO] Strict stride={args.stride} "
            f"-> {len(frame_indices)} frames selected"
        )

    else:

        frame_indices = list(range(total_frames))

        print(
            f"[INFO] Using all "
            f"{len(frame_indices)} frames"
        )

    print(f"[INFO] Parsing {args.ndx1}")
    groups1 = parse_ndx_keep_order(args.ndx1)

    print(f"[INFO] Parsing {args.ndx2}")
    groups2 = parse_ndx_keep_order(args.ndx2)

    if not groups1 or not groups2:

        print(
            "[ERROR] One or more "
            "index files contain no groups"
        )

        sys.exit(1)

    compute_order_for_groups(
        universe,
        groups1,
        args.out1,
        frame_indices,
        verbose=args.verbose
    )

    universe.trajectory.rewind()

    compute_order_for_groups(
        universe,
        groups2,
        args.out2,
        frame_indices,
        verbose=args.verbose
    )

    print("[INFO] DONE")


if __name__ == "__main__":
    main()
