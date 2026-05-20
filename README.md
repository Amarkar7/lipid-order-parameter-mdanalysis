# Lipid Order Parameter Calculation using MDAnalysis

A Python-based lipid tail order parameter (\(S_{CD}\)) calculation tool using MDAnalysis with strict frame skipping behavior comparable to GROMACS `gmx order`.

This script computes membrane lipid tail ordering from molecular dynamics trajectories by analyzing the orientation of C–H bonds relative to the membrane normal (Z-axis).

---

# Features

- Strict stride frame processing
- GROMACS-like frame skipping behavior
- Supports separate sn1 and sn2 lipid tails
- Reads GROMACS `.tpr` and `.xtc`
- Uses `.ndx` index groups
- Computes averaged \(S_{CD}\) order parameters
- Automatic bonded hydrogen detection
- Distance-based fallback hydrogen search
- Produces `.xvg` output compatible with plotting tools

---

# Scientific Background

The deuterium order parameter is defined as:

\[
S_{CD} = \frac{1}{2}(3\cos^2\theta -1)
\]

Where:

- \(\theta\) = angle between the C–H bond vector and membrane normal
- \(S_{CD}\) = lipid tail order parameter

Interpretation:

| SCD Value | Meaning |
|---|---|
| High positive | Ordered lipid tails |
| Near zero | Disordered/random |
| Negative | Perpendicular orientation |

This quantity is commonly compared with NMR experimental measurements.

---

# Requirements

## Python

- Python >= 3.8

## Python Packages

- MDAnalysis
- NumPy

Install dependencies:

```bash
pip install -r requirements.txt
```

---

# Installation

Clone the repository:

```bash
git clone https://github.com/YOUR_USERNAME/lipid-order-parameter-mdanalysis.git
```

Enter the directory:

```bash
cd lipid-order-parameter-mdanalysis
```

Install requirements:

```bash
pip install -r requirements.txt
```

---

# Input Files

## Topology

```text
step7_production.tpr
```

## Trajectory

```text
md_noPBC_1.xtc
```

## Index files

```text
sn1.ndx
sn2.ndx
```

Each index group should contain equivalent carbon atoms across lipid residues.

Example:

```text
[ C2 ]
12 45 78 90

[ C3 ]
13 46 79 91
```

---

# Usage

Basic run:

```bash
python3 calc_order_mdanalysis_strict_stride.py
```

Run with stride:

```bash
python3 calc_order_mdanalysis_strict_stride.py --stride 100
```

Verbose mode:

```bash
python3 calc_order_mdanalysis_strict_stride.py --stride 100 --verbose
```

---

# Command Line Options

| Argument | Description |
|---|---|
| `--top` | Topology file |
| `--traj` | Trajectory file |
| `--ndx1` | sn1 index file |
| `--ndx2` | sn2 index file |
| `--out1` | sn1 output |
| `--out2` | sn2 output |
| `--stride` | Frame skipping interval |
| `--verbose` | Detailed logging |

---

# Output

The script generates:

```text
lipids_sn1_mda.xvg
lipids_sn2_mda.xvg
```

Format:

```text
1    0.212
2    0.198
3    0.175
```

Where:

- Column 1 = carbon group index
- Column 2 = averaged order parameter

---

# Workflow

1. Load trajectory and topology
2. Parse lipid tail index groups
3. Identify hydrogens bonded to carbons
4. Compute C–H bond vectors
5. Calculate angle relative to membrane normal
6. Compute \(S_{CD}\)
7. Average across:
   - frames
   - lipids
   - hydrogens
8. Write `.xvg` output

---

# Important Notes

- Membrane normal is assumed along the Z-axis
- Trajectory should be:
  - centered
  - rotationally aligned
  - PBC corrected

Recommended preprocessing using GROMACS:

```bash
gmx trjconv -pbc mol -center
```

---

# Example Use Cases

- Lipid bilayer simulations
- Membrane phase analysis
- Cholesterol ordering studies
- Protein-membrane interaction studies
- Comparing ordered vs disordered membranes

---

# Performance

The calculation scales approximately with:

```text
Frames × Lipids × Tail Carbons × Hydrogens
```

Large trajectories may require significant computational time.

---

# Citation

If you use this script in published work, please cite:

## MDAnalysis

Michaud-Agrawal et al., J Comput Chem. 2011.

## GROMACS

Abraham et al., SoftwareX. 2015.

---

# Author

Amar Kar  
Bioinformatics Centre  
Indian Institute of Technology Kharagpur

---

# License

MIT License
