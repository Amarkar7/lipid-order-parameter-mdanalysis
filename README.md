# Lipid Order Parameter Calculation using MDAnalysis

Python-based lipid tail order parameter (SCD) calculation tool using MDAnalysis.

This script computes membrane lipid tail ordering from molecular dynamics trajectories by analyzing the orientation of C–H bonds relative to the membrane normal (Z-axis).

---

# Features

- Strict frame skipping
- GROMACS-like behavior
- Supports sn1 and sn2 lipid tails
- Reads topology and trajectory files
- Uses index groups
- Computes averaged lipid order parameters
- Generates XVG output

---

# Installation

Clone repository:

```bash
git clone https://github.com/YOUR_USERNAME/lipid-order-parameter-mdanalysis.git
```

Install dependencies:

```bash
pip install -r requirements.txt
```

---

# Required Inputs

The script requires:

- topology file
- trajectory file
- sn1 index file
- sn2 index file

Supported formats depend on MDAnalysis.

---

# Usage

```bash
python3 calc_order_mdanalysis_strict_stride.py \
    --top topology.tpr \
    --traj trajectory.xtc \
    --ndx1 sn1.ndx \
    --ndx2 sn2.ndx \
    --stride 100
```

---

# Output

The script generates:

- lipids_sn1_mda.xvg
- lipids_sn2_mda.xvg

Example output:

```text
1    0.212
2    0.198
3    0.175
```

---

# Scientific Background

The order parameter is calculated as:

SCD = 0.5 × (3cos²θ − 1)

where θ is the angle between the C–H bond vector and membrane normal.

Higher values indicate more ordered lipid tails.

---

# Notes

- Membrane normal is assumed along Z-axis
- Trajectory should ideally be PBC corrected
- Large trajectory files are intentionally excluded from this repository

---

# Author

Amar Kar  
Bioinformatics Centre  
Indian Institute of Technology Kharagpur

---

# License

MIT License
