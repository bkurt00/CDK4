# Chemical Space Representativeness Analysis

## Purpose

This analysis addresses whether selecting the first 50,000 molecules from ChEMBL 33 in CHEMBL_ID order introduces systematic chemotype bias compared to random sampling. The strategy splits the 50k subset into two non-overlapping halves by ID order and compares their molecular descriptor distributions using two-sample Kolmogorov–Smirnov (KS) tests. If ID ordering introduces bias, the two halves will show distributional divergence; if not, they will be statistically indistinguishable at practical effect sizes.

## Method

1. Read 50,000 SMILES from `50k_library.txt`
2. Parse each molecule with RDKit and compute six drug-relevant descriptors:
   - Molecular weight (Da)
   - Wildman–Crippen LogP
   - Number of hydrogen-bond acceptors (HBA)
   - Number of hydrogen-bond donors (HBD)
   - Topological polar surface area (TPSA, Å²)
   - Heavy-atom count
3. Compute Murcko generic scaffolds for all molecules (structural diversity metric)
4. Split the dataset at the midpoint: 1st half (ID 1–25,000) vs 2nd half (ID 25,001–50,000)
5. For each descriptor, run a two-sample KS test between the two halves
6. Generate a 6-panel overlay histogram figure and a summary statistics table

## Key Result

All KS statistics remained below 0.06 (range: 0.014–0.059). Despite formal statistical significance at N = 25,000 per group, the effect sizes are negligible, confirming that CHEMBL_ID ordering does not impose systematic clustering in molecular-property space. Murcko generic decomposition identified 12,423 unique scaffolds across 49,998 parseable molecules (diversity ratio = 0.249).

## Usage

```
python chembl_representativeness.py 50k_library.txt
```

Requirements: Python 3.8+, RDKit, NumPy, SciPy, Matplotlib

## Output

- `chembl_50k_representativeness.png` — 6-panel descriptor distribution figure (300 dpi)
- `chembl_50k_statistics.csv` — summary statistics with KS test results

## Contents

| File | Description |
|------|-------------|
| `50k_library.txt` | Raw SMILES input (50,000 molecules, one per line) |
| `chembl_representativeness.py` | Analysis script |
| `chembl_50k_representativeness.png` | Distribution figure |
| `chembl_50k_statistics.csv` | Descriptor statistics and KS test results |
| `README.md` | This file |

## Reference

This analysis accompanies the manuscript: Kurt, B. "AI-Assisted Identification of a Putative Allosteric Ligand Targeting the CDK4/Cyclin D1 Protein–Protein Interface." The dataset corresponds to Section 2.1 of the manuscript.
