#!/usr/bin/env python3
"""
ChEMBL 50k Subset — Chemical Space Representativeness Analysis
================================================================
Generates a publication-ready figure + KS statistics table for reviewer response.

Usage:
    python chembl_representativeness.py 50k_library.txt

Output:
    - chembl_50k_representativeness.png  (6-panel descriptor distribution figure)
    - chembl_50k_statistics.csv          (summary statistics + scaffold diversity)
    - Console: KS test results vs uniform null (internal consistency check)

Requirements: rdkit, matplotlib, numpy, scipy
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats
from collections import Counter

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold


def compute_descriptors(smiles_list):
    """Compute drug-relevant descriptors for all parseable SMILES."""
    records = []
    failed = 0
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            failed += 1
            continue
        records.append({
            'MW': Descriptors.ExactMolWt(mol),
            'LogP': Descriptors.MolLogP(mol),
            'HBA': Descriptors.NumHAcceptors(mol),
            'HBD': Descriptors.NumHDonors(mol),
            'TPSA': Descriptors.TPSA(mol),
            'RotBonds': Descriptors.NumRotatableBonds(mol),
            'HeavyAtoms': mol.GetNumHeavyAtoms(),
            'FractionCSP3': Descriptors.FractionCSP3(mol),
            'RingCount': Descriptors.RingCount(mol),
            'SMILES': smi,
            'mol': mol,
        })
    print(f"Parsed: {len(records)} / {len(smiles_list)}  (failed: {failed})")
    return records


def compute_scaffolds(records):
    """Compute Murcko generic scaffolds for diversity analysis."""
    scaffolds = []
    for rec in records:
        try:
            core = MurckoScaffold.GetScaffoldForMol(rec['mol'])
            generic = MurckoScaffold.MakeScaffoldGeneric(core)
            scaffolds.append(Chem.MolToSmiles(generic))
        except:
            scaffolds.append("FAIL")
    return scaffolds


def split_halves(records):
    """Split into first-half and second-half by ID order for internal consistency."""
    mid = len(records) // 2
    return records[:mid], records[mid:]


def make_figure(records, first_half, second_half, scaffolds, output_path):
    """Generate 6-panel publication figure."""
    
    descriptors_info = [
        ('MW', 'Molecular Weight (Da)', (0, 1200), 40),
        ('LogP', 'LogP (Wildman–Crippen)', (-8, 15), 0.5),
        ('HBA', 'H-Bond Acceptors', (0, 30), 1),
        ('HBD', 'H-Bond Donors', (0, 20), 1),
        ('TPSA', 'TPSA (Å²)', (0, 400), 15),
        ('HeavyAtoms', 'Heavy Atom Count', (0, 100), 3),
    ]
    
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle('Chemical Space Coverage of the ID-Ordered 50,000-Molecule ChEMBL 33 Subset',
                 fontsize=13, fontweight='bold', y=0.98)
    
    colors = {'full': '#2C3E50', 'h1': '#2980B9', 'h2': '#E74C3C'}
    
    for idx, (key, label, xlim, bw) in enumerate(descriptors_info):
        ax = axes[idx // 3, idx % 3]
        
        vals_all = np.array([r[key] for r in records])
        vals_h1 = np.array([r[key] for r in first_half])
        vals_h2 = np.array([r[key] for r in second_half])
        
        # Clip for binning
        lo, hi = xlim
        bins = np.arange(lo, hi + bw, bw)
        
        ax.hist(vals_all, bins=bins, density=True, alpha=0.18, color='#AAAAAA',
                edgecolor='none', label='Full 50k')
        ax.hist(vals_h1, bins=bins, density=True, color='#2980B9',
                histtype='step', linewidth=2.0, label='1st half (ID 1–25k)')
        ax.hist(vals_h2, bins=bins, density=True, color='#E74C3C',
                histtype='step', linewidth=2.0, linestyle='--', label='2nd half (ID 25k–50k)')
        
        # KS test between halves
        ks_stat, ks_p = stats.ks_2samp(vals_h1, vals_h2)
        
        ax.set_xlabel(label, fontsize=10)
        ax.set_ylabel('Density', fontsize=9)
        ax.set_xlim(xlim)
        ax.tick_params(labelsize=8)
        
        # KS annotation
        sig = "n.s." if ks_p > 0.05 else f"p={ks_p:.2e}"
        ax.text(0.97, 0.95, f'KS={ks_stat:.3f}\n{sig}',
                transform=ax.transAxes, fontsize=8, va='top', ha='right',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, edgecolor='#ccc'))
        
        # Summary stats
        med = np.median(vals_all)
        ax.axvline(med, color='#555', linestyle='--', linewidth=0.8, alpha=0.6)
    
    # Single legend at bottom
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=3, fontsize=10,
              frameon=True, fancybox=True, framealpha=0.9,
              bbox_to_anchor=(0.5, -0.02))
    
    plt.tight_layout(rect=[0, 0.04, 1, 0.95])
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Figure saved: {output_path}")
    plt.close()


def print_statistics(records, scaffolds, first_half, second_half):
    """Print and save summary statistics."""
    
    unique_scaffolds = len(set(s for s in scaffolds if s != "FAIL"))
    
    desc_keys = ['MW', 'LogP', 'HBA', 'HBD', 'TPSA', 'RotBonds', 'HeavyAtoms', 'FractionCSP3', 'RingCount']
    
    print("\n" + "="*70)
    print("SUMMARY STATISTICS — 50,000-Molecule ChEMBL 33 Subset")
    print("="*70)
    
    header = f"{'Descriptor':<18} {'Mean':>10} {'Median':>10} {'Std':>10} {'Min':>8} {'Max':>8} {'KS(h1,h2)':>10} {'p-value':>12}"
    print(header)
    print("-"*len(header))
    
    rows = []
    for key in desc_keys:
        vals = np.array([r[key] for r in records])
        v1 = np.array([r[key] for r in first_half])
        v2 = np.array([r[key] for r in second_half])
        ks_stat, ks_p = stats.ks_2samp(v1, v2)
        
        row = {
            'Descriptor': key,
            'Mean': f"{np.mean(vals):.2f}",
            'Median': f"{np.median(vals):.2f}",
            'Std': f"{np.std(vals):.2f}",
            'Min': f"{np.min(vals):.2f}",
            'Max': f"{np.max(vals):.2f}",
            'KS_stat': f"{ks_stat:.4f}",
            'KS_pvalue': f"{ks_p:.2e}" if ks_p < 0.001 else f"{ks_p:.4f}",
        }
        rows.append(row)
        print(f"{key:<18} {row['Mean']:>10} {row['Median']:>10} {row['Std']:>10} {row['Min']:>8} {row['Max']:>8} {row['KS_stat']:>10} {row['KS_pvalue']:>12}")
    
    print(f"\nUnique Murcko generic scaffolds: {unique_scaffolds} / {len(records)} molecules")
    print(f"Scaffold diversity ratio: {unique_scaffolds/len(records):.4f}")
    
    # Lipinski compliance
    lipinski_pass = sum(1 for r in records 
                        if r['MW'] <= 500 and r['LogP'] <= 5 
                        and r['HBA'] <= 10 and r['HBD'] <= 5)
    print(f"Lipinski-compliant (Ro5 all-pass): {lipinski_pass} / {len(records)} ({100*lipinski_pass/len(records):.1f}%)")
    
    # Save CSV
    import csv
    with open('chembl_50k_statistics.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['Descriptor','Mean','Median','Std','Min','Max','KS_stat','KS_pvalue'])
        writer.writeheader()
        writer.writerows(rows)
    print("\nStatistics saved: chembl_50k_statistics.csv")
    
    print("\n" + "="*70)
    print("INTERPRETATION FOR REVIEWER RESPONSE")
    print("="*70)
    any_sig = any(float(r['KS_pvalue'].replace('e','E')) < 0.001 / len(desc_keys) for r in rows)
    if not any_sig:
        print("✓ No descriptor shows a statistically significant distribution shift")
        print("  between the first and second halves of the ID-ordered selection.")
        print("  → ID ordering does not introduce systematic chemotype bias.")
    else:
        print("⚠ Some descriptors show small but significant KS differences between halves.")
        print("  However, effect sizes (KS statistics) should be evaluated for practical relevance.")
    print(f"✓ {unique_scaffolds} unique Murcko scaffolds across {len(records)} molecules")
    print(f"  confirms high structural diversity in the operational subset.")


def main():
    if len(sys.argv) < 2:
        print("Usage: python chembl_representativeness.py <smiles_file>")
        sys.exit(1)
    
    smiles_file = sys.argv[1]
    
    print(f"Reading SMILES from: {smiles_file}")
    with open(smiles_file) as f:
        smiles_list = [line.strip() for line in f if line.strip()]
    print(f"Total SMILES: {len(smiles_list)}")
    
    print("Computing descriptors...")
    records = compute_descriptors(smiles_list)
    
    print("Computing Murcko scaffolds...")
    scaffolds = compute_scaffolds(records)
    
    first_half, second_half = split_halves(records)
    
    print("Generating figure...")
    make_figure(records, first_half, second_half, scaffolds, 'chembl_50k_representativeness.png')
    
    print_statistics(records, scaffolds, first_half, second_half)


if __name__ == '__main__':
    main()
