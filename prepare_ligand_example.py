from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

def prepare_ligand_for_docking(input_sdf, output_pdbqt):
    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
    mol = supplier[0]

    if mol is None:
        raise ValueError("Invalid ligand structure.")

    mol = Chem.AddHs(mol)

    if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass

    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    pdbqt_string, ok, _ = PDBQTWriterLegacy.write_string(setups[0])

    if not ok:
        raise RuntimeError("PDBQT conversion failed.")

    with open(output_pdbqt, "w") as f:
        f.write(pdbqt_string)
