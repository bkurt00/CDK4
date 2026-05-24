import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def prepare_3d_ligands(input_csv="cdk46_admet_passed.csv", output_dir="Ligands_3D"):
    print(f">> 3D Dönüşüm Operasyonu Başlıyor: {input_csv}...")
    
    # Çıktı klasörünü oluştur
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    df = pd.read_csv(input_csv)
    success_count = 0
    
    for idx, row in df.iterrows():
        smi = row['SMILES']
        mol_name = f"Ligand_{idx:03d}"  # Örn: Ligand_000, Ligand_001
        
        # 1. SMILES'tan Molekül Yarat
        mol = Chem.MolFromSmiles(smi)
        if not mol: continue
        
        # 2. Hidrojenleri Ekle (Doking için KRİTİKTİR!)
        mol = Chem.AddHs(mol)
        
        # 3. 3D Koordinatları Üret (ETKDGv3 Algoritması)
        # RDKit bazen zor moleküllerde başarısız olabilir, kontrol ediyoruz.
        embed_status = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        
        if embed_status == -1:
            print(f"   [HATA] {mol_name} 3D'ye çevrilemedi (Sterik gerilim olabilir). Atlanıyor.")
            continue
            
        # 4. Geometri Optimizasyonu (Kuvvet Alanı ile rahatlatma)
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception as e:
            print(f"   [UYARI] {mol_name} MMFF94 optimizasyonu yapılamadı. Ham 3D kullanılacak.")
        
        # 5. SDF Olarak Kaydet (Smina'ya doğrudan verilecek)
        mol.SetProp("_Name", mol_name)  # Dosya içine ismini göm
        sdf_path = os.path.join(output_dir, f"{mol_name}.sdf")
        
        writer = Chem.SDWriter(sdf_path)
        writer.write(mol)
        writer.close()
        
        success_count += 1
        
    print(f">> İŞLEM TAMAM! {success_count} molekül 3 boyutlu olarak '{output_dir}' klasörüne dizildi.")
    print(">> Smina/Vina ordusu savaşa hazır.")

if __name__ == "__main__":
    prepare_3d_ligands()