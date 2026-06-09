import pandas as pd
from chembl_webresource_client.new_client import new_client

# Düzeltilen kısım burası: molecule objesini new_client üzerinden tanımlıyoruz
molecule = new_client.molecule

def read_smiles_from_file(file_path):
    """Dosyadaki SMILES'ları okur ve bir liste olarak döndürür."""
    try:
        df = pd.read_csv(file_path, header=None)
        smiles_list = [str(s).strip() for s in df.iloc[:, 0].dropna() if len(str(s).strip()) > 5]
        return set(smiles_list)
    except Exception as e:
        print(f"Dosya okunurken hata ({file_path}): {e}")
        return set()

print("Dosyalar okunuyor ve SMILES'lar çıkarılıyor...")
smiles_43 = read_smiles_from_file('first_filter.csv')
smiles_25 = read_smiles_from_file('cdk4_admet_passed.csv')

smiles_18_filtered = smiles_43 - smiles_25

print("-" * 40)
print(f"İlk Filtre (43 aday): {len(smiles_43)} molekül")
print(f"ADMET'i Geçen (25 aday): {len(smiles_25)} molekül")
print(f"ADMET'te Elenen: {len(smiles_18_filtered)} molekül")
print("-" * 40)

def get_phases_from_smiles(smiles_set):
    results = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 'Bulunamayan': 0}
    smiles_list = list(smiles_set)
    
    for i, smi in enumerate(smiles_list):
        try:
            records = molecule.filter(molecule_structures__canonical_smiles=smi).only(['max_phase'])
            
            if len(records) > 0:
                phase = records[0].get('max_phase')
                if phase in results:
                    results[phase] += 1
                elif phase is None:
                    results[0] += 1
                else:
                    results['Bulunamayan'] += 1
            else:
                results['Bulunamayan'] += 1
        except Exception as e:
            results['Bulunamayan'] += 1
            
        print(f"İşlenen: {i+1}/{len(smiles_list)}", end='\r')
        
    print()
    return results

if len(smiles_18_filtered) > 0:
    print("\n[1/2] ADMET Aşamasında Elenen Moleküller Analiz Ediliyor...")
    phases_18 = get_phases_from_smiles(smiles_18_filtered)
    print(f"Elenen Molekül Faz Dağılımı: {phases_18}")

if len(smiles_25) > 0:
    print("\n[2/2] ADMET'i Geçen Nihai 25 Molekül Analiz Ediliyor...")
    phases_25 = get_phases_from_smiles(smiles_25)
    print(f"Kalan 25 Molekül Faz Dağılımı: {phases_25}")