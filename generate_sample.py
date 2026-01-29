from rdkit import Chem
from rdkit.Chem import AllChem
import os

# Create directory
os.makedirs("data/raw", exist_ok=True)

# List of diverse compounds: (Name, SMILES, Category)
# Categories: Drug, PAINS, Fragment, Failed_Prop, Research
datasets = [
    # --- KNOWN DRUGS (Mostly Safe) ---
    ("Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O", "NSAID"),
    ("Paracetamol", "CC(=O)NC1=CC=C(O)C=C1", "Analgesic"),
    ("Ibuprofen", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "NSAID"),
    ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Stimulant"),
    ("Diazepam", "CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3", "Benzo"),
    ("Metformin", "CN(C)C(=N)NC(=N)N", "Diabetes"),
    ("Warfarin", "CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O", "Blood"),
    ("Sildenafil", "CCCC1=NN(C2=C1NC(=NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C", "Vasodilator"),
    ("Atorvastatin", "CC(C)C1=C(C(=C(N1CC[C@H](O)C[C@H](O)CC(=O)O)C2=CC=CC=C2)C3=CC=C(C=C3)F)C(=O)NC4=CC=CC=C4", "Statin"),
    ("Omeprazole", "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC", "PPI"),
    ("Imatinib", "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5", "Kinase"),
    ("Ranitidine", "CNCC1=CC=C(O1)CSCC2=C(NC(=C2[N+](=O)[O-])NC)C", "H2Block"),
    ("Loratadine", "CCOC(=O)N1CCC(=C2C3=C(CCC2)C=C(C=C3)Cl)C4=C1C=CC=N4", "Antihistamine"),
    ("Fluoxetine", "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F", "SSRI"),
    ("Ciprofloxacin", "C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O", "Antibiotic"),

    # --- FAILURES / BORDERLINE (Large, Greasy, etc) ---
    ("Vancomycin_Core", "CC1C(C2C(C(C(OC2OC3=C(OC4=C(C(C(C(C(C(C5=CC(=CC(=C54)O)Cl)O)C(=O)N6C(C(=O)NC(C7=C(C(=C(C=C7)O)Cl)O)C(=O)N1)CC8=CC(=C(C=C8)O)OC6O)N)C(=O)O)Cl)CO)O)O)(C)N)O", "Antibiotic_Huge"),
    ("Cyclosporine", "CCC1C(=O)N(CC(=O)N(C(C(=O)NC(C(=O)N(C(C(=O)NC(C(=O)NC(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)NC(C(=O)O1)C(C)C)C)C(C)C)C)C(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C", "Immunosuppressant"),
    ("Very_Long_Alkane", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "Grease"),
    ("PEG_Chain", "COCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO", "Polymer"),
    ("Large_TPSA_Fail", "OC(O)C(O)C(O)C(O)C(O)C(=O)NC(CO)(CO)CO", "Polar"),

    # --- PAINS (Pan-Assay Interference) & TOXIC ---
    ("Rhodanine_PAINS", "O=C1NC(=S)S/C1=C\\c2ccc(O)cc2", "PAINS"),
    ("Quinone_PAINS", "O=C1C=CC(=O)C=C1", "Quinone"),
    ("Curcumin_Aggregator", "COC1=C(C=CC(=C1)/C=C/C(=O)CC(=O)/C=C/C2=CC(=C(C=C2)O)OC)O", "Aggregator"),
    ("Michael_Acceptor", "C=CC(=O)C1=CC=CC=C1", "Warhead"),
    
    # --- CHEMICAL SERIES (Scaffold Analysis) ---
    # Series A: Indoles (Generic scaffold)
    ("Indole_1", "C1=CC=C2C(=C1)C=CN2C", "Indole_Series"),
    ("Indole_2", "C1=CC=C2C(=C1)C=CN2CC", "Indole_Series"),
    ("Indole_3", "C1=CC=C2C(=C1)C=CN2CC(=O)O", "Indole_Series"),
    ("Indole_4", "C1=CC=C2C(=C1)C=CN2CCC3=CC=CC=C3", "Indole_Series"),
    ("Indole_5", "FC1=CC=C2C(=C1)C=CN2C", "Indole_Series"),
    ("Indole_6", "ClC1=CC=C2C(=C1)C=CN2C", "Indole_Series"),
    
    # Series B: Sulfonamides (Often more polar)
    ("Sulfa_1", "NS(=O)(=O)C1=CC=CC=C1", "Sulfa_Series"),
    ("Sulfa_2", "NS(=O)(=O)C1=CC=C(C)C=C1", "Sulfa_Series"),
    ("Sulfa_3", "NS(=O)(=O)C1=CC=C(Cl)C=C1", "Sulfa_Series"),
    ("Sulfa_4", "NS(=O)(=O)C1=CC=C(NO2)C=C1", "Sulfa_Series"),
    ("Sulfa_5", "CC(=O)NS(=O)(=O)C1=CC=CC=C1", "Sulfa_Series"),
    
    # --- FRAGMENTS (Small MW) ---
    ("Frag_Benzene", "C1=CC=CC=C1", "Fragment"),
    ("Frag_Pyridine", "C1=CC=NC=C1", "Fragment"),
    ("Frag_Aniline", "NC1=CC=CC=C1", "Fragment"),
    ("Frag_Phenol", "OC1=CC=CC=C1", "Fragment"),
    ("Frag_Acid", "C1=CC=C(C=C1)C(=O)O", "Fragment"),
    
    # --- RANDOM/TEST ---
    ("Test_Halogen_1", "FC(F)(F)C1=CC=CC=C1", "Test"),
    ("Test_Halogen_2", "ClC1=CC=C(Br)C=C1I", "Test"),
    ("Test_Macrocycle_Safe", "C1CCCCCCCC1", "Test"), # Simple
]

# Generate variations to reach 50+
# Let's simple add duplicates with slight mods if needed, or just more diversity
extra_compounds = [
    (f"Var_Alkane_{i}", "C" * (i+5), "Alkane") for i in range(10)
]

full_list = datasets + extra_compounds

w = Chem.SDWriter("data/raw/sample.sdf")

count = 0
for name, smiles, cat in full_list:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mol.SetProp("_Name", name)
        mol.SetProp("Category", cat)
        try:
            AllChem.Compute2DCoords(mol)
        except:
            pass
        w.write(mol)
        count += 1

w.close()
print(f"Created data/raw/sample.sdf with {count} compounds.")
