from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import random
import os

# Define a diverse set of seed SMILES (Drugs, Fragments, Toxins)
SEEDS = [
    ("Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"),
    ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
    ("Atorvastatin", "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=CC=C2)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4"),
    ("Vancomycin", "CC1C(C(CC(O1)OCC2C(C(C(O2)OC3C(OC(C(C3O)N)O)C4=C(C=C(C(=C4)C(C(C(=O)NC5C(C(C=C(C5=O)O)OC6C(C(C(C(O6)CO)O)O)OC7CC(C(C(O7)C)O)N(C)C)Cl)NC(=O)C(C8=CC(=C(C=C8)O)C9=C(C=C(C(=C9)O)O)C(NC(=O)C(C(C1=CC(=C(O)C=C1)Cl)O)NC(=O)C(C(C1=CC(=C(O)C=C1)Cl)O)NC(=O)C(CC(C)C)NC)C(=O)O)CC(=O)N)O)Cl)O)O)C)O)N)O"), 
    ("Rhodanine", "C1CS(=O)(=O)N1"), # Simple replacement
    ("Benzene", "c1ccccc1"),
    ("Fragment1", "N#Cc1ccccc1"),
    ("Fragment2", "OC(=O)c1ccccc1"),
    ("Sulfonamide", "NS(=O)(=O)c1ccc(C)cc1"),
    ("Indole", "c1ccc2[nH]ccc2c1"),
    ("Quinone", "O=C1C=CC(=O)C=C1"),
    ("Penicillin", "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C"),
    ("Diazepam", "CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3"),
]

def mutate_smiles(smi):
    """
    Simulates finding an analog by slight random changes (Virtual Library approach).
    This doesn't use true chemical reactions, but uses RDKit to randomize and potentially add atoms
    by simple string/mol manipulation or enumeration.
    
    Actually, to be robust, let's use the 'Enumeration' logic we built!
    We will treat each seed as a scaffold (if possible) or just add R-groups to it.
    """
    try:
        mol = Chem.MolFromSmiles(smi)
        if not mol: return []
        
        analogs = []
        
        # Strategy 1: Replace random H with methyl, chloro, hydroxyl (Simple Mutation)
        # We find a random atom that acts as a connector?
        # Simpler: Add simple R-groups to the existing SMILES string if it has a sensible place? No.
        
        # Better Strategy: Just use the seed.
        # But we need 1000 different ones.
        
        # Correct Strategy: Reaction based mutation.
        # RXN: [C:1] >> [C:1]C (Analyze Methylation)
        # RXN: [c:1] >> [c:1]Cl (Chlorination)
        
        rxns = [
            AllChem.ReactionFromSmarts('[CH:1]>>[C:1](C)'), # Methylate aliphatic
            AllChem.ReactionFromSmarts('[cH:1]>>[c:1](Cl)'), # Chlorinate aromatic
            AllChem.ReactionFromSmarts('[cH:1]>>[c:1](O)'), # Hydroxylate aromatic
            AllChem.ReactionFromSmarts('[OH:1]>>[O:1](C)'), # Ether formation
            AllChem.ReactionFromSmarts('[NH:1]>>[N:1](C(=O)C)'), # Acetylation
        ]
        
        for rxn in rxns:
            ps = rxn.RunReactants((mol,))
            for p in ps:
                try:
                    p0 = p[0]
                    Chem.SanitizeMol(p0)
                    analogs.append(Chem.MolToSmiles(p0))
                except:
                    pass
        
        return list(set(analogs))
    except:
        return []

def generate_large_dataset(target_count=1000):
    print(f"Generating ~{target_count} compounds from {len(SEEDS)} seeds...")
    
    data = []
    
    # 1. Add Seeds
    for name, smi in SEEDS:
        data.append({"Compound_ID": name, "SMILES": smi})
        
    # 2. Mutate Seeds until we hit target
    generated_count = len(data)
    
    # Keep list of source miles to mutate
    pool = [s for n, s in SEEDS]
    
    import random
    
    iterations = 0
    while generated_count < target_count and iterations < 50: # Avoid infinite loop
        iterations += 1
        new_pool = []
        for smi in pool:
            # Generate analogs
            analogs = mutate_smiles(smi)
            for ana in analogs:
                if generated_count >= target_count: break
                
                # Check uniqueness roughly (in this run)
                # (Skipping strict unique check for speed, relying on ID)
                
                # Add to DB
                cid = f"ZINC_Sim_{generated_count+1:05d}"
                data.append({"Compound_ID": cid, "SMILES": ana})
                new_pool.append(ana)
                generated_count += 1
            
            if generated_count >= target_count: break
        
        # Expand the pool for next generation (iterative growth)
        if new_pool:
            pool = new_pool
            # Shuffle to avoid drilling down just one line
            random.shuffle(pool) 
            # Limit pool size to keep diversity?
            pool = pool[:50] 
        else:
            break
            
    # Create DataFrame
    df = pd.DataFrame(data)
    output_path = "data/raw/large_scale.csv"
    df.to_csv(output_path, index=False)
    print(f"✅ Generated {len(df)} compounds. Saved to {output_path}")

if __name__ == "__main__":
    generate_large_dataset(1500) # Aim high to account for failures
