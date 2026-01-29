import pandas as pd
from rdkit import Chem
from .descriptors import compute_descriptors
from .rules import analyze_risk
from .alerts import check_structure_alerts

def run_pipeline(sdf_path):
    """
    Reads an SDF file, processes each molecule through the ADME engine,
    and returns a DataFrame of results.
    
    Args:
        sdf_path (str): Path to the input SDF file.
        
    Returns:
        pd.DataFrame: Processed results including IDs, descriptors, risk analysis, and ranks.
    """
    try:
        suppl = Chem.SDMolSupplier(sdf_path)
    except OSError:
        print(f"Error: Could not open SDF file at {sdf_path}")
        return pd.DataFrame()
        
    results = []
    
    print(f"Processing molecules from {sdf_path}...")
    
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
            
        mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Cmpd_{i+1}"
        
        # 1. Compute Descriptors
        desc = compute_descriptors(mol)
        if desc is None:
            continue
            
        # 2. Analyze Risk
        risk = analyze_risk(desc)
        if risk is None:
            continue
            
        # 3. Check Alerts
        alerts = check_structure_alerts(mol)
        alert_str = ",".join(alerts) if alerts else "None"
            
        # 4. Combine Data
        # We assume explicit SMILES isn't strictly needed for the dashboard if we have name/ID, 
        # but storing it is useful.
        smiles = Chem.MolToSmiles(mol)
        
        entry = {
            "Compound_ID": mol_name,
            "SMILES": smiles,
            "Alerts": alert_str
        }
        entry.update(desc)
        entry.update(risk)
        
        results.append(entry)
        
    df = pd.DataFrame(results)
    
    if df.empty:
        print("No valid molecules processed.")
        return df
        
    # 4. Rank Compounds
    # Logic: Prioritize SAFE-ROBUST, then SAFE-FRAGILE, then FAIL.
    # Within each tier, rank by Min_Margin (higher is better).
    # Since we want Rank 1 to be best, we sort descending on a combined score or sort key.
    
    # Map Tier to a numeric priority for sorting
    tier_priority = {
        "SAFE-ROBUST": 3,
        "SAFE-FRAGILE": 2,
        "FAIL": 1
    }
    
    df["Tier_Score"] = df["Tier"].map(tier_priority)
    
    # Sort: First by Tier (desc), then by Min_Margin (desc)
    df = df.sort_values(by=["Tier_Score", "Min_Margin"], ascending=[False, False])
    
    # Assign Rank
    df["Rank"] = range(1, len(df) + 1)
    
    # Clean up aux columns
    df = df.drop(columns=["Tier_Score"])
    
    print(f"Successfully processed {len(df)} compounds.")
    return df
