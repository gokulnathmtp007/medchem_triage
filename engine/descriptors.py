from rdkit import Chem
from rdkit.Chem import Descriptors, Scaffolds
from rdkit.Chem import MolSurf
from rdkit.Chem import Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold

def compute_descriptors(mol):
    """
    Computes physicochemical descriptors for a given RDKit molecule.
    
    Args:
        mol (rdkit.Chem.Mol): The molecule to analyze.
        
    Returns:
        dict: A dictionary of descriptors or None if calculation fails.
    """
    if mol is None:
        return None
    
    try:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        rotb = Lipinski.NumRotatableBonds(mol)
        
        # Determine Murcko Scaffold
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            scaffold_smiles = Chem.MolToSmiles(scaffold)
        except:
            scaffold_smiles = "Unknown"
        
        return {
            "MW": mw,
            "LogP": logp,
            "TPSA": tpsa,
            "HBD": hbd,
            "HBA": hba,
            "RotB": rotb,
            "Scaffold": scaffold_smiles
        }
    except Exception as e:
        # In a production context we might log this, for now just print or ignore
        print(f"Error computing descriptors: {e}")
        return None
