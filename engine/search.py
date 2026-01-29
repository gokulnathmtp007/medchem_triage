from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

def compute_fingerprint(mol):
    """Computes Morgan Fingerprint (Radius 2, 2048 bits)."""
    if mol is None: return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

def search_similarity(query_mol, target_mols, threshold=0.7):
    """
    Finds compounds with Tanimoto similarity > threshold.
    
    Args:
        query_mol: RDKit molecule
        target_mols: List of (id, mol) tuples or just mols
        threshold: float (0.0 to 1.0)
        
    Returns:
        List of results: [(id, similarity_score, mol), ...]
    """
    if query_mol is None: return []
    
    query_fp = compute_fingerprint(query_mol)
    results = []
    
    for item in target_mols:
        # Handle input formats: (id, mol) or just mol
        if isinstance(item, tuple):
            mid, mol = item
        else:
            mid, mol = "Unknown", item
            
        if mol is None: continue
            
        fp = compute_fingerprint(mol)
        if fp:
            sim = DataStructs.TanimotoSimilarity(query_fp, fp)
            if sim >= threshold:
                results.append((mid, sim, mol))
                
    # Sort by similarity desc
    results.sort(key=lambda x: x[1], reverse=True)
    return results

def search_substructure(query_mol, target_mols):
    """
    Finds compounds containing the query substructure.
    """
    if query_mol is None: return []
    
    results = []
    
    for item in target_mols:
        if isinstance(item, tuple):
            mid, mol = item
        else:
            mid, mol = "Unknown", item
            
        if mol and mol.HasSubstructMatch(query_mol):
            results.append((mid, 1.0, mol)) # Score 1.0 for match
            
    return results
