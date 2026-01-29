from rdkit import Chem
from rdkit.Chem import AllChem

def enumerate_library(core_smiles, r_groups_smiles):
    """
    Generates a virtual library by attaching R-groups to a Core.
    
    Assumption:
    - Core has one or more dummy atoms `*` or `[*]` indicating attachment points.
    - R-groups are fragments that will replace the dummy atoms.
    - For simplicity in this demo: We replace the FIRST dummy atom found in Core with the R-group.
      (Combinatorial complexity for multiple sites is higher, we stick to 1-site or sequential replacement).
    """
    core = Chem.MolFromSmiles(core_smiles)
    if not core: return []
    
    products = []
    
    rxn_smarts = "[*:1]-[*-:2]>>[*:1]-[*:3]" 
    # This is too generic. 
    # Better approach for Robustness without defining strict reactions:
    # Use ReplaceSubstructs.
    
    # Simple approach: Replace '*' in core with the R-group structure
    # This usually requires the R-group to NOT have a dummy, or have one that we remove.
    # Let's assume standard input: "c1ccccc1[*]" (Core) + "CC" (R-group) -> "c1ccccc1C"
    
    for r_smi in r_groups_smiles:
        r_mol = Chem.MolFromSmiles(r_smi)
        if not r_mol: continue
        
        # Method: Replace * with R-group
        # We use Chem.ReplaceSubstructs
        # Pattern to replace: "*"
        dummy = Chem.MolFromSmarts("[*]")
        
        try:
            # RDKit ReplaceSubstructs returns a tuple of mols
            res = Chem.ReplaceSubstructs(core, dummy, r_mol, replaceAll=False)
            
            for prod in res:
                try:
                    Chem.SanitizeMol(prod)
                    prod_smi = Chem.MolToSmiles(prod)
                    products.append(prod_smi)
                except:
                    pass
        except:
            pass
            
    # DBSCAN or unique check
    return list(set(products))
