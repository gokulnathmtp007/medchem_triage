from rdkit import Chem
from rdkit.Chem import FilterCatalog

def get_pains_catalog():
    """Initializes and returns the PAINS filter catalog."""
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog.FilterCatalog(params)

# Global instance to avoid reloading overhead
_pains_catalog = get_pains_catalog()

def check_structure_alerts(mol):
    """
    Checks if a molecule matches any PAINS (Pan-Assay Interference Compounds) filters.
    
    Args:
        mol (rdkit.Chem.Mol): Molecule to check.
        
    Returns:
        list: A list of alert descriptions found (e.g., ["PAINS_A", "PAINS_B"]). 
              Returns empty list if clean.
    """
    if mol is None:
        return []

    alerts = []
    
    # Check PAINS
    if _pains_catalog.HasMatch(mol):
        matches = _pains_catalog.GetMatches(mol)
        for match in matches:
            description = match.GetDescription()
            alerts.append(description)
            
    # We could add more custom alerts here (e.g. Michael Acceptors, Aldehydes)
    
    return alerts
