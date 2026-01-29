from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina

def cluster_compounds(mol_list, cutoff=0.7):
    """
    Clusters compounds using Butina algorithm based on Tanimoto distance.
    Distance = 1 - Similarity.
    
    Args:
        mol_list: List of (id, mol) tuples.
        cutoff: Distance cutoff (default 0.7 means similarity > 0.3 are clustered? 
                Wait, Butina needs dist. If we want similarity > 0.7, dist < 0.3.
                Usually cutoff represents DISTANCE. So 0.3 cutoff = >0.7 sim.
                Let's simplify: User usually thinks in Similarity. 
                If input cutoff is 0.7 (high sim), distance is 0.3.
    
    Returns:
        dict: {cluster_id: [List of compound_ids], ...} 
        list: sorted clusters by size
    """
    # 1. Generate FPs
    fps = []
    ids = []
    valid_mols = []
    
    for mid, mol in mol_list:
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fps.append(fp)
            ids.append(mid)
            valid_mols.append(mol)
            
    if not fps: return {}, []

    # 2. Compute Distance Matrix (Linear lower triangular)
    dists = []
    npts = len(fps)
    for i in range(1, npts):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        for sim in sims:
            dists.append(1 - sim)

    # 3. Cluster (Butina)
    # cutoff: elements within this range of centroid are in cluster
    # If we want Similarity > 0.7, then Dist < 0.3. 
    # Let's interpret the arg 'cutoff' as Similarity Threshold (e.g. 0.7).
    dist_cutoff = 1.0 - cutoff
    
    clusters = Butina.ClusterData(dists, npts, dist_cutoff, isDistData=True)
    
    # 4. Format Results
    # clusters is a tuple of tuples: ((idx1, idx2...), (idx3...), ...)
    # The first element of each tuple is the centroid.
    
    cluster_map = {}
    sorted_clusters = []
    
    for i, cluster in enumerate(clusters):
        c_ids = [ids[idx] for idx in cluster]
        cluster_map[i+1] = c_ids
        sorted_clusters.append({
            "Cluster_ID": i+1,
            "Centroid_ID": ids[cluster[0]],
            "Size": len(cluster),
            "Members": c_ids
        })
        
    return cluster_map, sorted_clusters
