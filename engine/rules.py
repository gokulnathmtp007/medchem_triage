
# Standard Lipinski-like thresholds
# MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10. 
# Added TPSA <= 140 and RotB <= 10 for completeness in many drug-like definitions.

RULES = {
    "MW": {"threshold": 500.0, "type": "upper"},
    "LogP": {"threshold": 5.0, "type": "upper"},
    "HBD": {"threshold": 5, "type": "upper"},
    "HBA": {"threshold": 10, "type": "upper"},
    "TPSA": {"threshold": 140.0, "type": "upper"},
    "RotB": {"threshold": 10, "type": "upper"}
}

def analyze_risk(descriptors):
    """
    Applies ADME rules to descriptors, calculates margins, and determines risk.
    
    Args:
        descriptors (dict): Dictionary of calculated descriptors.
        
    Returns:
        dict: Risk analysis results including status, tier, margins, and limiting rule.
    """
    if not descriptors:
        return None
        
    margins = {}
    violations = []
    
    # Calculate margins and check violations
    # Margin = Threshold - Value (for upper limit)
    # Positive margin = Safe
    # Negative margin = Violation
    
    for prop, rule in RULES.items():
        val = descriptors.get(prop)
        if val is None:
            continue
            
        threshold = rule["threshold"]
        margin = threshold - val
        
        margins[f"Margin_{prop}"] = margin
        
        if margin < 0:
            violations.append(prop)
            
    # Determine Status
    if len(violations) > 0:
        status = "FAIL"
    else:
        status = "SAFE"
        
    # Find limiting rule (smallest margin)
    # We consider all margins to find the 'closest to failure' even if safe
    min_margin = float('inf')
    limiting_rule = None
    
    for prop in RULES.keys():
        m_key = f"Margin_{prop}"
        if m_key in margins:
            if margins[m_key] < min_margin:
                min_margin = margins[m_key]
                limiting_rule = prop
                
    # Determine Tier
    # FAIL
    # SAFE-FRAGILE (min margin < 2.0) - Arbitrary but reasonable heuristic for "close to edge" (esp for LogP/MW scales)
    # SAFE-ROBUST  (min margin >= 2.0)
    # Note: Margin scales are different (MW is hundreds, LogP is units). 
    # A normalized margin would be better scientifically, but user asked for "Threshold - Value".
    # We will stick to the requested simple definition but might need to adjust the "2.0" heuristic 
    # if it doesn't make sense for MW. 
    # For MW: 500 - 499 = 1. This is "Fragile". 500 - 300 = 200. Robust.
    # For LogP: 5 - 4 = 1. Fragile. 5 - 2 = 3. Robust. 
    # The heuristic "2" works surprisingly okay for LogP, HBD, HBA. For MW/TPSA it's very strict (rarely matches).
    # However, user explicitly said: "SAFE-FRAGILE (min margin < 2)". We follow instructions.
    
    if status == "FAIL":
        tier = "FAIL"
    elif min_margin < 2.0:
        tier = "SAFE-FRAGILE"
    else:
        tier = "SAFE-ROBUST"
        
    result = {
        "Status": status,
        "Tier": tier,
        "Min_Margin": min_margin,
        "Limiting_Rule": limiting_rule,
        "Violations_Count": len(violations),
        "Violations_List": ",".join(violations) if violations else "None"
    }
    
    # Merge margins into result
    result.update(margins)
    
    # Calculate MPO Score
    result["MPO_Score"] = calculate_mpo(descriptors)
    
    return result

def calculate_mpo(desc):
    """
    Calculates a simple 6-point MPO score based on Pfizer attributes (approximated).
    Score ranges from 0 (poor) to 6 (excellent).
    Simple discrete scoring for this version:
    - LogP < 3 (1.0), < 5 (0.5)
    - MW < 360 (1.0), < 500 (0.5)
    - TPSA 40-90 (1.0), < 140 (0.5)
    - HBD < 1 (1.0), < 3 (0.5)
    - HBA < 3 (1.0), < 8 (0.5)
    - RotB < 5 (1.0), < 10 (0.5)
    """
    score = 0
    
    # LogP
    logp = desc.get("LogP", 10)
    if logp <= 3: score += 1
    elif logp <= 5: score += 0.5
    
    # MW
    mw = desc.get("MW", 1000)
    if mw <= 360: score += 1
    elif mw <= 500: score += 0.5
    
    # TPSA
    tpsa = desc.get("TPSA", 200)
    if 40 <= tpsa <= 90: score += 1
    elif tpsa <= 140: score += 0.5
    
    # HBD
    hbd = desc.get("HBD", 10)
    if hbd <= 1: score += 1
    elif hbd <= 3: score += 0.5
    
    # HBA
    hba = desc.get("HBA", 20)
    if hba <= 3: score += 1
    elif hba <= 8: score += 0.5
    
    # RotB
    rotb = desc.get("RotB", 20)
    if rotb <= 5: score += 1
    elif rotb <= 10: score += 0.5
    
    return score
