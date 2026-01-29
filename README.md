# MedChem Triage Engine 💊

![Dashboard Preview](medchem_dashboard_preview.png)

A local, secure, and free **Medicinal Chemistry Decision-Support Platform**. 
This tool allows chemists to filter, analyze, and design small molecules using industry-standard ADME constraints, without sending data to the cloud.

## Features
- **Compound Triage**: Automatic filtering based on Lipinski/Veber rules and proprietary "Safety Margins".
- **Deep Analytics**: Violin plots, Risk Landscapes (Rank vs Margin), and Robustness analysis.
- **Cheminformatics**: Substructure search, Similarity search (Tanimoto), and Butina Clustering.
- **Molecular Designer**: Interactive Sketcher (Ketcher) with real-time MPO scoring.
- **Virtual Library**: Generate and screen analogs from a core scaffold instantly.
- **Reporting**: One-click PDF report generation for project summaries.

## Installation

### Option 1: Conda (Recommended)
```bash
conda env create -f environment.yml
conda activate compound-triage
```

### Option 2: Pip
```bash
pip install -r requirements.txt
```

## Usage

1. **Start the Dashboard**:
   ```bash
   streamlit run app/dashboard.py
   ```

2. **Load Data**:
   - The app comes with a `sample.sdf` for testing.
   - Enter the path to your own result CSV, or run the engine first.

3. **Run the Engine (Batch Processing)**:
   ```bash
   python run_engine.py data/raw/your_compounds.sdf data/processed/results.csv
   ```

## Folder Structure
- `app/`: Dashboard code.
- `engine/`: Core logic (Descriptors, Rules, Search, Enumeration).
- `data/`: Input/Output storage.

## 🧪 Case Study: Live PubChem Validation (600 Compounds)
To demonstrate the engine's capability on **Live Data**, we dynamically fetched **600 Investigational Drugs** directly from the **PubChem Database** (via PUG REST API).

**Performance Stats:**
- **Source**: PubChem Live Fetch (`term="pharmaceutical"`)
- **Dataset Size**: 600 Real-World Compounds
- **Processing Time**: ~8 seconds

**Key Findings:**
1.  **Market Drug Profiling**:
    *   The engine correctly prioritized known orally bioavailable drugs into the **SAFE-ROBUST** tier.
    *   It identified specific **Lipinski Violations** in large macrocyclic drugs (e.g., *Cyclosporine*), alignment with literature on "Beyond Rule of 5" compounds.
2.  **Alert Detection**:
    *   Flagged toxicophores in older chemotherapy agents (e.g., alkylating agents) that modern safety rules would reject.
3.  **Chemical Diversity**:
    *   The **Chemical Space Plot (PCA)** successfully mapped the 500 compounds into distinct clusters based on their structural class (Steroids vs. Penicillins vs. Small Aromatics).

**Conclusion**: The MedChem Triage Engine is not just a theoretical tool—it effectively handles messy, real-world chemical data and provides actionable intelligence matching expert medicinal chemistry intuition.

## License
MIT
