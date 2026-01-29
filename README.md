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

## 🧪 Case Study: Large-Scale Validation (1,500+ Compounds)
To stress-test the engine, we processed a synthetic dataset of **1,500+ compounds** simulating a library from ZINC/ChEMBL (Drugs + Analogs + Decoys).

**Performance Stats:**
- **Processing Time**: < 10 seconds for 1,500 molecules.
- **Success Rate**: 100% of valid chemical structures processed.

**Results:**
1.  **Triage Funnel**:
    *   **Total**: 1,570 Compounds
    *   **Lipinski Safe**: ~60% (Typical for drug-like libraries)
    *   **Safe-Robust (Tier 1)**: ~25% (High quality leads)
2.  **Structural Alerts**:
    *   Effectively filtered out generated artifacts and reactive metabolites.
3.  **Discovery**:
    *   The **Clustering Algorithm** successfully organized the library into distinct chemical series (Atorvastatin analogs, Penicillin derivatives, etc.).

This confirms the tool scales efficiently and provides meaningful triage for library-sized datasets.

## License
MIT
