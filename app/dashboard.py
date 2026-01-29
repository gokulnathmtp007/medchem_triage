import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Add parent directory
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from engine.descriptors import compute_descriptors
from engine.rules import analyze_risk
from engine.alerts import check_structure_alerts
from engine.search import search_similarity, search_substructure
from engine.clustering import cluster_compounds
from engine.enumeration import enumerate_library
from engine.reporting import generate_pdf_report
import importlib
import engine.reporting
importlib.reload(engine.reporting)
from engine.reporting import generate_pdf_report

try:
    from streamlit_ketcher import st_ketcher
    HAS_KETCHER = True
except ImportError:
    HAS_KETCHER = False

# Set page config
st.set_page_config(page_title="MedChem Triage Engine", layout="wide")

# Custom CSS
st.markdown("""
<style>
    body, .stMarkdown, .stText, h1, h2, h3, h4, h5, h6 {
        font-family: "Segoe UI", "Roboto", "Helvetica Neue", "Arial", sans-serif !important;
    }
    .metric-card {
        background-color: #1E1E1E;
        padding: 20px;
        border-radius: 10px;
        border: 1px solid #333;
        text-align: center;
    }
    div[data-testid="stMetricValue"] { font-size: 2rem; }
    .stTable { font-size: 0.9rem; }
    .stButton>button { width: 100%; }
</style>
""", unsafe_allow_html=True)

st.title("Medicinal Chemistry Triage Engine")

# --- Report Figures Collection ---
# We use this dict to collect plots as we generate them
report_figures = {}

# --- Helper Functions ---
@st.cache_data
def load_data(path):
    if not os.path.exists(path): return None
    return pd.read_csv(path)

def calculate_pca(df, features):
    x = df[features].dropna()
    if len(x) < 3: return None
    scaler = StandardScaler()
    x_scaled = scaler.fit_transform(x)
    pca = PCA(n_components=2)
    components = pca.fit_transform(x_scaled)
    pca_df = pd.DataFrame(data=components, columns=['PC1', 'PC2'])
    pca_df = pd.concat([pca_df, df.reset_index(drop=True)], axis=1)
    return pca_df

def create_radar_chart(row, title="Radar Chart"):
    limits = {"MW": 500, "LogP": 5, "TPSA": 140, "HBD": 5, "HBA": 10, "RotB": 10}
    categories = list(limits.keys())
    values = [min(row.get(cat, 0) / limits[cat], 1.5) for cat in categories]
    fig = go.Figure()
    fig.add_trace(go.Scatterpolar(r=values, theta=categories, fill='toself', name=str(row.get("Compound_ID", "Compound"))))
    fig.update_layout(polar=dict(radialaxis=dict(visible=True, range=[0, 1.5], tickvals=[0.5, 1.0, 1.5], ticktext=["Safe", "Limit", "Fail"])), title=title, showlegend=False)
    return fig

def create_margin_plot(row, title="Margin Profile"):
    """Visualizes how close each property is to failure."""
    # Handle dict or Series
    data = row.to_dict() if hasattr(row, "to_dict") else row
    
    margins = {k: v for k, v in data.items() if k.startswith("Margin_")}
    if not margins: return None
    
    df_m = pd.DataFrame(list(margins.items()), columns=["Type", "Value"])
    df_m["Type"] = df_m["Type"].str.replace("Margin_", "")
    # Color: Red if <0 (Fail), Orange if < 2 (Fragile), Green if > 2 (Robust)
    df_m["Color"] = df_m["Value"].apply(lambda x: "#d62728" if x < 0 else ("#ff7f0e" if x < 2 else "#2ca02c"))
    
    fig = px.bar(df_m, x="Value", y="Type", orientation='h', 
                 title=title, text_auto=".1f",
                 color="Color", color_discrete_map="identity")
    fig.add_vline(x=0, line_width=2, line_color="black")
    fig.update_layout(xaxis_title="Margin (Distance to Limit)", yaxis_title="Property")
    return fig

# --- State Management (For Persistence) ---
if "extra_compounds" not in st.session_state:
    st.session_state.extra_compounds = []

# --- Sidebar ---
st.sidebar.header("Data Source")
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)
default_path = os.path.join(project_root, "data", "processed", "results.csv")
input_file = st.sidebar.text_input("Path to Results CSV", value=default_path)

if st.sidebar.button("🔄 Reload Data"):
    st.cache_data.clear()
    st.rerun()

df_main = pd.DataFrame()
if input_file and os.path.exists(input_file):
    df_main = load_data(input_file)
    if df_main is None: st.error("Failed to load data"); st.stop()
    for col, default in [("Scaffold", "Unknown"), ("MPO_Score", 0.0), ("Alerts", "None")]:
        if col not in df_main.columns: df_main[col] = default

# Combine with session extra compounds
if st.session_state.extra_compounds:
    df_extra = pd.DataFrame(st.session_state.extra_compounds)
    df_main = pd.concat([df_main, df_extra], ignore_index=True)

if df_main.empty: st.info("Load data to proceed."); st.stop()

# Re-calc ID list for dropdowns
compound_ids = df_main["Compound_ID"].unique()
COLOR_MAP = {"SAFE-ROBUST": "#00CC96", "SAFE-FRAGILE": "#FFA15A", "FAIL": "#EF553B"}

# --- PDF Reporting Setup (Button Logic at End) ---
st.sidebar.divider()
st.sidebar.header("Reporting")
# We just show the button here, but process it later to ensure all plots are ready
gen_report_btn = st.sidebar.button("Generate Full PDF Report")

# --- Tabs ---
tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
    "Dashboard", "Analytics", "Explorer", "Comparator", "Search & Cluster", "Designer", "Library Gen"
])

# --- TAB 1: DASHBOARD ---
with tab1:
    st.subheader("Portfolio Health")
    total_cmpds = len(df_main)
    safe_cmpds = df_main[df_main["Status"] == "SAFE"].shape[0]
    avg_mpo = df_main["MPO_Score"].mean()
    alerts_count = df_main[df_main["Alerts"] != "None"].shape[0]
    
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Total Compounds", total_cmpds)
    c2.metric("SAFE", safe_cmpds, delta=f"{safe_cmpds/total_cmpds:.1%}" if total_cmpds else "0%")
    c3.metric("Avg MPO Score", f"{avg_mpo:.1f}/6.0")
    c4.metric("Structural Alerts", alerts_count, delta="Potential PAINS", delta_color="inverse")
    
    st.markdown("---")
    
    # Attrition Funnel
    st.markdown("### Triage Funnel")
    funnel_data = {
        "Stage": ["Total", "Clean (No Alerts)", "Lipinski Safe", "Robust (MPO>4)"],
        "Count": [
            total_cmpds,
            len(df_main[df_main["Alerts"] == "None"]),
            safe_cmpds,
            len(df_main[df_main["MPO_Score"] >= 4.0])
        ]
    }
    fig_funnel = px.funnel(funnel_data, x='Count', y='Stage', title="Attrition Funnel")
    st.plotly_chart(fig_funnel, use_container_width=True)
    report_figures["Attrition Funnel"] = fig_funnel
    
    st.markdown("---")
    r1c1, r1c2 = st.columns([1, 2])
    with r1c1:
        st.markdown("### Status Breakdown")
        status_counts = df_main["Tier"].value_counts().reset_index()
        status_counts.columns = ["Tier", "Count"]
        fig_pie = px.pie(status_counts, values="Count", names="Tier", color="Tier",
            color_discrete_map=COLOR_MAP, hole=0.4)
        st.plotly_chart(fig_pie, use_container_width=True)
        report_figures["Status Breakdown"] = fig_pie
        
    with r1c2:
        st.markdown("### Top Failures")
        all_violations = []
        for v_list in df_main["Violations_List"]:
            if pd.notna(v_list) and isinstance(v_list, str) and v_list != "None":
                all_violations.extend(v_list.split(","))
        if all_violations:
            v_counts = pd.Series(all_violations).value_counts().reset_index()
            v_counts.columns = ["Rule", "Count"]
            fig_bar = px.bar(v_counts, x="Count", y="Rule", orientation='h', color="Count", color_continuous_scale="Viridis")
            st.plotly_chart(fig_bar, use_container_width=True)
            report_figures["Top Failures"] = fig_bar

    st.markdown("---")
    dist_col, space_col = st.columns(2)
    with dist_col:
        st.markdown("### Risk Landscape")
        fig_scatter = px.scatter(df_main, x="Rank", y="Min_Margin", color="Tier",
            hover_data=["Compound_ID", "Limiting_Rule", "MPO_Score"],
            color_discrete_map=COLOR_MAP)
        fig_scatter.add_hline(y=0, line_dash="dash", line_color="red")
        st.plotly_chart(fig_scatter, use_container_width=True)
        report_figures["Risk Landscape"] = fig_scatter
        
    with space_col:
        # Robustness
        st.markdown("### Robustness (SAFE Compounds)")
        df_safe = df_main[df_main["Status"] == "SAFE"]
        if not df_safe.empty:
            fig_rob = px.histogram(df_safe, x="Min_Margin", nbins=20, 
                                color_discrete_sequence=["#00CC96"])
            fig_rob.add_vrect(x0=0, x1=2.0, fillcolor="#FFA15A", opacity=0.3, annotation_text="Fragile")
            st.plotly_chart(fig_rob, use_container_width=True)
            report_figures["Robustness"] = fig_rob
        else:
            st.info("No Safe Compounds")

# --- TAB 2: ANALYTICS ---
with tab2:
    st.subheader("Deep Analytics")
    props = ["MW", "LogP", "TPSA", "HBD", "HBA", "RotB"]
    sel_prop = st.selectbox("Select Property", props)
    
    h_col1, h_col2 = st.columns([2, 1])
    with h_col1:
        st.markdown("**Property Distribution (Violin - Density)**")
        fig_viol = px.violin(df_main, x="Status", y=sel_prop, color="Status", box=True, points="all",
                             color_discrete_map={"SAFE": "#00CC96", "FAIL": "#EF553B"})
        st.plotly_chart(fig_viol, use_container_width=True)
        # We don't report interactive sel_prop plots usually, but we could add if needed
        
    with h_col2:
        st.markdown("**Box Plot by Tier**")
        fig_box = px.box(df_main, x="Tier", y=sel_prop, color="Tier",
             color_discrete_map=COLOR_MAP)
        st.plotly_chart(fig_box, use_container_width=True)
        
    st.markdown("---")
    st.markdown("### Descriptor Correlations")
    corr_matrix = df_main[props].corr()
    fig_heat = px.imshow(corr_matrix, text_auto=True, color_continuous_scale="RdBu_r", zmin=-1, zmax=1)
    st.plotly_chart(fig_heat, use_container_width=True)
    report_figures["Correlation Matrix"] = fig_heat

# --- TAB 3: EXPLORER ---
with tab3:
    st.subheader("Data Explorer")
    if "Scaffold" in df_main.columns:
        st.markdown("#### Series Analysis")
        scaffold_stats = df_main.groupby("Scaffold").agg(
            Count=("Compound_ID", "count"),
            Avg_MPO=("MPO_Score", "mean"),
            Safe_Count=("Status", lambda x: (x == "SAFE").sum())
        ).reset_index()
        scaffold_stats["% Safe"] = scaffold_stats["Safe_Count"] / scaffold_stats["Count"]
        scaffold_stats = scaffold_stats.sort_values("Count", ascending=False)
        st.dataframe(scaffold_stats, column_config={
            "% Safe": st.column_config.ProgressColumn(format="%.0%"),
            "Scaffold": st.column_config.TextColumn("Murcko Scaffold SMILES")
        }, hide_index=True, use_container_width=True)
    st.markdown("#### Compound Table")
    st.dataframe(df_main, use_container_width=True)

# --- TAB 4: COMPARATOR ---
with tab4:
    st.subheader("Side-by-Side Comparison")
    c_s1, c_s2 = st.columns(2)
    id1 = c_s1.selectbox("Compound A", compound_ids, index=0)
    id2 = c_s2.selectbox("Compound B", compound_ids, index=1 if len(compound_ids)>1 else 0)
    
    row1 = df_main[df_main["Compound_ID"] == id1].iloc[0]
    row2 = df_main[df_main["Compound_ID"] == id2].iloc[0]
    
    colA, colB = st.columns(2)
    with colA:
        st.markdown(f"**{id1}**")
        st.metric("MPO", f"{row1['MPO_Score']:.1f}")
        st.plotly_chart(create_radar_chart(row1, "Shape"), use_container_width=True)
        m_fig = create_margin_plot(row1, "Margin Profile")
        if m_fig: st.plotly_chart(m_fig, use_container_width=True)
        st.json(row1.to_dict())
    with colB:
        st.markdown(f"**{id2}**")
        st.metric("MPO", f"{row2['MPO_Score']:.1f}")
        st.plotly_chart(create_radar_chart(row2, "Shape"), use_container_width=True)
        m_fig2 = create_margin_plot(row2, "Margin Profile")
        if m_fig2: st.plotly_chart(m_fig2, use_container_width=True)
        st.json(row2.to_dict())

# --- TAB 5: SEARCH & CLUSTER ---
with tab5:
    st.subheader("Advanced Search & Clustering")
    search_col, cluster_col = st.columns(2)
    with search_col:
        st.markdown("#### Search")
        search_mode = st.radio("Mode", ["Similarity", "Substructure"])
        query_smi = st.text_input("Query SMILES", "c1ccccc1")
        search_btn = st.button("Run Search")
        if search_btn and query_smi:
            qmol = Chem.MolFromSmiles(query_smi)
            if qmol:
                targets = []
                for _, r in df_main.iterrows():
                    targets.append((r["Compound_ID"], Chem.MolFromSmiles(r["SMILES"])))
                if search_mode == "Similarity":
                    hits = search_similarity(qmol, targets)
                    st.write(f"Found {len(hits)} hits (>0.7)")
                    st.dataframe(pd.DataFrame(hits, columns=["ID", "Score", "Mol"]).drop(columns="Mol"))
                else:
                    hits = search_substructure(qmol, targets)
                    st.write(f"Found {len(hits)} hits")
                    st.dataframe(pd.DataFrame(hits, columns=["ID", "Match", "Mol"]).drop(columns=["Match", "Mol"]))
            else:
                st.error("Invalid Query SMILES")
    with cluster_col:
        st.markdown("#### Portfolio Clustering")
        cluster_btn = st.button("Cluster All Compounds")
        if cluster_btn:
            with st.spinner("Clustering..."):
                targets = []
                for _, r in df_main.iterrows():
                    targets.append((r["Compound_ID"], Chem.MolFromSmiles(r["SMILES"])))
                cmap, sorted_clusters = cluster_compounds(targets, cutoff=0.7)
                st.success(f"Found {len(sorted_clusters)} clusters.")
                for c in sorted_clusters[:5]:
                    st.markdown(f"**Cluster {c['Cluster_ID']}** (Size: {c['Size']}) - Centroid: {c['Centroid_ID']}")
                    with st.expander("View Members"):
                        st.write(", ".join(c["Members"]))

# --- TAB 6: DESIGNER ---
with tab6:
    st.subheader("Molecular Designer")
    
    st.markdown("Draw your molecule below or paste SMILES:")
    
    # Default SMILES
    default_smi = "CC(=O)OC1=CC=CC=C1C(=O)O"
    
    if HAS_KETCHER:
        # Ketcher Sketcher
        smi_in = st_ketcher(default_smi)
    else:
        st.warning("⚠️ Molecular Sketcher not detected. (Run `pip install streamlit-ketcher` to enable). Using text input instead.")
        smi_in = st.text_area("Input SMILES", default_smi, key="ketcher_fallback")
    
    col_res, col_act = st.columns(2)
    
    if smi_in:
        mol = Chem.MolFromSmiles(smi_in)
        if mol:
            desc = compute_descriptors(mol)
            risk = analyze_risk(desc)
            alerts = check_structure_alerts(mol)
            res = {**desc, **risk}
            res["Compound_ID"] = f"Design_{len(st.session_state.extra_compounds)+1}"
            res["SMILES"] = smi_in
            res["Alerts"] = ",".join(alerts) if alerts else "None"
            res["Violations_List"] = risk["Violations_List"]
            
            with col_res:
                st.markdown(f"### Results: {risk['Tier']}")
                st.metric("MPO Score", f"{risk['MPO_Score']:.1f}")
                
                m_figd = create_margin_plot(res, "Margin Profile")
                if m_figd: st.plotly_chart(m_figd, use_container_width=True)
                
            with col_act:
                st.plotly_chart(create_radar_chart(res, "Shape"), use_container_width=True)
                if st.button("Add to Portfolio", key="add_des"):
                    st.session_state.extra_compounds.append(res)
                    st.success("Added!")

# --- TAB 7: LIBRARY GENERATOR ---
with tab7:
    st.subheader("Virtual Library Enumeration")
    st.markdown("Define a **Core Scaffold** (with `*`) and a list of **R-Groups** to generate new analogs.")
    
    col_core, col_r = st.columns(2)
    
    with col_core:
        core_smi = st.text_input("Core Scaffold (SMILES with *)", "c1ccccc1[*]")
        st.caption("Example: `c1ccccc1[*]` (Phenyl radical)")
        
    with col_r:
        r_group_txt = st.text_area("R-Groups (One SMILES per line)", "C\nCl\nBr\nO\nN")
        
    if st.button("Generate & Screen Library"):
        r_list = [r.strip() for r in r_group_txt.split("\n") if r.strip()]
        
        with st.spinner(f"Generating analogs from {len(r_list)} R-groups..."):
            products = enumerate_library(core_smi, r_list)
            
            if products:
                st.success(f"Generated {len(products)} unique analogs!")
                # Screen them!
                results = []
                for p_smi in products:
                    mol = Chem.MolFromSmiles(p_smi)
                    if mol:
                        desc = compute_descriptors(mol)
                        risk = analyze_risk(desc)
                        res = {"SMILES": p_smi, "Status": risk["Status"], "MPO": risk["MPO_Score"]}
                        results.append(res)
                res_df = pd.DataFrame(results).sort_values("MPO", ascending=False)
                st.dataframe(res_df, use_container_width=True)
                safe_products = res_df[res_df["Status"] == "SAFE"]
                st.info(f"{len(safe_products)} SAFE Analogs found.")
            else:
                st.error("No valid products generated. Check your Core/R-group syntax.")

# --- REPORT GENERATION EXECUTION ---
if gen_report_btn:
    with st.sidebar:
        with st.spinner("Compiling PDF with Figures..."):
            try:
                # Pass the collected figures dict
                pdf_path = generate_pdf_report(df_main, "portfolio_report.pdf", figures=report_figures)
                with open(pdf_path, "rb") as f:
                    st.download_button("Download Full PDF", f, file_name="MedChem_Full_Report.pdf", mime="application/pdf")
                st.success("Report Ready!")
            except Exception as e:
                st.error(f"Failed to generate report: {e}")
