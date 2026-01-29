from fpdf import FPDF
import pandas as pd
import os
from datetime import datetime

class ReportGenerator(FPDF):
    def header(self):
        self.set_font('Arial', 'B', 15)
        self.cell(0, 10, 'MedChem Triage Report', 0, 1, 'C')
        self.ln(5)

    def footer(self):
        self.set_y(-15)
        self.set_font('Arial', 'I', 8)
        self.cell(0, 10, f'Page {self.page_no()}', 0, 0, 'C')

def generate_pdf_report(df, filename="report.pdf", figures=None):
    pdf = ReportGenerator()
    pdf.add_page()
    
    # 1. Overview
    pdf.set_font('Arial', 'B', 12)
    pdf.cell(0, 10, f'Project Summary - {datetime.now().strftime("%Y-%m-%d")}', 0, 1)
    
    total = len(df)
    safe = len(df[df["Status"] == "SAFE"])
    fail = total - safe
    
    pdf.set_font('Arial', '', 10)
    pdf.cell(0, 8, f"Total Compounds Processed: {total}", 0, 1)
    pdf.cell(0, 8, f"Safe Candidates: {safe} ({safe/total:.1%})", 0, 1)
    pdf.cell(0, 8, f"Failed Candidates: {fail}", 0, 1)
    pdf.ln(5)

    # 2. Visualizations
    if figures:
        pdf.set_font('Arial', 'B', 12)
        pdf.cell(0, 10, 'Key Visualizations', 0, 1)
        
        import plotly.io as pio
        import tempfile
        
        for title, fig in figures.items():
            try:
                pdf.set_font('Arial', 'B', 10)
                pdf.cell(0, 10, title, 0, 1)
                
                # Create temp file
                with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
                    # Write image (requires kaleido)
                    fig.write_image(tmp.name, width=600, height=400, scale=2)
                    # Embed
                    pdf.image(tmp.name, w=170)
                    tmp_path = tmp.name
                
                # Cleanup
                try:
                    os.unlink(tmp_path)
                except:
                    pass
                
                pdf.ln(5)
            except Exception as e:
                pdf.set_font('Arial', 'I', 8)
                pdf.cell(0, 10, f"Could not render plot {title}: {str(e)}", 0, 1)

    pdf.add_page()
    
    # 3. Top Candidates Table
    pdf.set_font('Arial', 'B', 12)
    pdf.cell(0, 10, 'Top 10 SAFE Candidates', 0, 1)
    
    # Columns to show
    cols = ["Compound_ID", "Tier", "MPO_Score", "MW", "LogP"]
    
    # Header
    pdf.set_font('Arial', 'B', 9)
    for col in cols:
        pdf.cell(38, 8, col, 1)
    pdf.ln()
    
    # Rows
    pdf.set_font('Arial', '', 9)
    safe_df = df[df["Status"] == "SAFE"].sort_values("MPO_Score", ascending=False).head(10)
    
    for _, row in safe_df.iterrows():
        for col in cols:
            val = str(round(row[col], 2)) if isinstance(row[col], float) else str(row[col])
            pdf.cell(38, 8, val, 1)
        pdf.ln()
        
    # 4. Failures
    pdf.ln(10)
    pdf.set_font('Arial', 'B', 12)
    pdf.cell(0, 10, 'Top Failure Rules', 0, 1)
    
    pdf.set_font('Arial', '', 10)
    if "Violations_List" in df.columns:
        violations = []
        for v in df["Violations_List"]:
            if isinstance(v, str) and v != "None": 
                violations.extend(v.split(","))
        
        if violations:
            from collections import Counter
            counts = Counter(violations).most_common(5)
            for rule, count in counts:
                pdf.cell(0, 8, f"- {rule}: {count} violations", 0, 1)
    
    pdf.output(filename)
    return filename
