import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go

# --- Page Configuration ---
st.set_page_config(page_title="VAF-TC Precision Analyzer", layout="wide")

st.title("🧬 VAF-TC Clinical Genetics Analyzer")
st.markdown("""
This decision-support tool applies **Knudson’s two-hit model** to differentiate between somatic and germline variants by analyzing the relationship between **Pathological Tumor Content (TC)** and **Variant Allele Frequency (VAF)**.
""")

# --- Sidebar Inputs ---
st.sidebar.header("📋 Patient Data Input")
tc_input = st.sidebar.slider("Pathological Tumor Content (%)", 10, 100, 50, help="Assessment by a pathologist is the gold standard.")
vaf_input = st.sidebar.slider("Observed VAF (%)", 1, 100, 25)
st.sidebar.info("Note: 'Pathological TC' is prioritized to maintain independent biological validation.")

# --- Mathematical Models (f = Tumor Fraction) ---
f = tc_input / 100
models = {
    "Somatic Heterozygous": (f / 2) * 100,
    "Somatic LOH (with Del)": (f / (2 - f)) * 100,
    "Germline Heterozygous": 50.0,
    "Germline LOH (with Del)": (1 / (2 - f)) * 100
}

# --- 1. Dynamic Clinical Alerts ---
st.subheader("🚨 Real-time Clinical Alerts")

# Point 1: The 50% VAF Trap (TC 60-75%)
if 60 <= tc_input <= 75:
    v_loh_trap = models["Somatic LOH (with Del)"]
    st.warning(f"""
    **Alert: The 50% VAF Trap**
    At TC {tc_input}%, the theoretical VAF for Somatic LOH is {v_loh_trap:.1f}%. Since this crosses the 50% threshold at TC ≈ 66.7%, 
    somatic drivers can mimic heterozygous germline variants. Use caution before assuming germline status.
    """)

# Point 2: Convergence Alert (TC >= 70%) - Corrected KeyError (Removed 'Marc')
if tc_input >= 70 and vaf_input >= models["Somatic LOH (with Del)"]:
    val_threshold = models["Somatic LOH (with Del)"]
    st.error(f"""
    **⚠️ LOH Convergence Alert (TC ≥ 70%):**
    Potential convergence of Germline LOH and Somatic LOH detected. 
    When TC is ≥ 70% and the VAF is at or above the theoretical line for Somatic LOH with deletion ({val_threshold:.1f}%), 
    these two biological events become increasingly difficult to distinguish based on VAF alone.
    """)

# Point 3: Mathematical Limit (TC >= 90%)
if tc_input >= 90:
    st.info("""
    **💡 Mathematical Convergence Zone (TC ≥ 90%):**
    In high-purity samples, Somatic and Germline LOH models mathematically converge toward 100%. 
    In this range, VAF is insufficient to distinguish the variant's origin (e.g., *TP53* somatic LOH vs. germline LOH). 
    Family history and clinical correlation are essential.
    """)

# --- 2. Compatible Models (±10% Range) ---
st.subheader("🔍 Compatible Theoretical Models")
compatible_data = []
for name, theory_vaf in models.items():
    diff = abs(vaf_input - theory_vaf)
    if diff <= 10.0:
        compatible_data.append({"Model Name": name, "Theoretical VAF": f"{theory_vaf:.2f}%", "Difference": f"{diff:.2f}%"})

if compatible_data:
    st.table(pd.DataFrame(compatible_data))
else:
    st.error("No standard models match within ±10%. Consider clonal heterogeneity or aneuploidy.")

# --- 3. Interpretation & Clinical Guidance ---
with st.expander("📝 Clinical Significance & Therapeutic Implications"):
    st.markdown("""
    ### Hereditary Cancer Inference
    Syndromes driven by tumor suppressor genes (**HBOC, Lynch Syndrome, FAP**) can be inferred when VAFs align with theoretical "two-hit" lines based on **Knudson’s two-hit model**.
    
    ### Therapeutic Guidance
    - **BRCA1/2-associated tumors:** May indicate sensitivity to **PARP inhibitors**.
    - **Lynch Syndrome (MMR-d):** May predict responsiveness to **Immune Checkpoint Inhibitors (ICIs)**. In Lynch syndrome, curative potential with ICIs is a significant clinical observation.
    
    | Organ | gBRCA/sBRCA Indication | Drug Class |
    | :--- | :--- | :--- |
    | Ovarian/Prostate | Both gBRCA and sBRCA | PARPi |
    | Breast/Pancreas | **gBRCA Only** | PARPi |
    """)

# --- 4. Visualization ---
st.subheader("📈 VAF-TC Projection")
tr = np.linspace(10, 100, 100)
fr = tr / 100
fig = go.Figure()
fig.add_trace(go.Scatter(x=tr, y=(fr/2)*100, name="Somatic Het", line=dict(color='blue', width=1)))
fig.add_trace(go.Scatter(x=tr, y=(fr/(2-fr))*100, name="Somatic LOH (Del)", line=dict(color='blue', dash='dash')))
fig.add_trace(go.Scatter(x=tr, y=[50]*100, name="Germline Het", line=dict(color='green', width=1)))
fig.add_trace(go.Scatter(x=tr, y=(1/(2-fr))*100, name="Germline LOH (Del)", line=dict(color='red', width=2)))
fig.add_trace(go.Scatter(x=[tc_input], y=[vaf_input], mode='markers+text', name="CASE", text=["CASE"], marker=dict(color='black', size=15, symbol='x')))

fig.update_layout(xaxis_title="Pathological Tumor Content (%)", yaxis_title="Theoretical VAF (%)",
                  legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))
st.plotly_chart(fig, use_container_width=True)

st.divider()
st.caption("Developed by Clinical Genetics Suite (Maintainer: Sawai1960). Version 2.4 (Hotfix: KeyError Resolved).")
