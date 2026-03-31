import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go

# --- Page Configuration ---
st.set_page_config(page_title="VAF-TC Precision Analyzer", layout="wide")

st.title("🧬 VAF-TC Clinical Genetics Analyzer")
st.markdown("""
This application analyzes the relationship between **Pathological Tumor Content (TC)** and **Variant Allele Frequency (VAF)** to distinguish between Somatic and Germline variants, considering LOH (Loss of Heterozygosity).
""")

# --- Sidebar Inputs ---
st.sidebar.header("📋 Patient Data Input")
tc_input = st.sidebar.slider("Pathological Tumor Content (%)", 10, 100, 50, help="Assessment by a pathologist is highly recommended.")
vaf_input = st.sidebar.slider("Observed VAF (%)", 1, 100, 25)
st.sidebar.info("Note: 'Pathological TC' is used as the gold standard for this model.")

# --- Mathematical Models ---
# f = TC / 100
f = tc_input / 100

models = {
    "Somatic Heterozygous": (f / 2) * 100,
    "Somatic LOH (with Del)": (f / (2 - f)) * 100,
    "Germline Heterozygous": 50.0,
    "Germline LOH (with Del)": (1 / (2 - f)) * 100
}

# --- 1. Dynamic Alert Logic (The 50% VAF Trap) ---
st.subheader("🚨 Real-time Clinical Alerts")

if 60 <= tc_input <= 75:
    st.warning(f"""
    **Alert: The 50% VAF Trap (Grey Zone)**
    At TC {tc_input}%, the theoretical VAF for **Somatic LOH (with Del)** is {models['Somatic LOH (with Del) Marc']:.1f}%.
    Crucially, at TC ≈ 66.7%, this somatic model crosses the **50% threshold**.
    In this range, a somatic mutation can perfectly mimic a heterozygous germline variant. 
    **Recommendation:** Do not assume VAF ≈ 50% is Germline. Consider paired-normal testing.
    """)
elif tc_input >= 90:
    st.info("""
    **High TC Context (≥90%):**
    At very high tumor purity, Somatic LOH and Germline variants converge toward similar VAF ranges. 
    While germline LOH is statistically frequent in certain contexts, somatic LOH (e.g., in *TP53*) 
    is equally plausible. High TC requires careful integration with clinical history.
    """)
else:
    st.success("No critical mathematical intersections detected for the current TC range.")

# --- 2. Compatible Models (±10% Range) ---
st.subheader("🔍 Compatible Theoretical Models")
st.write(f"Listing all models within a ±10% margin of the observed VAF ({vaf_input}%):")

compatible_data = []
for name, theory_vaf in models.items():
    diff = abs(vaf_input - theory_vaf)
    if diff <= 10.0:
        compatible_data.append({
            "Model Name": name,
            "Theoretical VAF (%)": round(theory_vaf, 2),
            "Difference (%)": round(diff, 2)
        })

if compatible_data:
    st.table(pd.DataFrame(compatible_data))
else:
    st.error("No standard models match the observed VAF within a ±10% margin. Consider clonal heterogeneity or aneuploidy.")

# --- 3. Interpretation & Factors ---
with st.expander("📝 Clinical Interpretation Notes"):
    st.markdown(f"""
    **Interpretation Factors:**
    - **NGS Variance:** Sequencing depth and library prep can cause ±5-10% fluctuations.
    - **Aneuploidy / Copy Number Changes:** Deviations from these models often suggest large-scale genomic gains or losses.
    - **Clonal Heterogeneity:** Subclonal mutations will present with lower-than-expected VAF.
    
    **PARP Inhibitor (PARPi) Indications:**
    - **Ovarian & Prostate Cancer:** Both **Germline (gBRCA)** and **Somatic (sBRCA)** variants with LOH (Biallelic inactivation) may confer sensitivity.
    - **Breast & Pancreatic Cancer:** Currently, regulatory approval is primarily restricted to **Germline (gBRCA)** carriers. 
    *Note: Always refer to the latest regional clinical guidelines.*
    """)

# --- 4. Visualization ---
st.subheader("📈 VAF-TC Theoretical Projection")
tc_range = np.linspace(10, 100, 100)
f_range = tc_range / 100

fig = go.Figure()
fig.add_trace(go.Scatter(x=tc_range, y=(f_range/2)*100, name="Somatic Het"))
fig.add_trace(go.Scatter(x=tc_range, y=(f_range/(2-f_range))*100, name="Somatic LOH (Del)", line=dict(dash='dash')))
fig.add_trace(go.Scatter(x=tc_range, y=[50]*100, name="Germline Het", line=dict(color='green')))
fig.add_trace(go.Scatter(x=tc_range, y=(1/(2-f_range))*100, name="Germline LOH (Del)", line=dict(color='red')))

# Highlight User Data
fig.add_trace(go.Scatter(x=[tc_input], y=[vaf_input], mode='markers+text', 
                         name="Current Case", text=["Case"], textposition="top center",
                         marker=dict(color='black', size=12, symbol='x')))

fig.update_layout(xaxis_title="Pathological Tumor Content (%)", yaxis_title="VAF (%)",
                  legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))
st.plotly_chart(fig, use_container_width=True)

st.divider()
st.caption("Developed for Clinical Genetics Suite. Version 2.0 (Internal Preview).")
