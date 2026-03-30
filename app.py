import streamlit as st
import plotly.graph_objects as go
import numpy as np

# Page Configuration
st.set_page_config(page_title="VAF-TC Visualizer", layout="wide")

st.title("🧬 VAF-TC Relationship Visualizer")
st.write("An interactive tool to visualize the theoretical relationship between Pathological Tumor Content (TC) and Variant Allele Fraction (VAF).")

# --- Sidebar Parameter Inputs ---
st.sidebar.header("📊 Input Parameters")
gene_name = st.sidebar.text_input("Gene Name", value="BRCA2")
tc_input = st.sidebar.slider("Pathological Tumor Content (TC %)", 0, 100, 50)
vaf_input = st.sidebar.slider("Variant Allele Fraction (VAF %)", 0, 100, 50)

# --- Mathematical Model Calculations ---
def calculate_curves():
    tc_range = np.linspace(0, 1, 101)
    # Germline
    g_hetero = np.full_like(tc_range, 0.5)
    g_loh_del = 1 / (2 - tc_range)
    g_cnloh = 0.5 * (1 + tc_range)
    # Somatic
    s_hetero = 0.5 * tc_range
    s_loh_del = tc_range / (2 - tc_range)
    s_cnloh = tc_range
    return tc_range * 100, g_hetero * 100, g_loh_del * 100, g_cnloh * 100, s_hetero * 100, s_loh_del * 100, s_cnloh * 100

tc_plot, g_het, g_del, g_cn, s_het, s_del, s_cn = calculate_curves()

# --- Visualization ---
fig = go.Figure()

fig.add_trace(go.Scatter(x=tc_plot, y=g_cn, name="Germline + cnLOH", line=dict(color='#D4AF37', width=2.5)))
fig.add_trace(go.Scatter(x=tc_plot, y=g_del, name="Germline + LOH (Del)", line=dict(color='red', width=2.5)))
fig.add_trace(go.Scatter(x=tc_plot, y=g_het, name="Germline (Hetero)", line=dict(color='brown', width=2.5)))
fig.add_trace(go.Scatter(x=tc_plot, y=s_cn, name="Somatic + cnLOH", line=dict(color='green', dash='dash')))
fig.add_trace(go.Scatter(x=tc_plot, y=s_del, name="Somatic + LOH (Del)", line=dict(color='#666', dash='dot')))

fig.add_vrect(x0=0, x1=30, fillcolor="rgba(200, 200, 200, 0.2)", layer="below", line_width=0, 
              annotation_text="Low Confidence Zone", annotation_position="top left")

fig.add_trace(go.Scatter(x=[tc_input], y=[vaf_input], name=f"{gene_name} Sample", mode='markers+text',
                         marker=dict(color='black', size=14, symbol='circle'),
                         text=[f"{gene_name}<br>TC:{tc_input}%<br>VAF:{vaf_input}%"], textposition="top right"))

fig.update_layout(
    xaxis=dict(title="Pathological Tumor Content (%)", range=[0, 100], dtick=10, gridcolor='#eee'),
    yaxis=dict(title="Variant Allele Fraction (%)", range=[0, 100], dtick=25, gridcolor='#eee'),
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    plot_bgcolor='white', margin=dict(l=40, r=40, t=80, b=40), height=650
)

st.plotly_chart(fig, use_container_width=True)

# --- Calculations for Logic ---
tc_val = tc_input / 100
# Calculate the threshold: Somatic LOH (del) curve value at current TC
s_loh_del_threshold = (tc_val / (2 - tc_val)) * 100 if tc_val < 2 else 100

# --- Alert & Interpretation Section ---
st.markdown("---")

# 1. Refined Convergence Alert (TC >= 70% AND VAF >= Somatic LOH line)
if tc_input >= 70 and vaf_input >= s_loh_del_threshold:
    st.warning(f"""
    ⚠️ **Convergence Risk (Gray Zone)**:  
    At TC **{tc_input}%** and VAF **{vaf_input}%**, theoretical curves for **Germline LOH** and **Somatic LOH** converge significantly. 
    Distinguishing between these events based on VAF alone is difficult in this range. 
    Clinical correlation (e.g., family history, drug response) is strongly recommended.
    """)
elif tc_input < 30:
    st.info("ℹ️ **Low Confidence Zone**: Interpretation reliability may be limited when TC is below 30%.")

# 2. Automated Interpretation matches
theoretical_vafs = {
    "Germline (Hetero)": 50.0,
    "Somatic (Hetero)": (0.5 * tc_val) * 100,
    "Somatic + cnLOH": (tc_val) * 100,
    "Somatic + LOH (Del)": s_loh_del_threshold,
    "Germline + LOH (Del)": (1 / (2 - tc_val)) * 100 if tc_val < 2 else 100,
    "Germline + cnLOH": (0.5 * (1 + tc_val)) * 100
}

matched_models = [m for m, v in theoretical_vafs.items() if abs(vaf_input - v) <= 2.0]

if matched_models:
    st.success(f"**Interpretation for {gene_name}:** \n" + 
               f"Observed VAF aligns with: {', '.join(matched_models)}")
else:
    st.info(f"**Interpretation for {gene_name}:** \n" +
            "VAF does not closely align with standard models. Consider clonal heterogeneity or complex CNAs.")

# --- Static Clinical Notes ---
st.markdown("---")
st.subheader("📝 Clinical Interpretation Notes")
st.write(f"""
- **Convergence Zone**: The "Gray Zone" occurs because the mathematical difference between somatic and germline VAFs narrows as tumor purity increases.
- **High-TC Case Insights**: As demonstrated in cases with TC $\ge$ 90%, variants with high VAFs can be misidentified as somatic. In our study, these were confirmed as **Germline LOH**, supported by clinical outcomes such as favorable responses to **PARP inhibitors**.
""")
