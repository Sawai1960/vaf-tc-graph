import streamlit as st
import plotly.graph_objects as go
import numpy as np

# Page Configuration for Academic Presentation
st.set_page_config(page_title="VAF-TC Visualizer", layout="wide")

st.title("🧬 VAF-TC Relationship Visualizer")
st.write("An interactive tool to visualize the theoretical relationship between Pathological Tumor Content (TC) and Variant Allele Fraction (VAF).")

# --- Sidebar Parameter Inputs ---
st.sidebar.header("📊 Input Parameters")

# Feature 1: Gene Name Input (Restored)
gene_name = st.sidebar.text_input("Gene Name", value="BRCA2", help="e.g., BRCA2, BRCA1, MSH6, APC")

# Sliders for TC and VAF
tc_input = st.sidebar.slider("Pathological Tumor Content (TC %)", 0, 100, 50)
vaf_input = st.sidebar.slider("Variant Allele Fraction (VAF %)", 0, 100, 50)

# --- Mathematical Model Calculations ---
def calculate_curves():
    tc_range = np.linspace(0, 1, 101)
    
    # Germline Models (Knudson's Two-Hit Theory)
    g_hetero = np.full_like(tc_range, 0.5)
    g_loh_del = 1 / (2 - tc_range)
    g_cnloh = 0.5 * (1 + tc_range)
    
    # Somatic Models
    s_hetero = 0.5 * tc_range
    s_loh_del = tc_range / (2 - tc_range)
    s_cnloh = tc_range
    
    return tc_range * 100, g_hetero * 100, g_loh_del * 100, g_cnloh * 100, s_hetero * 100, s_loh_del * 100, s_cnloh * 100

tc_plot, g_het, g_del, g_cn, s_het, s_del, s_cn = calculate_curves()

# --- Visualization Construction ---
fig = go.Figure()

# Plotting Theoretical Trajectories
fig.add_trace(go.Scatter(x=tc_plot, y=g_cn, name="Germline + cnLOH", line=dict(color='#D4AF37', width=2.5)))
fig.add_trace(go.Scatter(x=tc_plot, y=g_del, name="Germline + LOH (Del)", line=dict(color='red', width=2.5)))
fig.add_trace(go.Scatter(x=tc_plot, y=g_het, name="Germline (Hetero)", line=dict(color='brown', width=2.5)))
fig.add_trace(go.Scatter(x=tc_plot, y=s_cn, name="Somatic + cnLOH", line=dict(color='green', dash='dash')))
fig.add_trace(go.Scatter(x=tc_plot, y=s_del, name="Somatic + LOH (Del)", line=dict(color='#666', dash='dot')))

# Highlighting the Low Confidence Zone (TC < 30%)
fig.add_vrect(x0=0, x1=30, fillcolor="rgba(200, 200, 200, 0.2)", layer="below", line_width=0, 
              annotation_text="Low Confidence Zone", annotation_position="top left")

# Feature 2: User-Defined Data Point (Updated with Gene Name)
fig.add_trace(go.Scatter(x=[tc_input], y=[vaf_input], name=f"{gene_name} Sample", mode='markers+text',
                         marker=dict(color='black', size=14, symbol='circle'),
                         text=[f"{gene_name}<br>TC:{tc_input}%<br>VAF:{vaf_input}%"], textposition="top right"))

# Axis Formatting (Fixed 0-100% range)
fig.update_layout(
    xaxis=dict(title="Pathological Tumor Content (%)", range=[0, 100], dtick=10, gridcolor='#eee'),
    yaxis=dict(title="Variant Allele Fraction (%)", range=[0, 100], dtick=25, gridcolor='#eee'),
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    plot_bgcolor='white',
    margin=dict(l=40, r=40, t=80, b=40),
    height=650
)

st.plotly_chart(fig, use_container_width=True)

# --- Dynamic Automated Interpretation (Restored & English) ---
st.markdown("---")
st.subheader("🔍 Automated Interpretation (Trial)")

# Convert input percentages to decimal for calculation
tc_val = tc_input / 100
vaf_val = vaf_input / 100
tolerance = 2.0  # Allowed deviation in VAF % for matching

# Theoretical VAFs at current TC
theoretical_vafs = {
    "Germline (Hetero)": 50.0,
    "Somatic (Hetero)": (0.5 * tc_val) * 100,
    "Somatic + cnLOH": (tc_val) * 100
}

# Avoid division by zero
if tc_val < 1.0:
    theoretical_vafs["Somatic + LOH (Del)"] = (tc_val / (2 - tc_val)) * 100
    theoretical_vafs["Germline + LOH (Del)"] = (1 / (2 - tc_val)) * 100
    theoretical_vafs["Germline + cnLOH"] = (0.5 * (1 + tc_val)) * 100
elif tc_val == 1.0:
    theoretical_vafs["Somatic + LOH (Del)"] = 100.0
    theoretical_vafs["Germline + LOH (Del)"] = 100.0
    theoretical_vafs["Germline + cnLOH"] = 100.0

# Identify matches within tolerance
matched_models = []
for model_name, theoretical_vaf in theoretical_vafs.items():
    if abs(vaf_input - theoretical_vaf) <= tolerance:
        matched_models.append(model_name)

# Display Interpretation Comments
if matched_models:
    st.success(f"**Automated Comment for {gene_name} Sample:**")
    st.write(f"The observed VAF ({vaf_input}%) closely aligns with the following theoretical model(s) (within ±{tolerance}% tolerance):")
    for model in matched_models:
        st.write(f"- {model}")
    
    if "Germline + LOH (Del)" in matched_models or "Somatic + LOH (Del)" in matched_models:
        st.write("**Note**: Since LOH (Deletion) models are suggested, please carefully differentiate between germline and somatic origin, especially if TC is high.")
else:
    st.info(f"**Automated Comment for {gene_name} Sample:**")
    st.write(f"The observed VAF ({vaf_input}%) does not closely align with standard theoretical models (tolerance ±{tolerance}%). This may suggest complex copy number alterations or clonal heterogeneity.")


# --- Alert System (Retained & Improved) ---
st.markdown("---")
st.subheader("⚠️ Quality & Interpretation Alerts")

# 1. Convergence & Mimicry Zone (Gray Zone) Alert (TC >= 60%)
if tc_input >= 60:
    st.warning(f"""
    **Convergence & Mimicry Zone (Gray Zone)**:  
    The current Tumor Content is **{tc_input}%**. Interpretation in this range requires extra caution:
    1. **Mimicry Risk (60-70% TC)**: Somatic LOH events can result in a VAF near 50%, mimicking a standard germline heterozygous state (e.g., as noted by reviewers).
    2. **Convergence Risk (>70% TC)**: Theoretical curves for Germline LOH and Somatic LOH converge significantly, making them difficult to distinguish by VAF alone.
    Clinical correlation (e.g., family history, drug response) is strongly recommended.
    """)

# 2. Low Confidence Zone Alert (TC < 30%)
elif tc_input < 30:
    st.info(f"**Low Confidence Zone**: Please note that interpretation reliability may be limited when Pathological Tumor Content is below 30%.")

# --- Clinical Context Section (Static but Organized) ---
st.markdown("---")
st.subheader("📝 Additional Clinical Interpretation Notes")
st.write(f"""
- **Fixed Scaling**: The graph utilizes fixed 0–100% ranges for both TC and VAF axes, ensuring a standardized visual perspective for clinical use, as requested.
- **Figure 5A Insight**: In samples with high TC ($\ge$ 90%), elevated VAFs are sometimes misidentified as somatic events. In our study, these were confirmed as **Germline LOH**, supported by clinical outcomes such as favorable responses to **PARP inhibitors** (e.g., SEC cases). Accurate distinction in this zone is crucial for identifying HBOC and Lynch syndrome and informing therapeutic decisions.
""")
