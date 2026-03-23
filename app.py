import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# --- Professional Page Configuration ---
st.set_page_config(page_title="VAF-TC Visualizer | Clinical Decision Support", layout="wide")
st.title("🧪 VAF–Tumor Content Graph Visualizer")

# Custom CSS for the Advice Box
st.markdown("""
<style>
.advice-box { padding: 15px; border-radius: 5px; border: 1px solid #d4af37; background-color: #f9f9f9; color: #333; }
.advice-title { font-weight: bold; font-size: 1.1em; color: #b71c1c; margin-bottom: 5px; }
</style>
""", unsafe_allow_html=True)

st.markdown("""
### Interactive interpretation of variant allele fraction in tumor-only sequencing data.
---
**Overview:** This application visualizes the relationship between **Variant Allele Fraction (VAF %)** and **Pathological Tumor Content (TC %)** based on the **Knudson two-hit hypothesis**. 
""")

# --- Sidebar for Data Input (Scale: 0-100) ---
st.sidebar.header("Patient Data Input (%)")
st.sidebar.markdown("---")
tc = st.sidebar.slider("Pathological Tumor Content (TC %)", 0, 100, 50, 1)
vaf = st.sidebar.slider("Observed VAF (%)", 0, 100, 50, 1)
gene = st.sidebar.text_input("Variant Identifier", "Variant X")

# --- Interpretation Logic (Calculation uses decimal, results converted to %) ---
tc_dec = tc / 100.0
vaf_dec = vaf / 100.0
tol = 0.05

germline_hetero = 0.5
germline_loh_del = 1 / (2 - tc_dec) if tc_dec < 1 else 1.0
germline_loh_cn = 0.5 * tc_dec + 0.5
somatic_hetero = 0.5 * tc_dec
somatic_loh_del = tc_dec / (2 - tc_dec) if tc_dec < 2 else 1.0
somatic_loh_cn = tc_dec

advice_text = ""
if tc < 30:
    advice_text = "**Low Confidence Zone (<30% TC):** Data reliability is limited at low tumor content. Clinical interpretation should be extremely cautious. These values may not reflect the true clonal architecture."
elif 60 <= tc <= 70:
    advice_text = "**TC Gray Zone Alert (60-70%):** At this tumor content, theoretical lines for germline and somatic LOH overlap significantly. Interpret with caution."
elif abs(vaf_dec - germline_loh_cn) < tol or abs(vaf_dec - germline_loh_del) < tol:
    advice_text = "**Likely Germline with LOH:** The high VAF relative to TC suggests a hereditary variant that has undergone biallelic inactivation (LOH) within the tumor."
elif abs(vaf_dec - germline_hetero) < tol:
    advice_text = "**Likely Heterozygous Germline:** Consistent with a constitutional heterozygous state regardless of tumor content."
elif abs(vaf_dec - somatic_loh_cn) < tol or abs(vaf_dec - somatic_loh_del) < tol:
    advice_text = "**Likely Somatic with LOH:** The VAF correlates with tumor content, suggesting an acquired variant that has undergone LOH."
elif abs(vaf_dec - somatic_hetero) < tol:
    advice_text = "**Likely Heterozygous Somatic:** The VAF is approximately half of the tumor content, indicating a typical acquired event without LOH."
else:
    advice_text = "**Atypical Distribution:** The observed VAF does not align with standard diploid models."

# --- Main Visualization Area (Scale: 0-100) ---
col1, col2 = st.columns([3, 1])

with col1:
    fig, ax = plt.subplots(figsize=(10, 7))
    x_range_dec = np.linspace(0.001, 1.0, 100)
    x_range = x_range_dec * 100  # Convert to % for plotting
    
    # Shading for Low Confidence Area (0 to 30%)
    ax.axvspan(0, 30, color='gray', alpha=0.15, label="Low Confidence Area (<30%)")
    ax.axvline(x=30, color='gray', linestyle='--', alpha=0.3)

    # Plotting Lines (Values converted to 0-100 scale)
    ax.plot(x_range, (0.5 * x_range_dec + 0.5) * 100, color='#D4AF37', label="Germline + cnLOH")
    ax.plot(x_range, (1 / (2 - x_range_dec)) * 100, color='red', label="Germline + LOH (Del)")
    ax.axhline(50, color='brown', linewidth=2, label="Germline (Hetero)")
    ax.plot(x_range, (x_range_dec) * 100, color='green', linestyle='--', alpha=0.5, label="Somatic + cnLOH")
    ax.plot(x_range, (x_range_dec / (2 - x_range_dec)) * 100, color='gray', linestyle=':', label="Somatic + LOH (Del)")
    ax.plot(x_range, (0.5 * x_range_dec) * 100, color='gray', linestyle='--', alpha=0.5, label="Somatic (Hetero)")
    
    # Patient point
    ax.scatter(tc, vaf, color='black', s=200, zorder=5, label=f"Patient: {gene}")
    
    # Formatting (Axis labels and limits in %)
    ax.set_xlabel("Tumor Content (%)", fontsize=12)
    ax.set_ylabel("Variant Allele Fraction (%)", fontsize=12)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax.grid(True, linestyle='--', alpha=0.2)
    st.pyplot(fig)

with col2:
    st.header("Summary")
    st.write(f"**Target:** {gene}")
    st.write(f"**TC:** {tc}% | **VAF:** {vaf}%")
    st.write("---")
    
    st.header("Clinical Guidance")
    st.markdown(f"<div class='advice-box'><div class='advice-title'>Automated Talking Points</div>{advice_text}</div>", unsafe_allow_html=True)
    
    st.write("---")
    st.header("Reference Guide")
    st.markdown("""
    - **Shaded Area:** Low Confidence (<30% TC).
    - **Solid Lines:** Possible germline origin.
    - **Dashed Lines:** Possible somatic origin.
    - **Gray Zone:** TC 60-70% overlap.
    """)

st.divider()
st.caption("Clinical Disclaimer: This tool is for supportive visual communication and educational purposes only.")
