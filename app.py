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
**Overview:** This application visualizes the relationship between **Variant Allele Fraction (VAF)** and **Pathological Tumor Content (TC)** based on the **Knudson two-hit hypothesis**. 
""")

# --- Sidebar for Data Input ---
st.sidebar.header("Patient Data Input")
st.sidebar.markdown("---")
tc = st.sidebar.slider("Pathological Tumor Content (TC)", 0.0, 1.0, 0.5, 0.01)
vaf = st.sidebar.slider("Observed VAF", 0.0, 1.0, 0.5, 0.01)
gene = st.sidebar.text_input("Variant Identifier", "Variant X")

# --- Interpretation Logic ---
tol = 0.05
germline_hetero = 0.5
germline_loh_del = 1 / (2 - tc) if tc < 1 else 1.0
germline_loh_cn = 0.5 * tc + 0.5
somatic_hetero = 0.5 * tc
somatic_loh_del = tc / (2 - tc) if tc < 2 else 1.0
somatic_loh_cn = tc

advice_text = ""
# Added advice for Low TC
if tc < 0.30:
    advice_text = "**Low Confidence Zone (<30% TC):** Data reliability is limited at low tumor content. Clinical interpretation should be extremely cautious. These values may not reflect the true clonal architecture."
elif abs(tc - 0.5) < 0.05 and abs(vaf - 0.5) < 0.05:
    advice_text = "**Inconclusive Intersection:** This data point lies where 'Heterozygous Germline' and 'Somatic with Copy-neutral LOH' converge. VAF analysis alone cannot distinguish the origin."
elif 0.60 <= tc <= 0.70:
    advice_text = "**TC Gray Zone Alert (60-70%):** At this tumor content, theoretical lines for germline and somatic LOH overlap significantly. Interpret with caution."
elif abs(vaf - germline_loh_cn) < tol or abs(vaf - germline_loh_del) < tol:
    advice_text = "**Likely Germline with LOH:** The high VAF relative to TC suggests a hereditary variant that has undergone biallelic inactivation (LOH) within the tumor."
elif abs(vaf - germline_hetero) < tol:
    advice_text = "**Likely Heterozygous Germline:** Consistent with a constitutional heterozygous state regardless of tumor content."
elif abs(vaf - somatic_loh_cn) < tol or abs(vaf - somatic_loh_del) < tol:
    advice_text = "**Likely Somatic with LOH:** The VAF correlates with tumor content, suggesting an acquired variant that has undergone LOH."
elif abs(vaf - somatic_hetero) < tol:
    advice_text = "**Likely Heterozygous Somatic:** The VAF is approximately half of the tumor content, indicating a typical acquired event without LOH."
else:
    advice_text = "**Atypical Distribution:** The observed VAF does not align with standard diploid models."

# --- Main Visualization Area ---
col1, col2 = st.columns([3, 1])

with col1:
    fig, ax = plt.subplots(figsize=(10, 7))
    x_range = np.linspace(0.001, 1.0, 100) # Start from slightly above 0 to avoid division by zero
    
    # 【Added】Shading for Low Confidence Area (0 to 0.3)
    ax.axvspan(0, 0.3, color='gray', alpha=0.15, label="Low Confidence Area (<30%)")
    ax.axvline(x=0.3, color='gray', linestyle='--', alpha=0.3)

    # Plotting Lines
    ax.plot(x_range, 0.5 * x_range + 0.5, color='#D4AF37', label="Germline + cnLOH")
    ax.plot(x_range, 1 / (2 - x_range), color='red', label="Germline + LOH (Del)")
    ax.axhline(0.5, color='brown', linewidth=2, label="Germline (Hetero)")
    ax.plot(x_range, x_range, color='green', linestyle='--', alpha=0.5, label="Somatic + cnLOH")
    ax.plot(x_range, x_range / (2 - x_range), color='gray', linestyle=':', label="Somatic + LOH (Del)")
    ax.plot(x_range, 0.5 * x_range, color='gray', linestyle='--', alpha=0.5, label="Somatic (Hetero)")
    
    # Patient point
    ax.scatter(tc, vaf, color='black', s=200, zorder=5, label=f"Patient: {gene}")
    
    # Formatting
    ax.set_xlabel("Tumor Content (TC)")
    ax.set_ylabel("Variant Allele Fraction (VAF)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax.grid(True, linestyle='--', alpha=0.2)
    st.pyplot(fig)

with col2:
    st.header("Summary")
    st.write(f"**Target:** {gene}")
    st.write(f"**TC:** {tc*100:.0f}% | **VAF:** {vaf*100:.0f}%")
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
