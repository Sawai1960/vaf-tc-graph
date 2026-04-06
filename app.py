import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import io

# 1. Page Configuration
st.set_page_config(page_title="VAF-TC Precision Analyzer", layout="wide")

# 2. Title
st.title("🧬 VAF-TC Precision Analyzer")
st.markdown("Interactive visual tool for germline/somatic variant differentiation in tumor-only sequencing.")
st.caption("⚠️ This tool is intended as a supportive aid for genetic counseling. It does not replace confirmatory germline testing or established clinical guidelines.")
st.caption("⚠️ Gene Reference System is based on the **Guidelines for GPV/PGPV Handling Procedures in Cancer Gene Panel Testing (2025 Edition)** (MHLW Research Grant) — T-only PGPV disclosure recommended gene list.")

# Important Note (startup)
with st.expander("⚠️ Important Note — Mathematical model assumptions and limitations", expanded=True):
    st.markdown("""
**The five theoretical VAF-TC curves visualized in this tool are derived from *Knudson's two-hit hypothesis* under a strict diploid model. Please consider the following before clinical interpretation:**

1. **Diploid assumption (Knudson's two-hit model)**: Each theoretical line represents a specific biallelic inactivation scenario assuming a diploid (2-copy) baseline. The five models are idealized mathematical references, not an exhaustive catalog of all possible mechanisms.

2. **Aneuploidy is not accounted for**: Real tumors frequently exhibit aneuploidy, whole-chromosome gains/losses, and subclonal heterogeneity. In such cases, observed VAF may deviate substantially from the theoretical lines, and the model matching must be interpreted with caution.

3. **Tumor content estimation carries ±10–20% error**: Pathological TC estimation is subject to **±10–20% variability** due to histological heterogeneity, sampling region, and inter-observer agreement. This tool applies a ±10% matching margin, but clinical interpretation should consider the full uncertainty range.

4. **Not a diagnostic tool**: This is a visual aid for genetic counseling. Confirmatory germline testing remains the standard for any clinical decision.
""")

# 3. Gene Reference Data — based on Guidelines for GPV/PGPV Handling Procedures in Cancer Gene Panel Testing (2025 Edition)
GENE_INFO = {
    # 🔴 Low VAF threshold (VAF >= 10%) — special handling
    "BRCA1":  ("low_threshold", "🔴 BRCA1 [GPV/PGPV Guidelines 2025]: **VAF ≥ 10%** threshold (lower than standard). Even low-VAF variants may be GPV. Expert panel review and confirmatory germline testing recommended. HBOC."),
    "BRCA2":  ("low_threshold", "🔴 BRCA2 [GPV/PGPV Guidelines 2025]: **VAF ≥ 10%** threshold (lower than standard). Even low-VAF variants may be GPV. Expert panel review and confirmatory germline testing recommended. HBOC."),
    # 🟠 Age-conditional (disclosure if SNV VAF ≥ 30% / indel ≥ 20% AND age of onset < 30 y)
    "APC":    ("age_cond", "🟠 APC [GPV/PGPV Guidelines 2025]: PGPV disclosure if **SNV VAF ≥ 30% (indel ≥ 20%) AND colorectal polyposis onset < 30 y**. FAP. Phenotype evaluation required (Box_E)."),
    "CDKN2A": ("age_cond", "🟠 CDKN2A [GPV/PGPV Guidelines 2025]: PGPV disclosure if **SNV VAF ≥ 30% (indel ≥ 20%) AND onset < 30 y**. Hereditary melanoma-pancreatic cancer syndrome."),
    "PTEN":   ("age_cond", "🟠 PTEN [GPV/PGPV Guidelines 2025]: PGPV disclosure if **SNV VAF ≥ 30% (indel ≥ 20%) AND onset < 30 y**. Cowden syndrome. Phenotype evaluation required (Box_E)."),
    "RB1":    ("age_cond", "🟠 RB1 [GPV/PGPV Guidelines 2025]: PGPV disclosure if **SNV VAF ≥ 30% (indel ≥ 20%) AND onset < 30 y**. Hereditary retinoblastoma. Phenotype evaluation required (Box_E)."),
    "TP53":   ("age_cond", "🟠 TP53 [GPV/PGPV Guidelines 2025]: PGPV disclosure if **SNV VAF ≥ 30% (indel ≥ 20%) AND onset < 30 y**. Li-Fraumeni syndrome. Phenotype evaluation required (Box_E). Note: clonal hematopoiesis possible at high VAF."),
    # 🟡 Standard (SNV VAF ≥ 30% / indel ≥ 20%) — 24 genes
    "ATM":    ("standard", "🟡 ATM [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. HBOC, ataxia-telangiectasia."),
    "BAP1":   ("standard", "🟡 BAP1 [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. BAP1 tumor predisposition syndrome."),
    "BARD1":  ("standard", "🟡 BARD1 [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. HBOC."),
    "BRIP1":  ("standard", "🟡 BRIP1 [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. HBOC."),
    "CHEK2":  ("standard", "🟡 CHEK2 [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. HBOC."),
    "DICER1": ("standard", "🟡 DICER1 [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. DICER1 syndrome (pleuropulmonary blastoma, etc.)."),
    "FH":     ("standard", "🟡 FH [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. HLRCC (hereditary leiomyomatosis and renal cell cancer)."),
    "FLCN":   ("standard", "🟡 FLCN [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. Birt-Hogg-Dubé syndrome."),
    "MLH1":   ("standard", "🟡 MLH1 [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. Lynch syndrome."),
    "MSH2":   ("standard", "🟡 MSH2 [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. Lynch syndrome."),
    "MSH6":   ("standard", "🟡 MSH6 [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. Lynch syndrome."),
    "MUTYH":  ("standard", "🟡 MUTYH [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**, **bi-allelic only**. MUTYH-associated polyposis."),
    "NF1":    ("standard", "🟡 NF1 [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. Neurofibromatosis type 1. Phenotype evaluation required (Box_E)."),
    "PALB2":  ("standard", "🟡 PALB2 [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. HBOC."),
    "PMS2":   ("standard", "🟡 PMS2 [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. Lynch syndrome."),
    "POLD1":  ("standard", "🟡 POLD1 [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. Hereditary colorectal cancer."),
    "POLE":   ("standard", "🟡 POLE [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. Hereditary colorectal cancer."),
    "RAD51C": ("standard", "🟡 RAD51C [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. HBOC."),
    "RAD51D": ("standard", "🟡 RAD51D [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. HBOC."),
    "RET":    ("standard", "🟡 RET [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. MEN2 (multiple endocrine neoplasia type 2)."),
    "SDHA":   ("standard", "🟡 SDHA [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. Hereditary paraganglioma-pheochromocytoma syndrome."),
    "SDHB":   ("standard", "🟡 SDHB [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. Hereditary paraganglioma-pheochromocytoma syndrome."),
    "TSC2":   ("standard", "🟡 TSC2 [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. Tuberous sclerosis complex."),
    "VHL":    ("standard", "🟡 VHL [GPV/PGPV Guidelines 2025]: Disclosure at **SNV VAF ≥ 30% (indel ≥ 20%)**. von Hippel-Lindau syndrome."),
}

def get_gene_message(gene):
    key = gene.upper().strip()
    if key in GENE_INFO:
        return GENE_INFO[key]
    return ("not_listed", f"⬜ {gene}: Not on the 2025 T-only PGPV disclosure recommended gene list (MHLW GPV/PGPV Guidelines). Germline disclosure priority is lower per the guidelines. Consult clinical guidelines and family history for comprehensive assessment.")

# 4. Sidebar Input Parameters
st.sidebar.header("📋 Patient Data Input")
st.sidebar.markdown("👉 **Please enter Gene Name, TC, and VAF.**")

gene_name = st.sidebar.text_input("Gene Name", value="BRCA2")
tc_input = st.sidebar.slider("Pathological Tumor Content (TC %)", 0, 100, 50)
vaf_input = st.sidebar.slider("Variant Allele Fraction (VAF %)", 0, 100, 50)

st.sidebar.markdown("---")
st.sidebar.caption("⚠️ When a CSV is uploaded below, the inputs above (Gene Name, TC, VAF) are disabled.")

# Multi-variant CSV Upload
st.sidebar.subheader("📂 Multi-variant Upload")
st.sidebar.caption("CSV format: Gene, TC, VAF")
uploaded_file = st.sidebar.file_uploader("Upload CSV", type=["csv"])

multi_df = None
if uploaded_file is not None:
    try:
        multi_df = pd.read_csv(uploaded_file)
        required_cols = {"Gene", "TC", "VAF"}
        if not required_cols.issubset(multi_df.columns):
            st.sidebar.error("CSV must have columns: Gene, TC, VAF")
            multi_df = None
        else:
            st.sidebar.success(f"{len(multi_df)} variants loaded.")
    except Exception as e:
        st.sidebar.error(f"Error reading CSV: {e}")
        multi_df = None

st.sidebar.markdown("---")
# CSV Template Download (sidebar)
st.sidebar.subheader("📊 Multi-variant Workflow")
st.sidebar.caption("💡 Download the template below, replace the sample genes with your own data, then upload using **Multi-variant Upload** above.")
template_df = pd.DataFrame({
    "Gene": [gene_name, "TP53", "MSH2"],
    "TC":   [tc_input,  tc_input, tc_input],
    "VAF":  [vaf_input, 0.0,      0.0]
})
csv_string = template_df.to_csv(index=False)
st.sidebar.download_button("📥 Download CSV Template", csv_string.encode("utf-8"), "VAF_TC_Template.csv", "text/csv")

# Theoretical Model Data Downloads (sidebar)
st.sidebar.subheader("📂 Theoretical Model Data")
try:
    with open("VAF_TC_theoretical_model.csv", "rb") as f:
        st.sidebar.download_button("📥 Download Theoretical Model (CSV)", f.read(), "VAF_TC_theoretical_model.csv", "text/csv")
except FileNotFoundError:
    st.sidebar.caption("VAF_TC_theoretical_model.csv not found.")
try:
    with open("VAF-TC theoretical_model.xlsx", "rb") as f:
        st.sidebar.download_button("📥 Download Theoretical Model (Excel)", f.read(), "VAF-TC_theoretical_model.xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
except FileNotFoundError:
    st.sidebar.caption("VAF-TC theoretical_model.xlsx not found.")

tc = tc_input / 100.0
vaf = vaf_input / 100.0

# 5. Mathematical Foundation (diploid model) — 5 theoretical models
x_range = np.linspace(0.01, 1.0, 100)
y_germ_cnloh  = (1 + x_range) / 2         # germline (cnLOH)
y_germ_del    = 1 / (2 - x_range)         # germline (LOH with Del)
y_germ_hetero = np.full_like(x_range, 0.5) # germline (Hetero)
y_som_del     = x_range / (2 - x_range)   # somatic (LOH with Del)
y_som_hetero  = x_range / 2               # somatic (Hetero)

# 6. Main Layout (Left 1 : Right 2)
col_alerts, col_graph = st.columns([1, 2])

# --- LEFT COLUMN: Clinical Interpretation & Alerts ---
with col_alerts:
    st.info(f"💡 **Analysis Mode:** {gene_name}")
    st.subheader("📋 Interpretation & Alerts")

    error_margin = 0.10

    def get_compatible_models(tc_val, vaf_val):
        f = tc_val / 100.0
        v = vaf_val / 100.0
        checks = {
            "germline (cnLOH)":        (1 + f) / 2,
            "germline (LOH with Del)": 1 / (2 - f),
            "germline (Hetero)":       0.5,
            "somatic (LOH with Del)":  f / (2 - f),
            "somatic (Hetero)":        f / 2,
        }
        return [(name, val) for name, val in checks.items() if abs(val - v) <= error_margin]

    def get_interpretation(compatible):
        names = [name for name, _ in compatible]
        if not names:
            return "warning", "VAF does not align with any standard model. Consider clonal heterogeneity, aneuploidy, or complex copy number changes."
        germ_cnloh  = "germline (cnLOH)"        in names
        germ_del    = "germline (LOH with Del)" in names
        germ_hetero = "germline (Hetero)"       in names
        som_del     = "somatic (LOH with Del)"  in names
        som_hetero  = "somatic (Hetero)"        in names

        has_germ = germ_cnloh or germ_del or germ_hetero
        has_som  = som_del or som_hetero

        # Ambiguous: both germline and somatic compatible
        if has_germ and has_som:
            return "error", "VAF is compatible with **both germline and somatic** models at this TC. VAF alone **cannot determine the origin**. **Pair-normal germline testing is essential.**"
        # Pure germline — single model
        if germ_cnloh and not germ_del and not germ_hetero:
            return "success", "Pattern is consistent with a **germline variant that has undergone copy-neutral LOH (UPD)**. Biallelic inactivation via germline + cnLOH."
        if germ_del and not germ_cnloh and not germ_hetero:
            return "success", "Pattern is consistent with a **germline variant with LOH by deletion**. Biallelic inactivation via germline + deletion."
        if germ_hetero and not germ_cnloh and not germ_del:
            return "success", "Pattern is consistent with a **heterozygous germline variant without LOH**. Only one allele is affected."
        # Multiple germline models
        if has_germ:
            return "info", "Multiple germline models are compatible. Germline origin is likely, but the LOH mechanism cannot be determined from VAF alone."
        # Pure somatic — single model
        if som_del and not som_hetero:
            return "info", "Pattern is consistent with a **somatic variant with LOH by deletion**. Germline origin is unlikely at this TC."
        if som_hetero and not som_del:
            return "info", "Pattern is consistent with a **somatic heterozygous variant (without LOH)**. Germline origin is unlikely at this TC."
        # Multiple somatic models
        if has_som:
            return "info", "Multiple somatic models are compatible. Somatic origin is likely at this TC."
        return "info", "Clinical correlation and pair-normal testing are recommended."

    def show_variant_interpretation(g, t, v):
        compatible = get_compatible_models(t, v)
        st.markdown(f"**{g}** (TC {t:.0f}%, VAF {v:.0f}%)")
        if compatible:
            for name, val in compatible:
                st.markdown(f"- **{name}** — theoretical VAF {val*100:.1f}%")
        level, msg = get_interpretation(compatible)
        if level == "success":
            st.success(f"➡️ {msg}")
        elif level == "error":
            st.error(f"➡️ {msg}")
        elif level == "warning":
            st.warning(f"➡️ {msg}")
        else:
            st.info(f"➡️ {msg}")
        # TC-based alerts
        if t <= 20:
            st.warning("⚠️ **Low TC (≤ 20%):** At low tumor content, theoretical lines are compressed into a narrow VAF range and model matching is less reliable. Subclonal variants, admixture with normal tissue, or technical noise may dominate.")
        if t >= 60:
            st.warning("⚠️ **High TC (≥ 60%):** At high tumor content, germline and somatic LOH lines begin to converge. Origin determination by VAF alone becomes increasingly difficult.")
        # Gene-specific message (GPV/PGPV Guidelines 2025)
        _, gene_msg = get_gene_message(g)
        st.info(gene_msg)
        st.caption("💡 Note: VAF–TC interpretation is based on mathematical models only. Gene-specific germline likelihood is provided separately per the GPV/PGPV Guidelines (2025 Edition).")

    if multi_df is not None:
        for _, row in multi_df.iterrows():
            show_variant_interpretation(str(row["Gene"]), float(row["TC"]), float(row["VAF"]))
            st.divider()
    else:
        show_variant_interpretation(gene_name, tc_input, vaf_input)

    # --- TC-based Clinical Alerts ---
    som_del_vaf  = tc / (2 - tc) * 100 if tc < 2 else 0
    germ_del_vaf = 1 / (2 - tc) * 100 if tc < 2 else 0

    if 61 <= tc_input <= 66:
        st.warning(
            f"⚠️ **Gray Zone (Somatic LOH Del):** At TC {tc_input}%, "
            f"somatic (LOH with Del) produces VAF = {som_del_vaf:.1f}%, "
            f"approaching germline (Hetero) at 50%. "
            f"Confirmation testing is recommended."
        )
    elif tc_input >= 67:
        if vaf_input >= tc / (2 - tc) * 100:
            st.error(
                f"🔴 **LOH Convergence Alert:** At TC {tc_input}% and VAF {vaf_input}%, "
                f"the variant falls at or above the somatic (LOH with Del) line "
                f"({som_del_vaf:.1f}%). In this region, germline (LOH with Del) = "
                f"{germ_del_vaf:.1f}% and somatic (LOH with Del) = {som_del_vaf:.1f}% "
                f"converge — origin cannot be determined by VAF alone. "
                f"Germline confirmation is essential."
            )
        if tc_input >= 90:
            st.warning(
                f"⚠️ **Extreme Tumor Purity:** At TC {tc_input}%, all theoretical "
                f"models compress into a narrow VAF range. Variants may still be "
                f"of somatic origin even at high VAF. "
                f"Family history review and germline testing are essential."
            )

    st.divider()

    # Multi-variant table display
    if multi_df is not None:
        st.subheader("📋 Uploaded Variants")
        st.dataframe(multi_df, use_container_width=True)
        st.divider()


# --- RIGHT COLUMN: Visualization ---
with col_graph:
    st.subheader("📈 VAF-TC Projection")
    fig = go.Figure()

    fig.add_trace(go.Scatter(x=x_range*100, y=y_germ_cnloh*100,  name="germline (cnLOH)",        line=dict(color='#d4af37', width=2)))
    fig.add_trace(go.Scatter(x=x_range*100, y=y_germ_del*100,    name="germline (LOH with Del)", line=dict(color='#e41a1c', width=2)))
    fig.add_trace(go.Scatter(x=x_range*100, y=y_germ_hetero*100, name="germline (Hetero)",       line=dict(color='#a65628', width=2)))
    fig.add_trace(go.Scatter(x=x_range*100, y=y_som_del*100,     name="somatic (LOH with Del)",  line=dict(color='#377eb8', dash='dot')))
    fig.add_trace(go.Scatter(x=x_range*100, y=y_som_hetero*100,  name="somatic (Hetero)",        line=dict(color='#4daf4a', dash='dash')))

    if multi_df is not None:
        colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00',
                  '#a65628','#f781bf','#999999','#66c2a5','#fc8d62']
        for i, row in multi_df.iterrows():
            fig.add_trace(go.Scatter(
                x=[row["TC"]], y=[row["VAF"]],
                mode='markers+text',
                name=str(row["Gene"]),
                text=[f"{row['Gene']}<br>TC:{row['TC']}%<br>VAF:{row['VAF']}%"],
                textposition="top right",
                marker=dict(color=colors[i % len(colors)], size=14, symbol='circle'),
                showlegend=True
            ))
    else:
        fig.add_trace(go.Scatter(
            x=[tc_input], y=[vaf_input],
            mode='markers+text',
            name=f"Current: {gene_name}",
            text=[f"{gene_name}<br>TC:{tc_input}%<br>VAF:{vaf_input}%"],
            textposition="top right",
            marker=dict(color='black', size=14, symbol='circle')
        ))

    fig.add_vrect(x0=0, x1=20, fillcolor="gray", opacity=0.1, layer="below", line_width=0,
                  annotation_text="Low Confidence Zone", annotation_position="top left")

    fig.update_layout(
        xaxis_title="Pathological Tumor Content (%)", yaxis_title="Variant Allele Fraction (%)",
        yaxis=dict(range=[0, 105]), xaxis=dict(range=[0, 105]),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        template="simple_white", height=600
    )
    st.plotly_chart(fig, use_container_width=True)

    # Gene Reference Table (GPV/PGPV Guidelines 2025)
    st.subheader("📖 Gene Reference (GPV/PGPV Guidelines 2025)")
    st.markdown("**🔴 Low VAF threshold (VAF ≥ 10%):**")
    st.caption("BRCA1, BRCA2")
    st.markdown("**🟠 Age-conditional (onset < 30 y):**")
    st.caption("APC (colorectal polyposis), CDKN2A, PTEN, RB1, TP53")
    st.markdown("**🟡 Standard (SNV VAF ≥ 30%, indel ≥ 20%):**")
    st.caption("ATM, BAP1, BARD1, BRIP1, CHEK2, DICER1, FH, FLCN, MLH1, MSH2, MSH6, MUTYH(bi), NF1, PALB2, PMS2, POLD1, POLE, RAD51C, RAD51D, RET, SDHA, SDHB, TSC2, VHL")
    st.caption("⬜ Genes not listed: not on the 2025 T-only PGPV list.")
    st.caption("Reference: MHLW Research Grant — Guidelines for GPV/PGPV Handling Procedures in Cancer Gene Panel Testing (2025 Edition)")

# 7. Footer
st.divider()
st.caption("VAF-TC Precision Analyzer | Clinical Genetics Suite | ver 3.4 ✅")
