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

# 3. Sidebar Input Parameters
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
st.sidebar.info(f"💡 **Analysis Mode:** {gene_name}")

tc = tc_input / 100.0
vaf = vaf_input / 100.0

# 4. Mathematical Foundation (diploid model)
x_range = np.linspace(0.01, 1.0, 100)
y_germ_cnloh = (1 + x_range) / 2
y_germ_del = 1 / (2 - x_range)
y_germ_hetero = np.full_like(x_range, 0.5)
y_som_cnloh = x_range
y_som_del = x_range / (2 - x_range)

# 5. Main Layout (Left 1 : Right 2)
col_alerts, col_graph = st.columns([1, 2])

# --- LEFT COLUMN: Clinical Interpretation & Alerts ---
with col_alerts:
    st.subheader("📋 Interpretation & Alerts")

    # Mathematical Match Analysis (±10% error margin)
    error_margin = 0.10

    def get_compatible_models(tc_val, vaf_val):
        f = tc_val / 100.0
        v = vaf_val / 100.0
        checks = {
            "Germline + cnLOH": (1 + f) / 2,
            "Germline + LOH (Del)": 1 / (2 - f),
            "Germline (Hetero)": 0.5,
            "Somatic + cnLOH": f,
            "Somatic + LOH (Del)": f / (2 - f)
        }
        return [(name, val) for name, val in checks.items() if abs(val - v) <= error_margin]

    if multi_df is not None:
        # Multi-variant interpretation
        for _, row in multi_df.iterrows():
            g, t, v = str(row["Gene"]), float(row["TC"]), float(row["VAF"])
            compatible = get_compatible_models(t, v)
            if compatible:
                st.success(f"**{g}** (TC {t:.0f}%, VAF {v:.0f}%)")
                for name, val in compatible:
                    st.markdown(f"- **{name}** — theoretical VAF {val*100:.1f}%")
            else:
                st.info(f"**{g}** (TC {t:.0f}%, VAF {v:.0f}%) — No model match within ±10%.")
    else:
        # Single-variant interpretation
        compatible_models = get_compatible_models(tc_input, vaf_input)
        if compatible_models:
            st.success(f"**Compatible Models for {gene_name} (±10%):**")
            for name, val in compatible_models:
                st.markdown(f"- **{name}** — theoretical VAF {val*100:.1f}%")
        else:
            st.info(f"**Insight:** VAF {vaf_input}% does not align with any standard model at TC {tc_input}% (±10%).")

    # --- Clinical Alerts ---

    # Pre-compute key thresholds
    som_cnloh_vaf = tc * 100
    som_del_vaf = tc / (2 - tc) * 100
    germ_del_vaf = 1 / (2 - tc) * 100

    # Alert 1: Somatic cnLOH Trap (TC 40–60%)
    if 40 <= tc_input <= 60:
        st.warning(
            f"⚠️ **Somatic cnLOH Trap:** At TC {tc_input}%, Somatic cnLOH (UPD) "
            f"produces VAF = {som_cnloh_vaf:.0f}%, which falls within ±10% of "
            f"Germline Heterozygous (50%). A somatic variant with cnLOH can "
            f"masquerade as a germline heterozygous variant. "
            f"Pair-normal testing is essential."
        )

    # Alert 2: Gray Zone (TC 61–66%)
    elif 61 <= tc_input <= 66:
        st.warning(
            f"⚠️ **Gray Zone (Somatic LOH Del):** At TC {tc_input}%, "
            f"Somatic LOH (deletion) produces VAF = {som_del_vaf:.1f}%, "
            f"approaching Germline Heterozygous (50%). "
            f"Confirmation testing is recommended."
        )

    # Alert 3: LOH Convergence Zone (TC ≥ 67%)
    elif tc_input >= 67:
        if vaf_input >= tc / (2 - tc) * 100:
            st.error(
                f"🔴 **LOH Convergence Alert:** At TC {tc_input}% and VAF {vaf_input}%, "
                f"the variant falls at or above the Somatic LOH (deletion) line "
                f"({som_del_vaf:.1f}%). In this region, Germline LOH (Del) = "
                f"{germ_del_vaf:.1f}% and Somatic LOH (Del) = {som_del_vaf:.1f}% "
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

    # CSV Template Download
    st.subheader("📊 Multi-variant Workflow")
    st.caption("💡 Download the template below, replace the sample genes with your own data, then upload the file using **Multi-variant Upload** in the left sidebar.")
    template_df = pd.DataFrame({
        "Gene": [gene_name, "TP53", "MSH2"],
        "TC":   [tc_input,  tc_input, tc_input],
        "VAF":  [vaf_input, 0.0,      0.0]
    })
    csv_string = template_df.to_csv(index=False)
    st.download_button("📥 Download CSV Template", csv_string.encode("utf-8"), "VAF_TC_Template.csv", "text/csv")

    # Theoretical Model Data Downloads
    st.subheader("📂 Theoretical Model Data")
    try:
        with open("VAF_TC_theoretical_model.csv", "rb") as f:
            st.download_button("📥 Download Theoretical Model (CSV)", f.read(), "VAF_TC_theoretical_model.csv", "text/csv")
    except FileNotFoundError:
        st.caption("VAF_TC_theoretical_model.csv not found.")
    try:
        with open("VAF-TC theoretical_model.xlsx", "rb") as f:
            st.download_button("📥 Download Theoretical Model (Excel)", f.read(), "VAF-TC_theoretical_model.xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    except FileNotFoundError:
        st.caption("VAF-TC theoretical_model.xlsx not found.")


# --- RIGHT COLUMN: Visualization ---
with col_graph:
    st.subheader("📈 VAF-TC Projection")
    fig = go.Figure()

    # Theoretical lines
    fig.add_trace(go.Scatter(x=x_range*100, y=y_germ_cnloh*100, name="Germline + cnLOH", line=dict(color='#d4af37', width=2)))
    fig.add_trace(go.Scatter(x=x_range*100, y=y_germ_del*100, name="Germline + LOH (Del)", line=dict(color='#e41a1c', width=2)))
    fig.add_trace(go.Scatter(x=x_range*100, y=y_germ_hetero*100, name="Germline (Hetero)", line=dict(color='#a65628', width=2)))
    fig.add_trace(go.Scatter(x=x_range*100, y=y_som_cnloh*100, name="Somatic + cnLOH", line=dict(color='#4daf4a', dash='dash')))
    fig.add_trace(go.Scatter(x=x_range*100, y=y_som_del*100, name="Somatic + LOH (Del)", line=dict(color='#377eb8', dash='dot')))

    # Multi-variant plot (CSV upload)
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
        # Single-variant plot
        fig.add_trace(go.Scatter(
            x=[tc_input], y=[vaf_input],
            mode='markers+text',
            name=f"Current: {gene_name}",
            text=[f"{gene_name}<br>TC:{tc_input}%<br>VAF:{vaf_input}%"],
            textposition="top right",
            marker=dict(color='black', size=14, symbol='circle')
        ))

    # Low Confidence Zone (TC < 30%)
    fig.add_vrect(x0=0, x1=30, fillcolor="gray", opacity=0.1, layer="below", line_width=0,
                  annotation_text="Low Confidence Zone", annotation_position="top left")

    fig.update_layout(
        xaxis_title="Pathological Tumor Content (%)", yaxis_title="Variant Allele Fraction (%)",
        yaxis=dict(range=[0, 105]), xaxis=dict(range=[0, 105]),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        template="simple_white", height=600
    )
    st.plotly_chart(fig, use_container_width=True)

# 6. Footer
st.divider()
st.caption("VAF-TC Precision Analyzer | Clinical Genetics Suite | ver 3.1 ✅")
