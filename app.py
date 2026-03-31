import streamlit as st
import plotly.graph_objects as go
import numpy as np

# 1. ページ設定
st.set_page_config(page_title="VAF-TC Relationship Visualizer", layout="wide")

# 2. タイトル
st.title("🧬 VAF-TC Relationship Visualizer")
st.markdown("Interactive visualization of theoretical Pathological Tumor Content (TC) and Variant Allele Fraction (VAF) relationships.")

# 3. サイドバー設定
st.sidebar.header("📊 Input Parameters")
gene_name = st.sidebar.text_input("Gene Name", value="BRCA2")
tc_input = st.sidebar.slider("Pathological Tumor Content (TC %)", 0, 100, 70)
vaf_input = st.sidebar.slider("Variant Allele Fraction (VAF %)", 0, 100, 57)

# 【追加】入力ガイド
st.sidebar.markdown("---")
st.sidebar.info("""
💡 **How to use:**
Please enter the **Gene Name**, **Pathological TC (%)**, and observed **VAF (%)** in the fields above. 
The plot and interpretation will update automatically.
""")

tc = tc_input / 100.0
vaf = vaf_input / 100.0

# 4. 理論曲線の計算
x = np.linspace(0.01, 1.0, 100)
y_germ_cnloh = (1 + x) / 2
y_germ_del = 1 / (2 - x)
y_germ_hetero = np.full_like(x, 0.5)
y_som_cnloh = x
y_som_del = x / (2 - x)

# 5. メインレイアウト (左 2 : 右 1)
main_col_left, main_col_right = st.columns([2, 1])

# --- 左カラム：グラフ表示 ---
with main_col_left:
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x*100, y=y_germ_cnloh*100, name="Germline + cnLOH", line=dict(color='#d4af37', width=2)))
    fig.add_trace(go.Scatter(x=x*100, y=y_germ_del*100, name="Germline + LOH (Del)", line=dict(color='#e41a1c', width=2)))
    fig.add_trace(go.Scatter(x=x*100, y=y_germ_hetero*100, name="Germline (Hetero)", line=dict(color='#a65628', width=2)))
    fig.add_trace(go.Scatter(x=x*100, y=y_som_cnloh*100, name="Somatic + cnLOH", line=dict(color='#4daf4a', dash='dash')))
    fig.add_trace(go.Scatter(x=x*100, y=y_som_del*100, name="Somatic + LOH (Del)", line=dict(color='#377eb8', dash='dot')))

    fig.add_trace(go.Scatter(
        x=[tc_input], y=[vaf_input],
        mode='markers+text',
        name=f"{gene_name} Sample",
        text=[f"{gene_name}<br>TC:{tc_input}%<br>VAF:{vaf_input}%"],
        textposition="top right",
        marker=dict(color='black', size=12)
    ))

    fig.add_vrect(x0=0, x1=30, fillcolor="gray", opacity=0.1, layer="below", line_width=0, annotation_text="Low Confidence Zone", annotation_position="top left")

    fig.update_layout(
        xaxis_title="Pathological Tumor Content (%)",
        yaxis_title="Variant Allele Fraction (%)",
        yaxis=dict(range=[0, 105]),
        xaxis=dict(range=[0, 105]),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1, font=dict(size=10)),
        template="simple_white",
        height=550,
        margin=dict(l=20, r=20, t=50, b=20)
    )
    st.plotly_chart(fig, use_container_width=True)

# --- 右カラム：解釈と臨床ノート ---
with main_col_right:
    st.subheader("📋 Interpretation")
    
    error_margin = 0.10
    models = {
        "Germline + cnLOH": (1 + tc) / 2,
        "Germline + LOH (Del)": 1 / (2 - tc),
        "Germline (Hetero)": 0.5,
        "Somatic + cnLOH": tc,
        "Somatic + LOH (Del)": tc / (2 - tc)
    }
    
    compatible_models = [name for name, val in models.items() if abs(val - vaf) <= error_margin]

    if compatible_models:
        st.success(f"**Compatible Models for {gene_name}:**")
        st.markdown("Considering a **±10% measurement error**, the observed VAF aligns with the following theoretical model(s):")
        for m in compatible_models:
            st.markdown(f"- **{m}**")
        st.caption("Factors such as NGS variance, aneuploidy, or copy number changes should be considered.")
    else:
        st.info(f"**{gene_name} Insight:** VAF does not closely align with any standard models (deviation > 10%). Consider complex genomic alterations or significant clonal heterogeneity.")

    # 収束リスク警告 (アラート文)
    if 60 <= tc_input <= 75:
        alert_text = f"⚠️ **Convergence Risk (Gray Zone):** At TC {tc_input}% and VAF {vaf_input}%, theoretical curves for **Germline LOH** and **Somatic LOH** converge significantly. Distinguishing between these events based on VAF alone is difficult in this range. Clinical correlation (e.g., family history, drug response) is strongly recommended."
        st.warning(alert_text)

    st.divider()

    # 臨床ノート (科学的根拠に基づく注釈)
    st.subheader("📝 Clinical Notes")
    notes = """
* **Measurement Tolerance:** In clinical NGS analysis, a variance of approximately 10% in VAF is common due to technical limitations and biological factors such as aneuploidy.
* **Tumor Purity:** To ensure accuracy, tumor content (TC) should be determined via **pathological assessment** by a specialist, as NGS-based estimations can carry higher uncertainty.
* **High TC Context:** In samples with high tumor content (TC ≥ 90%), variants with high VAFs are statistically more likely to represent **Germline LOH** rather than somatic events.
* **Therapeutic Implication:** **Biallelic inactivation (LOH)** is the critical indicator for PARP inhibitor sensitivity, regardless of whether the initial variant is germline or somatic in origin.
    """
    st.markdown(notes)
