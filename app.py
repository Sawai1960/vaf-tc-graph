import streamlit as st
import plotly.graph_objects as go
import numpy as np

# ページ設定
st.set_page_config(page_title="VAF-TC Relationship Visualizer", layout="wide")

# タイトルと説明
st.title("🧬 VAF-TC Relationship Visualizer")
st.markdown("An interactive tool to visualize the theoretical relationship between Pathological Tumor Content (TC) and Variant Allele Fraction (VAF).")

# サイドバー設定
st.sidebar.header("📊 Input Parameters")
gene_name = st.sidebar.text_input("Gene Name", value="BRCA2")
tc_input = st.sidebar.slider("Pathological Tumor Content (TC %)", 0, 100, 70)
vaf_input = st.sidebar.slider("Variant Allele Fraction (VAF %)", 0, 100, 57)

# 数値変換 (0-100 -> 0.0-1.0)
tc = tc_input / 100.0
vaf = vaf_input / 100.0

# 理論曲線の計算
x = np.linspace(0.01, 1.0, 100)
y_germ_cnloh = (1 + x) / 2
y_germ_del = 1 / (2 - x)
y_germ_hetero = np.full_like(x, 0.5)
y_som_cnloh = x
y_som_del = x / (2 - x)

# プロット作成
fig = go.Figure()

# 理論線を追加
fig.add_trace(go.Scatter(x=x*100, y=y_germ_cnloh*100, name="Germline + cnLOH", line=dict(color='#d4af37', width=2)))
fig.add_trace(go.Scatter(x=x*100, y=y_germ_del*100, name="Germline + LOH (Del)", line=dict(color='#e41a1c', width=2)))
fig.add_trace(go.Scatter(x=x*100, y=y_germ_hetero*100, name="Germline (Hetero)", line=dict(color='#a65628', width=2)))
fig.add_trace(go.Scatter(x=x*100, y=y_som_cnloh*100, name="Somatic + cnLOH", line=dict(color='#4daf4a', dash='dash')))
fig.add_trace(go.Scatter(x=x*100, y=y_som_del*100, name="Somatic + LOH (Del)", line=dict(color='#377eb8', dash='dot')))

# 入力データを追加
fig.add_trace(go.Scatter(
    x=[tc_input], y=[vaf_input],
    mode='markers+text',
    name=f"{gene_name} Sample",
    text=[f"{gene_name}<br>TC:{tc_input}%<br>VAF:{vaf_input}%"],
    textposition="top right",
    marker=dict(color='black', size=12)
))

# 低信頼領域（TC < 30%）のシェーディング
fig.add_vrect(x0=0, x1=30, fillcolor="gray", opacity=0.1, layer="below", line_width=0, annotation_text="Low Confidence Zone", annotation_position="top left")

# レイアウト調整
fig.update_layout(
    xaxis_title="Pathological Tumor Content (%)",
    yaxis_title="Variant Allele Fraction (%)",
    yaxis=dict(range=[0, 105]),
    xaxis=dict(range=[0, 105]),
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    template="simple_white",
    height=600
)

st.plotly_chart(fig, use_container_width=True)

# --- 動的な解釈テキストの生成 ---

# 各モデルとの誤差を計算
models = {
    "Germline + cnLOH": (1 + tc) / 2,
    "Germline + LOH (Del)": 1 / (2 - tc),
    "Germline (Hetero)": 0.5,
    "Somatic + cnLOH": tc,
    "Somatic + LOH (Del)": tc / (2 - tc)
}

closest_model = min(models, key=lambda k: abs(models[k] - vaf))
diff = abs(models[closest_model] - vaf) * 100

# 解釈セクションの表示
st.subheader("📋 Clinical Interpretation")

# 10%の誤差範囲内かどうかでメッセージを分岐
if diff <= 10.0:
    st.success(f"**Interpretation for {gene_name}:** The observed VAF is within the expected range of **measurement error (approx. ±10%)** for standard models. It aligns most closely with **{closest_model}**, but potential factors including **NGS variance, aneuploidy, or copy number changes** should be considered.")
else:
    st.info(f"**Interpretation for {gene_name}:** VAF does not closely align with standard models (deviation > 10%). Consider complex genomic alterations, significant clonal heterogeneity, or multiple copy number changes.")

# 60-70%の収束リスク警告
if 60 <= tc_input <= 75:
    st.warning("⚠️ **Convergence Risk (Gray Zone):** At this TC range, theoretical curves for Germline LOH and Somatic LOH converge significantly. Distinguishing between these events based on VAF alone is clinically challenging.")

# --- 臨床ノートセクション (第一著者の指摘を反映) ---
st.divider()
st.subheader("📝 Clinical Interpretation Notes")

col1, col2 = st.columns(2)

with col1:
    st.markdown("""
    * **Measurement Tolerance:** Based on our study, a variance of approximately 10% in VAF is common due to NGS technical limitations and biological factors such as aneuploidy.
    * **Tumor Purity:** To ensure clinical reliability, tumor content (TC) must be determined via **pathological assessment** by a pathologist, as NGS-based estimation often carries higher uncertainty.
    """)

with col2:
    st.markdown("""
    * **Clinical Context (High TC):** In samples with TC ≥ 90%, variants with high VAFs are statistically more likely to represent **Germline LOH** rather than somatic events, particularly in contexts like ovarian cancer.
    * **Biallelic Loss and Therapy:** Both germline and somatic variants can lead to biallelic inactivation (LOH), which is a key indicator for **PARP inhibitor sensitivity**, regardless of the variant's origin.
    """)
