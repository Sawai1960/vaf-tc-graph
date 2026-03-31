import streamlit as st
import plotly.graph_objects as go
import numpy as np

# ページ設定
st.set_page_config(page_title="VAF-TC Relationship Visualizer", layout="wide")

# タイトル
st.title("🧬 VAF-TC Relationship Visualizer")
st.markdown("Interactive visualization of theoretical Pathological Tumor Content (TC) and Variant Allele Fraction (VAF) relationships.")

# --- サイドバー設定 ---
st.sidebar.header("📊 Input Parameters")
gene_name = st.sidebar.text_input("Gene Name", value="BRCA2")
tc_input = st.sidebar.slider("Pathological Tumor Content (TC %)", 0, 100, 70)
vaf_input = st.sidebar.slider("Variant Allele Fraction (VAF %)", 0, 100, 57)

# 数値変換 (0-100 -> 0.0-1.0)
tc = tc_input / 100.0
vaf = vaf_input / 100.0

# --- 理論曲線の計算 ---
x = np.linspace(0.01, 1.0, 100)
y_germ_cnloh = (1 + x) / 2
y_germ_del = 1 / (2 - x)
y_germ_hetero = np.full_like(x, 0.5)
y_som_cnloh = x
y_som_del = x / (2 - x)

# --- レイアウト作成 (左 2 : 右 1) ---
main_col_left, main_col_right = st.columns([2, 1])

with main_col_left:
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
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1, font=dict(size=10)),
        template="simple_white",
        height=550,
        margin=dict(l=20, r=20, t=50, b=20)
    )

    st.plotly_chart(fig, use_container_width=True)

with main_col_right:
    # --- 動的な解釈エリア ---
    st.subheader("📋 Interpretation")
    
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

    # 10%の誤差範囲（測定誤差）に基づくメッセージ
    if diff <= 10.0:
        st.success(f"**{gene_name} Insight:** The observed VAF is within the expected range of **measurement error (approx. ±10%)** for standard models. It aligns most closely with **{closest_model}**, but potential factors including **NGS variance, aneuploidy, or copy number changes** should be considered.")
    else:
        st.info(f"**{gene_name} Insight:** VAF does not closely align with standard models (deviation > 10%). Consider complex genomic alterations or significant clonal heterogeneity.")

    # 60-70%の収束リスク警告
    if 60 <= tc_input <= 75:
        st.warning("⚠️ **Convergence Risk (Gray Zone):** Theoretical curves for Germline LOH and Somatic LOH converge significantly in this TC range. Clinical correlation is strongly recommended.")

    st.divider()

    # --- 臨床ノートエリア (図番号を削除し、内容を一般化) ---
    st.subheader("📝 Clinical Notes")
    st.markdown(f"""
    * **Measurement Tolerance:** In clinical NGS analysis, a variance of approximately 10% in VAF is commonly observed due to technical limitations and biological factors such as aneuploidy. 
    * **Tumor Purity:** To ensure accuracy, tumor content (TC) should be determined via **pathological assessment**, as NGS-based estimations can carry higher uncertainty. [cite: 13]
    * **High TC Context:** In samples with high tumor content (TC ≥ 90%), variants with high VAFs are statistically more likely to be **Germline LOH** rather than somatic events. [cite: 15]
    * **Therapeutic Implication:** Biallelic inactivation (LOH) is a critical indicator for **PARP inhibitor sensitivity**, regardless of whether the initial variant is germline or somatic in origin. [cite: 18]
    """)
