import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import io

# --- ページ設定 ---
st.set_page_config(page_title="VAF-TC Precision Analyzer", layout="wide")

st.title("🧬 VAF-TC Precision Analyzer")
st.markdown("Clinical Decision-Support Tool for Germline/Somatic Differentiation")

# --- サイドバー入力 ---
st.sidebar.header("📋 Patient Data Input")
tc_input = st.sidebar.slider("Pathological Tumor Content (%)", 10, 100, 50)
vaf_input = st.sidebar.slider("Observed VAF (%)", 1, 100, 25)

# --- 数学モデル（5本の曲線） ---
f = tc_input / 100
models = {
    "Somatic Heterozygous": (f / 2) * 100,
    "Somatic LOH (deletion)": (f / (2 - f)) * 100,
    "Somatic cnLOH (UPD)": f * 100,
    "Germline Heterozygous": 50.0,
    "Germline LOH (deletion)": (1 / (2 - f)) * 100
}

# --- レイアウト構築 ---
col_alerts, col_graph = st.columns([1, 2])

with col_alerts:
    st.subheader("🚨 Clinical Alerts")
    
    # 1. The 50% VAF Trap (TC 60-75%)
    if 60 <= tc_input <= 75:
        st.warning("**Alert: The 50% VAF Trap.** Somatic LOH mimics germline heterozygous (VAF ≈ 50%).")

    # 2. Convergence Alert (TC >= 70%)
    if tc_input >= 70 and vaf_input >= models["Somatic LOH (deletion)"]:
        st.error("**⚠️ LOH Convergence Alert.** Somatic and Germline LOH are indistinguishable.")

    # 3. Mathematical Limit (TC >= 90%)
    if tc_input >= 90:
        st.info("**💡 Mathematical Convergence Zone.** VAF alone is insufficient for origin classification.")

    # 適合モデル表
    compatible_data = []
    for name, theory_vaf in models.items():
        if abs(vaf_input - theory_vaf) <= 10.0:
            compatible_data.append({"Model": name, "Theory": f"{theory_vaf:.1f}%"})
    if compatible_data:
        st.markdown("### 🔍 Compatible Models")
        st.table(pd.DataFrame(compatible_data))

    # --- Feature 3: Excel Workflow ---
    st.subheader("📊 Multi-variant Workflow")
    df_template = pd.DataFrame({"Variant": ["BRCA1", "TP53"], "TC": [tc_input, tc_input], "VAF": [vaf_input, 0.0]})
    buffer = io.BytesIO()
    df_template.to_csv(buffer, index=False)
    st.download_button(
        label="📥 Download Excel/CSV Template",
        data=buffer.getvalue(),
        file_name="VAF_TC_Case_Study.csv",
        mime="text/csv",
        help="Use this template for hypermutated tumors (Lynch/POLE)."
    )

with col_graph:
    st.subheader("📈 VAF-TC Projection")
    tr = np.linspace(10, 100, 100)
    fr = tr / 100
    fig = go.Figure()

    # 5本の線（README + cnLOH）
    fig.add_trace(go.Scatter(x=tr, y=(fr/2)*100, name="Somatic Het", line=dict(color='blue', width=1)))
    fig.add_trace(go.Scatter(x=tr, y=(fr/(2-fr))*100, name="Somatic LOH (Del)", line=dict(color='blue', dash='dash')))
    fig.add_trace(go.Scatter(x=tr, y=fr*100, name="Somatic cnLOH", line=dict(color='blue', dash='dot')))
    fig.add_trace(go.Scatter(x=tr, y=[50]*100, name="Germline Het", line=dict(color='green', width=1)))
    fig.add_trace(go.Scatter(x=tr, y=(1/(2-fr))*100, name="Germline LOH (Del)", line=dict(color='red', width=2)))
    
    fig.add_trace(go.Scatter(x=[tc_input], y=[vaf_input], mode='markers+text', name="CASE", text=["CASE"], marker=dict(color='black', size=15, symbol='x')))
    
    fig.update_layout(xaxis_title="Tumor Content (%)", yaxis_title="VAF (%)", legend=dict(orientation="h", y=-0.2))
    st.plotly_chart(fig, use_container_width=True)

# --- Therapeutic Implications (Full Text) ---
st.divider()
st.subheader("🩺 Therapeutic Implications & Clinical Notes")
c1, c2 = st.columns(2)
with c1:
    st.markdown("""
    **Hereditary Cancer Inference:**
    Syndromes like **HBOC, Lynch, and FAP** can be inferred when VAF aligns with two-hit models.
    
    **BRCA1/2-associated tumors:**
    - Ovarian/Prostate: Both gBRCA and sBRCA sensitivity to **PARPi**.
    - Breast/Pancreas: Generally **gBRCA Only** (includes Talazoparib for Breast).
    """)
with c2:
    st.markdown("""
    **Lynch Syndrome (MMR-d):**
    High responsiveness to **ICIs**. The curative potential of ICIs in Lynch syndrome is a critical clinical differentiator from epigenetic dMMR.
    """)

st.caption("Version 5.0 - README Full Integration. ✅ Clinical Genetics Suite")
