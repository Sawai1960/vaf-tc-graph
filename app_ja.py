import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import io

# 1. ページ設定
st.set_page_config(page_title="VAF-TC 精密解析ツール", layout="wide")

# 2. タイトル
st.title("🧬 VAF-TC 精密解析ツール")
st.markdown("腫瘍単独シーケンシングにおける生殖細胞系列・体細胞変異の鑑別を支援するインタラクティブ可視化ツール。")
st.caption("⚠️ 本ツールは遺伝カウンセリングの補助ツールです。確認的な生殖細胞系列検査や確立された臨床ガイドラインの代替とはなりません。")
st.caption("⚠️ 遺伝子リファレンスシステムは **厚生労働科学研究費補助金 がん遺伝子パネル検査におけるGPV/PGPV対応手順に関する指針(2025版)** のT-only検査PGPV開示推奨遺伝子リストに準拠しています。")

# 重要な注意事項(起動時)
with st.expander("⚠️ 重要な注意事項 ― 数理モデルの前提と限界", expanded=True):
    st.markdown("""
**本ツールで可視化される5つの理論的VAF-TC曲線は、二倍体モデルにおける *Knudsonの二段階発癌説* から導出されています。臨床解釈の前に以下の前提を必ず確認してください：**

1. **二倍体仮定(Knudsonの二段階発癌説)**：各理論線は、二倍体(2コピー)をベースラインとした特定のバイアレリック不活化シナリオを表しています。5つのモデルは理想化された数学的リファレンスであり、すべての可能な機構を網羅するものではありません。

2. **異数性は考慮されない**：実腫瘍ではしばしば異数性、染色体全体のゲイン/ロス、サブクローン異質性を示します。これらの場合、観測されるVAFは理論線から大きく乖離する可能性があり、モデルマッチングは慎重に解釈する必要があります。

3. **TC推定には ±10〜20% の誤差がある**：病理学的腫瘍含有率推定は、組織学的異質性・サンプリング部位・観察者間一致度により、通常 **±10〜20% の変動** を伴います。本ツールのモデルマッチングには ±10% マージンを適用していますが、臨床解釈では全体の不確実性を考慮する必要があります。

4. **診断ツールではありません**：本ツールは遺伝カウンセリングの視覚的補助ツールです。臨床判断には確認的な生殖細胞系列検査が標準です。
""")

# 3. 遺伝子リファレンスデータ ― がん遺伝子パネル検査におけるGPV/PGPV対応手順に関する指針(2025版)準拠
GENE_INFO = {
    # 🔴 低VAF閾値グループ(VAF ≥ 10%) ― 特別扱い
    "BRCA1":  ("low_threshold", "🔴 BRCA1 [GPV/PGPV対応指針 2025版]：**VAF ≥ 10%** の低閾値（標準より低い）。低VAFでもGPVの可能性があります。エキスパートパネル検討と生殖細胞系列確認検査を推奨。HBOC。"),
    "BRCA2":  ("low_threshold", "🔴 BRCA2 [GPV/PGPV対応指針 2025版]：**VAF ≥ 10%** の低閾値（標準より低い）。低VAFでもGPVの可能性があります。エキスパートパネル検討と生殖細胞系列確認検査を推奨。HBOC。"),
    # 🟠 年齢条件付きグループ(SNV VAF ≥ 30% / Indel ≥ 20% かつ 発症年齢 < 30歳)
    "APC":    ("age_cond", "🟠 APC [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)かつ大腸ポリポーシス発症 < 30歳** で開示推奨。家族性大腸腺腫症(FAP)。Box_E：表現型評価が必要。"),
    "CDKN2A": ("age_cond", "🟠 CDKN2A [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)かつ発症年齢 < 30歳** で開示推奨。遺伝性黒色腫・膵がん症候群。"),
    "PTEN":   ("age_cond", "🟠 PTEN [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)かつ発症年齢 < 30歳** で開示推奨。Cowden症候群。Box_E：表現型評価が必要。"),
    "RB1":    ("age_cond", "🟠 RB1 [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)かつ発症年齢 < 30歳** で開示推奨。遺伝性網膜芽細胞腫。Box_E：表現型評価が必要。"),
    "TP53":   ("age_cond", "🟠 TP53 [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)かつ発症年齢 < 30歳** で開示推奨。Li-Fraumeni症候群。Box_E：表現型評価が必要。注：高VAFではクローン性造血の可能性も考慮。"),
    # 🟡 標準グループ(SNV VAF ≥ 30% / Indel ≥ 20%) ― 24遺伝子
    "ATM":    ("standard", "🟡 ATM [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。HBOC・毛細血管拡張性運動失調症。"),
    "BAP1":   ("standard", "🟡 BAP1 [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。BAP1腫瘍素因症候群。"),
    "BARD1":  ("standard", "🟡 BARD1 [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。HBOC。"),
    "BRIP1":  ("standard", "🟡 BRIP1 [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。HBOC。"),
    "CHEK2":  ("standard", "🟡 CHEK2 [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。HBOC。"),
    "DICER1": ("standard", "🟡 DICER1 [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。DICER1症候群(胸膜肺芽腫など)。"),
    "FH":     ("standard", "🟡 FH [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。HLRCC(遺伝性平滑筋腫症・腎細胞がん)。"),
    "FLCN":   ("standard", "🟡 FLCN [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。Birt-Hogg-Dubé症候群。"),
    "MLH1":   ("standard", "🟡 MLH1 [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。Lynch症候群。"),
    "MSH2":   ("standard", "🟡 MSH2 [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。Lynch症候群。"),
    "MSH6":   ("standard", "🟡 MSH6 [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。Lynch症候群。"),
    "MUTYH":  ("standard", "🟡 MUTYH [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)・両アレル病的バリアントのみ** で開示推奨。MUTYH関連ポリポーシス。"),
    "NF1":    ("standard", "🟡 NF1 [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。神経線維腫症1型。Box_E：表現型評価が必要。"),
    "PALB2":  ("standard", "🟡 PALB2 [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。HBOC。"),
    "PMS2":   ("standard", "🟡 PMS2 [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。Lynch症候群。"),
    "POLD1":  ("standard", "🟡 POLD1 [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。遺伝性大腸がん。"),
    "POLE":   ("standard", "🟡 POLE [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。遺伝性大腸がん。"),
    "RAD51C": ("standard", "🟡 RAD51C [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。HBOC。"),
    "RAD51D": ("standard", "🟡 RAD51D [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。HBOC。"),
    "RET":    ("standard", "🟡 RET [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。多発性内分泌腫瘍症2型(MEN2)。"),
    "SDHA":   ("standard", "🟡 SDHA [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。遺伝性パラガングリオーマ・褐色細胞腫症候群。"),
    "SDHB":   ("standard", "🟡 SDHB [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。遺伝性パラガングリオーマ・褐色細胞腫症候群。"),
    "TSC2":   ("standard", "🟡 TSC2 [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。結節性硬化症。"),
    "VHL":    ("standard", "🟡 VHL [GPV/PGPV対応指針 2025版]：**SNV VAF ≥ 30%(Indel ≥ 20%)** で開示推奨。フォン・ヒッペル・リンダウ症候群。"),
}

def get_gene_message(gene):
    key = gene.upper().strip()
    if key in GENE_INFO:
        return GENE_INFO[key]
    return ("not_listed", f"⬜ {gene}：GPV/PGPV対応手順指針(2025版)のT-only検査PGPV開示推奨遺伝子リストに含まれていません。本指針における生殖細胞系列開示の優先度は低く、体細胞変異の可能性が高いと考えられます。臨床ガイドラインと家族歴による総合的な評価が必要です。")

# 4. サイドバー入力
st.sidebar.header("📋 患者データ入力")
st.sidebar.markdown("👉 **遺伝子名・TC・VAFを入力してください。**")

gene_name = st.sidebar.text_input("遺伝子名", value="BRCA2")
tc_input = st.sidebar.slider("病理学的腫瘍含有率(TC %)", 0, 100, 50)
vaf_input = st.sidebar.slider("変異アレル頻度(VAF %)", 0, 100, 50)

st.sidebar.markdown("---")
st.sidebar.caption("⚠️ 下記でCSVをアップロードすると、上記の入力(遺伝子名・TC・VAF)は無効になります。")

# 複数変異CSVアップロード
st.sidebar.subheader("📂 複数変異アップロード")
st.sidebar.caption("CSV形式：Gene, TC, VAF")
uploaded_file = st.sidebar.file_uploader("CSVをアップロード", type=["csv"])

multi_df = None
if uploaded_file is not None:
    try:
        multi_df = pd.read_csv(uploaded_file)
        required_cols = {"Gene", "TC", "VAF"}
        if not required_cols.issubset(multi_df.columns):
            st.sidebar.error("CSVにはGene・TC・VAF列が必要です。")
            multi_df = None
        else:
            st.sidebar.success(f"{len(multi_df)} 件の変異を読み込みました。")
    except Exception as e:
        st.sidebar.error(f"CSV読み込みエラー：{e}")
        multi_df = None

st.sidebar.markdown("---")
# CSVテンプレートダウンロード(サイドバー)
st.sidebar.subheader("📊 複数変異ワークフロー")
st.sidebar.caption("💡 下記のテンプレートをダウンロードし、サンプル遺伝子を自分のデータに書き換えてから、上の「複数変異アップロード」でアップロードしてください。")
template_df = pd.DataFrame({
    "Gene": [gene_name, "TP53", "MSH2"],
    "TC":   [tc_input,  tc_input, tc_input],
    "VAF":  [vaf_input, 0.0,      0.0]
})
csv_string = template_df.to_csv(index=False)
st.sidebar.download_button("📥 CSVテンプレートをダウンロード", csv_string.encode("utf-8"), "VAF_TC_Template.csv", "text/csv")

# 理論モデルデータダウンロード(サイドバー)
st.sidebar.subheader("📂 理論モデルデータ")
try:
    with open("VAF_TC_theoretical_model.csv", "rb") as f:
        st.sidebar.download_button("📥 理論モデルをダウンロード(CSV)", f.read(), "VAF_TC_theoretical_model.csv", "text/csv")
except FileNotFoundError:
    st.sidebar.caption("VAF_TC_theoretical_model.csv が見つかりません。")
try:
    with open("VAF-TC theoretical_model.xlsx", "rb") as f:
        st.sidebar.download_button("📥 理論モデルをダウンロード(Excel)", f.read(), "VAF-TC_theoretical_model.xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
except FileNotFoundError:
    st.sidebar.caption("VAF-TC theoretical_model.xlsx が見つかりません。")

tc = tc_input / 100.0
vaf = vaf_input / 100.0

# 5. 数理モデル(二倍体モデル) ― 5つの理論モデル
x_range = np.linspace(0.01, 1.0, 100)
y_germ_cnloh  = (1 + x_range) / 2         # germline (cnLOH)
y_germ_del    = 1 / (2 - x_range)         # germline (LOH with Del)
y_germ_hetero = np.full_like(x_range, 0.5) # germline (Hetero)
y_som_del     = x_range / (2 - x_range)   # somatic (LOH with Del)
y_som_hetero  = x_range / 2               # somatic (Hetero)

# 6. メインレイアウト(左1：右2)
col_alerts, col_graph = st.columns([1, 2])

# --- 左カラム：解釈とアラート ---
with col_alerts:
    st.info(f"💡 **解析モード：** {gene_name}")
    st.subheader("📋 解釈とアラート")

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
            return "warning", "VAFはいずれの標準モデルとも一致しません。クローン異質性・異数性・複雑なコピー数変化を考慮してください。"
        germ_cnloh  = "germline (cnLOH)"        in names
        germ_del    = "germline (LOH with Del)" in names
        germ_hetero = "germline (Hetero)"       in names
        som_del     = "somatic (LOH with Del)"  in names
        som_hetero  = "somatic (Hetero)"        in names

        has_germ = germ_cnloh or germ_del or germ_hetero
        has_som  = som_del or som_hetero

        # 両立：生殖細胞系列・体細胞の両方が適合
        if has_germ and has_som:
            return "error", "VAFは **生殖細胞系列・体細胞の両方** のモデルに適合します。このTCではVAFのみで起源を判別できません。**ペア正常検体による生殖細胞系列検査が必須です。**"
        # 生殖細胞系列 ― 単一モデル
        if germ_cnloh and not germ_del and not germ_hetero:
            return "success", "**コピー数中立LOH(UPD)を伴う生殖細胞系列変異** のパターンと一致します。生殖細胞系列 + cnLOHによるバイアレリック不活化。"
        if germ_del and not germ_cnloh and not germ_hetero:
            return "success", "**欠失によるLOHを伴う生殖細胞系列変異** のパターンと一致します。生殖細胞系列 + 欠失によるバイアレリック不活化。"
        if germ_hetero and not germ_cnloh and not germ_del:
            return "success", "**LOHを伴わないヘテロ接合性生殖細胞系列変異** のパターンと一致します。1つのアレルのみが影響を受けています。"
        # 複数の生殖細胞系列モデル
        if has_germ:
            return "info", "複数の生殖細胞系列モデルが該当します。生殖細胞系列起源の可能性が高いですが、LOH機構はVAFのみでは判別できません。"
        # 体細胞 ― 単一モデル
        if som_del and not som_hetero:
            return "info", "**欠失によるLOHを伴う体細胞変異** のパターンと一致します。このTCでは生殖細胞系列の可能性は低いです。"
        if som_hetero and not som_del:
            return "info", "**LOHを伴わないヘテロ接合性体細胞変異** のパターンと一致します。このTCでは生殖細胞系列の可能性は低いです。"
        # 複数の体細胞モデル
        if has_som:
            return "info", "複数の体細胞モデルが該当します。このTCでは体細胞起源の可能性が高いです。"
        return "info", "臨床的文脈との照合とペア正常検査を推奨します。"

    def show_variant_interpretation(g, t, v):
        compatible = get_compatible_models(t, v)
        st.markdown(f"**{g}**(TC {t:.0f}%、VAF {v:.0f}%)")
        if compatible:
            for name, val in compatible:
                st.markdown(f"- **{name}** ― 理論VAF {val*100:.1f}%")
        level, msg = get_interpretation(compatible)
        if level == "success":
            st.success(f"➡️ {msg}")
        elif level == "error":
            st.error(f"➡️ {msg}")
        elif level == "warning":
            st.warning(f"➡️ {msg}")
        else:
            st.info(f"➡️ {msg}")
        # TCベースアラート
        if t <= 20:
            st.warning("⚠️ **低TC(≤ 20%)：** 低TCでは理論線が狭いVAF範囲に圧縮され、モデルマッチングの信頼性が低下します。サブクローン変異・正常組織の混入・技術的ノイズが支配的になる可能性があります。")
        if t >= 60:
            st.warning("⚠️ **高TC(≥ 60%)：** 高TCでは生殖細胞系列と体細胞のLOH線が収束し始めます。VAFのみでの起源判別が困難になります。")
        # 遺伝子別メッセージ(GPV/PGPV対応指針 2025版)
        _, gene_msg = get_gene_message(g)
        st.info(gene_msg)
        st.caption("💡 注：この解釈はVAF-TCの数理モデルのみに基づいています。遺伝子別の生殖細胞系列確率はGPV/PGPV対応手順指針(2025版)に準拠して別途提示しています。")

    if multi_df is not None:
        for _, row in multi_df.iterrows():
            show_variant_interpretation(str(row["Gene"]), float(row["TC"]), float(row["VAF"]))
            st.divider()
    else:
        show_variant_interpretation(gene_name, tc_input, vaf_input)

    # --- TCベースの臨床アラート ---
    som_del_vaf  = tc / (2 - tc) * 100 if tc < 2 else 0
    germ_del_vaf = 1 / (2 - tc) * 100 if tc < 2 else 0

    if 61 <= tc_input <= 66:
        st.warning(
            f"⚠️ **グレーゾーン(Somatic LOH Del)：** TC {tc_input}%では、"
            f"somatic (LOH with Del) の理論VAF = {som_del_vaf:.1f}%となり、"
            f"germline (Hetero) の50%に接近します。確認検査を推奨します。"
        )
    elif tc_input >= 67:
        if vaf_input >= tc / (2 - tc) * 100:
            st.error(
                f"🔴 **LOH収束アラート：** TC {tc_input}%、VAF {vaf_input}%では、"
                f"変異が somatic (LOH with Del) ライン({som_del_vaf:.1f}%)以上に位置します。"
                f"この領域では germline (LOH with Del) = {germ_del_vaf:.1f}% と "
                f"somatic (LOH with Del) = {som_del_vaf:.1f}% が収束し、"
                f"VAFのみでは起源を判別できません。生殖細胞系列確認検査が必須です。"
            )
        if tc_input >= 90:
            st.warning(
                f"⚠️ **極高腫瘍純度：** TC {tc_input}%では、すべての理論モデルが"
                f"狭いVAF範囲に圧縮されます。高VAFでも体細胞起源の可能性があります。"
                f"家族歴の確認と生殖細胞系列検査が必須です。"
            )

    st.divider()

    # 複数変異一覧表示
    if multi_df is not None:
        st.subheader("📋 アップロードされた変異")
        st.dataframe(multi_df, use_container_width=True)
        st.divider()


# --- 右カラム：グラフ ---
with col_graph:
    st.subheader("📈 VAF-TC 投影グラフ")
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
            name=f"現在：{gene_name}",
            text=[f"{gene_name}<br>TC:{tc_input}%<br>VAF:{vaf_input}%"],
            textposition="top right",
            marker=dict(color='black', size=14, symbol='circle')
        ))

    fig.add_vrect(x0=0, x1=20, fillcolor="gray", opacity=0.1, layer="below", line_width=0,
                  annotation_text="低信頼ゾーン", annotation_position="top left")

    fig.update_layout(
        xaxis_title="病理学的腫瘍含有率(%)", yaxis_title="変異アレル頻度(%)",
        yaxis=dict(range=[0, 105]), xaxis=dict(range=[0, 105]),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        template="simple_white", height=600
    )
    st.plotly_chart(fig, use_container_width=True)

    # 遺伝子リファレンス表(GPV/PGPV対応指針 2025版)
    st.subheader("📖 遺伝子リファレンス(GPV/PGPV対応指針 2025版)")
    st.markdown("**🔴 低VAF閾値(VAF ≥ 10%)：**")
    st.caption("BRCA1, BRCA2")
    st.markdown("**🟠 年齢条件付き(発症年齢 < 30歳)：**")
    st.caption("APC(大腸ポリポーシス), CDKN2A, PTEN, RB1, TP53")
    st.markdown("**🟡 標準(SNV VAF ≥ 30%, Indel ≥ 20%)：**")
    st.caption("ATM, BAP1, BARD1, BRIP1, CHEK2, DICER1, FH, FLCN, MLH1, MSH2, MSH6, MUTYH(bi), NF1, PALB2, PMS2, POLD1, POLE, RAD51C, RAD51D, RET, SDHA, SDHB, TSC2, VHL")
    st.caption("⬜ リスト外の遺伝子：2025年版のT-only PGPVリストに含まれません。")
    st.caption("出典：厚生労働科学研究費補助金 がん遺伝子パネル検査におけるGPV/PGPV対応手順に関する指針(2025版)")

# 7. フッター
st.divider()
st.caption("VAF-TC 精密解析ツール | Clinical Genetics Suite | ver 3.4 ✅")
