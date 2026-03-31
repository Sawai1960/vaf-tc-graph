# VAF-TC Precision Analyzer: Clinical Genetics Support Tool

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://vaf-tc-app.streamlit.app/)

## 🧬 Overview
**VAF-TC Precision Analyzer** is a clinical decision-support tool designed to differentiate between somatic and germline variants by modeling the mathematical relationship between **Pathological Tumor Content (TC)** and **Variant Allele Frequency (VAF)**. 

By applying theoretical lines based on **Knudson’s two-hit model**, this tool facilitates the identification of hereditary cancer syndromes and helps interpret the clinical significance of variants in high-purity or hypermutated tumor samples.

---

## 🚀 Key Clinical Features

### 1. The "50% VAF Trap" Alert
In clinical NGS, VAF ≈ 50% is often incorrectly assumed to be a germline variant. 
- **The Intersection:** At TC ≈ 66.7%, a **Somatic LOH (deletion)** event yields a theoretical VAF of **50%**.
- **Dynamic Alert:** The app warns users when TC is in the 60-75% range, preventing misidentifying somatic drivers as hereditary findings.

### 2. Convergence Alert (TC ≥ 70%)
We implemented an alert to indicate the potential convergence of germline LOH and somatic LOH when tumor content (TC) is ≥ 70% and the variant allele frequency (VAF) is at or above the theoretical line for somatic LOH with deletion. This specific alert addresses the clinical "Grey Zone" where somatic events can be indistinguishable from germline findings.

### 3. Support for Hypermutated Tumors (Practical Excel Workflow)
Particularly useful for interpreting cases with high mutational burdens, such as **MMR-deficient (Lynch syndrome)** and **POLE-mutant** tumors. 
- **Practical Marking:** Users can download a dedicated Excel template to plot and manually mark multiple variants for a single patient case. 
- **Clinical Records:** This facilitates the creation of a visual, case-by-case record for complex clinical documentation and paper-based reporting.

### 4. Mathematical Limit Analysis (TC ≥ 90%)
In high-purity samples (TC ≥ 90%), the app flags a "Mathematical Convergence Zone." As TC approaches 100%, the theoretical difference between somatic and germline LOH models falls within the margin of sequencing error. The tool reminds users that definitive classification in this range requires clinical correlation rather than VAF alone.

---

## 🩺 Clinical Significance

### Inferring Hereditary Cancer Syndromes
Hereditary cancers driven by tumor suppressor gene alterations (e.g., **HBOC, Lynch Syndrome, and FAP**) can be efficiently inferred using theoretical lines based on Knudson’s two-hit model. The ability to visualize whether an observed VAF aligns with a "two-hit" theoretical line provides strong supportive evidence for the variant's clinical relevance.

### Therapeutic Implications
- **BRCA1/2-associated tumors:** May indicate sensitivity to **PARP inhibitors**.
- **Lynch Syndrome (MMR-d):** May predict responsiveness to **Immune Checkpoint Inhibitors (ICIs)**.
*Note: In Lynch syndrome, the curative potential with ICIs is a significant clinical observation, differentiating these cases from epigenetic dMMR tumors.*

---

## 🌐 Live Application
👉 **[https://vaf-tc-app.streamlit.app/](https://vaf-tc-app.streamlit.app/)**

---

## 📊 Mathematical Foundation
The app utilizes the following frameworks ($f$ = Tumor Fraction):
- **Somatic Heterozygous:** $VAF = f / 2$
- **Somatic LOH (Deletion):** $VAF = f / (2 - f)$
- **Germline Heterozygous:** $VAF = 0.5$
- **Germline LOH (Deletion):** $VAF = 1 / (2 - f)$

---

## 👥 Authors & Contribution
- **Organization:** Clinical Genetics Suite
- **Maintainer:** Sawai1960
