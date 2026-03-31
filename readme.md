# VAF-TC Relationship Visualizer 🧬

## Overview
The **VAF-TC Relationship Visualizer** is an interactive clinical tool designed to assist in the interpretation of genetic variants by modeling the mathematical relationship between **Pathological Tumor Content (TC)** and **Variant Allele Fraction (VAF)**. 

This tool helps clinicians and researchers evaluate the likelihood of germline vs. somatic events, providing a theoretical framework based on Knudson's Two-Hit Theory and various copy number alteration models.

## 🚀 Live Application
Access the interactive web tool here:
**[https://vaf-tc-app.streamlit.app/](https://vaf-tc-app.streamlit.app/)**

## 📂 Downloadable Data Resources
For users who prefer offline analysis, batch processing of multiple variants, or custom editing, the following raw data files are available in this repository:
- **[VAF_TC_theoretical_model.xlsx](./VAF_TC_theoretical_model.xlsx):** Excel format for easy editing and paper-based use.
- **[VAF_TC_theoretical_model.csv](./VAF_TC_theoretical_model.csv):** Raw data in CSV format.
- **[data_dictionary.txt](./data_dictionary.txt):** Detailed descriptions of all data columns.

## Key Features
* **Automated Interpretation:** Dynamically identifies theoretical models within a **±10% measurement error threshold**.
* **TMB-high Tumor Application:** Highly effective for interpreting hypermutated cases, including **MMRd** and **POLEm** tumors, as discussed in the associated study.
* **Convergence Zone (Gray Zone) Alert:** A warning system for samples where germline LOH and somatic LOH curves converge (typically TC 60–75%).
* **Pathological Integration:** Designed for use with **Pathological TC (%)** to ensure clinical reliability.

## Clinical Significance
Distinguishing between germline and somatic variants is critical in tumor-only sequencing. High-TC samples (TC ≥ 90%) with elevated VAFs are often misidentified as somatic events when they actually represent **Germline LOH**. 

Accurate identification of **Biallelic inactivation (LOH)** is therapeutically significant for identifying HBOC/Lynch Syndrome and determining sensitivity to **PARP inhibitors**, regardless of the variant's origin.

## How to Use
1. **Online:** Use the [Web App](https://vaf-tc-app.streamlit.app/) for quick, interactive visualization.
2. **Offline/Batch:** Download the `.xlsx` file to process multiple mutant genes in a single clinical case or for paper-based records.

## Citation
If you use this tool or data in your research, please cite:
*Integrated VAF-TC graph: a novel tool for differentiating germline from somatic variants in tumor-only sequencing* (In preparation/submission).

---
© 2024 Clinical-Genetics-Suite
