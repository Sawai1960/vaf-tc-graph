# VAF-TC Precision Analyzer

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://vaf-tc-app.streamlit.app/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An interactive visual tool for differentiating germline and somatic variants in **tumor-only sequencing**, based on the mathematical relationship between pathological **Tumor Content (TC)** and **Variant Allele Fraction (VAF)**.

> **Disclaimer:** This tool is intended as a supportive aid for genetic counseling. It does not replace confirmatory germline testing or established clinical guidelines. Further prospective validation is required.
> Gene Reference System is based on the **Guidelines for GPV/PGPV Handling Procedures in Cancer Gene Panel Testing (2025 Edition)** (MHLW Research Grant) — T-only PGPV disclosure recommended gene list.

## Live Application

**https://vaf-tc-app.streamlit.app/**

## Background

In tumor-only comprehensive genomic profiling (CGP), distinguishing germline from somatic variants is a fundamental challenge. While VAF of approximately 50% is often assumed to indicate a germline heterozygous variant, **somatic variants with LOH can produce the same VAF** depending on tumor content — a diagnostic trap.

This tool visualizes five theoretical VAF-TC models derived from **Knudson's two-hit hypothesis** (diploid model) and provides automated clinical alerts for known ambiguity zones.

## Important Note — Mathematical Model Assumptions

1. **Diploid assumption (Knudson's two-hit model):** Each theoretical line represents a specific biallelic inactivation scenario assuming a diploid (2-copy) baseline. The five models are idealized mathematical references, not an exhaustive catalog of all possible mechanisms.
2. **Aneuploidy is not accounted for:** Real tumors frequently exhibit aneuploidy, whole-chromosome gains/losses, and subclonal heterogeneity. Observed VAF may deviate substantially from the theoretical lines.
3. **TC estimation carries +/-10-20% error:** Pathological TC estimation is subject to +/-10-20% variability due to histological heterogeneity, sampling region, and inter-observer agreement. A +/-10% matching margin is applied, but clinical interpretation should consider the full uncertainty range.
4. **Not a diagnostic tool:** This is a visual aid for genetic counseling. Confirmatory germline testing remains the standard for any clinical decision.

## Mathematical Models

Given tumor content *f* (0-1):

| Model | Formula | Description |
|-------|---------|-------------|
| germline (cnLOH) | VAF = (1 + f) / 2 | Germline variant with copy-neutral LOH (UPD) |
| germline (LOH with Del) | VAF = 1 / (2 - f) | Germline variant with LOH by deletion |
| germline (Hetero) | VAF = 0.5 | Germline heterozygous variant without LOH |
| somatic (LOH with Del) | VAF = f / (2 - f) | Somatic variant with LOH by deletion |
| somatic (Hetero) | VAF = f / 2 | Somatic heterozygous variant without LOH |

A **+/-10% error margin** is applied for model matching to account for variability in pathological TC estimation.

## Clinical Alert System

The app generates five context-dependent alerts based on TC and VAF:

### Alert 1 — Gray Zone (TC 61-66%)

As TC increases toward 66.7%, somatic (LOH with Del) = f/(2-f) approaches 50% from below. In this range, the somatic LOH deletion line is close enough to 50% to create ambiguity with germline (Hetero).

### Alert 2 — LOH Convergence Zone (TC >= 67%)

At TC = 2/3 (approximately 66.7%), somatic (LOH with Del) **exactly equals** germline (Hetero) at 50%. Above this TC, the somatic and germline LOH lines converge. This alert fires when:

- **TC >= 67%**, AND
- **VAF >= somatic (LOH with Del) line** at current TC

Both conditions must be met. The alert displays the actual theoretical values for clinical reference.

### Alert 3 — Extreme Tumor Purity (TC >= 90%)

At very high purity, all five theoretical models compress into a narrow VAF range. Variants may still be of somatic origin even at high VAF. This alert fires regardless of VAF, as germline testing becomes essential in all cases.

### Alert 4 — Low TC (TC <= 20%)

At low tumor content, theoretical lines are compressed into a narrow VAF range and model matching is less reliable. Subclonal variants, admixture with normal tissue, or technical noise may dominate.

### Alert 5 — High TC (TC >= 60%)

At high tumor content, germline and somatic LOH lines begin to converge. Origin determination by VAF alone becomes increasingly difficult.

## Gene Reference System (GPV/PGPV Guidelines 2025 Edition)

The app provides gene-specific contextual messages based on the **Guidelines for GPV/PGPV Handling Procedures in Cancer Gene Panel Testing (2025 Edition)** (MHLW Research Grant), T-only PGPV disclosure recommended gene list (31 genes).

| Category | VAF Threshold | Genes | Notes |
|---|---|---|---|
| 🔴 Low threshold | VAF >= 10% | BRCA1, BRCA2 | Even low-VAF variants may be GPV. Expert panel review recommended. |
| 🟠 Age-conditional | SNV >= 30%, indel >= 20% AND onset < 30 y | APC, CDKN2A, PTEN, RB1, TP53 | Box_E: phenotype evaluation required for APC, PTEN, RB1, TP53 |
| 🟡 Standard | SNV >= 30%, indel >= 20% | ATM, BAP1, BARD1, BRIP1, CHEK2, DICER1, FH, FLCN, MLH1, MSH2, MSH6, MUTYH(bi), NF1, PALB2, PMS2, POLD1, POLE, RAD51C, RAD51D, RET, SDHA, SDHB, TSC2, VHL | MUTYH: bi-allelic only. NF1: Box_E phenotype evaluation. |
| ⬜ Not listed | — | All other genes | Not on 2025 T-only PGPV list. Consult clinical guidelines and family history. |

## Features

- **Interactive graph** with five theoretical VAF-TC curves (Plotly)
- **Model matching** with +/-10% error margin and theoretical VAF display
- **Automated interpretation** based on compatible model combinations
- **Gene-specific messages** for 31 genes per GPV/PGPV Guidelines (2025 Edition)
- **Five clinical alerts** based on TC values
- **Low Confidence Zone** shading for TC < 20%
- **Important Note** on startup with model assumptions and limitations
- **Multi-variant CSV upload** to plot multiple variants simultaneously on the graph
- **CSV template download** for multi-variant workflows
- **Theoretical model data download** (CSV and Excel) directly from the app

## Multi-variant Upload

Multiple variants from a single patient can be uploaded as a CSV file and plotted simultaneously on the graph. This is particularly useful for cases with high mutational burden (e.g., Lynch syndrome, POLE-mutant tumors).

**CSV format:**

```
Gene,TC,VAF
BRCA2,70,57
TP53,70,35
MSH2,70,68
```

Each variant is plotted with a distinct color and gene label. Interpretation and gene-specific messages are shown for each variant. A template CSV can be downloaded from within the app.

## Getting Started

### Requirements

- Python 3.9+
- Dependencies: streamlit, plotly, numpy, pandas

### Installation

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Repository Contents

| File | Description |
|------|-------------|
| app.py | Main Streamlit application (ver 3.4) |
| requirements.txt | Python dependencies |
| VAF-TC theoretical_model.xlsx | Excel file for generating theoretical VAF-TC curves |
| VAF_TC_theoretical_model.csv | CSV version of the theoretical model data |
| data_dictionary.txt | Variable definitions for the theoretical model |

## Changelog (ver 3.4)

- **Moved** Multi-variant Workflow download and Theoretical Model Data download to sidebar (vertically aligned with upload)
- **Moved** Gene Reference from sidebar to right column below graph (wider display, no longer collapsible)
- **Moved** Analysis Mode indicator from sidebar to top of left column (prominent banner)

## Changelog (ver 3.3)

- **Removed** Somatic + cnLOH model; **added** somatic (Hetero) = f/2
- **Renamed** all model labels to lowercase format: germline (cnLOH), germline (LOH with Del), germline (Hetero), somatic (LOH with Del), somatic (Hetero)
- **Deleted** Alert 1 (Somatic cnLOH Trap); alerts renumbered (5 total)
- **Changed** Low VAF / High VAF alerts to **Low TC / High TC** alerts
- **Changed** Low Confidence Zone from TC < 30% to **TC < 20%**
- **Added** Important Note on startup (Knudson assumptions, aneuploidy, TC estimation error)
- **Rebuilt** Gene Reference System per **MHLW GPV/PGPV Guidelines (2025 Edition)** (31 genes, 3 tiers)

## Citation

If you use this tool in your research, please cite:

> Kashima M, Tsubamoto H, et al. "VAF-Tumor Content Graph: A Simple Visual Tool for Discriminating Germline and Somatic Variants in Tumor-Only Sequencing." *Journal of Human Genetics* (submitted).

## Authors

**Clinical Genetics Suite** - Hyogo Medical University

## License

MIT License
