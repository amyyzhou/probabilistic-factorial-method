# Probabilistic Factorial Design Analysis of Drug Effects on Gene Expression

This repository contains the code and report for a project investigating the impact of various drug treatments and their combinations on gene expression in the A375 human malignant melanoma cell line, utilizing a probabilistic factorial design approach.

## Table of Contents

* [Introduction](#introduction)
* [Key Features](#key-features)
* [Methodology Highlights](#methodology-highlights)
* [Results Summary](#results-summary)
* [Repository Structure](#repository-structure)
* [Setup and Running the Project](#setup-and-running-the-project)
    * [Prerequisites](#prerequisites)
    * [Data](#data)
    * [Installation](#installation)
    * [Compiling the LaTeX Report](#compiling-the-latex-report)
* [Citation](#citation)

## Introduction

Understanding how pharmaceutical compounds influence gene expression is crucial for drug discovery and personalized medicine. Traditional high-throughput screening often focuses on single-drug effects, but biological systems involve complex interactions, especially with combination therapies.

This project applies a "probabilistic factorial experimental design" method, as described by Shyamal et al. (2025), to systematically identify both individual (main) drug effects and their pairwise interaction effects on gene expression. The approach involves constructing a Fourier basis design matrix to represent drug treatments and then fitting Ridge regression models to quantify their influence on target genes.

## Key Features

* **Data Preprocessing:** Handles a subset of the Broad LINCS L1000 dataset for the A375 cell line.
* **Fourier Basis Design Matrix Construction:** Generates a design matrix capturing main and second-order (pairwise) drug interaction effects.
* **Ridge Regression Modeling:** Fits regularized linear regression models to predict gene expression based on drug treatments and combinations.
* **Coefficient Analysis:** Identifies the most influential individual drugs and drug combinations for specific gene expression modulation.
* **Comprehensive Report:** A detailed LaTeX report outlining the methodology, results, and discussion.

## Methodology Highlights

The analysis was performed on a subset of the GSE70138 Broad LINCS L1000 dataset, focusing on the A375 human malignant melanoma cell line.

* **Data Subset:** A pre-processed HDF5 file (`GSE70138_A375_subset_expression.h5`) containing expression data for 978 landmark genes across 1000 randomly selected A375 active compound signatures was used.
* **Drug Selection:** 120 unique drug treatments were randomly selected for modeling.
* **Interaction Order:** The model included up to second-order interactions (`K_ORDER_INTERACTION = 2`), resulting in 7261 Fourier basis terms (features).
* **Gene Selection:** 10 landmark genes were randomly selected for detailed modeling.
* **Regression Model:** Ridge Regression (with `alpha=1.0`) was used to fit the models.

### Limitation
The R-squared and MSE metrics presented here rely on in-sample evaluation, a decision driven by two factors. 
1. The primary goal was to validate the implementation of the probabilistic factorial design method with its Fourier basis construction in a real-world context.
2. An initial misunderstanding of data compatibility between LINCS L1000 Phase I and II (due to differences in perturbation structure, preprocessing like MODZ vs. COMPZ, and treatment metadata) precluded straightforward external validation. This realization, alongside computational constraints, led us to prioritize internal consistency for this initial phase. Future work will crucially involve implementing robust cross-validation strategies to assess model generalizability.

## Results Summary

The Ridge regression models consistently achieved high R-squared values (ranging from 0.6698 to 0.9804), indicating strong predictive power of the drug features on gene expression changes. A key finding is the **predominance of pairwise drug-drug interaction effects** among the top estimated coefficients, highlighting the critical role of drug combinations.

Specific examples of influential drugs and their interactions include:

* **Tioguanine:** Frequently involved in interactions, both upregulating (e.g., NUP133, KIAA0196, BNIP3) and downregulating (e.g., NFE2L2, RNH1) gene expression.
* **TAS-103:** A common partner in strong interactions, often leading to downregulation (e.g., NFE2L2, CAMSAP2) but also involved in upregulation (e.g., BNIP3, SPTAN1).
* **LY-2090314:** Shows strong involvement in both upregulating (e.g., BNIP3) and downregulating (e.g., SPTAN1, CAMSAP2) effects.
* **Epirizole, Docetaxel, Lorglumide, Withaferin-a:** Also identified as key players in various synergistic and antagonistic interactions.

These findings provide data-driven hypotheses for further investigation into combination therapies for melanoma.

## Repository Structure
```bash
.
├── README.md                     # This README file
├── main.tex                  # Main LaTeX report document
├── references.bib            # Bibliography file for the LaTeX report
├── analysis_script.ipynb     # Python script for data processing and modeling
├── cut_gctx_data.py          # Python script to generate the subsetted HDF5 data
├── GSE70138_A375_subset_expression.h5 # Subsetted gene expression data (managed by Git LFS)
├── GSE70138_Broad_LINCS_sig_info_2017-03-06.txt # Metadata: Signature information
├── GSE70138_Broad_LINCS_pert_info_2017-03-06.txt # Metadata: Perturbagen (drug) information
├── GSE70138_Broad_LINCS_gene_info_2017-03-06.txt # Metadata: Gene information
├── GSE70138_Broad_LINCS_cell_info_2017-04-28.txt # Metadata: Cell line information
├── .gitignore                # Specifies files/folders to ignore (e.g., original .gctx)
└── .gitattributes            # Git LFS configuration file
```

## Setup and Running the Project

### Prerequisites

* Python 3.x
* Git
* Git Large File Storage (Git LFS)
* A LaTeX distribution (e.g., TeX Live, MiKTeX) or an online LaTeX editor (e.g., Overleaf)

### Data

The primary dataset used in this project is the **Broad LINCS L1000 dataset (GSE70138)**. The full dataset can be accessed and downloaded from the Gene Expression Omnibus (GEO) or the LINCS Data Portal:

* **GEO Accession:** [GSE70138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138)

**Important**: Due to its large size (~5.8 GB), the original raw gene expression data (`GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx`) is not included in this repository. While a smaller subset (`GSE70138_A375_subset_expression.h5`) is committed via Git LFS, for reliable operation, it is highly recommended to generate this subsetted file locally.

To do this:
1. **Download the original GCTX file:**
Obtain the ```GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz``` file from one of these sources:
* **GEO Accession:** [GSE70138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138) -> Look for the GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz file under "SRA Run Selector" or "FTP link"

2. **Unzip the downloaded file:**
The file is gzipped (```.gz```). Unzip it to get ```GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx```.

3. **Place the unzipped ```.gctx``` file:**
Move the unzipped ```GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx``` file into the directory of your cloned repository.

### Installation

Follow these steps to set up the project locally:

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/amyyzhou/probabilistic-factorial-method.git
    ```
2.  **Navigate into the project directory:**
    ```bash
    cd "probabilistic-factorial-method"
    ```
3.  **Install Git LFS** (if you haven't already):
    Follow instructions at [https://git-lfs.com/](https://git-lfs.com/). Then run:
    ```bash
    git lfs install
    ```
4.  **Pull Git LFS files:**
    Ensure the large data file (`.h5`) is downloaded.
    ```bash
    git lfs pull
    ```
5.  **Set up Python virtual environment and install dependencies:**
    It's highly recommended to use a virtual environment to manage project dependencies.
    ```bash
    python3 -m venv venv
    source venv/bin/activate  # On Windows, use: .\venv\Scripts\activate
    pip install pandas numpy scikit-learn
    ```

### Compiling the LaTeX Report
If you wish to view report on web, please visit [Overleaf](https://www.overleaf.com/read/jqsmhnhrgfyj#ff08aa).

The project includes a detailed report in LaTeX format (main.tex). Make sure you are in the Probabilistic Factorial Design Project directory.
1. Ensure references.bib is present: This file contains the bibliography entry and should be in the same directory as main.tex.
2. Compile the LaTeX document: For proper compilation of citations and bibliography using biblatex and biber, you need to run a specific sequence of commands.
- If using a local LaTeX distribution (like TeX Live or MiKTeX), navigate to the Probabilistic Factorial Design Project directory in your terminal and run:
pdflatex main.tex
biber main
pdflatex main.tex
pdflatex main.tex
- If using Overleaf (recommended), simply click the "Recompile" or "Build & View" button 2-3 times after making any changes to ensure all references are resolved.
The compiled PDF report (main.pdf) will be generated in the same directory.

## Citation
If you find this work useful, please consider citing the original paper:
```bash
@misc{shyamal2025probabilisticfactorialexperimentaldesign,
      title={Probabilistic Factorial Experimental Design for Combinatorial Interventions}, 
      author={Divya Shyamal and Jiaqi Zhang and Caroline Uhler},
      year={2025},
      eprint={2506.03363},
      archivePrefix={arXiv},
      primaryClass={cs.LG},
      url={https://arxiv.org/abs/2506.03363}, 
}
```
