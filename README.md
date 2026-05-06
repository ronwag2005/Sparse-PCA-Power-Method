# Sparse-PCA-Power-Method: Dimensionality Reduction for Sparse Movie Ratings Data

## Overview
This repository contains a manual implementation of the **Simultaneous Power Method** used to perform Principal Component Analysis (PCA) on the **MovieLens 100K** dataset. The project focuses on recovering the latent structures of movie preferences while navigating the technical constraints of a **93.7% sparse ratings matrix**. Code for the efficient **Simultaneous Power Method** has been made accessible and open-source for implementing further projects.

 **[View the Interactive PCA Loadings Plot](http://rpubs.com/ronwag2005/1429898)**

## Key Technical Achievements
* **Implicit Centering Algorithm:** Developed a memory-efficient two-step "Forward-Backward" pass logic that performs PCA on centered data without explicitly densifying the sparse matrix.
* **Subspace Stability:** Manually utilized the **Gram-Schmidt** process within the power iteration loop to prevent subspace collapse, successfully recovering the first 10 principal components.
* **Mathematical Optimization:** By avoiding explicitly centering the data matrix ($X - \mathbf{1}\mu^T$), the implementation maintains $O(\text{nnz})$ computational complexity, allowing for analysis of large, sparse data matrices on standard consumer hardware.
* **Gender-Stratified Subspace Analysis:** Conducted independent PCA on demographic submatrices, verifying the structural stability of the spectral gap across isolated male and female cohorts.

## Analytical Insights & Demographic Inferences
* **Dimensionality Reduction:** Condensed 1,682 movie dimensions into a low-rank space that explains over **35% of total market variation**.
* **Latent Market Axes:** Identified the fundamental drivers of cinematic preference: a primary "Market Engagement" axis (PC1) capturing mainstream saturation, and a secondary "Genre Pivot" (PC2) separating visceral blockbusters from critical classics.
* **Age & Preference Evolution:** Identified a significant **0.369 correlation** between user age and the second principal component (PC2). This suggests a systematic shift toward "Classic/Critically Acclaimed" titles as users age.
* **Gender-Based Divergence:** Utilizing **two-tailed t-tests**, the analysis identified a highly significant divergence ($p < 0.001$) in user engagement volume (PC1) between genders, while confirming that taste distributions along the stylistic Genre-based axis (PC2) remain statistically indistinguishable.

## Methodology: The Power Method
The first k (here k=10) Singular Vectors and Singular Values are recovered by approximating the centered ratings matrix $X$ using SVD: $X \approx U \Sigma V^{T}$.

| Step | Operation | Purpose |
| :--- | :--- | :--- |
| **Forward Pass** | $U_{tmp} = X Q - \mathbf{1}(\mu^T Q)$ | Projecting subspace while implicitly centering |
| **Backward Pass** | $V_{raw} = X^T U_{tmp} - \mu(\mathbf{1}^T U_{tmp})$ | Recovering movie loadings |
| **Orthonormalization** | Gram-Schmidt Process | Preventing subspace collapse into a single vector |

## Repository Structure
```text
├── README.md                   # Project overview and methodology
├── main_analysis.R             # Core R script containing the Power Method and visualizations
├── ml-100k/                    # MovieLens 100K dataset directory (requires download)
│   ├── u.data                  # Sparse ratings data
│   ├── u.user                  # User demographics
│   ├── u.item                  # Movie metadata
│   └── u.genre                 # Genre mappings
├── outputs/                    # Exported visualizations and interactive HTML
│   ├── PCA_Interactive_Loadings.html
│   ├── Female_Loadings_Plot.png
│   ├── Male_Loadings_Plot.png
│   └── Gender_Diagnostics_Plot.png
└── report/                     # LaTeX documentation
    └── Final_Report.tex        # Full technical analysis, equations, and interpretations
 ```
## Getting Started

### Dependencies

This project requires **R** and the following libraries:

* `Matrix` (for sparse matrix operations)[cite: 1]

* `ggplot2` & `plotly` (for static and interactive visualization)[cite: 1]

* `dplyr` & `stringr` (for data wrangling)[cite: 1]



### Data Setup

1. Download the MovieLens 100K dataset from [GroupLens](https://grouplens.org/datasets/movielens/100k/).

2. Ensure the `ml-100k` folder is placed in your working directory.

***

**Author:** Rohan Wagle  

**Date:** April 2026  

**Institution:** Ashoka University

**Acknowledgement:** This project has been done as the final project for the Applied Statistics 1 undergraduate course at Ashoka University in the Spring 2026 semester.
AI tools were utilized during this project to assist with code structuring, debugging syntax errors, and formatting LaTeX outputs. All mathematical logic, analytical interpretations, and dataset subsets were independently directed and verified.
