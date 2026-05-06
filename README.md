
# Sparse-PCA-Power-Method: Dimensionality Reduction for Sparse Movie Ratings Data

## Overview
This repository contains a manual implementation of the **Simultaneous Power Method** used to perform Principal Component Analysis (PCA) on the **MovieLens 100K** dataset[cite: 1]. The project focuses on recovering the latent structures of movie preferences while navigating the technical constraints of a **93.7% sparse ratings matrix**[cite: 1]. Code for the efficient **Simultaneous Power Method** has been made accessible and open-source for implementing further projects.

 **[View the Interactive PCA Loadings Plot](YOUR_HOSTED_LINK_HERE)**

## Key Technical Achievements
* **Implicit Centering Algorithm:** Developed a memory-efficient two-step "Forward-Backward" pass logic that performs PCA on centered data without explicitly densifying the sparse matrix[cite: 1].
* **Subspace Stability:** Manually utilized the **Gram-Schmidt** process within the power iteration loop to prevent subspace collapse, successfully recovering the first 10 principal components[cite: 1].
* **Mathematical Optimization:** By avoiding explicitly centering the data matrix ($X - \mathbf{1}\mu^T$), the implementation maintains $O(\text{nnz})$ computational complexity, allowing for analysis of large, sparse data matrices on standard consumer hardware[cite: 1].

## Analytical Insights & Demographic Inferences
* **Dimensionality Reduction:** Condensed 1,682 movie dimensions into a low-rank space that explains over **35% of total market variation**[cite: 1].
* **Age & Preference Evolution:** Identified a significant **0.369 correlation** between user age and the second principal component (PC2)[cite: 1]. This suggests a systematic shift toward "Classic/Critically Acclaimed" titles as users age[cite: 1].
* **Gender-Based Divergence:** Utilizing **Welch’s T-test**, the analysis identified a distinct divergence in user centroids along the latent preference axes, providing predictive power for gender-based movie recommendations[cite: 1].

## Methodology: The Power Method
The first k (here k=10) Singular Vectors and Singular Values are recovered by approximating the centered ratings matrix $X$ using SVD: $X \approx U \Sigma V^{T}$[cite: 1].

| Step | Operation | Purpose |
| :--- | :--- | :--- |
| **Forward Pass** | $U_{tmp} = X Q - \mathbf{1}(\mu^T Q)$ | Projecting subspace while implicitly centering[cite: 1] |
| **Backward Pass** | $V_{raw} = X^T U_{tmp} - \mu(\mathbf{1}^T U_{tmp})$ | Recovering movie loadings[cite: 1] |
| **Orthonormalization** | Gram-Schmidt Process | Preventing subspace collapse into a single vector[cite: 1] |

## Getting Started
### Dependencies
This project requires **R** and the following libraries:
* `Matrix` (for sparse matrix operations)[cite: 1]
* `ggplot2` & `plotly` (for static and interactive visualization)[cite: 1]
* `dplyr` & `stringr` (for data wrangling)[cite: 1]

### Data Setup
1. Download the MovieLens 100K dataset from [GroupLens](https://grouplens.org/datasets/movielens/100k/).
2. Ensure the `ml-100k` folder is placed in your working directory.

## Known Limitations
* **Categorical Heuristic:** The current genre-matching strategy uses a "first-seen" rule from the dataset[cite: 1]. While this serves as a computational proxy for IMDb’s foremost genre, it may mask the nuances of multi-label hybrid films[cite: 1].

## References
* [1] Wagle, R. (2026). *PCA by Efficient Power Method Implementation on Sparse Movie Ratings Data*. Final Project, Ashoka University[cite: 1].
* [2] Harper, F. M., & Konstan, J. A. (2015). *The MovieLens Datasets: History and Context*. ACM Transactions on Interactive Intelligent Systems.

***
**Author:** Rohan Wagle  
**Date:** April 2026  
**Institution:** Ashoka University
