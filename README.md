# Sparse-PCA-Power-Method: Dimensionality Reduction for Sparse Movie Ratings Data

## Overview
This repository contains a manual implementation of the **Simultaneous Power Method** used to perform Principal Component Analysis (PCA) on the **MovieLens 100K** dataset[cite: 1]. The project focuses on recovering the latent structures of movie preferences while navigating the technical constraints of a **93.7% sparse ratings matrix**[cite: 1].

## Key Technical Achievements
* **Implicit Centering Algorithm:** Developed a memory-efficient "Forward-Backward" pass logic that performs PCA on centered data without explicitly densifying the sparse matrix, maintaining $O(\text{nnz})$ computational complexity[cite: 1].
* **Subspace Stability:** Manually utilized the **Gram-Schmidt** process within the power iteration loop to prevent subspace collapse, successfully recovering the first 10 principal components[cite: 1].
* **Convergence Rigor:** Demonstrated log-linear convergence of eigenvectors over ~600 iterations, providing a mathematically sound basis for the resulting latent space[cite: 1].

## Analytical Insights
* **Dimensionality Reduction:** Condensed 1,682 movie dimensions into a low-rank space that explains over **35% of total data variation**[cite: 1].
* **Demographic Correlation:** Identified a significant **0.369 correlation** between user age and the second principal component (PC2), which acts as a "Classic vs. Blockbuster" preference pivot[cite: 1].
* **Interactive Exploration:** Includes a Plotly-based interactive loadings plot that allows users to audit the "Full Genre Profile" and release year of movies, validating the model's mathematical accuracy against qualitative metadata[cite: 1].

## Repository Structure
* `/code`: R implementation featuring the manual Power Method loop and implicit centering logic[cite: 1].
* `/report`: Final LaTeX report detailing the mathematical proofs for convergence and demographic findings[cite: 1].
* `/output`: Static and interactive visualizations, including Scree plots, log-error convergence charts, and the interactive Subspace Explorer[cite: 1].

## Future Work
* **Categorical Refinement:** Integration of the IMDb API to replace the current "first-seen" genre heuristic with primary editorial labels for sharper cluster interpretation[cite: 1].
* **Statistical Validation:** Implementation of **Parallel Analysis** to further validate the selection of $k$ components against random noise.

---
**Author:** Rohan Wagle  
**Date:** April 2026  
**Institution:** Ashoka University

## References
[1] Wagle, R. (2026). PCA by Efficient Power Method Implementation on Sparse Movie Ratings Data. Final Project Presentation, Ashoka University.
