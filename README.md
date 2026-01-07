# MAC-explorer

Reproducible pipeline for genome classification based on **k-mer frequency features** and **machine learning**.  
It extracts sequence features, trains Random Forest models, and analyzes hyperparameter behavior.

---

## 📂 Repository Structure
```r
Scripts/
├── 01_feature_extraction.R # Compute k-mer frequency matrices from DNA sequences
├── 02_model_training.R # Train Random Forest classifiers using caret::train
├── 03_hyperparam_analysis.R # Aggregate and visualize hyperparameters
README.md
```
---

## ⚙️ Requirements

**Language:** R (≥ 4.2)

Analyses were run in R v4.3.2. Visualization and data handling used dplyr v1.1.4 , ggplot2 v3.4.4 , and rstatix v0.7.2.

**Core packages:**
```r
install.packages(c(
  "Biostrings", "caret", "ranger", "data.table",
  "ggplot2", "progress", "dplyr"
))
Optional: reshape2, gridExtra
```
Pipeline Overview

Step 1 — Feature Extraction
```r
#Script: Scripts/01_feature_extraction.R
#Reads .fna / .fa genome files or simulates sequences.
#Computes k-mer frequency matrices for k ∈ {1, 2, 3, 4, 6, 8}.
#Output: kmers_features_by_batch{k}.RData
```
Step 2 — Model Training
```r
#Script: Scripts/02_model_training.R
#Loads k-mer matrices and taxonomy labels.
#Runs 50 independent 70/30 stratified train/test splits.
#Performs 3-fold CV with random hyperparameter search.
#Saves models and results as models_results*.RData.
```
Step 3 — Hyperparameter Analysis
```r
#Script: Scripts/03_hyperparam_analysis.R
#Scans all .RData files containing caret::train objects.
#Extracts grids (mtry, min.node.size, splitrule, num.trees) and metrics (Accuracy, Kappa).
```
Outputs:
```r
#hyperparams_agg.csv — aggregated hyperparameter-performance table
#Supplementary_Figures_S2_S4.pdf — figures for publication (S2–S4)
```

Output Summary
File and	Description
```r
#hyperparams_agg.csv	-> Combined hyperparameter grid and metrics
#Supplementary_Figures_S2_S4.pdf -> 	Figures: Kappa distributions, num.trees effect, heatmap
#models_results*.RData	->  Trained models and confusion matrices
#kmers_features_by_batch{k}.RData	-> Computed k-mer frequency matrices
```
Reproducibility Example
```r
# Step 1: Feature extraction
Rscript Scripts/01_feature_extraction.R

# Step 2: Model training
Rscript Scripts/02_model_training.R

# Step 3: Hyperparameter aggregation & visualization
Rscript Scripts/03_hyperparam_analysis.R
```
R session information and package references
```r
R version: >= 4.2.0

The following R packages were used for data analysis, statistical testing, clustering,
machine learning, and visualization:

Core statistical analysis and multivariate methods
- stats (base R): PCA, hierarchical clustering (hclust), Wilcoxon tests, Spearman correlations
- FactoMineR (v2.9): principal component analysis
- minpack.lm (v1.2.4): nonlinear least-squares fitting (double-sigmoid models)
- effsize (>= 0.8.1): effect size estimation (Cohen’s d)

Clustering comparison and concordance metrics
- mclust (v6.1.1): Adjusted Rand Index (ARI)
- aricode (v1.0.3): normalized mutual information (NMI)
- clue (v0.3.66): clustering comparison utilities

Phylogenetics and tree handling
- ape (v5.8.1): tree manipulation and export (Newick format)

Distance-based ecological and genomic analyses
- vegan (>= 2.6-4): PERMANOVA (adonis), Mantel tests
- proxy (>= 0.4-27): Jaccard distance computation

Pan-genome and accessory genome analyses
- tidyverse (>= 2.0.0): data manipulation and visualization
  - dplyr
  - tidyr
  - ggplot2
  - readr
  - purrr

Machine learning and classification
- Biostrings (>= 2.66.0): k-mer frequency computation
- ranger (v0.17.0): random forest implementation
- caret (v6.0-94): model training, resampling, and evaluation
- e1071 (>= 1.7-14): auxiliary ML utilities (caret dependency)

Reproducibility and utilities
- here (>= 1.0.1): project-relative file paths
- cowplot (>= 1.1.1): figure assembly
```

# Citation
If you use this repository, please cite as: 
Poignon, C., Matondo, M., Giai Gianetto, Q., Douché, T., Saffarian, A., Veziris, N., Aubry, A. & Godmer, A.
MAC-Explorer: bridging genome-based taxonomy and an identification tool for the Mycobacterium avium complex.
