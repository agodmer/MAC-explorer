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

**Core packages:**
```r
install.packages(c(
  "Biostrings", "caret", "ranger", "data.table",
  "ggplot2", "progress", "dplyr"
))
Optional: reshape2, gridExtra

🚀 Pipeline Overview
Step 1 — Feature Extraction
Script: Scripts/01_feature_extraction.R

Reads .fna / .fa genome files or simulates sequences.

Computes k-mer frequency matrices for k ∈ {1, 2, 3, 4, 6, 8}.

Output: kmers_features_by_batch{k}.RData

Step 2 — Model Training
Script: Scripts/02_model_training.R

Loads k-mer matrices and taxonomy labels.

Runs 50 independent 70/30 stratified train/test splits.

Performs 3-fold CV with random hyperparameter search.

Saves models and results as models_results*.RData.

Step 3 — Hyperparameter Analysis
Script: Scripts/03_hyperparam_analysis.R

Scans all .RData files containing caret::train objects.

Extracts grids (mtry, min.node.size, splitrule, num.trees) and metrics (Accuracy, Kappa).

Outputs:

hyperparams_agg.csv — aggregated hyperparameter-performance table

Supplementary_Figures_S2_S4.pdf — figures for publication (S2–S4)

📊 Output Summary
File	Description
hyperparams_agg.csv	Combined hyperparameter grid and metrics
Supplementary_Figures_S2_S4.pdf	Figures: Kappa distributions, num.trees effect, heatmap
models_results*.RData	Trained models and confusion matrices
kmers_features_by_batch{k}.RData	Computed k-mer frequency matrices

🧪 Reproducibility Example
bash
Copier le code
# Step 1: Feature extraction
Rscript Scripts/01_feature_extraction.R

# Step 2: Model training
Rscript Scripts/02_model_training.R

# Step 3: Hyperparameter aggregation & visualization
Rscript Scripts/03_hyperparam_analysis.R
🧬 Citation
If you use this repository, please cite as:
