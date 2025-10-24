############################################################
# Example reproducible pipeline for k-mer + Random Forest
# Author: [Your Name]
# Description:
#   - Generate artificial DNA sequences
#   - Compute k-mer frequency matrices (k = {1,2,3,4,6,8})
#   - Train Random Forest classifiers using caret::train
#   - Store results and metrics
############################################################

suppressPackageStartupMessages({
  library(Biostrings)
  library(caret)
  library(data.table)
  library(ranger)
  library(progress)
  library(dplyr)
})

set.seed(42)

###########################
### 1. Generate example data
###########################

# Parameters
n_genomes <- 60         # number of genomes
seq_length <- 5000      # bp per genome
classes <- c("A","B","C")

# Random DNA sequences
generate_random_dna <- function(n, len) {
  alphabet <- c("A","C","G","T")
  replicate(n, paste0(sample(alphabet, len, replace = TRUE), collapse = ""))
}

dna_sequences <- generate_random_dna(n_genomes, seq_length)
labels <- sample(classes, n_genomes, replace = TRUE)
meta <- data.frame(Genome = paste0("Genome_", seq_len(n_genomes)),
                   Class = labels, stringsAsFactors = FALSE)

###########################
### 2. Compute k-mer features
###########################

k_set <- c(1,2)
df_kmers_list <- list()
pb <- progress_bar$new(total = length(k_set), format = "⏳ [:bar] k=:current ETA=:eta")

extract_kmers <- function(seq, k) {
  oligonucleotideFrequency(DNAString(seq), width = k)
}

for (k in k_set) {
  pb$tick()
  km <- t(sapply(dna_sequences, extract_kmers, k = k))
  rownames(km) <- meta$Genome
  df_kmers_list[[as.character(k)]] <- as.data.frame(km)
}
cat("\n✅ K-mer matrices generated.\n")

###########################
### 3. Machine learning classification
###########################

n_splits <- 30
results_all <- list()
pb2 <- txtProgressBar(min = 0, max = n_splits * length(k_set), style = 3)

for (k in k_set) {
  X <- df_kmers_list[[as.character(k)]]
  y <- make.names(meta$Class)
  
  # Upsample minority classes to balance
  df_bal <- upSample(x = X, y = as.factor(y))
  X <- as.matrix(df_bal[, setdiff(names(df_bal), "Class")])
  y <- df_bal$Class
  
  for (s in seq_len(n_splits)) {
    set.seed(1000 + s + k)
    idx <- createDataPartition(y, p = 0.7, list = FALSE)
    x_train <- X[idx, ]
    y_train <- y[idx]
    x_test  <- X[-idx, ]
    y_test  <- y[-idx]
    
    ctrl <- trainControl(
      method = "cv",
      number = 3,
      classProbs = TRUE,
      savePredictions = TRUE
    )
    
    ntree_val <- sample(seq(250, 1000, by = 50), 1)
    model <- train(
      x = x_train, y = y_train,
      method = "ranger",
      metric = "Kappa",
      trControl = ctrl,
      preProcess = c("center","scale"),
      tuneLength = 5,
      importance = "impurity",
      num.trees = ntree_val
    )
    
    pred <- predict(model, newdata = x_test)
    cm <- confusionMatrix(pred, y_test)
    res <- data.table(
      k = k,
      split = s,
      num.trees = ntree_val,
      Accuracy = cm$overall["Accuracy"],
      Kappa = cm$overall["Kappa"]
    )
    
    results_all[[paste(k, s, sep = "_")]] <- res
    setTxtProgressBar(pb2, (which(k_set == k) - 1) * n_splits + s)
  }
}

close(pb2)
results_df <- rbindlist(results_all)
fwrite(results_df, "example_results_hyperparams.csv")
cat("\n✅ Training complete, results saved to example_results_hyperparams.csv\n")

###########################
### 4. Visualization (optional)
###########################
library(ggplot2)

pdf("example_hyperparam_figures.pdf", width = 7, height = 5)

# Distribution of Kappa per k
ggplot(results_df, aes(x = factor(k), y = Kappa)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2) +
  labs(title = "Distribution of Cohen's Kappa across k", x = "k-mer length", y = "Kappa") -> p1
print(p1)

# Effect of num.trees
ggplot(results_df, aes(x = num.trees, y = Kappa, color = factor(k))) +
  geom_point(alpha = 0.6) +
  geom_smooth(se = FALSE) +
  labs(title = "Effect of num.trees on model performance", x = "num.trees", y = "Kappa") -> p2
print(p2)

dev.off()
cat("✅ Figures saved to example_hyperparam_figures.pdf\n")

###########################
### 5. Reproducibility info
###########################
sink("sessionInfo_example.txt")
cat("## System and session info ##\n")
print(Sys.info())
print(sessionInfo())
sink()
cat("✅ Session info saved.\n")
