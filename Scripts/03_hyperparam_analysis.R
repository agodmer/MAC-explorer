# Build supplementary figures S2–S5 from aggregated results.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dir.create("results", showWarnings = FALSE)

agg_file <- "example_results_hyperparams.csv"
if(!file.exists(agg_file)) stop("Missing results/hyperparams_agg.csv. Run 02_model_training.R first.")
A <- fread(agg_file)
setnames(A, tolower(names(A)))
if("k" %in% names(A)) A[, k := as.factor(k)]

pdf("results/Supplementary_Figures_S2_S5.pdf", width = 7, height = 5)

# S2. Distribution of Cohen's kappa by k
if(all(c("k","kappa") %in% names(A))){
  p2 <- ggplot(A[!is.na(kappa)], aes(x = k, y = kappa)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.2) +
    labs(title = "S2. Distribution of Cohen's kappa by k", x = "k-mer length (k)", y = "Kappa")
  print(p2)
}

# S3. Effect of num.trees (if present; with tuneLength we fixed num.trees, so this may be flat)
if(all(c("num.trees","kappa") %in% names(A))){
  p3 <- ggplot(A[!is.na(num.trees) & !is.na(kappa)],
               aes(x = num.trees, y = kappa, group = k, linetype = k, shape = k)) +
    stat_summary(fun = mean, geom = "line") +
    stat_summary(fun = mean, geom = "point") +
    labs(title = "S3. Effect of num.trees", x = "num.trees", y = "Mean kappa") +
    guides(linetype = guide_legend(title = "k"), shape = guide_legend(title = "k"))
  print(p3)
}

# S4. Heatmap mtry x min.node.size
if(all(c("mtry","min.node.size","kappa") %in% names(A))){
  H <- A[!is.na(mtry) & !is.na(min.node.size) & !is.na(kappa),
         .(kappa = mean(kappa, na.rm = TRUE)), by = .(k, mtry, min.node.size)]
  p4 <- ggplot(H, aes(x = mtry, y = min.node.size, fill = kappa)) +
    geom_tile() +
    facet_wrap(~ k, scales = "free") +
    labs(title = "S4. Kappa heatmap by mtry and min.node.size",
         x = "mtry", y = "min.node.size", fill = "kappa")
  print(p4)
}

# S5. Top-20 features by permutation importance (if exported)
imp_files <- list.files("results", pattern = "^importance_k[0-9]+_split[0-9]+\\.csv$", full.names = TRUE)
if(length(imp_files) > 0){
  IMP <- rbindlist(lapply(imp_files, fread), fill = TRUE)
  setnames(IMP, tolower(names(IMP)))
  if(all(c("feature","importance","k") %in% names(IMP))){
    TOP <- IMP[, .SD[order(-importance)][1:20], by = k]
    p5 <- ggplot(TOP, aes(x = reorder(feature, importance), y = importance)) +
      geom_col() +
      coord_flip() +
      facet_wrap(~ k, scales = "free_y") +
      labs(title = "S5. Top-20 informative k-mers",
           x = "k-mer", y = "Permutation importance")
    print(p5)
  }
}

dev.off()
message("Wrote results/Supplementary_Figures_S2_S5.pdf")
