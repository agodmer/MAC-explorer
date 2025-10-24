# Build k-mer feature matrices from simple input.
# If no FASTA/TSV is found, create a small synthetic dataset with >=2 samples per class.

suppressPackageStartupMessages({
  library(data.table)
  # library(Biostrings) # Uncomment if you want to parse FASTA with Biostrings
})

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

k_set <- c(1,2,3,4,6,8)

# --- Helper: naive k-mer counting for a single string ---
count_kmer <- function(s, k){
  s <- toupper(gsub("[^ACGT]", "N", s))
  n <- nchar(s)
  if(n < k) return(integer(0))
  kmers <- substring(s, 1:(n - k + 1), k:n)
  kmers <- kmers[!grepl("N", kmers, fixed = TRUE)]
  if(length(kmers) == 0) return(integer(0))
  as.integer(table(factor(kmers, levels = unique(kmers))))
  # Returned as counts aligned to unique(kmers) is not needed for matrix build below.
}

# --- Build k-mer matrix for a data.frame with columns: id, seq, class ---
build_matrix <- function(dt, k){
  # Collect per-sample k-mer named integer vectors
  vecs <- lapply(dt$seq, function(s){
    z <- substring(s, 1:(nchar(s)-k+1), k:nchar(s))
    z <- toupper(z)
    z <- z[!grepl("N", z, fixed = TRUE)]
    table(z)
  })
  # Determine the full vocabulary
  vocab <- sort(unique(unlist(lapply(vecs, names))))
  if(length(vocab) == 0){
    mat <- matrix(0, nrow = nrow(dt), ncol = 0)
    return(list(X = mat, y = dt$class))
  }
  # Fill matrix
  M <- matrix(0, nrow = nrow(dt), ncol = length(vocab), dimnames = list(dt$id, vocab))
  for(i in seq_along(vecs)){
    v <- vecs[[i]]
    if(length(v)) M[i, names(v)] <- as.numeric(v)
  }
  # Row-normalize to relative frequencies
  rs <- rowSums(M)
  rs[rs == 0] <- 1
  M <- M / rs
  list(X = M, y = dt$class)
}

# --- Load input or create a synthetic demo dataset with >=2 per class ---
tsv <- file.path("data", "sim_sequences.tsv")
if(!file.exists(tsv)){
  # Synthetic minimal dataset ensuring at least 2 samples per class
  dna <- c(
    "ATGCATGCAT", "ATGCGGATAT", "TTTTGGCCAATG", "GCCATGTTAT",
    "ATGATGCCAT", "CCGGTTAAAT", "ATGCCCATTT", "GGGATTAACC"
  )
  classes <- c("A","A","A","A","B","B","B","B")
  fwrite(data.table(id = paste0("sim", seq_along(dna)), seq = dna, class = classes), tsv)
  message("Wrote synthetic data to data/sim_sequences.tsv")
}
D <- fread(tsv)  # expects columns: id, seq, class
stopifnot(all(c("id","seq","class") %in% names(D)))

# --- Save feature matrices for all k values ---
for(k in k_set){
  fx <- build_matrix(D, k)
  saveRDS(fx, file = sprintf("results/features_k%d.rds", k))
}
message("Features saved to results/features_k{1,2,3,4,6,8}.rds")
