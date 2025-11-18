#!/bin/Rscript
if(!require(clue)){install.packages("clue");require(clue)}

# # Definition of the function to align cluster
align_Q_down <- function(QK, QKminus1, fix_thr = 0.99) {
  # QK:        K × N (larger)
  # QKminus1: (K−1) × N (smaller)
  
  K  <- nrow(QK)
  Km <- nrow(QKminus1)
  
  # Identify near-fixed individuals
  fixed_K  <- QK > fix_thr
  fixed_Km <- QKminus1 > fix_thr
  
  # Overlap matrix K × (K−1)
  overlap <- matrix(0, K, Km)
  for (i in 1:K) {
    for (j in 1:Km) {
      #overlap[i, j] <- sum(fixed_K[i, ] & fixed_Km[j, ])
      overlap[i, j] <- mean(sqrt((QK[i, ] - QKminus1[j, ])^2))
    }
  }
  
  # LSAP requires rows <= columns → transpose (Km × K)
  overlap_t <- t(overlap)
  
  # Hungarian assignment
  #assignment <- solve_LSAP(overlap_t, maximum = TRUE)
  assignment <- solve_LSAP(overlap_t, maximum = FALSE)
  
  # assignment[j] = row index in QK best matching row j in QKminus1
  matched_rows <- as.integer(assignment)   # ensure numeric
  
  # Reorder Q(K−1) to match cluster order at K
  QKminus1_reordered <- QKminus1[order(matched_rows), , drop = FALSE]
  
  return(list(
    Q = QKminus1_reordered,     # aligned Q(K−1)
    matched = matched_rows      # numeric row indices of Q(K)
  ))
}
