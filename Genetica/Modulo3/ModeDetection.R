#!/bin/Rscript
library(clue)

# Función para alinear las columnas de diferentes iteraciones (runs) de sNMF
# para un mismo valor de K usando correlación máxima entre componentes
# Función para alinear las columnas de diferentes iteraciones (runs) de sNMF
# para un mismo valor de K usando correlación máxima entre componentes
align_Q_runs <- function(Qlist) {
  # Qlist: lista de matrices Q, cada elemento es un run para el mismo K
  # Devuelve una lista de matrices Q con columnas reordenadas
  
  if(length(Qlist) < 2) return(Qlist)

  Q_ref <- Qlist[[1]]  # Primera corrida como referencia
  Q_aligned <- list(Q_ref)
  numK=ncol(Q_ref)  
  for(i in 2:length(Qlist)) {
    Q_curr <- Qlist[[i]]
    
    # Matriz de correlación entre columnas
    dist_mat <- matrix(NA,numK,numK)
    for(c1 in c(1:numK)){
      for(c2 in c(1:numK)){
        dist_mat[c1,c2]<-dist_mat[c2,c1]<-mean(sqrt((Q_ref[,c1] - Q_ref[,c2])^2))
      }
    }
    # Encontrar el mejor emparejamiento usando el algoritmo húngaro
    # install.packages("clue") si no está
    library(clue)
    match_idx <- solve_LSAP(dist_mat, maximum = FALSE)
    
    # Reordenar columnas de la corrida actual
    Q_curr_aligned <- Q_curr[, match_idx]
    
    Q_aligned[[i]] <- Q_curr_aligned
  }
  
  return(Q_aligned)
}


###euclidian distance among different Q 
dist_among_runs<-function(Qlist){
  numK=ncol(Qlist[[1]])  
  numR=length(Qlist)
  dist_mat <- matrix(NA,numR,numR)
  for(r1 in c(1:numR)){
    Q1=Qlist[[r1]]
    for(r2 in c(1:numR)){
      Q2=Qlist[[r2]]
      tmpDist <- 0
      for(c1 in c(1:numK)){
          tmpDist<-tmpDist+(Q1[,c1] - Q2[,c1])^2
      }
      tmpDist<-sqrt(tmpDist)
      dist_mat[r1,r2]<-dist_mat[r2,r1]<-mean(tmpDist)
  
    }
  }
  return(dist_mat)
}



# Function to get the best representive runs (the most similar to all runs listed)
get_representative_run <- function(alignedList,listnames) {
  dist_pairruns <- dist_among_runs(alignedList)
  avg_dist_per_run<-c()
  for(i in c(1:nrow(dist_pairruns))){
    avg_dist_per_run[i]<-mean(dist_pairruns[i,-i])
  }
  rep_run <- listnames[which.min(avg_dist_per_run)]
  return(rep_run)
}


# Function to get mismatch score across runs (only euclidien distance)
compute_mismatch <- function(alignedList) {
  nRuns <- length(alignedList)
  if(nRuns == 1) return(0)  # only one run, no mismatch
  # compute all pairwise mean absolute differences
  pairwise_dist <- dist_among_runs(alignedList)
  # mismatch index: average across all pairwise distances
  mismatch <- sum(pairwise_dist)/(nRuns^2-3)
  return(mismatch)
}

compute_Qmean <- function(alignedList) {
  nRuns <- length(alignedList)
  if(nRuns == 1) return(alignedList[[1]])
  Qout<-alignedList[[1]]
  for(i in c(2:nRuns)){
    Qout<-Qout+alignedList[[i]]
  }
  return(Qout/nRuns)
}


##################
# Definición de la función de detección de modo para un K
##################

process_K <- function(k,proj,dirRes) {
  message("Processing K = ", k)
  
  ce <- cross.entropy(proj,K=k)
  nruns<-nrow(ce)
  message("Number of runs detected: ", nruns)
  
  # ----------------------------------------------------
  # Extract Q matrices from all runs
  # ----------------------------------------------------
  Qs <- lapply(1:nruns, function(r) Q(proj, K = k, run = r))
  
  # ----------------------------------------------------
  # Align Q matrices (fix label switching)
  # ----------------------------------------------------
  aligned <- align_Q_runs(Qs)
  
  # Convert the aligned list into a single matrix
  # ----------------------------------------------------
  # Compute distances among runs
  # ----------------------------------------------------
  D <- computeDist(aligned)
  
  # ----------------------------------------------------
  # Cluster runs into modes (hierarchical clustering)
  # ----------------------------------------------------
  hc <- hclust(D, method = "average")
  
  # Automatic mode detection:
  # rule: cut at 10% of maximum height (adjust if needed)
  threshold <- 0.10 * max(hc$height)
  modes <- cutree(hc, h = threshold)
  
  nmodes <- length(unique(modes))
  message("Identified modes for K=", k, ": ", nmodes)
  
  # ----------------------------------------------------
  # Compute average Q matrix for each mode
  # ----------------------------------------------------
  
  listOutput<-list()
  
  for(mode in unique(sort(modes))){
    
    selected <- aligned[modes == mode]
    listRuns <- which(modes==mode)
    
    
    listOutput[[paste("mode "),mode]]<-list(
      listRuns=listRuns,
      representativeRun=get_representative_run(selected,paste0("Run",listRuns)),
      mismatch=compute_mismatch(selected),
      meanQ=compute_Qmean(selected)
    )
    outdir <- paste0(dirRes,"/K", k,"Mode",mode)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
    # Save mode runs
    write.table(
      data.frame(Run = 1:nruns, Mode = modes),
      file = file.path(outdir, "run_modes.txt"),
      row.names = FALSE, quote = FALSE
    )
  
  # Save mode-average Q matrices
  for (i in seq_along(mode_Q)) {
    write.table(
      mode_Q[[i]],
      file = file.path(outdir, paste0("mode_", i, "_Qmatrix.txt")),
      row.names = FALSE, col.names = FALSE, quote = FALSE
    )
  }
  
  # ----------------------------------------------------
  # Plot modes (simple barplot for each mode)
  # ----------------------------------------------------
  for (i in seq_along(mode_Q)) {
    
    df <- data.frame(
      Individual = 1:nrow(mode_Q[[i]]),
      mode_Q[[i]]
    )
    df_long <- reshape2::melt(df, id.vars = "Individual")
    
    p <- ggplot(df_long, aes(x = Individual, y = value, fill = variable)) +
      geom_bar(stat = "identity", width = 1) +
      theme_minimal() +
      ggtitle(paste("K =", K, " | Mode", i)) +
      ylab("Ancestry proportion")
    
    ggsave(file.path(outdir, paste0("mode_", i, "_plot.png")),
           p, width = 10, height = 4)
  }
  
  # Return summary
  list(K = K, nruns = nruns, nmodes = nmodes)
}
