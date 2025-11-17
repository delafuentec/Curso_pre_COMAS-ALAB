#!/bin/Rscript
if(! require(clue)){install.packages("proxy");require(proxy)}
if(! require(proxy)){install.packages("proxy");require(proxy)}
if(! require(stringr)){install.packages("stringr");require(stringr)}
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
        dist_mat[c1,c2]<-dist_mat[c2,c1]<-sqrt(mean((Q_ref[,c1] - Q_ref[,c2])^2))
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
          tmpDist<-tmpDist+mean((Q1[,c1] - Q2[,c1])^2)
      }
      tmpDist<-sqrt(tmpDist)
      dist_mat[r1,r2]<-dist_mat[r2,r1]<-mean(tmpDist)
  
    }
  }
  colnames(dist_mat)<-row.names(dist_mat)<-paste0("Run",c(1:numR))
  return(dist_mat)
}



# Function to get the best representive runs (the most similar to all runs listed)
get_representative_run <- function(alignedList,listnames) {
  dist_pairruns <- dist_among_runs(alignedList)
  avg_dist_per_run<-c()
  for(i in c(1:nrow(dist_pairruns))){
    avg_dist_per_run[i]<-mean(dist_pairruns[i,-i,drop=FALSE])
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
  mismatch <- mean(dist(pairwise_dist))
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


##definition de grandes en un hclust
groupsHC <- function(hc, alpha = 0.05) {
  # root merge height

  Hroot <- max(hc$height)-min(hc$height)
  
  # splits to keep must exceed threshold
  
  threshold <- alpha * Hroot
  # find all internal nodes whose height ≥ threshold
  good_nodes <- which(hc$height >= threshold)
  
  # No significant splits → only 1 mode
  if (length(good_nodes) == 0)
    return(list(modes = list(unlist(as.list(hc$labels))), n = 1))
  
  # Each such node defines a subtree → mode
  extract_tips <- function(node, merge, labels) {
    if (node < 0)
      return(labels[-node])
    children <- merge[node, ]
    c(extract_tips(children[1], merge, labels),
      extract_tips(children[2], merge, labels))
  }
  
  modes <- c()
  used_tips <- character()
  
  for (node in good_nodes) {
    tips <- extract_tips(node, hc$merge, hc$labels)
    # Avoid duplicates from nested nodes
    if (!any(tips %in% used_tips)) {
      modes[[length(modes)+1]] <- tips
      used_tips <- c(used_tips, tips)
    }#else{
    #  modes[[length(modes)+1]] <- tips[ ! tips %in% used_tips]
    #  used_tips <- unique(c(used_tips, tips))
    #}
  }
  
  
  return(list(modes = modes, n = length(modes)))
}




##################
# Definición de la función de detección de modo para un K
##################

detectModes <- function(k,proj,alpha=0.25) {
  ##get number of runs
  ce <- cross.entropy(proj,K=k)
  nruns<-nrow(ce)

  ## Extract Q matrices from all runs
  Qs<-list()
  for(run in c(1:nruns)){
    Qs[[run]]<-Q(proj,K=k,run=run)
  }
  ## Align Q matrices (match components across runs)
  aligned <- align_Q_runs(Qs)
  
  ## Compute distances among runs
  D <- dist_among_runs(aligned)
  
  #print(D)
  # Cluster runs into modes (hierarchical clustering)
  hc <- hclust(dist(D), method = "average")
  plot(hc,main=paste0("hclust from euclidian distance among runs\nK = ",k),hang=-1,ylim=c(4,0))
  # Automatic mode detection:
  # rule: cut at 10% of maximum height (can be adjusted if needed)
  groupsF <- groupsHC(hc,alpha)
  nmodes<-groupsF[[2]]
  Listmodes<-groupsF[[1]]
  
  
  message("Identified modes for K=", k, ": ", nmodes)
  
  # Keep results
  listOutput<-list()
  for(mode in c(1:length(Listmodes))){
  
    modes<-as.numeric(str_remove(Listmodes[[mode]],"Run"))
    selected <- aligned[modes]
    mismatch=compute_mismatch(selected)
    representativeRun=get_representative_run(selected,paste0("Run",listRuns))
    meanQ=compute_Qmean(selected)
    
    listOutput[[paste0("mode",mode)]]<-list(
      listRuns=modes,
      representativeRun=representativeRun,
      mismatch=mismatch,
      meanQ=meanQ
      )
    message(paste("K=", k, ", mode ",mode,
                  ". listRuns:",paste(listOutput[[paste0("mode",mode)]][1],collapse=","),
                  ", Most representative run: ",listOutput[[paste0("mode",mode)]][2],
                  ", Mismatch: ",listOutput[[paste0("mode",mode)]][3],sep=""))
                                                     
                  
            
  }  
  # Return summary
  return(listOutput)
}
