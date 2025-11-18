#!/bin/Rscript
if(! require(clue)){install.packages("proxy");require(proxy)}
if(! require(proxy)){install.packages("proxy");require(proxy)}
if(! require(stringr)){install.packages("stringr");require(stringr)}
if(! require(mclust)){install.packages("mclust");require(mclust)}
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
        dist_mat[c1,c2]<-dist_mat[c2,c1]<-sqrt(sum((Q_ref[,c1] - Q_curr[,c2])^2))
      }
    }
    
    # Encontrar el mejor emparejamiento usando el algoritmo húngaro
    # install.packages("clue") si no está
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
  numI<-nrow(Qlist[[1]])
  numR=length(Qlist)
  dist_mat <- matrix(0,numR,numR)
  for(r1 in c(1:numR)){
    Q1=Qlist[[r1]]
    
    for(r2 in c(1:numR)){
      Q2=Qlist[[r2]]
      tmpDist <- rep(0,numI)
      for(c1 in c(1:numK)){
          tmpDist<-tmpDist+(Q1[,c1] - Q2[,c1])^2
      }
      
      dist_mat[r1,r2]<-dist_mat[r2,r1]<-mean(sqrt(tmpDist))
  
    }
  }
  colnames(dist_mat)<-row.names(dist_mat)<-paste0("Run",c(1:numR))
  return(dist_mat)
}


max_dist_among_runs<-function(Qlist){
  numK=ncol(Qlist[[1]])
  numI<-nrow(Qlist[[1]])
  numR=length(Qlist)
  dist_mat <- matrix(0,numR,numR)
  for(r1 in c(1:numR)){
    Q1=Qlist[[r1]]
    
    for(r2 in c(1:numR)){
      Q2=Qlist[[r2]]
      tmpDist <- -Inf
      for(c1 in c(1:numK)){
        tmpDist<-max(c(tmpDist,abs(Q1[,c1] - Q2[,c1])))
      }
      
      dist_mat[r1,r2]<-dist_mat[r2,r1]<-tmpDist
      
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


library(mclust)


# Function for automatic thresholds (for mean and max dist)
auto_thresholds <- function(D) {
  dm <- as.vector(D[lower.tri(D)])
  
  # Fit 2-component Gaussian mixture
  mc <- Mclust(dm, G = 2,modelNames="E")
  
  mu  <- mc$parameters$mean
  sdv <- sqrt(mc$parameters$variance$sigmasq)
  # identify low and high components
  low  <- which.min(mu)
  high <- setdiff(1:2, low)
  
  # Threshold  (mean between high and low)
  thr <- ((mu[low]+2*sdv) + (mu[high]-2*sdv)) / 2
  
  
  
  return(list(thr,
       mu = mu,
       sd = sdv))
}

#Function of Mode inference from distance thresholds
find_modes <- function(Dmean,Dmax, thr_mean, thr_max) {
  if(!is.matrix(Dmean)!="matrix"){Dmean=as.matrix(Dmean)}
  if(!is.matrix(Dmax)!="matrix"){Dmax=as.matrix(Dmax)}
  R <- nrow(Dmean)
  remaining <- c(1:R)
  modes <- list()
  # Precompute distances
  
  avg_dist <- sapply(c(1:R), function(i) mean(Dmean[i, -i]))

  while (length(remaining) > 0) {
    
    # start with run that is most "central"
    seed <- remaining[ which.min(avg_dist[remaining]) ]
    group <- seed
    
    added <- TRUE
    while (added) {
      added <- FALSE
      candidates <- setdiff(remaining, group)
      
      for (c in candidates) {
        
        Meandists <- Dmean[c, group]
        Maxdists <- Dmax[c, group]
        if (mean(Meandists) <= thr_mean && max(Maxdists) <= thr_max) {
          group <- c(group, c)
          added <- TRUE
        }
      }
    }
    
    modes[[length(modes) + 1]] <- group
    remaining <- setdiff(remaining, group)
  }
  
  return(list(modes=modes,nmodes=length(modes),thr_mean=thr_mean,thr_max=thr_max))
}



##################
# Definición de la función de detección de modo para un K
##################

detectModes <- function(k,proj,th_mean="auto",th_max=0.1) {
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
  meanD <- dist_among_runs(aligned)
  maxD <- max_dist_among_runs(aligned)
  if(th_mean=="auto"){
    th_mean <- auto_thresholds(meanD)[[1]]
  }

  groupsF <- find_modes(meanD,maxD,th_mean,th_max)
  nmodes<-groupsF[[2]]
  Listmodes<-groupsF[[1]]

  message("Identified modes for K=", k, ": ", nmodes)
  
  # Keep results
  listOutput<-list()
  for(mode in c(1:length(Listmodes))){
  
    modes<-as.numeric(str_remove(Listmodes[[mode]],"Run"))
    selected <- aligned[modes]
    mismatch=compute_mismatch(selected)
    representativeRun=get_representative_run(selected,Listmodes[[mode]])
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
