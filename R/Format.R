#' Read the parameters of the function
#'
#' @param z
#'
#'
#' @keywords internal
read.Xarg <- function(z){
  type <- NULL
  issou <- rep(NA, length(z))
  for (i in 1:length(z)){
    issou[i] <- is.null(z[i])
    if (issou[i]==FALSE) type <- c(type,z[i]$type)
  }
  return(type)
}



#' Factor partitions finder
#'
#' This function is used to find all the unique partitions of k factors into 2 groups
#'
#' @param Factor
#' @param id
#'
#' @keywords internal
Fact.partitions <- function(Factor, id){

  U <- unique(Factor)
  P <- Part.facts[[length(U)]]
  L <- list()
  for (k in 1:nrow(P)){
    w <- which(P[k,]==0)
    U_courant <- U[w]
    W <- NULL
    for (m in U_courant){
      W <- c(W,which(Factor==m))
    }
    L[[k]] <- id[W]
  }
  return(L)
}


#' Ordonne
#'
#' @param X
#' @param time
#' @param id
#'
#'
#' @keywords internal
ordonne <- function(X , time , id){
  mat  <- matrix(NA, length(unique(id)), length(unique(time)))
  for( i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    t_w <- time[w]
    w_time <- NULL
    for (j in 1:length(w)){
      w_time <- c(w_time, which(unique(time)==t_w[j]))
    }
    mat[i,w_time] <- X[w]
  }
  return(mat)
}
