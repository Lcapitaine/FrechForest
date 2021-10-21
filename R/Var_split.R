#' Classical Variable Split function
#'
#' @param X
#' @param Y
#' @param timeScale
#'
#' @import kmlShape
#' @import Evomorph
#'
#' @keywords internal
var_split <- function(X ,Y,timeScale=0.1){
  # Pour le moment on se concentre sur le cas des variables courbes ::
  impur <- rep(0,dim(X$X)[length(dim(X$X))])
  toutes_imp <- list()
  split <- list()
  centers <- list() # On va stocker les centres associes aux kmeans
  Pure <- FALSE

  for (i in 1:dim(X$X)[length(dim(X$X))]){

    if (X$type=="factor"){
      if (length(unique(X$X[,i]))>1){
        L <- Fact.partitions(X$X[,i],X$id)
        split_courant <- list()
        impur_courant <- rep(NA,length(L))
        toutes_imp_courant <- list()
        # Il faut maintenant regarder quelles sont les meilleures combinaisons ::
        for (k in 1:length(L)){
          split_courant[[k]] <- rep(2,length(X$id))
          for (l in L[[k]]){
            split_courant[[k]][which(X$id==l)] <- 1
          }
          # Il faut maintenant regarder la qualite du decoupage ::
          impurete <- impurity_split(Y,split_courant[[k]])
          impur_courant[k] <- impurete$impur
          toutes_imp_courant[[k]] <- impurete$imp_list
        }
        select <- which.min(impur_courant)
        split[[i]] <- split_courant[[select]]
        impur[i] <- impur_courant[select]
        toutes_imp[[i]] <- toutes_imp_courant[[select]]
      }
      else {
        impur[i] <- Inf
        split[[i]] <- Inf
      }
    }

    if( X$type=="curve"){
      mclds <- kmlShape::cldsWide(ordonne(X$X[,i], X$time, X$id), unique(X$time), unique(X$id))
      crit <- kmlShape::kmlShape(mclds, nbClusters = 2, timeScale = timeScale, toPlot="none")
      att <- attributes(crit)
      split[[i]] <- att$clusters
      impurete <- impurity_split(Y,split[[i]], timeScale)
      impur[i] <- impurete$impur
      toutes_imp[[i]] <- impurete$imp_list
    }

    if(X$type=="scalar"){
      if (length(unique(X$X[,i]))>2){
        sp <- kmeans(X$X[,i], centers=2)
        split[[i]] <- sp$cluster
        impurete <- impurity_split(Y,split[[i]], timeScale)
        impur[i] <- impurete$impur
        toutes_imp[[i]] <- impurete$imp_list
      }

      if (length(unique(X$X[,i]))==2){
        split[[i]] <- rep(2,length(X$X[,i]))
        split[[i]][which(X$X[,i]==unique(X$X[,i])[1])] <- 1
        impurete <- impurity_split(Y,split[[i]], timeScale)
        impur[i] <- impurete$impur
        toutes_imp[[i]] <- impurete$imp_list
      }

      if (length(unique(X$X[,i]))==1) {
        impur[i] <- Inf
        split[[i]] <- Inf
      }
    }
  }

  if (length(unique(impur))==1 & is.element(Inf,impur)==TRUE){
    return(list(Pure=TRUE))
  }
  true_split <- which.min(impur)
  split <- split[[true_split]]
  return(list(split=split, impurete=min(impur),impur_list = toutes_imp[[true_split]], variable=which.min(impur), Pure=Pure))
}
