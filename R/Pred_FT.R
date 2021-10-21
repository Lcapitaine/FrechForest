#' Predict Frechet tree
#'
#' @param tree : Frechet tree.
#' @param Curve [list]: A list that contains the input curves.
#' @param Scalar [list]: A list that contains the input scalars.
#' @param Factor [list]: A list that contains the input factors.
#' @param Shape [list]: A list that contains the input shape.
#' @param Image [list]: A list that contains the input images.
#' @param aligned.shape [logical]: \code{TRUE} if the input shapes are aligned and normalized (\code{aligned.shape=FALSE} by default)
#' @param timeScale [numeric]: Time scale for the input and output curves (\code{timeScale=0.1} by default)
#'
#' @import stringr
#' @import geomorph
#' @import kmlShape
#' @import Evomorph
#'
#' @export
#'
pred.FT <- function(tree, Curve=NULL,Scalar=NULL,Factor=NULL,Shape=NULL,Image=NULL, aligned.shape=FALSE ,timeScale=0.1){

  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs

  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  id.pred <- unique(get(Inputs[1])$id)

  if (is.element("shape",inputs)==TRUE & aligned.shape==FALSE){
    for (k in 1:dim(Shape$X)[length(dim(Shape$X))]){
      Shape$X[,,,k] <- gpagen(Shape$X[,,,k],print.progress = FALSE)$coords
    }
  }


  if (tree$Y$type=="factor"){
    pred <- factor(rep(NA, length(id.pred)),levels=tree$Ylevels)
  }

  else{pred <- rep(NA,length(id.pred))}

  for (i in 1:length(id.pred)){

    if (is.element("curve",inputs)==TRUE) wCurve <- which(Curve$id==id.pred[i])
    if (is.element("scalar",inputs)==TRUE) wScalar <- which(Scalar$id==id.pred[i])
    if (is.element("factor",inputs)==TRUE) wFactor <- which(Factor$id==id.pred[i])
    if (is.element("shape",inputs)==TRUE) wShape <- which(Shape$id==id.pred[i])
    if (is.element("image",inputs)==TRUE) wImage <- which(Image$id==id.pred[i])

    noeud_courant <- 1

    while (is.element(noeud_courant, tree$feuilles)==FALSE){

      X <- get(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),1]))
      type <- str_to_lower(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),1]))
      var.split <- as.numeric(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),3]))

      # Maintenant il nous faut regarder la difference entre la moyenne à gauche et a droite et conclure :

      meanG <- tree$hist_nodes[[2*noeud_courant]]
      meanD <- tree$hist_nodes[[2*noeud_courant+1]]

      if (type=="curve"){
        distG <- distFrechet(meanG[,1], meanG[,2], X$time[wCurve], X$X[wCurve,var.split], timeScale = timeScale)
        distD <- distFrechet(meanD[,1], meanD[,2], X$time[wCurve], X$X[wCurve,var.split], timeScale = timeScale)
      }
      if (type=="scalar"){
        distG <- abs(meanG- X$X[wScalar,var.split])
        distD <- abs(meanD-X$X[wScalar,var.split])
      }

      if (type=="shape"){
        elementz <- array(X$X[,,wShape,var.split],dim = c(nrow(meanG),ncol(meanG),1))
        distG <- ShapeDist(elementz,meanG)
        distD <- ShapeDist(elementz, meanD)
      }

      if (type=="image"){
        distG <- mean((X$X[wImage,,var.split]-meanG)^2)
        distD <- mean((X$X[wImage,,var.split]-meanD)^2)
      }

      if (type=="factor"){
        distG <- -1*(is.element(X$X[wFactor,var.split],meanG))
        distD <- -1*(is.element(X$X[wFactor,var.split],meanD))
      }

      if (distG <= distD) { noeud_courant <- 2*noeud_courant}
      if (distD < distG) {noeud_courant <- 2*noeud_courant +1}


    }

    if(tree$Y$type=="curve" || tree$Y$type=="image" || tree$Y$type=="shape"){
      pred[i] <- noeud_courant
    }

    else{
      pred[i] <- tree$Y_pred[[noeud_courant]]
    }
  }
  return(pred)
}
