#' Predict with Frechet random forests
#'
#' @param object : Frechet random forest
#' @param Curve [list]:
#' @param Scalar [list]:
#' @param Factor [list]:
#' @param Shape [list]:
#' @param Image [list]:
#' @param timeScale [numeric]:
#' @param d_out [numeric]:
#' @param ... : optional parameters to be passed to the low level function
#'
#' @import kmlShape
#' @import stringr
#' @import Evomorph
#' @import geomorph
#'
#' @export
#'
predict.FrechForest <- function(object, Curve=NULL,Scalar=NULL,Factor=NULL,Shape=NULL, Image=NULL, timeScale=0.1, d_out=0.1,...){
  # La première etape est de toujours lire les predicteurs ::

  if (is.null(Curve)==FALSE){
    Curve <- list(type="curve",X=Curve$X,id=Curve$id,time=Curve$time)
  }
  if (is.null(Scalar)==FALSE){
    Scalar <- list(type="scalar",X=Scalar$X,id=Scalar$id)
  }
  if (is.null(Factor)==FALSE){
    Factor <- list(type="factor",X=Factor$X,id=Factor$id)
  }
  if (is.null(Shape)==FALSE){
    Shape <- list(type="shape",X=Shape$X,id=Shape$id)
  }
  if (is.null(Image)==FALSE){
    Image <- list(type="image",X=Image$X,id=Image$id)
  }

  ## Puis on prend les predicteurs:

  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs
  # On va les lires en mettant la maj sur les differents elements qui le constituent :

  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  Id.pred <- unique(get(Inputs[1])$id)
  pred.feuille <- matrix(0, ncol(object$rf), length(Id.pred))

  if (object$type=="factor"){
    pred.feuille <- as.data.frame(matrix(0, ncol(object$rf), length(Id.pred)))
  }

  for (t in 1:ncol(object$rf)){
    pred.feuille[t,] <- pred.FT(object$rf[,t], Curve = Curve,Scalar = Scalar,Factor=Factor,Shape=Shape,Image=Image, timeScale)
  }

  if (object$type=="scalar"){
    pred <- apply(pred.feuille, 2, "mean")
    return(pred)
  }

  if (object$type=="factor"){
    pred.all <- apply(pred.feuille, 2, "table")
    val <- factor(rep(NA, length(pred.all)), levels=object$levels)
    proba <- rep(NA, length(pred.all))
    for (k in 1:length(pred.all)){
      val[k] <- factor(attributes(which.max(pred.all[[k]])))
      proba[k] <- max(pred.all[[k]])/ncol(object$rf)
    }
    prediction <- data.frame(pred=val, prob=proba)
    return(prediction)
  }

  if (object$type=="shape"){
    pred <- array(0, dim=c(object$size[1], object$size[2],length(Id.pred)))
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- array(0,dim=c(object$size[1],object$size[2],ncol(object$rf)))
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant[,,k] <- object$rf[,k]$Y_pred[[pred.feuille[k,l]]]
      }
      Ms <- mshape(pred_courant)
      M <- matrix(0,dim(Ms)[1], dim(Ms)[2])
      M[,1] <- Ms[,1]
      M[,2] <- Ms[,2]
      pred[,,l] <- M
    }
  }

  if (object$type=="image"){
    pred <- matrix(0, length(Id.pred), object$size[2])
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- matrix(0,ncol(object$rf),object$size[2])
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant[,,k] <- object$rf[,k]$Y_pred[[pred.feuille[k,l]]]
      }
      pred[,,l] <-  apply(pred_courant, 2, "mean")
    }
  }

  if (object$type=="curve"){
    pred <- NULL
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- NULL
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant <- rbind(pred_courant, cbind(rep(k,dim(object$rf[,k]$Y_pred[[pred.feuille[k,l]]])[1]),object$rf[,k]$Y_pred[[pred.feuille[k,l]]]))
      }
      predouille <- kmlShape::meanFrechet(pred_courant, timeScale = timeScale)
      predouille <- cbind(predouille, rep(Id.pred[l],dim(predouille)[1]))
      pred <- rbind(pred, predouille)
    }
    names(pred) <- c("times", "traj", "ID")
  }
  return(pred)
}
