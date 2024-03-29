#' OOB for random Forest
#'
#' @param rf
#' @param Curve
#' @param Scalar
#' @param Factor
#' @param Shape
#' @param Image
#' @param Y
#' @param timeScale
#' @param d_out
#'
#' @import stringr
#' @import kmlShape
#' @import Evomorph
#' @import geomorph
#'
#' @keywords internal
OOB.rfshape <- function(rf, Curve=NULL, Scalar=NULL, Factor=NULL, Shape=NULL, Image=NULL, Y, timeScale=0.1, d_out=0.1){

  ### Pour optimiser le code il faudra virer cette ligne et ne le calculer qu'une seule fois !
  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs


  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }


  err <- rep(NA,length(unique(Y$id)))

  Curve_courant <- NULL
  Scalar_courant <- NULL
  Factor_courant <- NULL
  Shape_courant <- NULL
  Image_courant <- NULL

  if (Y$type=="curve"){
    oob.pred <- list()
    #errdp <- rep(NA,length(unique(id)))

    for (i in 1:length(unique(Y$id))){
      indiv <- unique(Y$id)[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- NULL
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        if (is.element(indiv, oob)== TRUE){

          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
          }

          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }

          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }

          if (is.element("shape",inputs)==TRUE){
            w_XShape <- which(Shape$id== indiv)
            Shape_courant <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
          }

          if (is.element("image",inputs)==TRUE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(type="image", X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }

          pred <- pred.FT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape=Shape_courant,Image=Image_courant, timeScale = timeScale)
          courbe <- rf$rf[,t]$Y_pred[[pred]]
          pred_courant <- rbind(cbind(rep(t,dim(courbe)[1]),courbe),pred_courant)
        }
      }
      mean_pred <- meanFrechet(pred_courant, timeScale = d_out)
      dp <- as.data.frame(Curve.reduc.times(mean_pred$times, mean_pred$traj, Y$time[w_y]))
      names(dp) <- c("x","y")
      oob.pred[[i]] <- dp
      err[i] <- distFrechet(dp$x, dp$y, Y$time[w_y], Y$Y[w_y], timeScale = d_out)^2
    }
    return(list(err=err,oob.pred=oob.pred))
  }

  if (Y$type=="scalar"){
    oob.pred <- rep(NA, length(unique(Y$id)))
    #errdp <- rep(NA,length(unique(id)))
    for (i in 1:length(Y$id)){
      indiv <- Y$id[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- NULL
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        if (is.element(indiv, oob)== TRUE){

          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
          }

          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }

          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }

          if (is.element("shape",inputs)==TRUE){
            w_XShape <- which(Shape$id== indiv)
            Shape_courant <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
          }

          if (is.element("image",inputs)==TRUE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(type="image", X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }

          pred <- pred.FT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape=Shape_courant,Image=Image_courant, timeScale = timeScale)
          pred_courant <- c(pred_courant, pred)
        }
      }
      oob.pred[i] <- mean(pred_courant)
      err[i] <- (oob.pred[i]-Y$Y[w_y])^2
    }
  }

  if (Y$type=="factor"){
    oob.pred <- factor(rep(NA, length(unique(Y$id))), levels=rf$levels)
    #errdp <- rep(NA,length(unique(id)))
    for (i in 1:length(Y$id)){
      indiv <- Y$id[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- factor(rep(NA, length(unique(Y$id))), levels=rf$levels)
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        if (is.element(indiv, oob)== TRUE){

          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
          }

          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }

          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }

          if (is.element("shape",inputs)==TRUE){
            w_XShape <- which(Shape$id== indiv)
            Shape_courant <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
          }

          if (is.element("image",inputs)==TRUE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(type="image", X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }

          pred <- pred.FT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape=Shape_courant,Image=Image_courant, timeScale = timeScale)
          pred_courant[t] <- pred
        }
      }
      pred_courant <- na.omit(pred_courant)
      oob.pred[i] <- as.factor(attributes(which.max(table(pred_courant))))
    }
    err <- 1*(oob.pred!=Y$Y)
  }

  if (Y$type=="shape"){
    oob.pred <- array(0,dim=dim(Y$Y))
    #errdp <- rep(NA,length(unique(id)))
    for (i in 1:length(Y$id)){
      indiv <- unique(Y$id)[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- array(0, dim=c(dim(Y$Y)[1],dim(Y$Y)[2],length(rf$rf)))
      selection <- NULL
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        if (is.element(indiv, oob)== TRUE){

          selection <- c(selection, t)

          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
          }

          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }

          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }

          if (is.element("shape",inputs)==TRUE){
            w_XShape <- which(Shape$id== indiv)
            Shape_courant <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
          }

          if (is.element("image",inputs)==TRUE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(type="image", X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }

          pred <- pred.FT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape=Shape_courant,Image=Image_courant, timeScale = timeScale)
          pred_courant[,,t] <- rf$rf[,t]$Y_pred[[pred]]
        }
      }
      pred_courant <- pred_courant[,,selection]
      mean_pred <- mshape(pred_courant)
      err[i] <- ShapeDist(Y$Y[,,w_y, drop=FALSE],mean_pred)^2
      M <- matrix(0,dim(mean_pred)[1], dim(mean_pred)[2])
      M[,1] <- mean_pred[,1]
      M[,2] <- mean_pred[,2]

      oob.pred[,,i] <- M
    }
    return(list(err=err,oob.pred=oob.pred))
  }

  if (Y$type=="image"){
    oob.pred <- array(0,dim=dim(Y$Y))
    err = array(NA, dim = dim(Y$Y))
    #errdp <- rep(NA,length(unique(id)))
    for (i in 1:length(Y$id)){
      indiv <- unique(Y$id)[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- array(NA, dim=c(ncol(rf$rf),ncol(Y$Y)))
      selection <- NULL
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        if (is.element(indiv, oob)== TRUE){
          selection <- c(selection, t)
          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
          }

          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }

          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }

          if (is.element("shape",inputs)==TRUE){
            w_XShape <- which(Shape$id== indiv)
            Shape_courant <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
          }

          if (is.element("image",inputs)==TRUE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(type="image", X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }

          pred <- pred.FT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape=Shape_courant,Image=Image_courant, timeScale = timeScale)
          pred_courant[t,] <- rf$rf[,t]$Y_pred[[pred]]
        }
      }
      oob.pred[i,] <-  apply(na.omit(pred_courant),2,"mean")
      err[i,] <- (oob.pred[i,]-Y$Y[w_y,])^2
    }
    return(list(err=err,oob.pred=oob.pred))
  }
  return(list(err=err,oob.pred=oob.pred))
}
