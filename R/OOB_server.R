#' OOB for random Forest
#'
#' @param rf [list]:
#' @param Curve [list]:
#' @param Scalar [list]:
#' @param Factor [list]:
#' @param Shape [list]:
#' @param Image [list]:
#' @param Y [list]:
#' @param timeScale [numeric]
#' @param d_out [numeric]
#' @param range [vector]:
#' @param ncores [numeric]:
#'
#' @import stringr
#' @import kmlShape
#' @import emdist
#'
#' @export
OOB.server <- function(Curve=NULL, Scalar=NULL, Factor=NULL, Shape=NULL, Image=NULL,
                               ncores=NULL,range=NULL, Y, timeScale=0.1, d_out=0.1){

  ### Pour optimiser le code il faudra virer cette ligne et ne le calculer qu'une seule fois !
  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs

  trees = list.files()
  ntree = length(trees)

  if(is.null(ncores)==TRUE){
    ncores <- detectCores()
  }

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

    for (i in 1:length(unique(Y$id)[range])){
      indiv <- unique(Y$id)[range[i]]
      w_y <- which(Y$id==indiv)

      pred_courant = NULL

      pred_courant <- foreach::foreach(t = 1:ntree,.packages = "kmlShape" ,.combine = "rbind") %dopar% {


        tree = get(load(trees[t]))
        BOOT <- tree$boot
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

          pred <- pred.FT(tree,Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape=Shape_courant,Image=Image_courant, timeScale = timeScale)
          res <- cbind(rep(t,dim(tree$Y_pred[[pred]])[1]),tree$Y_pred[[pred]])
        }
      }


      mean_pred <- meanFrechet(pred_courant, timeScale = d_out)
      dp <- as.data.frame(Curve.reduc.times(mean_pred$times, mean_pred$traj, Y$time[w_y]))
      names(dp) <- c("x","y")
      err[i] <- distFrechet(dp$x, dp$y, Y$time[w_y], Y$Y[w_y], timeScale = d_out)^2
    }
  }

  if (Y$type=="scalar"){

    for (i in 1:length(range)){
      indiv <- Y$id[range[i]]
      w_y <- which(Y$id==indiv)
      pred_courant <- NULL

      pred_courant <- foreach::foreach(t = 1:ntree,.packages = "kmlShape" ,.combine = "c") %dopar% {

        tree = get(load(trees[t]))
        BOOT <- tree$boot
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

          res <- pred.FT(tree,Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,
                          Shape=Shape_courant,Image=Image_courant, timeScale = timeScale)
        }
      }

      err[i] <- (mean(pred_courant)-Y$Y[w_y])^2

    }
  }


  if (Y$type=="image"){

    for (i in 1:length(range)){
      indiv <- Y$id[range[i]]
      w_y <- which(Y$id==indiv)
      pred_courant <- NULL

      pred_courant <- foreach::foreach(t = 1:ntree,.packages = "kmlShape" ,.combine = "c") %dopar% {

        tree = get(load(trees[t]))
        BOOT <- tree$boot
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

          res <- pred.FT(tree,Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,
                         Shape=Shape_courant,Image=Image_courant, timeScale = timeScale)
        }
      }

      err[i] <- (apply(pred_courant,2,'mean')-Y$Y[w_y,])^2
    }
  }

  return(err)

}
