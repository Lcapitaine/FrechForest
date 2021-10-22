#' Randomized Frechet tree
#'
#' @param Curve [list]: A list that contains the different input curves. It must contain the following elements (no choice): \code{X} the matrix of the different curves, each column code for a different curve variable; \code{id} is the vector of the identifiers for the different trajectories contained in \code{X}; \code{time} is the vector of the measurement times associated with the trajectories contained in \code{X}.
#' @param Scalar [list]: A list that contains the different input scalars. It must contain the following elements (no choice):  \code{X} the matrix of the scalars, each column code for a different variable; \code{id} is the vector of the identifiers for each individual.
#' @param Factor [list]: A list that contains the different input factors. It must contain the following elements (no choice):  \code{X} the matrix of the factors, each column code for a different variable; \code{id} is the vector of the identifiers for each individual.
#' @param Shape [list]: A list that contains the different input shapes. It must contain the following elements (no choice):  \code{X} the array of the shapes of dimension \code{n}x2x\code{l}x\code{p} where \code{n} is the number of points for composing each shape, \code{l} is the number of shapes and \code{p} is the number of shapes variables, \code{id} is the vector of the identifiers for each individual.
#' @param Image [list]: A list that contains the different input images. It must contain the following elements (no choice):  \code{X} the array of the images of dimension \code{n}x\code{m}x\code{l}x\code{p} where \code{n}*\code{m} is the size of each image, \code{l} is the number of images and \code{p} is the number of shapes variables; \code{id} is the vector of the identifiers for each individual.
#' @param Y [list]: A list that contains the output, It must contain the following elements (no choice): \code{type} defines the nature of the output, can be "\code{curve}", "\code{sclalar}", "\code{factor}", "\code{shape}", "\code{image}"; \code{Y} is the output variable; \code{id} is the vector of the identifiers for each individuals, they should be the same as the identifiers of the inputs.
#' @param mtry [integer]: Number of variables randomly sampled as candidates at each split. The default value \code{p/3}
#' @param ERT [logical]: Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.
#' @param timeScale [numeric]: Allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping. Only used when there are trajectories either in input or output.
#' @param ntry [numeric]: Only with \code{ERT=TRUE}, allows to manage with randomness of the trees.
#' @param ... :  optional parameters to be passed to the low level function
#'
#' @import kmlShape
#' @import stringr
#' @import Evomorph
#' @import geomorph
#'
#' @export
Rtmax <- function(Curve=NULL, Scalar=NULL, Factor=NULL, Shape=NULL, Image=NULL,Y,mtry,ERT=FALSE,ntry=3, timeScale=0.1, ...){


  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs

  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  impurity_feuilles <- NULL
  V_split <- NULL
  hist_nodes <- list()
  id_boot <- unique(sample(unique(Y$id), length(unique(Y$id)), replace=TRUE))
  boot <- id_boot
  decoupe <- 1

  wXCurve <- NULL
  wXScalar <- NULL
  wXFactor <- NULL
  wXShape <- NULL
  wXImage <- NULL
  wY <- NULL

  for (k in id_boot){
    wY <- c(wY, which(Y$id==k))
    if (is.element("curve",inputs)==TRUE) wXCurve <- c(wXCurve, which(Curve$id==k))
    if (is.element("scalar",inputs)==TRUE) wXScalar <- c(wXScalar, which(Scalar$id==k))
    if (is.element("factor",inputs)==TRUE) wXFactor <- c(wXFactor, which(Factor$id==k))
    if (is.element("shape",inputs)==TRUE) wXShape <- c(wXShape, which(Shape$id==k))
    if (is.element("image",inputs)==TRUE) wXImage <- c(wXImage, which(Image$id==k))
  }

  Y_pred <- list()
  Y_pred_surv  <- list()

  if (is.element("curve",inputs)==TRUE) Curve_boot <- list(type=Curve$type,   X=Curve$X[wXCurve,, drop=FALSE], id= Curve$id[wXCurve], time = Curve$time[wXCurve]) ### bootstrap pour les courbes
  if (is.element("scalar",inputs)==TRUE) Scalar_boot <- list(type=Scalar$type,   X=Scalar$X[wXScalar,, drop=FALSE], id= Scalar$id[wXScalar]) ### bootstrap pour les courbes
  if (is.element("factor",inputs)==TRUE) Factor_boot <- list(type=Factor$type,   X=Factor$X[wXFactor,, drop=FALSE], id= Factor$id[wXFactor])
  if (is.element("shape",inputs)==TRUE) Shape_boot <- list(type=Shape$type,   X=Shape$X[,,wXShape, , drop=FALSE], id= Shape$id[wXShape])
  if (is.element("image",inputs)==TRUE) Image_boot <- list(type=Image$type,   X=Image$X[wXImage,,, drop=FALSE], id= Image$id[wXImage])


  if (Y$type=="curve") {Y_boot <- list(type=Y$type,Y=Y$Y[wY], id=Y$id[wY], time=Y$time[wY])} ### idem pour Y
  if (Y$type=="shape") {Y_boot <- list(type=Y$type, Y=Y$Y[,,wY], id=Y$id[wY])}
  if (Y$type=="image") {Y_boot <- list(type=Y$type, Y=Y$Y[wY,], id=Y$id[wY])}
  if (Y$type=="factor" || Y$type=="scalar") {Y_boot <- list(type=Y$type,Y=Y$Y[wY], id=Y$id[wY])}


  imp_nodes <- list()
  imp_nodes[[1]] = Inf
  impurete = Inf

  id_feuille <- rep(1,length(Y_boot$id)) #### localisation des feuilles de l'arbre
  id_feuille_prime <- id_feuille

  for (p in 1:(length(unique(Y_boot$id))/2-1)){
    count_split <- 0
    for (i in 1:length(unique(id_feuille))){
      # Il faut que l'on regarde le tirage des variables de manière aleatoire :
      V <- NULL
      for (v in Inputs){
        V <- c(V, rep(get(v)$type,dim(get(v)$X)[length(dim(get(v)$X))]))
      }
      variables <- sample(V,mtry) # Maintenant on sait combien on doit en tirer dans chaque espace
      # On ne va regarder que les espaces tires :
      split.spaces <- unique(variables)

      # variables <- sample(c(1:dim(X_boot$X[,,drop=FALSE])[2]),mtry)
      w <- which(id_feuille==unique(id_feuille)[i])
      wXCurve <- NULL
      wXScalar <- NULL
      wXFactor <- NULL
      wXShape <- NULL
      wXImage <- NULL

      for (l in unique(Y_boot$id[w])){
        if (is.element("curve",inputs)==TRUE) wXCurve <- c(wXCurve, which(Curve_boot$id==l))
        if (is.element("scalar",inputs)==TRUE) wXScalar <- c(wXScalar, which(Scalar_boot$id==l))
        if (is.element("factor",inputs)==TRUE) wXFactor <- c(wXFactor, which(Factor_boot$id==l))
        if (is.element("shape",inputs)==TRUE) wXShape <- c(wXShape, which(Shape_boot$id==l))
        if (is.element("image",inputs)==TRUE) wXImage <- c(wXImage, which(Image_boot$id==l))
      }

      if (length(unique(Y_boot$id[w]))>1 & imp_nodes[[unique(id_feuille)[i]]] >0){

        # On est ici

        if (is.element("curve",split.spaces)==TRUE){
          tirageCurve <- sample(1:ncol(Curve$X),length(which(variables=="curve")))
          Curve_courant <- list(type = Curve_boot$type, X=Curve_boot$X[wXCurve,tirageCurve, drop=FALSE], id=Curve_boot$id[wXCurve, drop=FALSE], time=Curve_boot$time[wXCurve, drop=FALSE])
        }

        if (is.element("scalar",split.spaces)==TRUE){
          tirageScalar <- sample(1:ncol(Scalar$X),length(which(variables=="scalar")))
          Scalar_courant <- list(type = Scalar_boot$type, X=Scalar_boot$X[wXScalar,tirageScalar, drop=FALSE], id=Scalar_boot$id[wXScalar, drop=FALSE])
        }

        if (is.element("factor",split.spaces)==TRUE){
          tirageFactor <- sample(1:ncol(Factor$X),length(which(variables=="factor")))
          Factor_courant <- list(type = Factor_boot$type, X=Factor_boot$X[wXFactor,tirageFactor, drop=FALSE], id=Factor_boot$id[wXFactor, drop=FALSE])
        }

        if (is.element("shape",split.spaces)==TRUE){
          tirageShape <- sample(1:dim(Shape$X)[length(dim(Shape$X))],length(which(variables=="shape")))
          Shape_courant <- list(type = Shape_boot$type, X=Shape_boot$X[,,wXShape,tirageShape, drop=FALSE], id=Shape_boot$id[wXShape, drop=FALSE])
        }

        if (is.element("image",split.spaces)==TRUE){
          tirageImage <- sample(1:dim(Image$X)[3],length(which(variables=="image")))
          Image_courant <- list(type = Image_boot$type, X=Image_boot$X[wXImage,,tirageImage, drop=FALSE], id=Image_boot$id[wXImage])
        }

        if (Y_boot$type=="curve"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[w], id=Y_boot$id[w], time=Y_boot$time[w])
        }

        if (Y_boot$type=="shape"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[,,w, drop=FALSE], id=Y_boot$id[w, drop=FALSE])
        }

        if (Y_boot$type=="image"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[w, ,drop=FALSE], id=Y_boot$id[w, drop=FALSE])
        }


        if (Y_boot$type=="factor" || Y_boot$type=="scalar"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[w, drop=FALSE], id=Y_boot$id[w, drop=FALSE])
        }


        F_SPLIT <- data.frame(vari = as.character(split.spaces), hetero= rep(Inf,length(split.spaces)))

        if (is.element("factor",split.spaces)==TRUE){

          if( ERT==FALSE){
            feuille_split_factor = list(Pure = TRUE)
            tryCatch({
              feuille_split_factor <-  var_split(Factor_courant,Y_courant,timeScale)
            }, error = function(sp){feuille_split_factor = list(Pure = TRUE)})
          }

          else{
            feuille_split_factor = list(Pure = TRUE)
            tryCatch({feuille_split_factor <- ERvar_split(X=Factor_courant,Y=Y_courant,timeScale=timeScale,ntry = ntry)
          }, error = function(sp){feuille_split_factor = list(Pure = TRUE)})
          }

          if (feuille_split_factor$Pure==FALSE){
            F_SPLIT[which(F_SPLIT[,1]=="factor"),2] <- feuille_split_factor$impurete}
        }

        if (is.element("curve",split.spaces)==TRUE){

          if( ERT==FALSE){
            feuille_split_curve = list(Pure = TRUE)
            tryCatch({
              feuille_split_curve <-  var_split(Curve_courant,Y_courant,timeScale)
            }, error = function(sp){feuille_split_curve = list(Pure = TRUE)})
          }

          else{
            feuille_split_curve = list(Pure = TRUE)
            tryCatch({feuille_split_curve <- ERvar_split(X=Curve_courant,Y=Y_courant,timeScale=timeScale,ntry = ntry)
            }, error = function(sp){feuille_split_curve = list(Pure = TRUE)})
          }

          if (feuille_split_curve$Pure==FALSE){
            F_SPLIT[which(F_SPLIT[,1]=="curve"),2] <- feuille_split_curve$impurete}

        }

        if (is.element("scalar",split.spaces)==TRUE){

          if( ERT==FALSE){
            feuille_split_scalar = list(Pure = TRUE)
            tryCatch({
              feuille_split_scalar <-  var_split(Scalar_courant,Y_courant,timeScale)
            }, error = function(sp){feuille_split_scalar = list(Pure = TRUE)})
          }

          else{
            feuille_split_scalar = list(Pure = TRUE)
            tryCatch({feuille_split_scalar <- ERvar_split(X=Scalar_courant,Y=Y_courant,timeScale=timeScale,ntry = ntry)
            }, error = function(sp){feuille_split_scalar = list(Pure = TRUE)})
          }

          if (feuille_split_scalar$Pure==FALSE){
            F_SPLIT[which(F_SPLIT[,1]=="scalar"),2] <- feuille_split_scalar$impurete}

        }

        if (is.element("shape",split.spaces)==TRUE){

          if( ERT==FALSE){
            feuille_split_shape = list(Pure = TRUE)
            tryCatch({
              feuille_split_shape <-  var_split(Shape_courant,Y_courant,timeScale)
            }, error = function(sp){feuille_split_shape = list(Pure = TRUE)})
          }

          else{
            feuille_split_shape = list(Pure = TRUE)
            tryCatch({feuille_split_shape <- ERvar_split(X=Shape_courant,Y=Y_courant,timeScale=timeScale,ntry = ntry)
            }, error = function(sp){feuille_split_shape = list(Pure = TRUE)})
          }

          if (feuille_split_shape$Pure==FALSE){
            F_SPLIT[which(F_SPLIT[,1]=="shape"),2] <- feuille_split_shape$impurete}

        }


        if (is.element("image",split.spaces)==TRUE){

          if( ERT==FALSE){
            feuille_split_image = list(Pure = TRUE)
            tryCatch({
              feuille_split_image <-  var_split(Image_courant,Y_courant,timeScale)
            }, error = function(sp){feuille_split_image = list(Pure = TRUE)})
          }

          else{
            feuille_split_image = list(Pure = TRUE)
            tryCatch({feuille_split_image <- ERvar_split(X=Image_courant,Y=Y_courant,timeScale=timeScale,ntry = ntry)
            }, error = function(sp){feuille_split_image = list(Pure = TRUE)})
          }

          if (feuille_split_image$Pure==FALSE){
            F_SPLIT[which(F_SPLIT[,1]=="image"),2] <- feuille_split_image$impurete}

        }


        if (min(F_SPLIT[,2])<Inf){

          TYPE <- as.character(F_SPLIT[which.min(F_SPLIT[,2]),1])
          TYPE2 <- paste(str_to_upper(str_sub(TYPE,1,1)),str_sub(TYPE,2,nchar(TYPE)),sep="")
          X <- get(TYPE2)
          X_boot <- get(paste(TYPE2,"_boot",sep=""))

          feuille_split <- get(paste("feuille_split_",TYPE, sep=""))

          vsplit_space <- get(paste("tirage",TYPE2, sep=""))[feuille_split$variable]

          #if (imp_apres_split<imp_avant_split){

          gauche_id <- unique(Y_boot$id[w])[which(feuille_split$split==1)]
          droit_id <- unique(Y_boot$id[w])[which(feuille_split$split==2)]


          imp_nodes[[2*unique(id_feuille)[i]]] <- feuille_split$impur_list[[1]]
          imp_nodes[[2*unique(id_feuille)[i]+1]] <- feuille_split$impur_list[[2]]


          V_split <- rbind(V_split,c(TYPE2,unique(id_feuille)[i],vsplit_space))

          wY_gauche <- NULL
          wY_droit <- NULL
          w_gauche <- NULL
          w_droit <- NULL


          for (k in 1:length(gauche_id)){
            w_gauche <- c(w_gauche, which(X_boot$id==gauche_id[k]))
            wY_gauche <- c(wY_gauche, which(Y_boot$id==gauche_id[k]))
          }

          for (k in 1:length(droit_id)){
            w_droit <- c(w_droit, which(X_boot$id==droit_id[k]))
            wY_droit <- c(wY_droit, which(Y_boot$id==droit_id[k]))
          }


          id_feuille_prime[wY_gauche] <- 2*(unique(id_feuille)[i])
          id_feuille_prime[wY_droit] <- 2*(unique(id_feuille)[i])+1

          #print(paste("Split on the variable", vsplit_space, "on the space of ", paste(TYPE,"s",sep="")))

          if (X$type=="curve"){
            trajG <- as.data.frame(cbind(X_boot$id[w_gauche], X_boot$time[w_gauche], X_boot$X[w_gauche,vsplit_space]))
            trajD <- as.data.frame(cbind(X_boot$id[w_droit], X_boot$time[w_droit], X_boot$X[w_droit,vsplit_space]))
            meanFg <- as.matrix(kmlShape::meanFrechet(trajG))
            meanFd <- as.matrix(kmlShape::meanFrechet(trajD))
          }

          if (X$type=="shape"){
            meanFg <- mshape(X_boot$X[,,w_gauche,vsplit_space])
            meanFd <- mshape(X_boot$X[,,w_droit,vsplit_space])
          }

          if (X$type=="image"){
            meanFg <- apply(X_boot$X[w_gauche,,vsplit_space, drop=FALSE],2,"mean")
            meanFd <- apply(X_boot$X[w_droit,,vsplit_space, drop=FALSE],2,"mean")
          }

          if (X$type=="factor"){
            meanFg <- unique(X_boot$X[w_gauche, vsplit_space])
            meanFd <- unique(X_boot$X[w_droit,vsplit_space])
          }

          if (X$type=="scalar"){
            meanFg <- mean(X_boot$X[w_gauche,vsplit_space])
            meanFd <- mean(X_boot$X[w_droit,vsplit_space])
          }


          hist_nodes[[2*(unique(id_feuille)[i])]] <- meanFg
          hist_nodes[[2*(unique(id_feuille)[i])+1]] <- meanFd
          count_split <- count_split+1

          feuilles_courantes <- unique(id_feuille_prime)
        }


      }
    }

    id_feuille <- id_feuille_prime

    if (count_split ==0 ){

      V_split <- data.frame(V_split)
      names(V_split) <- c("type","num_noeud", "var_split")
      for (q in unique(id_feuille)){
        w <- which(id_feuille == q)
        if (Y$type=="curve"){
          Y_pred[[q]] <- kmlShape::meanFrechet(data.frame(Y_boot$id[w], Y_boot$time[w], Y_boot$Y[w]))
        }
        if(Y$type=="scalar"){
          Y_pred[[q]]<- mean(Y_boot$Y[w])
        }
        if(Y$type=="factor"){
          Table <- which.max(table(Y_boot$Y[w]))
          Y_pred[[q]] <-  as.factor(attributes(Table)$names)
        }

        if (Y$type=="shape"){
          Y_pred[[q]] <-  mshape(Y_boot$Y[,,w, drop=FALSE])
        }

        if (Y$type=="image"){
          Y_pred[[q]] <- apply(Y_boot$Y[w,,drop=FALSE], 2, "mean")
        }

      }
      if (Y$type=="factor"){
        Ylevels <- unique(Y_boot$Y)
        return(list(feuilles = id_feuille, idY=Y_boot$id,Ytype=Y_boot$type, V_split=V_split, hist_nodes=hist_nodes, Y_pred = Y_pred, time = time, Y=Y, boot=boot, Ylevels=Ylevels))
      }
      return(list(feuilles = id_feuille, idY=Y_boot$id,Ytype=Y_boot$type, V_split=V_split, hist_nodes=hist_nodes, Y_pred = Y_pred, time = time, Y=Y, boot=boot, Y_pred_surv=Y_pred_surv))
    }
  }


  V_split <- data.frame(V_split)
  names(V_split) <- c("type","num_noeud", "var_split")
  for (q in unique(id_feuille)){
    w <- which(id_feuille == q)
    if (Y$type=="curve"){
      Y_pred[[q]] <- kmlShape::meanFrechet(data.frame(Y_boot$id[w], Y_boot$time[w], Y_boot$Y[w]))
    }

    if (Y$type=="image"){
      Y_pred[[q]] <- apply(Y_boot$Y[w,,drop=FALSE], 2, "mean")
    }

    if(Y$type=="scalar"){
      Y_pred[[q]]<- mean(Y_boot$Y[w])
    }

    if(Y$type=="factor"){
      Table <- which.max(table(Y_boot$Y[w]))
      Y_pred[[q]] <-  as.factor(attributes(Table)$names)
    }

    if (Y$type=="shape"){
      Y_pred[[q]] <- mshape(Y_boot$Y[,,w, drop=FALSE])
    }

  }
  if (Y$type=="factor"){
    Ylevels <- unique(Y_boot$Y)
    return(list(feuilles = id_feuille, idY=Y_boot$id,Ytype=Y_boot$type, V_split=V_split, hist_nodes=hist_nodes, Y_pred = Y_pred, time = time, Y=Y, Ylevels=Ylevels, boot=boot))
  }
  return(list(feuilles = id_feuille,Ytype=Y_boot$type, idY=Y_boot$id, V_split=V_split, hist_nodes=hist_nodes, Y_pred= Y_pred, time=time, Y=Y, boot=boot, Y_pred_surv=Y_pred_surv))
}
