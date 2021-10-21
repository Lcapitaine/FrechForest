#' Maximal Frechet tree
#'
#' @param Curve [list]:
#' @param Scalar [list]:
#' @param Factor [list]:
#' @param Shape [list]:
#' @param Image [list]:
#' @param Y [list]:
#' @param timeScale [list]:
#' @param ... : Optional parameters
#'
#' @import kmlShape
#' @import stringr
#' @import Evomorph
#' @import geomorph
#'
#' @keywords internal
#'
Tmax <- function(Curve=NULL, Scalar=NULL, Factor=NULL, Shape=NULL, Image=NULL ,Y, timeScale=0.1, ...){

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
    for (k in 1:dim(Shape$X)[length(dim(Shape$X))]){
      Shape$X[,,,k] <- gpagen(Shape$X[,,,k], print.progress = FALSE)$coords
    }
  }

  if (is.null(Image)==FALSE){
    Image <- list(type="image",X=Image$X,id=Image$id)
  }
  # On commence par lire les arguments :
  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs
  # On va les lires en mettant la maj sur les differents elements qui le constituent :

  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }


  id_feuille <- rep(1,length(Y$id))
  id_feuille_prime <- id_feuille
  V_split <- NULL
  hist_nodes <- list()
  decoupe <- 1
  imp_nodes <- list()
  imp_nodes[[1]] = Inf
  impurete = Inf

  Y_pred <- list()
  Y_pred_imputation <- list()
  hist_imp_nodes <- as.matrix(cbind(1, impurete,length(unique(Y$id))))

  for (p in 1:(length(unique(Y$id))/2-1)){
    count_split <- 0
    for (i in 1:length(unique(id_feuille))){

      w_Y <- which(id_feuille==unique(id_feuille)[i])
      #On rÃ©cupÃ¨re les identifiants des individus dans ces feuilles
      which_feuilles <- unique(Y$id[w_Y])

      ### Il faut trouver les moyen d'avoir de manière automatique les indexes de chaque entree

      w_XCurve <- NULL
      w_XScalar <- NULL
      w_XFactor <- NULL
      for (l in which_feuilles){

        if (is.element("curve",inputs)==TRUE)  w_XCurve <- c(w_XCurve, which(Curve$id==l))
        if (is.element("scalar",inputs)==TRUE) w_XScalar <- c(w_XScalar, which(Scalar$id==l))
        if (is.element("factor",inputs)==TRUE) w_XFactor <- c(w_XFactor, which(Factor$id==l))
      }

      if (length(unique(Y$id[w_Y]))>1 & imp_nodes[[unique(id_feuille)[i]]] >0){

        #Il nous faut les entrees :

        if (is.element("curve",inputs)==TRUE) Curve_courant <- list(type=Curve$type, X=Curve$X[w_XCurve,],id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
        if (is.element("scalar",inputs)==TRUE) {Scalar_courant <- list(type=Scalar$type, X=Scalar$X[w_XScalar,], id=Scalar$id[w_XScalar])}
        if (is.element("factor",inputs)==TRUE) {Factor_courant <- list(type=Factor$type, X=Factor$X[w_XFactor,], id=Factor$id[w_XFactor])}

        #Il nous faut les sorties maintenant

        if (Y$type=="curve"){Y_courant <- list(type=Y$type, Y=Y$Y[w_Y], id=Y$id[w_Y], time=Y$time[w_Y])}
        if(Y$type=="factor" || Y$type=="scalar"){Y_courant <- list(type=Y$type, Y=Y$Y[w_Y], id=Y$id[w_Y])}
        if (Y$type=="shape") { Y_courant <- list(type="shape",Y=Y$Y[,,w_Y], id=Y$id[w_Y])}
        if (Y$type=="image"){Y_courant <- list(type="image",Y=Y$Y[w_Y,], id=Y$id[w_Y])}

        # Il nous faut maintenant faire le split sur toutes les differents types :

        F_SPLIT <- NULL
        decoupe <- 0

        if (is.element("curve",inputs)==TRUE){
          feuille_split_Curve <- var_split(Curve_courant,Y_courant, timeScale)
          F_SPLIT <- rbind(F_SPLIT,c("Curve",feuille_split_Curve$impurete))
          decoupe <- decoupe+1
        }

        if (is.element("scalar",inputs)==TRUE){
          feuille_split_Scalar <- var_split(Scalar_courant,Y_courant,timeScale)
          if (feuille_split_Scalar$Pure==FALSE){
            F_SPLIT <- rbind(F_SPLIT,c("Scalar",feuille_split_Scalar$impurete))
            decoupe <- decoupe +1
          }
        }

        if (is.element("factor",inputs)==TRUE){
          feuille_split_Factor <- var_split(Factor_courant,Y_courant,timeScale)
          if (feuille_split_Factor$Pure==FALSE){
            F_SPLIT <- rbind(F_SPLIT,c("Factor",feuille_split_Factor$impurete))
            decoupe <- decoupe +1
          }
        }


        if (decoupe >0){
          TYPE <- F_SPLIT[which.min(F_SPLIT[,2]),1]
          X <- get(TYPE)
          w_X <- get(paste("w_X",TYPE, sep=""))
          # Ensuite on va tout comparer

          feuille_split <- get(paste("feuille_split_",TYPE, sep="")) ####  renvoie l'impurete ainsi que le split ? gauche et a droite

          #if (imp_apres_split<imp_avant_split){

          if (Y$type=="curve"){
            Y_pred[[unique(id_feuille)[i]]] <- kmlShape::meanFrechet(as.data.frame(cbind(Y$id[w_Y], Y$time[w_Y], Y$Y[w_Y])))}

          if (Y$type=="scalar"){
            Y_pred[[unique(id_feuille)[i]]] <- mean(Y$Y[w_Y])
          }

          if (Y$type=="shape"){
            Y_pred[[unique(id_feuille)[i]]] <- mshape(Y$Y[,,w_Y, drop=FALSE])
          }

          if (Y$type=="factor"){ #On va renvoyer la classe dominante ::
            Table <- which.max(table(Y$Y[w_Y]))
            Y_pred[[unique(id_feuille)[i]]] <- as.factor(attributes(Table)$names)
          }

          if (Y$type=="image"){
            Y_pred[[unique(id_feuille)[i]]] <- apply(Y$Y[w_Y,, drop=FALSE], 2, "mean")
          }

          else {
            imp_nodes[[2*unique(id_feuille)[i]]] <- feuille_split$impur_list[[1]]
            imp_nodes[[2*unique(id_feuille)[i]+1]] <- feuille_split$impur_list[[2]]
          }


          hist_imp_nodes <- rbind(hist_imp_nodes, c(2*unique(id_feuille)[i],imp_nodes[[2*unique(id_feuille)[i]]], length(which(feuille_split$split==1))))
          hist_imp_nodes <- rbind(hist_imp_nodes, c(2*unique(id_feuille)[i]+1,imp_nodes[[2*unique(id_feuille)[i]+1]], length(which(feuille_split$split==2))))


          gauche_id <- unique(Y$id[w_Y])[which(feuille_split$split==1)]
          droit_id <- unique(Y$id[w_Y])[which(feuille_split$split==2)]

          V_split <- rbind(V_split,c(TYPE,unique(id_feuille)[i],feuille_split$variable))
          wY_gauche <- NULL
          wY_droit <- NULL
          w_gauche <- NULL
          w_droit <- NULL

          for (k in 1:length(gauche_id)){
            w_gauche <- c(w_gauche, which(X$id==gauche_id[k]))
            wY_gauche <- c(wY_gauche, which(Y$id==gauche_id[k]))
          }

          for (k in 1:length(droit_id)){
            w_droit <- c(w_droit, which(X$id==droit_id[k]))
            wY_droit <- c(wY_droit, which(Y$id==droit_id[k]))
          }

          id_feuille_prime[wY_gauche] <- 2*(unique(id_feuille)[i])
          id_feuille_prime[wY_droit] <- 2*(unique(id_feuille)[i])+1

          if (X$type=="curve"){
            trajG <- as.data.frame(cbind(X$id[w_gauche], X$time[w_gauche], X$X[w_gauche,feuille_split$variable]))
            trajD <- as.data.frame(cbind(X$id[w_droit], X$time[w_droit], X$X[w_droit,feuille_split$variable]))
            meanFg <- as.matrix(kmlShape::meanFrechet(trajG))
            meanFd <- as.matrix(kmlShape::meanFrechet(trajD))
          }

          if (X$type=="factor"){
            meanFg <- unique(X$X[w_gauche, feuille_split$variable])
            meanFd <- unique(X$X[w_droit,feuille_split$variable])
          }

          if (X$type=="scalar"){
            meanFg <- mean(X$X[w_gauche,feuille_split$variable])
            meanFd <- mean(X$X[w_droit,feuille_split$variable])
          }

          hist_nodes[[2*(unique(id_feuille)[i])]] <- meanFg
          hist_nodes[[2*(unique(id_feuille)[i])+1]] <- meanFd

          count_split <- count_split+1

          # QUelles sont les feuilles courantes
          feuilles_courantes <- unique(id_feuille_prime)
          info_feuilles <- hist_imp_nodes[is.element(hist_imp_nodes[,1], feuilles_courantes),]
          impurete <- c(impurete, sum(info_feuilles[,2]*info_feuilles[,3]/length(unique(Y$id))))
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
          Y_pred[[q]] <- kmlShape::meanFrechet(data.frame(Y$id[w], Y$time[w], Y$Y[w]))
        }
        if(Y$type=="scalar"){
          Y_pred[[q]]<- mean(Y$Y[w])
        }
        if(Y$type=="factor"){
          Table <- which.max(table(Y$Y[w]))
          Y_pred[[q]] <-  as.factor(attributes(Table)$names)
        }

        if (Y$type=="shape"){
          Y_pred[[q]] <- mshape(Y$Y[,,w, drop=FALSE])
        }

        if (Y$type=="image"){
          Y_pred[[q]] <- apply(Y$Y[w,, drop=FALSE], 2, "mean")
        }

      }
      if (Y$type=="factor"){
        Ylevels <- unique(Y$Y)
        return(list(feuilles = id_feuille,  V_split=V_split, impurity=impurete, hist_nodes=hist_nodes, Y_pred = Y_pred, Y=Y, hist_imp_nodes=hist_imp_nodes, Alpha =0, Ylevels=Ylevels))
      }
      return(list(feuilles = id_feuille,  V_split=V_split, impurity=impurete, hist_nodes=hist_nodes, Y_pred = Y_pred, Y=Y, hist_imp_nodes=hist_imp_nodes, Alpha =0))
    }
  }

  V_split <- data.frame(V_split)
  names(V_split) <- c("type","num_noeud", "var_split")
  for (q in unique(id_feuille)){
    w <- which(id_feuille == q)
    if (Y$type=="curve"){
      Y_pred[[q]] <- kmlShape::meanFrechet(data.frame(Y$id[w], Y$time[w], Y$Y[w]))
    }
    if(Y$type=="scalar"){
      Y_pred[[q]]<- mean(Y$Y[w])
    }
    if(Y$type=="factor"){
      Table <- which.max(table(Y$Y[w]))
      Y_pred[[q]] <-  as.factor(attributes(Table)$names)
    }
    if (Y$type=="shape"){
      Y_pred[[q]] <- mshape(Y$Y[,,w, drop=FALSE])
    }

    if (Y$type=="image"){
      Y_pred[[q]] <- apply(Y$Y[w_Y,, drop=FALSE], 2, "mean")
    }
  }
  if (Y$type=="factor"){
    Ylevels <- unique(Y$Y)
    return(list(feuilles = id_feuille, V_split=V_split, impurity=impurete, hist_nodes=hist_nodes, Y_pred = Y_pred, Y=Y, hist_imp_nodes=hist_imp_nodes, Alpha =0, Ylevels=Ylevels))
  }
  return(list(feuilles = id_feuille, V_split=V_split, impuity=impurete, hist_nodes=hist_nodes, Y_pred= Y_pred, Y=Y, hist_imp_nodes=hist_imp_nodes, Alpha=0))
}

