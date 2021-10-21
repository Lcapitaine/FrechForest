#' General Frechet Tree
#'
#' @param Curve [list]:
#' @param Scalar [list]:
#' @param Factor [list]:
#' @param Y [list]:
#' @param timeScale [numeric]:
#' @param ncores [numeric]:
#'
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import kmlShape
#' @import Evomorph
#' @import pbapply
#'
#' @export
#'
FrechetTree <- function(Curve=NULL,Scalar=NULL,Factor=NULL,Y,timeScale=0.1, ncores=NULL){

  ### Il faut normaliser les elements des formes ::

  if (is.null(ncores)==TRUE) ncores <- detectCores()-1

  if (Y$type=="shape"){
    Y$Y <- gpagen(Y$Y,print.progress = FALSE)$coords
  }

  TMAX <- Tmax(Curve=Curve,Scalar = Scalar,Factor=Factor, Shape=NULL,Image=NULL,Y,timeScale = timeScale)

  if (Y$type=="shape") dime <- dim(Y$Y)[1:2]


  elag_max <- elagage(TMAX)
  ALPHA <- rep(NA, length(elag_max))
  for (i in 1:length(ALPHA)){
    ALPHA[i] <- elag_max[[i]]$Alpha
  }
  #### on transforme le tout en beta
  beta <- rep(NA, length(ALPHA))
  beta[length(ALPHA)] <- ALPHA[length(ALPHA)]
  for (i in 1:(length(ALPHA)-1)){
    beta[i] <- sqrt(abs(ALPHA[i]*ALPHA[i+1]))
  }

  #### Il faut faire faire les sous ensemble de validation crois?e::::
  #ELAG <- list()
  n_folds <- 10
  VC <- sample(rep(1:n_folds, length.out= length(unique(Y$id))))
  #tmax <- list()
  #APP <- list()
  #err <- matrix(0, length(beta), n_folds)


  Scalar.app <- Scalar
  Curve.app <- Curve
  Factor.app <- Factor

  Scalar.val <- NULL
  Factor.val <- NULL
  Curve.val <- NULL

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  #p=1
  #err <- foreach(p = 1:n_folds, .packages='kmlShape', .combine="cbind") %dopar% {

  err <- pbsapply(1:n_folds, FUN=function(i){


    res <- rep(NA,length(beta))
    app <- unique(Y$id)[which(VC!=i)] ### on r?cup?re les identifiants
    w <- NULL
    wCurve <- NULL
    wScalar <- NULL
    wFactor <- NULL

    for (a in app){
      w <- c(w, which(Y$id==a))
      if (is.null(Scalar)!=TRUE) wScalar <- c(wScalar, which(Scalar$id==a))
      if (is.null(Factor)!=TRUE) wFactor <- c(wFactor, which(Factor$id==a))
      if (is.null(Curve)!=TRUE) wCurve <- c(wCurve, which(Curve$id==a))

    }
    APP <- w

    ### On prend les elements d'apprentissage maintenant :::

    if (is.null(Scalar)!=TRUE){
      Scalar.app <- list(X=Scalar$X[wScalar,,drop=FALSE], id=Scalar$id[wScalar])
      Scalar.val <- list(type="scalar",X=Scalar$X[-wScalar,,drop=FALSE], id=Scalar$id[-wScalar])
    }


    if (is.null(Factor)!=TRUE){
      Factor.app <- list(X=Factor$X[wFactor,,drop=FALSE], id=Factor$id[wFactor])
      Factor.val <- list(type="factor",X=Factor$X[-wFactor,,drop=FALSE], id=Factor$id[-wFactor])
    }

    if (is.null(Curve)!=TRUE){
      Curve.app <- list(X=Curve$X[wCurve,,drop=FALSE], id=Curve$id[wCurve], time=Curve$time[wCurve])
      Curve.val <- list(type="curve",X=Curve$X[-wCurve,,drop=FALSE], id=Curve$id[-wCurve], time=Curve$time[-wCurve])
    }

    if (Y$type=="curve"){
      Y.app <- list(type="curve",Y=Y$Y[w],id=Y$id[w],time=Y$time[w])
      Y.val <- list(type="curve",Y=Y$Y[-w],id=Y$id[-w],time=Y$time[-w])
    }

    if (Y$type=="factor" || Y$type=="scalar"){
      Y.app <- list(type=Y$type,Y=Y$Y[w],id=Y$id[w],time=Y$time[w])
      Y.val <- list(type=Y$type,Y=Y$Y[-w],id=Y$id[-w],time=Y$time[-w])
    }

    if (Y$type=="shape"){
      Y.app <- list(type=Y$type,Y=Y$Y[,,w],id=Y$id[w])
      Y.val <- list(type=Y$type,Y=Y$Y[,, -w],id=Y$id[-w])
    }

    if (Y$type=="image"){
      Y.app <- list(type=Y$type,Y=Y$Y[w,],id=Y$id[w])
      Y.val <- list(type=Y$type,Y=Y$Y[-w,],id=Y$id[-w])
    }


    tmax <- Tmax(Curve = Curve.app,Scalar = Scalar.app,Factor = Factor.app,Y=Y.app, timeScale = timeScale)

    ELAG <- elagage(tmax)

    pen <- rep(NA,length(ELAG))
    #pen[length(ELAG[[p]])] <- ELAG[[p]][[length(ELAG[[p]])]]$Alpha
    for (l in 1:length(pen)){
      pen[l] <- ELAG[[l]]$Alpha
    }

    for (k in 1:length(beta)){
      sous_arbre <- ELAG[[which.min(abs(pen-beta[k]))]]
      where <- pred.FT(sous_arbre,Curve = Curve.val,Scalar=Scalar.val,Factor = Factor.val,timeScale = timeScale) #### on doit trouver les feuilles de pr?diction :::
      ##### il nous faut maintenant pr?dire les diff?rentes courbes ::::
      err_courant <- rep(0, length(where))

      for (j in 1:length(where)){
        ww <- which(Y.val$id == unique(Y.val$id)[j])
        #mean_courant <- DouglasPeuckerEpsilon(sous_arbre$Y_curves[[where[j]]][,1],sous_arbre$Y_curves[[where[j]]][,2], 0.01)
        ## On regarde si on a bien une sortie qui est une courbe:
        if (Y$type=="curve") err_courant[j] <-  kmlShape::distFrechet(Y.val$time[ww], Y.val$Y[ww],sous_arbre$Y_pred[[where[j]]][,1], sous_arbre$Y_pred[[where[j]]][,2])
        if (Y$type=="scalar") err_courant[j] <- (Y.val$Y[ww]-where[j])^2
        if (Y$type=="factor") err_courant[j] <- 1*(Y.val$Y[ww]==where[j])
        if (Y$type=="image") err_courant[j] <- mean((Y.val[ww,]-sous_arbre$Y_pred[[where[j]]])^2)
        if (Y$type=="shape") err_courant[j] <- ShapeDist(array(sous_arbre$Y_pred[[where[j]]],dim=c(dime[1],dime[2],1)),Y.val$Y[,,ww])
      }
      res[k] <- mean(err_courant)
    }
    res
  },cl=cl)

  parallel::stopCluster(cl)

  SD <- apply(err, 1, "sd")
  err_M <- apply(err, 1, "mean")

  ### On prend le meilleur modèle:
  seuil <- min(err_M) + SD[which.min(err_M)]
  ### il faut s?lectionner le meilleur arbre
  optimal.tree <- max(which(err_M<=seuil))

  beta.opt <- beta[optimal.tree]
  final_tree <- elag_max[[optimal.tree]]

  #### On va s?lectionner l'arbre optimal pour chaque ensemble d'apprentissage puis calculer l'importance des variables sur ceux-ci::
  #Importance <- matrix(0, n_folds, dim(X)[2])
  #err_arbres_select <- rep(NA, n_folds)
  #for (k in 1:n_folds){
  ### on r?cup?re les ?l?ments de validation::::
  #X.val <- X[-APP[[k]],]
  #Y.val <- Y[-APP[[k]]]
  #time.val <- time[-APP[[k]]]
  #id.val <- id[-APP[[k]]]

  # pen <- rep(NA,length(ELAG[[k]]))
  #for (l in 1:length(pen)){
  #  pen[l] <- ELAG[[k]][[l]]$Alpha
  #}
  #}
  ## On va faire l'affichage de la selection de l'abre
  #plot(err_M)
  #lines(rep(seuil, length(err_M)),col=2)
  #points(optimal.tree,err_M[optimal.tree], col=2)

  m_leafs <- max(unique(final_tree$feuilles))
  return(list(feuilles = final_tree$feuilles, V_split=final_tree$V_split, hist_nodes=final_tree$hist_nodes[1:m_leafs], Y_pred=final_tree$Y_pred[1:m_leafs], err_elag = err_M, seuil=seuil, Y=Y))
}
