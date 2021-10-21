#' Extremely randomized split
#'
#' @param X
#' @param Y
#' @param timeScale
#' @param ntry
#'
#' @import kmlShape
#' @import Evomorph
#'
#' @keywords internal
ERvar_split <- function(X ,Y,ntry=3,timeScale=0.1){

  impur <- rep(0,dim(X$X)[length(dim(X$X))])
  toutes_imp <- list()
  split <- list()
  Pure <- FALSE

  for (i in 1:dim(X$X)[length(dim(X$X))]){

    if (X$type=="factor"){

      if (length(unique(X$X[,i]))>1){
        L <- Fact.partitions(X$X[,i],X$id)
        split_courant <- list()
        impur_courant <- rep(NA,length(L))
        toutes_imp_courant <- list()

        # On tire une partition au hasard
        tirage <- sample(1:length(L), 1)
        # Il faut maintenant regarder quelles sont les meilleures combinaisons ::

        split[[i]] <- rep(2,length(X$id))
        for (l in L[[tirage]]){
          split[[i]][which(X$id==l)] <- 1
        }
        # Il faut maintenant regarder la qualite du decoupage ::
        impurete <- impurity_split(Y,split[[i]])
        impur[i] <- impurete$impur
        toutes_imp[[i]] <- impurete$imp_list
      }
      else {
        impur[i] <- Inf
        split[[i]] <- Inf
      }
    }

    if( X$type=="curve"){

      # Il faut commencer par tirer les multiples centres ::

      id_centers <- matrix(NA,ntry,2)
      for (l in 1:ntry){
        id_centers[l,] <- sample(unique(X$id),2)
      }

      ### Il faut ensuite boucler sur le ntry
      split_prime <- matrix(2,ntry,length(unique(X$id)))
      u <- 0
      impurete2 <- list()
      qui <- NULL
      imp <- NULL

      for (c in 1:ntry){

        w_gauche <- which(X$id==id_centers[c,1])
        w_droit <- which(X$id==id_centers[c,2])

        for (l in 1:length(unique(X$id))){

          w <- which(X$id==unique(X$id)[l])
          dg <- distFrechet(X$time[w_gauche],X$X[w_gauche,i],X$time[w],X$X[w,i], timeScale = timeScale)
          dd <- distFrechet(X$time[w_droit],X$X[w_droit,i],X$time[w],X$X[w,i], timeScale = timeScale)
          if (dg<=dd) split_prime[c,l] <- 1
        }

        if (length(unique(split_prime[c,]))>1){
          u <- u+1
          qui <- c(qui, c)
          impurete2[[c]] <- impurity_split(Y,split_prime[c,], timeScale)
          imp <- c(imp,impurete2[[c]]$impur)
        }

      }

      if (u>0){
        gagnant <- qui[which.min(imp)]
        split[[i]] <- split_prime[gagnant,]
        impurete <- impurete2[[gagnant]]
        impur[i] <- impurete$impur
        toutes_imp[[i]] <- impurete$imp_list
      }

      else{
        impur[i] <- Inf
        split[[i]] <- Inf
      }

    }

    if (X$type=="shape"){
      if (dim(X$X[,,,i])[3]>2){

        id_centers <- matrix(NA,ntry,2)
        for (l in 1:ntry){
          id_centers[l,] <- sample(X$id,2)
        }

        split_prime <- matrix(2,ntry,length(X$id))

        u <- 0
        qui <- NULL
        impurete2 <- list()
        imp <- NULL

        for (c in 1:ntry){
          dg <- ShapeDist(X$X[,,,i],X$X[,,which(X$id==id_centers[c,1]),i])
          dd <- ShapeDist(X$X[,,,i],X$X[,,which(X$id==id_centers[c,2]),i])
          for (l in 1:length(unique(X$id))){
            if (dg[l]<=dd[l]) split_prime[c,l] <- 1
          }
          if (length(split_prime[c,])>1){
            u <- u+1
            qui <- c(qui,c)
            impurete2[[c]] <- impurity_split(Y,split_prime[c,], timeScale)
            imp <- c(imp, impurete2[[c]]$impur)
          }
        }
        ## Il nous faut calculer les distances à gauche et à droite pour chaque element

        if (u>0){
          gagnant <- qui[which.min(imp)]
          split[[i]] <- split_prime[gagnant,]

          impurete <- impurete2[[gagnant]]
          impur[i] <- impurete$impur
          toutes_imp[[i]] <- impurete$imp_list
        }

        else{
          impur[i] <- Inf
          split[[i]] <- Inf
        }

      }

      else{
        split[[i]] <- c(1,2)
        impurete <- impurity_split(Y,split[[i]], timeScale)
        impur[i] <- impurete$impur
        toutes_imp[[i]] <- impurete$imp_list
      }
    }

    if (X$type=="image"){
      if (nrow(X$X)>2){
        id_centers <- matrix(NA,ntry,2)
        for (l in 1:ntry){
          id_centers[l,] <- sample(X$id,2)
        }

        split_prime <- matrix(2,ntry,length(X$id))


        u <- 0
        qui <- NULL
        impurete2 <- list()
        imp <- NULL

        for (c in 1:ntry){

          w_g <- which(X$id==id_centers[c,1])
          w_d <- which(X$id==id_centers[c,2])
          ### Il nous faut calculer la distance :
          dg = apply(apply(X$X[,,i],1,"-",X$X[w_g,,i])^2,2,"mean")
          dd = apply(apply(X$X[,,i],1,"-",X$X[w_d,,i])^2,2,"mean")

          split_prime[c,which((dg<=dd)==TRUE)]=1
          if (length(unique(split_prime[c,]))>1){
            u <-u+1
            qui <- c(qui,c)
            impurete2[[c]] <- impurity_split(Y,split_prime[c,], timeScale)
            imp <- c(imp,impurete2[[c]]$impur)
          }

        }



        if (u>0){
          gagnant <- qui[which.min(imp)]
          split[[i]] <- split_prime[gagnant,]
          impurete <- impurete2[[gagnant]]
          impur[i] <- impurete$impur
          toutes_imp[[i]] <- impurete$imp_list
        }

        else{
          impur[i] <- Inf
          split[[i]] <- Inf
        }

      }

      else{
        split[[i]] <- c(1,2)
        impurete <- impurity_split(Y,split[[i]], timeScale)
        impur[i] <- impurete$impur
        toutes_imp[[i]] <- impurete$imp_list
      }
    }

    if(X$type=="scalar"){
      if (length(unique(X$X[,i]))>2){

        ### On doit tier les centres
        #centers <- sample(X$X[,i],2)

        centers <- matrix(NA,ntry,2)
        for (l in 1:ntry){
          centers[l,] <- sample(X$X[,i],2)
        }

        #split[[i]] <- rep(2,length(X$X[,i]))
        split_prime <- matrix(2,ntry,length(X$X[,i]))

        for (l in 1:length(X$X[,i])){
          for (k in 1:ntry){
            if (abs(centers[k,1]-X$X[l,i])<= abs(centers[k,2]-X$X[l,i])) split_prime[k,l] <- 1
          }
        }

        u <- 0
        qui <- NULL
        impurete2 <- list()
        imp <- NULL
        for (k in 1:ntry){
          if (length(unique(split_prime[k,]))>1){
            u <- u+1
            qui <- c(qui,k)
            impurete2[[k]] <- c(impurete2,impurity_split(Y,split_prime[k,], timeScale))
            imp <- c(imp, impurete2[[k]]$impur)
          }
        }

        if (u>0){
          gagnant <- qui[which.min(imp)]
          split[[i]] <- split_prime[gagnant,]
          impurete <- impurete2[[gagnant]]
          impur[i] <- impurete$impur
          toutes_imp[[i]] <- impurete$imp_list
        }

        else{
          impur[i] <- Inf
          split[[i]] <- Inf
        }
      }

      else {
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

