#' Optimized Importance variable calculation for servers
#'
#' @param Curve [list]:
#' @param Scalar [list]:
#' @param Factor [list]:
#' @param Shape [list]:
#' @param Image [list]:
#' @param Y [list]:
#' @param range [vector]:
#' @param ncores [numeric]:
#' @param timeScale [numeric]:
#'
#' @export
#'
Importance_Curve <- function(Curve=NULL,Scalar=NULL, Factor=NULL, Shape=NULL,
                              Image=NULL ,Y, range=NULL,ncores=NULL, timeScale=0.1){

  if(is.null(ncores)==TRUE){
    ncores <- detectCores()
  }

  trees = list.files()
  ntree = length(trees)

  imp = matrix(NA, length(range), ncol(Y$Y))

  Curve.err <- matrix(NA, ntree, length(range))


  for (p in 1:length(range)){

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    k=1

    Curve.err <- foreach::foreach(k = 1:ntree,.packages = "kmlShape" ,.combine = "cbind") %dopar% {

      tree <- get(load(trees[k]))
      BOOT <- tree$boot
      nboot <- length(unique(Y$id))- length(BOOT)

      id_boot_Curve <- NULL
      for (i in 1:length(BOOT)){
        id_boot_Curve <- c(id_boot_Curve, which(Curve$id==BOOT[i]))
      }

      Curve.perm <- Curve

      Curve.perm$X[-id_boot_Curve,range[p]] <- permutation_courbes(Curve$X[-id_boot_Curve,range[p]], Curve$id[-id_boot_Curve])


      res <- OOB.tree(tree, Curve=Curve.perm, Scalar = Scalar, Factor=Factor,Shape=Shape, Image=Image, Y, timeScale=timeScale)-
        OOB.tree(tree, Curve=Curve, Scalar = Scalar, Factor=Factor,Shape=Shape, Image=Image, Y, timeScale=timeScale)
    }
    parallel::stopCluster(cl)
    imp[p,] <- apply(Curve.err, 2, "mean")
  }

  return(imp)
}
