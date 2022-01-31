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
Importance_Factor <- function(Curve=NULL,Scalar=NULL, Factor=NULL, Shape=NULL,
                              Image=NULL ,Y, range=NULL,ncores=NULL, timeScale=0.1){

  if(is.null(ncores)==TRUE){
    ncores <- detectCores()
  }

  trees = list.files()
  ntree = length(trees)

  imp = rep(NA,length(range))

  Factor.err <- matrix(NA, ntree, length(range))


  for (p in 1:length(range)){

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    k=1

    Factor.err <- foreach::foreach(k = 1:ntree,.packages = "kmlShape" ,.combine = "cbind") %dopar% {

      tree <- get(load(trees[k]))
      BOOT <- tree$boot
      nboot <- length(unique(Y$id))- length(BOOT)

      id_boot_Factor <- NULL
      for (i in 1:length(BOOT)){
        id_boot_Factor <- c(id_boot_Factor, which(Factor$id==BOOT[i]))
      }

      Factor.perm <- Factor

      Factor.perm$X[-id_boot_Factor,range[p]] <- sample(Factor.perm$X[-id_boot_Factor,range[p]])


      res <- OOB.tree(tree, Curve=Curve, Scalar = Scalar, Factor=Factor.perm,Shape=Shape, Image=Image, Y, timeScale=timeScale)-
        OOB.tree(tree, Curve=Curve, Scalar = Scalar, Factor=Factor,Shape=Shape, Image=Image, Y, timeScale=timeScale)
    }
    parallel::stopCluster(cl)
    imp[p] <- mean(Factor.err)
  }

  return(imp)
}
