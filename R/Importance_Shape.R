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
Importance_Shape <- function(Curve=NULL,Scalar=NULL, Factor=NULL, Shape=NULL,
                              Image=NULL ,Y, range=NULL,ncores=NULL, timeScale=0.1){

  if(is.null(ncores)==TRUE){
    ncores <- detectCores()
  }

  trees = list.files()
  ntree = length(trees)

  imp = rep(NA,length(range))

  Shape.err <- matrix(NA, ntree, length(range))


  for (p in 1:length(range)){

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    k=1

    Shape.err <- foreach::foreach(k = 1:ntree,.packages = "kmlShape" ,.combine = "cbind") %dopar% {

      tree <- get(load(trees[k]))
      BOOT <- tree$boot
      nboot <- length(unique(Y$id))- length(BOOT)

      id_boot_Shape <- NULL
      for (i in 1:length(BOOT)){
        id_boot_Shape <- c(id_boot_Shape, which(Shape$id==BOOT[i]))
      }

      Shape.perm <- Shape

      Shape.perm$X[,,-id_boot_Shape,range[p]] <- permutation_shapes(Shape.perm$X[,,-id_boot_Shape, range[p]], Shape.perm$id[-id_boot_Shape])


      res <- OOB.tree(tree, Curve=Curve, Scalar = Scalar, Factor=Factor,Shape=Shape.perm, Image=Image, Y, timeScale=timeScale)-
        OOB.tree(tree, Curve=Curve, Scalar = Scalar, Factor=Factor,Shape=Shape, Image=Image, Y, timeScale=timeScale)
    }
    parallel::stopCluster(cl)
    imp[p] <- mean(Shape.err)
  }

  return(imp)
}
