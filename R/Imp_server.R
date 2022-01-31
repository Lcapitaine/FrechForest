#' Optimized Importance variable calculation for servers
#'
#' @param Trees [string]:
#' @param Curve [list]:
#' @param Scalar [list]:
#' @param Factor [list]:
#' @param Shape [list]:
#' @param Image [list]:
#' @param Y [list]:
#' @param type [string]:
#' @param range [vector]:
#' @param ncores [numeric]:
#' @param timeScale [numeric]:
#' @param xerror [vector]
#'
#' @keywords internal
#'
Importance_server <- function(Trees,Curve=NULL,Scalar=NULL, Factor=NULL, Shape=NULL,
                              Image=NULL ,Y ,type=NULL, range=NULL,ncores=NULL, timeScale=0.1, xerror){

  if(is.null(ncores)==TRUE){
    ncores <- detectCores()
  }

  trees = list.files(Trees)
  ntree = length(trees)

  p=1

  if (type=="curve"){
    imp = NULL
    Curve.err <- matrix(NA, ntree, length(range))
    Curve.perm <- Curve

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    imp <- foreach::foreach(p = 1:length(range),.packages = "kmlShape" ,.combine = "rbind") %dopar% {

      for (k in 1:ntree){

        tree <- get(load(trees[k]))
        BOOT <- tree$boot
        nboot <- length(unique(Y$id))- length(BOOT)

        id_boot_Curve <- NULL
        for (i in 1:length(BOOT)){
          id_boot_Curve <- c(id_boot_Curve, which(Curve$id==BOOT[i]))
        }

        #Permutation time

        Curve.perm$X[-id_boot_Curve,range[p]] <- permutation_courbes(Curve$X[-id_boot_Curve,range[p]], Curve$id[-id_boot_Curve])


        Curve.err[k,p] <- OOB.tree(tree, Curve=Curve.perm, Scalar = Scalar, Factor=Factor,Shape=Shape, Image=Image, Y, timeScale=timeScale)

      }
      Curve.perm$X[,range[p]] <- Curve$X[,range[p]]
      res <- data.frame(vari = range[p],imp = mean(Curve.err[,p]- xerror))
    }
    parallel::stopCluster(cl)
  }

  if (type=="scalar"){

    imp=NULL
    Scalar.err <- matrix(NA, length(trees), length(range))
    Scalar.perm <- Scalar

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    imp <- foreach::foreach(p =1:length(range),.packages = "kmlShape" ,.combine = "rbind") %dopar% {

      for (k in 1:ntree){

        tree <- get(load(trees[k]))

        BOOT <- tree$boot
        nboot <- length(unique(Y$id))- length(BOOT)

        id_boot_Scalar <- NULL
        for (i in 1:length(BOOT)){
          id_boot_Scalar <- c(id_boot_Scalar, which(Scalar$id==BOOT[i]))
        }


        Scalar.perm$X[-id_boot_Scalar,range[p]] <- sample(Scalar.perm$X[-id_boot_Scalar,range[p]])

        Scalar.err[k,p] <- OOB.tree(tree, Curve=Curve, Scalar = Scalar.perm, Factor=Factor,Shape=Shape, Image=Image, Y, timeScale=timeScale)

      }
      Scalar.perm$X[,range[p]] <- Scalar$X[,range[p]]
      res <- data.frame(vari=range[p],imp = mean(Scalar.err[,p]- xerror))
    }

    parallel::stopCluster(cl)
  }
  return(imp)
}

