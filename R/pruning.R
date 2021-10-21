#' General pruning function
#'
#' @param tree
#'
#'
#' @keywords internal
elagage <- function(tree){

  t_feuilles <- NULL
  t_id <- NULL
  t_time <- NULL
  t_split <- NULL
  t_hist <- NULL
  t_Y <- tree$Y


  tree_courant <- tree
  nb_feuilles <- length(unique(tree$feuilles))
  n_max <- nb_feuilles
  courant <- 2
  TREES <- list()
  TREES[[1]] <- tree
  ##### il faut aussi trouver les d?coupe superficielles :::: on garde un historique des d?coupes :::::
  while(nb_feuilles >1){
    deg <- noeuds_deg(tree_courant)
    if (dim(deg)[1]>1) deg <- apply(deg, 2, sort, decreasing=TRUE)
    t_feuilles_courant <- tree_courant$feuilles
    for (t in deg[,1]){
      b <- branche(tree_courant, t)
      feuilles_b <- unique(b$feuilles)
      w_feuilles <- NULL
      for (f in feuilles_b ){
        w_feuilles <- c(w_feuilles, which(tree_courant$feuilles==f))
      }
      t_feuilles_courant[w_feuilles] <- t
      #### il faut maintenant retirer toute la branche de t

      nodes <- as.numeric(as.character(b$V_split[,2]))
      w_nodes <- NULL
      for (node in nodes){
        w_nodes <- c(w_nodes, which(tree_courant$V_split[,2]==node))
      }

      t_split_courant <- tree_courant$V_split[-w_nodes,, drop = FALSE]
      ##### il faut alors recalculer l'importance dans les feuilles
      tree_courant <- list(feuilles=t_feuilles_courant, V_split = t_split_courant,hist_nodes=tree$hist_nodes, Y=tree$Y, hist_imp_nodes=tree$hist_imp_nodes, Y_pred = tree$Y_pred, Alpha=unique(deg[,2]))
      TREES[[courant]] <- tree_courant
    }
    courant <- courant+1
    nb_feuilles <- length(unique(tree_courant$feuilles))
  }

  return(TREES)
}
