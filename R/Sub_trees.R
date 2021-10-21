#' Sub trees  extractor
#'
#' @param tree
#' @param t
#'
#'
#' @keywords internal
branche <- function(tree, t){
  Y <- list()
  f <- unique(tree$feuilles)
  sous_split <- tree$V_split[which(tree$V_split[,2]==t),]
  N <- 2
  g <- which(tree$V_split[,2]==2*t)
  d <- which(tree$V_split[,2]==2*t+1)
  noeuds_courants <- as.numeric(as.character(tree$V_split[c(g,d),2]))
  noeuds_courants1 <- noeuds_courants
  sous_split <- rbind(sous_split, tree$V_split[c(g,d),])
  sous_feuilles <- NULL
  hist_nodes <- list()
  if (length(g)>0) {hist_nodes[[2*t]] <- tree$hist_nodes[[2*t]]}
  if (length(d)>0) {hist_nodes[[2*t+1]] <- tree$hist_nodes[[2*t+1]]}
  if (length(d)== 0) {sous_feuilles <- c(sous_feuilles, 2*t+1)
  Y[[2*t+1]] <- tree$Y_pred[[2*t+1]]}
  if (length(g)== 0) {sous_feuilles <- c(sous_feuilles, 2*t)
  Y[[2*t]] <- tree$Y_pred[[2*t]]}
  racine <- t
  if (length(noeuds_courants)>0) {
    while(N>0){
      p <- 0
      courant_prime <- NULL
      for (l in noeuds_courants){
        g <- which(tree$V_split[,2]==2*l)
        d <- which(tree$V_split[,2]==2*l+1)

        if (length(g)>0){ p <- p+2
        courant_prime <- c(courant_prime, as.numeric(as.character(tree$V_split[g,2])))
        sous_split <- rbind(sous_split, tree$V_split[g,])
        hist_nodes[[2*l]] <- tree$hist_nodes[[2*l]]}

        if (length(d)>0){ p <- p+2
        courant_prime <- c(courant_prime, as.numeric(as.character(tree$V_split[d,2])))
        sous_split <- rbind(sous_split, tree$V_split[d,])
        hist_nodes[[2*l+1]] <- tree$hist_nodes[[2*l+1]]}

        if(length(g)==0) {sous_feuilles <- c(sous_feuilles,2*l)
        Y[[2*l]] <- tree$Y_pred[[2*l]]}

        if (length(d)==0) { sous_feuilles <- c(sous_feuilles, 2*l+1)
        Y[[2*l+1]] <- tree$Y_pred[[2*l+1]]}
      }
      noeuds_courants <- courant_prime
      N <-p
    }
  }

  if (length(noeuds_courants1)==0) {sous_feuilles <- c(2*t, 2*t+1)}

  ## C'est maintenant que ca devient coton :::
  # Il faut recuperer les id des gens qui sont

  s_feuilles <- NULL
  s_id <- NULL
  s_time <- NULL
  s_Y <- NULL

  for(f in unique(sous_feuilles)){
    w <- which(tree$feuilles==f)
    s_feuilles <- c(s_feuilles, tree$feuilles[w])
    s_id <- c(s_id, tree$Y$id[w])
    if (tree$Y$type=="curve"){
      s_time <- c(s_time,tree$Y$time[w])
    }
    #s_time <- c(s_time, tree$time[w])
    if (tree$Y$type=="shape" || tree$Y$type=="image") s_Y <- c(s_Y,w)
    else s_Y <- c(s_Y, tree$Y$Y[w])
  }
  if (tree$Y$type=="shape" || tree$Y$type=="image") s_Y <- tree$Y$Y[,,s_Y,drop=FALSE]
  #### il faut maintenant calculer l'impurete de la branche ainsi que celle du noeud t
  #### impurete dans le noeud racine :::
  impurity_racine <- tree$hist_imp_nodes[which(tree$hist_imp_nodes[,1]==racine),2]
  n_racine <- tree$hist_imp_nodes[which(tree$hist_imp_nodes[,1]==racine),3]
  n_base <- tree$hist_imp_nodes[1,3]
  impurity_racine <- impurity_racine*(n_racine/n_base)

  impurity_T <- 0
  for (i in unique(s_feuilles)){
    w <- which(tree$hist_imp_nodes[,1]==i)
    prop <- tree$hist_imp_nodes[w,3]/n_base
    impurity_T <- impurity_T + tree$hist_imp_nodes[w,2]*prop
  }
  if (tree$Y$type=="curve"){
    sous_Y <- list(type=tree$Y$type, Y=s_Y, id = s_id, time=s_time)
  }
  else sous_Y <- list(type=tree$Y$type, Y=s_Y, id = s_id)
  return(list(feuilles=s_feuilles, V_split = sous_split, hist_nodes=hist_nodes, Y=sous_Y, impurity_T = impurity_T, impurity_racine = impurity_racine, n_racine=n_racine, Y_pred=Y))
}


#' Detect and destroy nodes
#'
#' @param tree
#'
#'
#' @keywords internal
noeuds_deg <- function(tree){
  noeuds <- as.numeric(as.character(tree$V_split$num_noeud))
  deg <- NULL
  alpha <- rep()
  mat_pen <- matrix(0, length(noeuds), 5)
  mat_pen[,1] <- noeuds
  for (t in noeuds){
    b <- branche(tree,t) ### on recupère la branche associee à t
    if (length(unique(b$feuilles))>1){
      mat_pen[which(noeuds==t), 2] <- b$impurity_racine
      mat_pen[which(noeuds==t), 3] <- b$impurity_T
      mat_pen[which(noeuds==t), 4] <- length(unique(b$feuilles))
      mat_pen[which(noeuds==t), 5] <- (b$impurity_racine-b$impurity_T)/(length(unique(b$feuilles))-1)}
    #pen <- mat_pen[which(noeuds==t), 5]
    #err <- b$impurity_T + pen*length(unique(b$feuilles)) - b$impurity_racine - pen
    #print(err)
  }
  alpha <- min(mat_pen[,5])
  err <- rep(0, length(noeuds))
  for (i in  1:dim(mat_pen)[1]){
    err[i] <- round(mat_pen[i,3] + alpha*mat_pen[i,4] - mat_pen[i,2] - alpha, 5)
    if (err[i]==0){
      deg <- rbind(deg, c(mat_pen[i,1], alpha))
    }
  }
  return(deg)
}
