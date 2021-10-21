#' Title
#'
#' @param Courbes
#' @param id
#'
#' @keywords internal
permutation_courbes <- function(Courbes,id){
  perm <- sample(unique(id), length(unique(id))) #### on change les identifiants ::
  new <- NULL
  for (i in perm){
    w <- which(id==i)
    new <- c(new,Courbes[w])
  }
  return(new)
}

#' Title
#'
#' @param Shapes
#' @param id
#'
#' @keywords internal
permutation_shapes <- function(Shapes, id){
  perm <- sample(id,length(id))
  new <- array(0,dim=dim(Shapes)[1:3])
  for (i in 1:length(id)){
    new[,,i] <- Shapes[,,which(id==perm[i])]
  }
  return(new)
}
