#' @title Twin/Untwin
#'
#' @description Put twin datasets into one object (a list of lists),
#'     or else combine a twin object into a merged list.
#'
#' @param spe1 dataset 1 species matrix
#' @param spe2 dataset 2 species matrix
#' @param id1 dataset 1 environmental/identifier matrix
#' @param id2 dataset 2 environmental/identifier matrix
#' @param x twin object
#' @param ... additional arguments passed to function
#'
#' @return
#' \code{twin} returns a list of two named items: `z1` and `z2`, each
#'     of which is a list with two named items: `spe` and `id` \cr
#' \code{untwin} returns a list with two named items: `spe` and `id`
#'
#' @examples
#' data(smoky)
#' spe <- smoky$spe
#' env <- smoky$env
#' tw  <- twin(spe, spe, env, env)
#' utw <- untwin(tw)
#'
#' # must have same number of rows per dataset:
#' \dontrun{twin(spe, spe[3:9,], env, env[3:9,])}
#'
#' @seealso \code{\link{nearestspecies}} to force same number of rows
#'
#' @export
#' @rdname twin
`twin` <- function(spe1, spe2, id1, id2, ... ){
     stopifnot(identical(row.names(spe1), row.names(id1)),
               identical(row.names(spe2), row.names(id2)))
     if(!(identical(dim(spe1)[[1]], dim(spe2)[[1]]))) {
          stop('Number of rows must match')
     }
     row.names(spe2) <- paste0('.', row.names(spe2))
     row.names(id2)  <- paste0('.', row.names(id2))
     out <- list( list(spe=spe1, id=id1),
                  list(spe=spe2, id=id2) )
     names(out) <- paste0('z', 1:2) # arbitrary names for each twin
     class(out) <- c('twin')
     out
}
#' @export
#' @rdname twin
`untwin` <- function(x, ...){
     if(class(x) != 'twin') stop('Input must be class `twin`')
     d1  <- x[[1]][['spe']]
     d2  <- x[[2]][['spe']]
     id1 <- x[[1]][['id']]
     id2 <- x[[2]][['id']]
     if(!identical(dimnames(d1)[2], dimnames(d2)[2]))
          cat('\nspecies not exact match, joining data anyway\n\n')
     spe <- ecole::mx_rbind_all(d1, d2)
     id  <- ecole::mx_rbind_all(id1, id2)
     list(spe, id)
}
