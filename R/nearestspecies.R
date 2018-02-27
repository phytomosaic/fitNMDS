#' @title Find nearest species compositional neighbors among two
#'     datasets
#
#' @description Given two species abundance matrices (of possibly
#'     unequal size), find nearest species compositional neighbors.
#'
#' @param spe1 species abundance matrix with greater number of sites
#'     (rows)
#'
#' @param spe2 species abundance matrix with lesser number of sites
#'     (rows)
#'
#' @param method dissimilarity index, per
#'     \code{\link[vegan]{vegdist}}; default uses Bray-Curtis
#'     dissimilarities with stepacross adjustment
#'
#' @param uniq logical, force unique values? (currently unimplemented)
#'
#' @param ... additional arguments passed to function
#'
#' @return
#' Numeric vector indexing sites in \code{spe1} which are nearest
#'     compositional neighbors of \code{spe2}.
#'
#' @details
#' Useful to subset the longer of two datasets for analyses that
#'     require compatible dimensions (e.g.,
#'     \code{\link[vegan]{mantel}}, \code{\link[vegan]{procrustes}}).
#'     Currently, any sample unit in \code{spe1} may be picked
#'     repeatedly (unique values not enforced).
#'
#' @examples
#' # two candidate datasets, full and partial
#' data(smoky)
#' full <- smoky$spe
#' part <- full[3:9,]
#'
#' # subset the full matrix, based on compositional nearest neighbors
#' keep <- nearestspecies(full, part)
#' keep
#' full <- full[keep,,]
#'
#' @export
#' @rdname nearestspecies
`nearestspecies` <- function(spe1, spe2, method, uniq=FALSE, ...){
     nr1 <- nrow(spe1)
     nr2 <- nrow(spe2)
     if(nr1 < nr2){
          stop('spe1 requires >= rows than spe2, try reverse order')
     }
     # calculate combined dissimilarity matrix
     both <- ecole::mx_rbind_all(spe1, spe2)
     if(missing(method)){
          D <- vegan::vegdist(both, 'bray')
          D <- vegan::stepacross(D, path='shortest', toolong=1)
     } else {
          D <- vegan::vegdist(both, method=method)
     }
     # index for each dataset
     i  <- c(rep('1',len=nr1),rep('2',len=nr2))
     # subset dissimilarity matrix where rows = spe1, cols = spe2
     D <- as.matrix(D)[i=='1',i=='2']
     dimnames(D)[[2]] <- paste0('_', dimnames(D)[[2]])

     ### TODO: handle duplicates
     if(uniq){
          ### TODO
     } else {
          # return index of lowest-ranked match
          #     (same SU may be picked >1x)
          keep <- apply(D, 2, order)[1,]
     }
     keep
}
