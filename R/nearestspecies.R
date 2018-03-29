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
#' @param ties logical, should ties (non-unique values) be allowed?
#'     default is \code{FALSE}
#'
#' @param ... additional arguments passed to function
#'
#' @return
#' Numeric vector indexing rows (sample units) in \code{spe1} which
#'     are nearest compositional neighbors of \code{spe2}. The order
#'     corresponds to row order in \code{spe2}.
#'
#' @details
#' Useful to subset the longer of two datasets for analyses that
#'     require compatible dimensions (e.g.,
#'     \code{\link[vegan]{mantel}}, \code{\link[vegan]{procrustes}}).
#'
#'     If \code{ties = FALSE}, then all returned values are unique;
#'     if \code{ties = TRUE}, then unique values are not enforced, so
#'     any sample unit in \code{spe1} may be picked repeatedly.
#'
#' @examples
#' # two candidate datasets, full and partial
#' set.seed(181)
#' data(smoky)
#' full <- smoky$spe
#' part <- (full + abs(rnorm(prod(dim(full)), 0, 25)))[1:6,]
#'
#' # subset the full matrix, based on compositional nearest neighbors
#' i <- nearestspecies(full, part, ties=FALSE)
#' i
#' full <- full[i,,]
#'
#' # default is to break ties
#' nearestspecies(full, part, ties=FALSE)
#' nearestspecies(full, part, ties=TRUE)
#'
#' @export
#' @rdname nearestspecies
`nearestspecies` <- function(spe1, spe2, method, ties = FALSE, ...){
     nr1 <- nrow(spe1)
     nr2 <- nrow(spe2)
     if (nr1 < nr2) {
          stop("spe1 requires >= rows than spe2, try reverse order")
     }
     both <- ecole::mx_rbind_all(spe1, spe2)
     if (missing(method)) {
          D <- vegan::vegdist(both, "bray")
          # need capture.output to silence unwanted message printing:
          dontprint <- capture.output(
               D <- vegan::stepacross(D, path = "shortest",
                                      toolong = 1), file=NULL
          )
     } else {
          D <- vegan::vegdist(both, method = method)
     }
     i <- c(rep("1", len = nr1), rep("2", len = nr2))
     D <- as.matrix(D)[i == "1", i == "2"]
     dimnames(D)[[2]] <- paste0("_", dimnames(D)[[2]])
     # o matrix gives indices of query SUs "most similar" to each
     #     target SU; ties may exist (same query SU may be picked >1x)
     o <- apply(D, 2, order)
     if ( !ties ) {
          ### disallow any ordering ties, keep most similar one
          `f` <- function(x){
               which(duplicated(v,MAR=0) |
                          duplicated(v,MAR=0,fromLast=T))
          }
          nro <- nrow(o)
          # iterate down rows of the ordering matrix 'o'
          for(i in 1:nro){
               # set NA if already in precedent row(s) above
               if( i > 1 ){
                    isabove <- which(na.exclude(o[i,])%in%o[1:(i-1),])
                    setNA      <- which(!is.na(o[i,]))[isabove]
                    o[i,setNA] <- NA
               }
               Dt    <- D[1:i,]         # temporary dissim matrix
               v     <- o[1:i,]         # temporary ordering matrix
               isna  <- which(is.na(v)) # find NA
               ties  <- f(v)            # o index that are ties
               ties  <- ties[!ties %in% isna] # dont compare NA
               utie  <- unique(v[ties]) # unique values of ties
               # handle multiple unique ties per o row:
               for( u in 1:length(utie) ){
                    ftie  <- ties[(v[ties]==utie[u])] # focal ties
                    minD  <- which.min(Dt[ftie])  # minD of focal ties
                    ftie  <- ftie[-minD]          # update, rm minD
                    o[1:i,][ftie] <- NA  # larger focal ties set NA
               }
               # NA fill subordinate columns to keep:
               colkeep <- !is.na(o[1:i,])
               if( is.vector(colkeep) ) colkeep <- t(colkeep)
               colkeep <- which(colSums(colkeep*1)>0)
               if(i+1 > nro) break
               o[(i+1):nro, colkeep] <- NA
               if( all(is.na(o[i+1,])) )  break
          }
          keep <- apply(o, 2, function(x) x[!is.na(x)] ) # NO ties
     } else {
          ### simply allow any ordering ties
          keep <- o[1, ]  # may contain ties
     }
     keep
}
