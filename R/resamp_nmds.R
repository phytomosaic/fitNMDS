#' @title Resampled NMDS
#'
#' @description Bootstrap or jackknife resampling of nonmetric
#'     multidimensional scaling (NMDS) configurations, to estimate
#'     internal sampling variability. Accompanied by summary and plot
#'     methods.
#'
#' @param spe species abundance matrix
#'
#' @param method dissimilarity index, per \code{\link[vegan]{vegdist}}
#'
#' @param zeroadj default \code{TRUE}: adjust dissimilarities to
#'      account for zero-sum rows?
#'
#' @param step default \code{TRUE}: use shortest-path stepacross
#'      dissimilarity adjustment?
#'
#' @param k number of dimensions sought in final NMDS solution
#'
#' @param B number of bootstrap replicates
#'
#' @param BS bootstrap size (currently constrained to not exceed the
#'     number of rows in original matrix). Default \code{BS=NULL}, or
#'     else \code{BS=1}, will perform jackknife NMDS (leave-one-out).
#'
#' @param x,object result from \code{resamp_nmds}
#'
#' @param noaxes default \code{TRUE}: cleanly plot without axes?
#'
#' @param col character vector of colors for spider lines
#'
#' @param type plot type, one of \code{'spider'} or \code{'points'}
#'     for ordination plots, or \code{'hist'} to plot the distribution
#'     of fit statistics across resampling replicates
#'
#' @param ... additional arguments passed to function
#'
#' @return
#' List of class \sQuote{resamp_nmds} with elements:\cr
#' \describe{
#'      \item{\code{rP_draw}:}{
#'           Vector of rP_internal values for each replicate.
#'      }
#'      \item{\code{stresses}:}{
#'           Vector of NMDS stress for each replicate.
#'      }
#'      \item{\code{SRV}:}{
#'           Scaled rank variance across all replicates.
#'      }
#'      \item{\code{rP_internal}:}{
#'           Internal sampling variability across all replicates.
#'      }
#'      \item{\code{points}:}{
#'           Scores of each replicate configuration.
#'      }
#' }
#'
#' @details
#' Bootstrap or jackknife resampling of NMDS with \code{resamp_nmds}
#' estimates internal sampling variability of a candidate dataset.
#' The premise of this is resampling with replacement from the
#' dissimilarities matrix, performing NMDS ordination (Kruskal 1964)
#' on each resampled replicate, then determining collective sampling
#' variability across all NMDS solutions as \code{rP_internal}.  Other
#' uses of resampled NMDS may include estimating confidence regions
#' for site scores, or testing the \sQuote{stability} of species
#' positions along ordination axes using Scaled Rank Variance of Knox
#' and Peet (1989).
#'
#' Users may select dissimilarity index per
#' \code{\link[vegan]{vegdist}}, as well as options for zero-adjusted
#' dissimilarities (Clarke et al. 2006) using \code{zeroadj=TRUE},
#' and/or shortest-path \code{\link[vegan]{stepacross}}
#' dissimilarities using \code{step=TRUE}.
#'
#' @examples
#' set.seed(231)
#' data(smoky)
#' x <- smoky$spe
#'
#' # try jackknife and bootstrap on a simplistic dataset
#' #    (vegan::metaMDS gives warning of near-zero stress)
#' j <- resamp_nmds(x, k=2)
#' b <- resamp_nmds(x, B=29, BS=NROW(x), k=2)
#' summary(j)
#' summary(b)
#' par(mfrow=c(1,2))
#' plot(j)
#' plot(b)
#'
#' @references
#' Clarke, K. R., P. J. Somerfield, and M. G. Chapman. 2006. On
#'     resemblance measures for ecological studies, including
#'     taxonomic dissimilarities and a zero-adjusted Bray–Curtis
#'     coefficient for denuded assemblages. J Exp Marine Biol and Ecol
#'     330:55–80.\cr
#' Knox, R. G., and R. K. Peet. 1989. Bootstrapped ordination: a
#'     method for estimating sampling effects in indirect gradient
#'     analysis. Vegetatio 80:153–165.\cr
#' Kruskal, J. B. 1964. Multidimensional scaling by optimizing
#'     goodness of fit to a nonmetric hypothesis. Psychometrika 29:
#'     1-27.\cr
#'
#' @family resampled NMDS functions
#' @seealso \code{\link{recip_nmds}} for reciprocal NMDS, and
#'     \code{\link[vegan]{metaMDS}} for the core NMDS algorithm
#'
#' @export
#' @rdname resamp_nmds
### resample NMDS core function
`resamp_nmds` <- function(spe, method='bray', zeroadj=TRUE,
                            step=TRUE, k=2, BS=NULL, B=9, ...){

     ### control
     time0 <- Sys.time()
     x     <- as.matrix(spe)
     M     <- NCOL(x)
     N     <- NROW(x)
     if (is.null(BS)) BS <- 1
     if (BS > N) BS <- N

     ### dissimilarities
     D <- Dx <- dissim(x, method, zeroadj, step, ...)

     ### reference ordination
     ord  <- ordfn(D, k, x, ...)

     ### resampled ordinations
     if (BS == 1) {
          ### jackknife
          Draws <- sapply(1:N, FUN=jackfn, x, Dx, k, N, old=ord)
          NC    <- N
     } else {
          ### bootstrap
          Draws <- replicate(n=B, expr=bootfn(x, Dx, k, BS, N,
                                              old=ord, ...))
          NC    <- B
     }

     ### summarize resampled values
     scr <- do.call(rbind,lapply(1:NC,function(nc)Draws[,nc]$points))
     rnk <- Draws['rnk', ] # species ranks for SRV
     xx  <- matrix(NA, nrow=prod(dim(Draws[,1]$rnk)), ncol=NC)
     for (r in 1:length(rnk)) {xx[,r] <- c(rnk[[r]])} # unfold
     rv_obs   <- mean(rowSums((xx - rowMeans(xx))^2)/(dim(xx)[2] - 1))
     rv_exp   <- (M^2 - 1)/12        # random expectation
     SRV      <- rv_obs / rv_exp     # Scaled Rank Variance
     rP_draw  <- unlist(Draws['rP', ], use.names=F) # rP fit all draws
     rP_medn  <- median(rP_draw)                    # rP fit median
     stresses <- unlist(Draws['stress', ], use.names=F)

     ### output list
     out <- list(rP_draw     = rP_draw,
                 stresses    = stresses,
                 SRV         = SRV,
                 rP_internal = rP_medn,
                 points      = scr)
     class(out) <- 'resamp_nmds'

     ### exit timing
     difftime <- as.numeric(Sys.time()-time0, units='mins')
     cat('\ntime elapsed for resamp_nmds:', difftime, 'minutes\n\n')
     return(out)
}

#' @export
#' @rdname resamp_nmds
### summary method for resampled statistics
`summary.resamp_nmds` <- function(object, ...){
     out <- c(object$SRV, object$rP_internal)
     out <- round(out, 3)
     names(out) <- c('SRV','rP_internal')
     out
}

#' @export
#' @rdname resamp_nmds
### plot resampled scores as a spider around originals (centroids)
`plot.resamp_nmds` <- function(x, noaxes=F, col='#00000040',
                             type='spider', ...){
     stopifnot(class(x)=='resamp_nmds')
     r     <- x$points
     # sloppy hack to handle punctuation in rownames:
     rn    <- row.names(r)
     rn    <- gsub('^[[:punct:]]', '', rn)
     left  <- gsub('\\.[[:alnum:]]+$', '', rn)
     right <- gsub('.*\\.', '', rn)
     if(all(left=='' | left=='.')) { rn <- right }else{ rn <- left }
     type <- match.arg(type, c('spider', 'points', 'hist'))
     if (type=='hist'){
          hist(x$rP_draw, col=col, main='',
               xlab=bquote(paste(R[P], ' (internal)')), ...)
     } else {
          if(noaxes){
               vegan::ordiplot(r,type='n',display='sites',
                               bty='n',axes=F, ylab='',xlab='')
          }else{
               vegan::ordiplot(r,type='n',display='sites',
                               bty='l',las=1,
                               xaxt='none',yaxt='none',
                               ylab='NMDS2',xlab='NMDS1')
          }
          if(type=='spider'){
               vegan::ordispider(r, groups=as.factor(rn), col=col,...)

          }
          if(type=='points'){
               b <- viridis(min(nlevels(as.factor(rn)), 99), alpha=.7)
               points(r, pch=16, col=b[as.factor(rn)])
          }
     }
}

######################################################################
###     UNEXPORTED functions:     #############################
######################################################################

### dissimilarity function (replaces defunct `step_bray0`)
###  zero-adjusted dissimilarities *iff* zero-sum rows exist
###  shortest-path stepacross *iff* no-share SU pairs exist
`dissim` <- function(...){
     `f` <- function(x, method, zeroadj, step, ...){
          x    <- as.matrix(x)
          zero <- any(rowSums(x, na.rm=TRUE) == 0)
          if (zero && zeroadj && method=='bray') {
               val <- min(x[x != 0], na.rm=TRUE)*0.5
               x   <- cbind(x, rep(val, nrow(x)))
          }
          D   <- vegan::vegdist(x, method=method, ...)
          nan <- is.nan(D)
          noshare <- any(vegan::no.shared(x))
          if ((any(nan)||noshare) && step) {
               D[nan] <- NA
               D <- vegan::stepacross(D, path='shortest',
                                      toolong=0, trace=FALSE)
          }
          D
     }
     f(...)
}
`subset_d` <- function(...){
     f <- function(Dx, ii, ...){
          as.dist(as.matrix(Dx)[ii,ii])
     }
     f(...)
}
`subset_x` <- function(...){
     f <- function(x, ii, ...){
          x[ii,]
     }
     f(...)
}
`ordfn` <- function(...){
     f <- function(D, k, x, ...){
          obj <- vegan::metaMDS(comm=D, k=k, trymax=99,
                                autotransform=FALSE, noshare=FALSE,
                                wascores=FALSE, trace=0,
                                plot=FALSE, weakties=TRUE)
          obj$rnk <- apply(vegan::wascores(obj$points, x),2,rank)
          snip    <- c('diss','dhat','dist','iidx','jidx')
          obj[names(obj)%in%snip] <- NULL # remove hefty items
          obj
     }
     f(...)
}
`jackfn` <- function(iii, x, Dx, k, N, old, ...){
     ii     <- !(1:N)%in%iii            # rows to keep
     x      <- subset_x(x, ii)          # downsample species matrix
     D      <- subset_d(Dx, ii, ...)    # downsample dissim matrix
     oldscr <- old$points[ii, ]         # downsample old scores
     neword <- ordfn(D, k, x, ...)
     newscr <- neword$points
     pro    <- vegan::protest(oldscr, newscr, symm=T, perm=0)
     neword$rP     <- pro$t0            # procrustes fit to reference
     neword$points <- pro$Yrot          # pass resampled scores
     neword
}
`bootfn` <- function(x, Dx, k, BS, N, old, ...){
     i <- 1
     while( i < 999 ){
          ii  <- sample(1:N, size=BS, replace=T)  # rows to keep
          x   <- subset_x(x, ii)        # downsample species matrix
          D   <- subset_d(Dx, ii, ...)  # downsample dissim matrix
          oldscr <- old$points[ii, ]    # downsample old scores
          zsum  <- (colSums(x, na.rm=T)==0)       # check empty rows
          if (ecole::mx_valid(x[,!zsum])) break   # check valid
          i <- i + 1
     }
     if (i>1) cat(sprintf('%d iters of rejection resampling\n',i))
     neword <- ordfn(D, k, x, ...)
     newscr <- neword$points
     pro    <- vegan::protest(oldscr, newscr, symm=T, perm=0)
     neword$rP     <- pro$t0            # procrustes fit to reference
     neword$points <- pro$Yrot          # pass resampled scores
     neword
}
