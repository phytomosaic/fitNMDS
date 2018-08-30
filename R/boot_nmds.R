#' @title Bootstrap or jackknife NMDS
#'
#' @description Bootstrap or jackknife nonmetric multidimensional
#'     scaling (NMDS) to estimate internal sampling variability.
#'     Accompanied by summary and plot methods.
#'
#' @param spe species abundance matrix
#'
#' @param B number of bootstrap replicates
#'
#' @param BS bootstrap size (currently constrained to not exceed the
#'     number of rows in original matrix). Default \code{BS=NULL}, or
#'     else \code{BS=1}, will perform jackknife NMDS (leave-one-out).
#'
#' @param k number of dimensions sought in final NMDS solution
#'
#' @param method dissimilarity index, per \code{\link[vegan]{vegdist}}
#'
#' @param zeroadj logical, default \code{TRUE}: adjust dissimilarities
#'     to account for zero-sum rows?
#'
#' @param step logical, default \code{TRUE}: use shortest-path
#'      stepacross dissimilarity adjustment?
#'
#' @param x,object result from \code{boot_nmds}
#'
#' @param noaxes logical, default \code{TRUE}: plot without axes?
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
#' List of class \sQuote{boot_nmds} with elements:\cr
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
#'      \item{\code{bootpt}:}{
#'           Scores of each replicate configuration.
#'      }
#' }
#'
#' @details
#' Bootstrap or jackknife resampling of NMDS with \code{boot_nmds}
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
#' x   <- smoky$spe
#'
#' # simplistic dataset yields a few poor solutions, with warning:
#' res <- boot_nmds(x, B=29, k=2)
#' summary(res)
#' plot(res)
#' plot(res, type='points')
#' plot(res, type='hist')
#'
#' # alternatively, plot replicates colored by each sample unit:
#' rn  <- gsub('\\.[[:alnum:]]+$', '', row.names(res$bootpt))
#' plot(res, col=ecole::colvec(as.factor(rn), alpha=0.5))
#'
#' @references
#' Clarke, K.R., P.J. Somerfield, and M.G. Chapman. 2006. On
#'     resemblance measures for ecological studies, including
#'     taxonomic dissim and a zero-adjusted Bray–Curtis coefficient
#'     for denuded assemblages. J Exp Marine Biol and Ecol 330:55–80.
#' Knox, R. G., and R. K. Peet. 1989. Bootstrapped ordination: a
#'     method for estimating sampling effects in indirect gradient
#'     analysis. Vegetatio 80:153–165.
#' Kruskal, J. B. 1964. Multidimensional scaling by optimizing
#'     goodness of fit to a nonmetric hypothesis. Psychometrika 29:
#'     1-27.
#'
#' @family bootstrap NMDS functions
#' @seealso \code{\link{recip_nmds}} for reciprocal NMDS, and
#'     \code{\link[vegan]{metaMDS}} for the core NMDS algorithm
#'
#' @export
#' @rdname boot_nmds
### bootstrap NMDS core function
`boot_nmds` <- function(spe, B=9, BS=NULL, k=2, method='bray',
                        zeroadj=TRUE, step=TRUE, ...){
     ### setup
     time0 <- Sys.time()
     M     <- NCOL(spe)
     N     <- NROW(spe)
     if (is.null(BS)) BS <- 1
     if (BS > N) BS <- N

     ### dissimilarity function (replaces defunct `step_bray0`)
     ###  zero-adjusted dissimilarities *iff* zero-sum rows exist
     ###  shortest-path stepacross *iff* no-share SU pairs exist
     `dissim` <- function(zeroadj, step, ...){
          x    <- as.matrix(spe)
          zero <- any(rowSums(x, na.rm=TRUE) == 0)
          if(zero && zeroadj && method=='bray'){
               val <- min(x[x != 0], na.rm=TRUE)*0.5
               x   <- cbind(x, rep(val, nrow(x)))
          }
          D   <- vegan::vegdist(x, method=method, ...)
          nan <- is.nan(D)
          noshare <- any(no.shared(x))
          if((any(nan)||noshare) && step){
               D[nan] <- NA
               # SILENT <- capture.output(
               D <- vegan::stepacross(D, path='shortest',
                                      toolong=0, trace=FALSE)
               # )
          }
          D
     }

     ### original (full) dissimilarities matrix Dx; D will change
     Dx <- D <- try(dissim(zeroadj, step, ...), TRUE)

     ###   general function to perform the ordinations
     `ordfn` <- function(...){
          ord <- vegan::metaMDS(D, k=k, trymax=99, autotransform=FALSE,
                                noshare=FALSE, wascores=FALSE, trace=0,
                                plot=FALSE, weakties=TRUE)
          # calc species ranks
          ord$rnk<-apply(vegan::wascores(ord$points,spp),2,rank)
          ord
     }
     ### original reference ordination
     spp <- spe
     old <- ordfn()

     ### for jackknife (leave out each SU once per ordn)
     if (BS == 1) {
          Draws <- sapply(1:N, FUN=function(sk){
               rj  <- (1:N)%in%sk            # rows to reject
               spp <<- spe[!rj, ]            # downsample species matrix
               oldscr <- old$points[!rj, ]   # downsample old scores
               D <<- as.dist(as.matrix(Dx)[!rj,!rj]) # downsample dissim mx
               neword <- ordfn()
               newscr <- neword$points
               pro    <- vegan::protest(oldscr, newscr, symm=T, perm=0)
               neword$rP     <- pro$t0   # procrustes fit to reference ordn
               neword$points <- pro$Yrot # pass bootstrapped scores
               neword
          })
          NC <- N
     } else {
          ### for bootstrap (randomly draw n SUs with replacement)
          `bootfn` <- function(...){
               i <- 1
               while( i < 999 ){
                    rj  <- sample(1:N, size=BS, replace=T) # rows to KEEP
                    # rj  <- (1:N)%in%sk            # rows to KEEP
                    spp <<- spe[rj, ]            # downsample species matrix
                    oldscr <- old$points[rj, ]   # downsample old scores
                    D <<- as.dist(as.matrix(Dx)[rj,rj]) # downsample dissim mx
                    zsum   <- (colSums(spp, na.rm=T)==0)
                    if(ecole::mx_valid(spp[,!zsum])) break
                    i <- i + 1
               }
               if(i>1) cat(sprintf('%d iters of rejection resampling\n',i))
               neword <- ordfn()
               newscr <- neword$points
               pro    <- vegan::protest(oldscr, newscr, symm=T, perm=0)
               neword$rP     <- pro$t0   # procrustes fit to reference ordn
               neword$points <- pro$Yrot # pass bootstrapped scores
               neword
          }
          Draws <- replicate(B, bootfn())
          NC <- B
     }

     ### summarize bootstrap/jackknife values
     # bootstrap scores (eval sample unit score variability)
     bscr <- do.call(rbind,lapply(1:NC,function(nc)Draws[,nc]$points))
     # SRV of replicate ordinations (eval species score variability)
     rnk <- Draws['rnk', ]
     xx  <- matrix(NA, nrow=prod(dim(Draws[,1]$rnk)), ncol=NC)
     for(r in 1:length(rnk)) { xx[,r] <- c(rnk[[r]]) } # unfold
     rv_obs  <- mean(rowSums((xx - rowMeans(xx))^2)/(dim(xx)[2] - 1))
     rv_exp  <- (M^2 - 1)/12        # random expectation
     SRV     <- rv_obs / rv_exp     # Scaled Rank Variance
     # median rP fit to original ordination, across all draws
     rP_draw <- unlist(Draws['rP', ], use.names=F)
     rP_medn <- median(rP_draw)
     stresses <- unlist(Draws['stress', ], use.names=F)

     ### output list
     out <- list(rP_draw=rP_draw,
                 stresses=stresses,
                 SRV=SRV,
                 rP_internal=rP_medn,
                 bootpt=bscr)
     class(out) <- 'boot_nmds'

     ### exit timing
     difftime <- as.numeric(Sys.time()-time0, units='mins')
     cat('\ntime elapsed for boot_nmds:', difftime, 'minutes\n\n')
     return(out)
}

#' @export
#' @rdname boot_nmds
### summary method for bootstrapped statistics
`summary.boot_nmds` <- function(object, ...){
     out <- c(object$SRV, object$rP_internal)
     out <- round(out, 3)
     names(out) <- c('SRV','rP_internal')
     out
}

#' @export
#' @rdname boot_nmds
### plot bootstraps as a spider around each centroid
`plot.boot_nmds` <- function(x, noaxes=F, col='#00000040',
                             type='spider', ...){
     stopifnot(class(x)=='boot_nmds')
     r     <- x$bootpt
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
