#' @title Bootstrap NMDS
#'
#' @description Bootstrap nonmetric multidimensional scaling (NMDS) to
#'     estimate internal sampling variability.  Accompanied by summary
#'     and plot methods.
#'
#' @param x list of species and environment matrices
#'
#' @param B number of bootstrap replicates
#'
#' @param BS bootstrap size (if not equal to number of rows in
#'     original matrix)
#'
#' @param k number of dimensions sought in final NMDS solution
#'
#' @param rot logical, should solution be rotated to similarity with
#'    an environmental gradient? Probably only sensible with known
#'    synthetic gradients.
#'
#' @param object result from \code{boot_nmds}
#'
#' @param noaxes logical, plot without axes?
#'
#' @param col character vector of colors for spider lines
#'
#' @param ... additional arguments passed to function
#'
#' @return
#' List of class \sQuote{boot_nmds} with elements:\cr
#' \describe{
#'   \item{rP_draw}{vector of rP_internal values for each replicate
#'   }
#'   \item{stresses}{vector of NMDS stress for each replicate
#'   }
#'   \item{sumtab}{summary across \emph{all} replicates (\strong{SRV}:
#'    scaled rank variance; \strong{rP_internal}: internal sampling
#'    variability)
#'   }
#'   \item{bootpt}{scores of each replicate configuration
#'   }
#'   }
#'
#' @details
#' Bootstrapped NMDS with \code{boot_nmds} estimates internal sampling
#' variability of a candidate dataset.  The premise of bootstrapped
#' ordination is resampling with replacement from the reference
#' dataset, performing NMDS ordination (Kruskal 1964) on each
#' bootstrap replicate, then determining collective sampling
#' variability across all NMDS solutions.  Other uses of bootstrap
#' NMDS may include estimating confidence regions for site scores, or
#' testing the stability of species positions along ordination axes.
#'
#' @examples
#' set.seed(231)
#' data(smoky)
#' x   <- list(spe=smoky$spe, id=smoky$env[,1:4])
#' res <- boot_nmds(x, B=29, k=2, rot=TRUE)
#' summary(res)
#' plot(res, col='#00000040') # 'spider' plot
#'
#' # alternative plotting:
#' rn  <- gsub('\\.[[:alnum:]]+$', '', row.names(res$bootpt))
#' plot(res, col=ecole::colvec(as.factor(rn), alpha=0.5))
#'
#' @references
#' Kruskal, J. B. 1964. Multidimensional scaling by optimizing
#'   goodness of fit to a nonmetric hypothesis. Psychometrika 29:
#'   1-27.
#'
#' @family bootstrap NMDS functions
#' @seealso \code{\link{recip_nmds}} for reciprocal NMDS, and
#'     \code{\link[vegan]{metaMDS}} for the core NMDS algorithm
#'
#' @export
#' @rdname boot_nmds
### bootstrap NMDS core function
`boot_nmds` <- function(x, B=9, BS, k=2, rot=TRUE, method='bray', ... ){
     time0 <- Sys.time()
     environment(ordfn)  <- environment()
     environment(bootfn) <- environment()
     environment(jackfn) <- environment()
     spe <- x[['spe']]
     id  <- x[['id']]
     M   <- ncol(spe)
     N   <- nrow(spe)
     if(missing(BS)) BS <- N
     if( BS > N ) BS <- N
     old <- ordfn(spe, id, k, rot, ...) # orig reference ordination
     # do B bootstrap or N jackknife draws
     if (BS==1) {
          Draws <- sapply(1:N, FUN=function(x){
               jackfn(spe, id, k, rot, old, indx=x, ...)})
          NC <- N
     } else {
          Draws <- replicate(B, bootfn(spe, id, k, rot, old, BS, ...))
          NC <- B
     }
     # get bootstrap scores
     bscr <- do.call(rbind, lapply(1:B, function(x)Draws[,x]$points))
     # get SRV of bootstrapped ordns (to assess species stability)
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
     out <- list(rP_draw=rP_draw,
                 stresses=stresses,
                 sumtab=round(data.frame(SRV=SRV,
                                         rP_internal=rP_medn),3),
                 bootpt=bscr)
     class(out) <- 'boot_nmds'
     difftime <- as.numeric(Sys.time()-time0, units='mins')
     cat('\ntime elapsed for bootstrap:', difftime, 'minutes\n\n')
     return(out)
}

#' @export
#' @rdname boot_nmds
### summary method for bootstrapped statistics
`summary.boot_nmds` <- function(object, ...){
     out <- c(t(object[['sumtab']]))
     names(out) <- c('SRV','rP_internal')
     out
}

#' @export
#' @rdname boot_nmds
### plot bootstraps as a spider around each centroid
###     (formerly `plot_boot_spider`)
`plot.boot_nmds` <- function(object, noaxes=F,
                             col=rgb(0,0,0,10,max=255), ...) {
     r     <- object$bootpt
     # sloppy hack to handle punctuation in rownames:
     rn    <- row.names(r)
     rn    <- gsub('^[[:punct:]]', '', rn)
     left  <- gsub('\\.[[:alnum:]]+$', '', rn)
     right <- gsub('.*\\.', '', rn)
     if(all(left=='' | left=='.')) { rn <- right }else{ rn <- left }
     if(noaxes){
          vegan::ordiplot(r,type='n',display='sites',bty='n',axes=F,
                          ylab='',xlab='')
     }else{
          vegan::ordiplot(r,type='n',display='sites',bty='l',las=1,
                          xaxt='none',yaxt='none',
                          ylab='NMDS2',xlab='NMDS1')
     }
     vegan::ordispider(r, groups=as.factor(rn), col=col)
}


### unexported functions:


# ### stepacross zero-adjusted Bray-Curtis dissimilarities
# ###     toolong=1 is fixed (data may be disconnected if <1)
# `step_bray0` <- function(x, ...){
#      # need capture.output to silence unwanted message printing:
#      dontprint <- capture.output(
#           out <- vegan::stepacross(ecole::bray0(x), path='shortest',
#                                    toolong=1),
#           file=NULL
#      )
#      out
# }

### dissimilarity convenience function (replaces defunct `step_bray0`)
###    use 'zero-adjusted' BC dissims only if zero-sum rows exist;
###    stepacross only if no-share sample unit pairs exist (and the
###    user allows).
`dissim` <- function(x, method='bray', zeroadj=TRUE,
                     allowstep=TRUE, ...){
     x    <- as.matrix(x)
     zero <- any(rowSums(x, na.rm = T) == 0)
     if(zero && zeroadj && method=='bray'){
          val <- min(x[x != 0], na.rm=TRUE)*0.5
          x   <- cbind(x, rep(val, nrow(x)))
     }
     D   <- vegan::vegdist(x, method=method, ...)
     nan <- is.nan(D)
     noshare <- any(no.shared(x))
     if( (any(nan) || noshare) && allowstep )   {
          D[nan] <- NA
          SILENT <- capture.output(
               D <- vegan::stepacross(D, path = 'shortest',
                                      toolong=0, ...)
          )
     }
     D
}

### general function to perform the ordinations
`ordfn` <- function(spe, id, k, rot, method=method, ......){
     D    <- try(dissim(spe, method=method, ...), TRUE)
     ord  <- vegan::metaMDS(D, k=k, trymax=99, autotransform=FALSE,
                            noshare=FALSE, wascores=FALSE, trace=0,
                            plot=FALSE, weakties=TRUE)
     scr  <- ord$points
     # rotate to similarity with an environmental gradient?
     if(rot){
          envscr  <- vegan::wcmdscale(vegan::vegdist(id,'euc'), eig=T)
          kmin    <- min(k, dim(envscr$points)[[2]])
          envscr  <- envscr$points[,1:kmin]
          scr     <- vegan::protest(envscr, scr, symm=T, perm=0)$Yrot
     }
     # species ranks:
     ord$rnk <- apply(vegan::wascores(scr, spe, expand=F), 2, rank)
     ord$points <- scr
     ord
}

### for jackknife (leave out each SU once per ordn)
`jackfn` <- function(spe, id, k, rot, old, indx, ... ){
     environment(ordfn) <- environment()
     oldscr <- old$points[-indx,]
     neword <- ordfn(spe[-indx,], id[-indx,], k, rot)
     newscr <- neword$points
     pro           <- vegan::protest(oldscr, newscr, symm=T, perm=0)
     neword$rP     <- pro$t0   # procrustes fit to ORIG reference ordn
     neword$points <- pro$Yrot # pass bootstrapped scores
     neword
}

### for bootstrap (randomly draw n SUs with replacement)
`bootfn` <- function(spe, id, k, rot, old, BS, ... ){
     environment(ordfn) <- environment()
     N <- nrow(spe)
     i <- 1
     while( i < 999 ){
          indx   <- sample(1:N, size=BS, replace=T)
          oldscr <- old$points[indx,]
          zsum   <- (colSums(spe[indx,], na.rm=T)==0)
          if(ecole::mx_valid( spe[indx,!zsum] )) break
          i <- i + 1
     }
     if(i>1) cat(sprintf('%d iters of rejection sampling\n',i))
     neword <- ordfn(spe[indx,], id[indx,], k, rot)
     newscr <- neword$points
     pro           <- vegan::protest(oldscr, newscr, symm=T, perm=0)
     neword$rP     <- pro$t0   # procrustes fit to ORIG reference ordn
     neword$points <- pro$Yrot # pass bootstrapped scores
     neword
}
