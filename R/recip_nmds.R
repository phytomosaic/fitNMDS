#' @title Reciprocal NMDS
#'
#' @description Reciprocal nonmetric multidimensional scaling (NMDS)
#'     to estimate external exchangeability of two datasets.
#'
#' @param x list of species and environment matrices, of class
#'     \sQuote{twin}
#'
#' @param k number of dimensions sought in final NMDS solution
#'
#' @param method dissimilarity index, per \code{\link[vegan]{vegdist}}
#'
#' @param object result from \code{recip_nmds}
#'
#' @param type for plotting, one of \code{c('points', 'text',
#'     'twinned')}
#'
#' @param leg logical, should legend be plotted?
#'
#' @param noaxes logical, should a clean configuration with no axes be
#'      plotted?
#'
#' @param ... additional arguments passed to function
#'
#' @return
#' List of class \sQuote{recip_nmds} with elements:\cr
#'
#' \describe{
#'   \item{m1}{reciprocal model 1, calibrated with dataset 1 then
#'          swapping in dataset 2)
#'   }
#'   \item{m2}{reciprocal model 2, calibrated with dataset 2 then
#'          swapping in dataset 1
#'   }
#'   \item{grp}{vector identifying group membership in dataset 1 or 2
#'   }
#'   \item{sumtab}{data.frame summary table, with elements:
#' \describe{
#'   \item{d1_rP, d2_rP}{\emph{partial} intermodel fit for group 1,
#'   group 2 data (i.e., how well each dataset fits the opposing
#'   calibration model)
#'   }
#'   \item{rP_external}{\emph{complete} intermodel fit (degree of
#'   overall external exchangebility among the two datasets, measured
#'   as Procrustean agreement of the two reciprocal models)
#'   }
#'   \item{stress1, stress2}{respective stress of model 1, model 2
#'   }
#'   \item{varexp1, varexp2}{respective variance explained of model 1,
#'   model 2
#'   }
#'   }
#'   }
#'   }
#'
#' @details
#' Reciprocal NMDS with \code{recip_nmds} estimates external
#' exchangeability of two candidate datasets.  The premise of
#' reciprocal NMDS is to calibrate each of two models with each
#' respective dataset, mutually exchanging datasets among the
#' alternative ordination models, then determining how well each
#' dataset fit the alternative model.  This is really just a special
#' case of 2-fold cross-validation as applied to NMDS.  The primary
#' use of reciprocal NMDS is to test the null hypothesis that two
#' multivariate datasets are exchangeable.
#'
#' @examples
#' # prepare two candidate datasets (here we just modify one)
#' set.seed(231)
#' data(smoky)
#' spe1 <- smoky$spe
#' env1 <- env2 <- smoky$env
#' spe2 <- spe1 + abs(rnorm(prod(dim(spe1)), 0, 2))  # add noise
#' tw   <- twin(spe1, spe2, env1, env2)
#'
#' # reciprocal NMDS
#' r <- recip_nmds(tw)
#' summary(r)
#' plot(r)
#' plot(r, 'text')
#' plot(r, 'twin')
#' plot_marg_grp(r)
#' plot_axismatch(r)
#' plot_axismatch(r, ann=FALSE)
#'
#' @references
#' Kruskal, J. B. 1964. Multidimensional scaling by optimizing
#'   goodness of fit to a nonmetric hypothesis. Psychometrika 29:
#'   1-27.
#'
#' @family reciprocal NMDS functions
#' @seealso \code{\link{boot_nmds}} for bootstrap NMDS, and
#'      \code{\link{nearestspecies}} to force same number of rows
#'      among candidate datasets
#'
#' @export
#' @rdname recip_nmds
### reciprocal NMS model fitting
`recip_nmds` <- function(x, k=2, method='bray', ...){
     if(class(x) != 'twin') stop('input must be class `twin`')
     d1 <- x[[1]][['spe']]    # dataset 1
     d2 <- x[[2]][['spe']]    # dataset 2
     i1 <- (1:nrow(d1))       # index d1
     i2 <- ((nrow(d1)+1):(nrow(d1)*2)) # index d2
     grp <- as.factor(c(rep('d1',nrow(d1)), rep('d2',nrow(d2))))
     # join both datasets for the joint model
     d12 <- ecole::mx_rbind_all(d1, d2)
     d21 <- ecole::mx_rbind_all(d2, d1)
     # distance matrices
     D1  <- dissim(d1, method, ...)
     D2  <- dissim(d2, method, ...)
     D12 <- try(dissim(d12, method, ...), TRUE)
     D21 <- try(dissim(d21, method, ...), TRUE)
     # calibration models, separate datasets
     o1_ <- vegan::metaMDS(D1, k=k, trymax=99, autotransform=FALSE,
                           noshare=FALSE, wascores=FALSE, trace=0,
                           plot=FALSE, weakties=TRUE)$points
     o2_ <- vegan::metaMDS(D2, k=k, trymax=99, autotransform=FALSE,
                           noshare=FALSE, wascores=FALSE, trace=0,
                           plot=FALSE, weakties=TRUE)$points
     # reciprocal models, predict scores on 'swapped' data
     o12 <- NMSpredict(scr=o1_, dis=D12, neighb=10, maxits=200)
     o21 <- NMSpredict(scr=o2_, dis=D21, neighb=10, maxits=200)
     stress1 <- o12$stress
     stress2 <- o21$stress
     o12     <- o12$newpoints     # dataset 2 scores in model 1
     o21     <- o21$newpoints     # dataset 1 scores in model 2
     # procrustes alignment
     pp <- vegan::protest(rbind(o1_,o12),rbind(o21,o2_),perm=0,symm=T)
     m1 <- pp$X               # M1 RECIPROCAL model
     m2 <- pp$Yrot            # M2 RECIPROCAL model
     colnames(m2) <- colnames(m1)
     # PARTIAL intermodel fit: compare ea dataset to ITSELF btwn mods
     d1_rP <- vegan::protest(m1[i1,], m2[i1,], perm=0, symm=T)$t0
     d2_rP <- vegan::protest(m2[i2,], m1[i2,], perm=0, symm=T)$t0
     # COMPLETE intermodel fit: compare each model to the other
     rP_external <- pp$t0  # M1 vs M2 fit
     # var expl by each configuration = R2 = coef of detn (as PCORD7)
     Ds  <- vegan::vegdist(d12, 'bray', diag=T, upper=T)
     Dz1 <- dist(m1)
     Dz2 <- dist(m2)
     ve1 <- cor(Ds,Dz1,meth='pear')^2
     ve2 <- cor(Ds,Dz2,meth='pear')^2
     # output
     out <- list(m1=m1, m2=m2, grp=grp,
                 sumtab=round(data.frame(
                      d1_rP   = d1_rP,      # partial intermodel fit
                      d2_rP   = d2_rP,      # partial intermodel fit
                      rP_external= rP_external,   # complete intermodel fit
                      stress1 = stress1,    # stress M1
                      stress2 = stress2,    # stress M2
                      varexp1 = ve1,        # var explained M1
                      varexp2 = ve2),3))    # var explained M2
     class(out) <- 'recip_nmds'
     out
}


#' @export
#' @rdname recip_nmds
### summary method for reciprocal fit statistics
`summary.recip_nmds` <- function(object, ...){
     out <- c(t(object[['sumtab']]))
     names(out) <- c('d1_rP','d2_rP','rP_external',
                     'stress1','stress2','varexp1','varexp2')
     out
}

#' @export
#' @rdname recip_nmds
### plot competing models as points, text, or side-by-side
`plot.recip_nmds` <- function(x, type='points', leg=FALSE,
                              noaxes=TRUE, ...){
     if(class(x) != 'recip_nmds') stop('input must be `recip_nmds`')
     type <- match.arg(type, c('points', 'text', 'twinned'))
     if(type=='points'){
          if(noaxes){
               vegan::ordiplot(x$m1,dis='sites',type='n',bty='n',
                               axes=F,xlab='',ylab='')
          }else{
               vegan::ordiplot(x$m1,dis='sites',type='n',bty='l',
                               xaxt='none',yaxt='none',
                               xlab='NMDS1',ylab='NMDS2')
          }
          points(x$m1, col=rgb(0,0,0,0.7), pch=16, cex=0.5)
          arrows(x$m1[,1], x$m1[,2], x$m2[,1], x$m2[,2],
                 col=rgb(0,0,0,0.3), len=0.02, angle=20)
          if(leg){
               legend('topright', c('m1','m2'), col=c(1,2),
                      border='transparent', pch=c(16,16), bty='o',
                      cex=0.7, title='Models', horiz=T)
          }
     }
     if(type=='text'){
          if(noaxes){
               vegan::ordiplot(x$m1,dis='sites',type='n',bty='n',
                               axes=F,xlab='',ylab='')
          }else{
               vegan::ordiplot(x$m1,dis='sites',type='n',bty='l',
                               xaxt='none',yaxt='none',
                               xlab='NMDS1',ylab='NMDS2')
          }
          text(x$m1, lab=row.names(x$m1), cex=0.7)
          segments(x$m1[,1],x$m1[,2],x$m2[,1],x$m2[,2],col=17)
          text(x$m2, lab=row.names(x$m2), cex=0.7, col=2)
          if(leg){
               legend('topright', c('m1','m2'), col=c(1,2),
                      border='transparent', pch=c(16,16), bty='o',
                      cex=0.7, title='Models', horiz=T)
          }
     }
     if(type=='twinned'){
          op <- par(mfrow=c(1,2), oma=rep(0,4)+.1, mar=c(2,2,2,0)+.1)
          vegan::ordiplot(x$m1,dis='sites',type='t',xlab='NMDS1',
                          ylab='NMDS2',main='Model 1: D1+D2')
          vegan::ordiplot(x$m2,dis='sites', type='t',xlab='NMDS1',
                          ylab='NMDS2',main='Model 2: D1+D2')
          par(op)
     }
}

#' @export
#' @rdname recip_nmds
### plot reciprocal models as grouped scatterplot with marg histogram
`plot_marg_grp` <- function(object, ...){
     x <- object
     if(class(x) != 'recip_nmds') stop('input must be `recip_nmds`')
     x1 <- x$m1[,1]
     y1 <- x$m1[,2]
     x2 <- x$m2[,1]
     y2 <- x$m2[,2]
     op <- par(no.readonly=TRUE)
     on.exit(par(op))
     x1den <- density(x1)
     x2den <- density(x2)
     y1den <- density(y1)
     y2den <- density(y2)
     xlimv <- c(min(x1den$x, x2den$x), max(x1den$x, x2den$x))*1.2
     ylimv <- c(0, max(x1den$y, x2den$y))
     xlimh <- c(0, max(y1den$y, y2den$y))
     ylimh <- c(min(y1den$x, y2den$x), max(y1den$x, y2den$x))*1.2
     layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
     par(mar=c(3,3,0,0))
     plot(x1, y1, xlim=xlimv, ylim=ylimh, las=1,
          xlab='NMDS 1', ylab='NMDS 2', col=rgb(0,0,0,0.5), pch=19)
     points(x2, y2, xlim=xlimv, ylim=ylimh, col=rgb(1,0,0,0.5),pch=19)
     par(mar=c(0,3,1,1))
     plot(x1den, axes=F, xlim=xlimv, ylim=ylimv,
          type='n', xlab='', ylab='', main='')
     polygon(x1den$x, x1den$y, col=rgb(0,0,0,0.5), xpd=NA)
     polygon(x2den$x, x2den$y, col=rgb(1,0,0,0.5), xpd=NA)
     par(mar=c(3,0,1,1))
     plot(y1den$y, y1den$x, axes=F, xlim=xlimh, ylim=ylimh,
          type='n', xlab='', ylab='', main='')
     polygon(y1den$y, y1den$x, col=rgb(0,0,0,0.5), xpd=NA)
     polygon(y2den$y, y2den$x, col=rgb(1,0,0,0.5), xpd=NA)
}

#' @export
#' @rdname recip_nmds
### plot axis-wise orthogonal comparison of NMS scores after rotation
`plot_axismatch` <- function(object, ...){
     x <- object
     if(class(x) != 'recip_nmds') stop('input must be `recip_nmds`')
     dm <- dim(x$m1)[2]
     w <- dimnames(x$m1)[[2]]
     w <- paste0('N',w)
     x  <- as.data.frame(cbind(x$m1, x$m2))
     if(dm==2){
          par(mfrow=c(1,2))
          ecole::plot_ortho(x[,1],x[,(dm+1)],xlab=w[1],ylab=w[1], ...)
          ecole::plot_ortho(x[,2],x[,(dm+2)],xlab=w[2],ylab=w[2], ...)
          par(mfrow=c(1,1))
     }
     if(dm==3){
          par(mfrow=c(1,3))
          ecole::plot_ortho(x[,1],x[,(dm+1)],xlab=w[1],ylab=w[1], ...)
          ecole::plot_ortho(x[,2],x[,(dm+2)],xlab=w[2],ylab=w[2], ...)
          ecole::plot_ortho(x[,3],x[,(dm+3)],xlab=w[3],ylab=w[3], ...)
          par(mfrow=c(1,1))
     }
}


### unexported functions:


############      #' @useDynLib vegan
### NMSpredict core algorithm lightly adapted from add.points()
`NMSpredict` <- function(scr, dis, neighb, maxits){
     if (!requireNamespace("vegan", quietly = FALSE)) {
          stop("Package package `vegan` required", call.=FALSE)
     }
     # convert original scores to local matrix, and get sizes
     points  <- list(points=scr)
     class(points) <- 'nmds'
     points <- points$points
     oldn <- nrow(points)
     ndim <- ncol(points)
     totn <- attr(dis,'Size')
     newn <- totn - oldn
     # TODO: allow 1 point instead of many - until then, force > 1
     if (newn < 1) stop('come on, try more than just 1 point!')
     # test correspondence
     if (!identical(dimnames(points)[[1]],attr(dis,'Labels')[1:oldn]))
          stop('ordination and dissimilarity matrix do not match')
     # decompose dissimilarity object to isolate new values
     diss <- as.matrix(dis)[1:oldn,(oldn+1):totn]
     # set up initial conditions
     ndis <- oldn * newn
     tmp <- matrix(rep(0,newn*ndim),ncol=ndim)
     for (i in 1:newn) {
          pnt <- seq(1,oldn)[order(diss[,i])][1:neighb]
          weight <- 1-diss[pnt,i]
          for (j in 1:ncol(points)) {
               tmp[i,j] <-stats::weighted.mean(points[pnt,j],w=weight)
          }
     }
     xinit <- rbind(points,tmp)
     dimnames(xinit)[[1]] <- attr(dis,'Labels')
     # set up indices
     iidx <- rep((1:oldn),newn)
     jidx <- NULL
     for (i in (oldn+1):totn) jidx <- c(jidx,rep(i,oldn))
     # set up ordination
     nfix <- oldn
     ngrp <- istart <- 1
     isform <- 2
     ities <- 1
     iregn <- 1
     iscal <- 0
     sratmx <- 0.99999
     strmin <- 1e-07
     sfgrmn <- 1e-05
     dist <- rep(0,ndis)
     dhat <- rep(0,ndis)
     x <- matrix(0,nrow=totn,ncol=ndim)
     stress <- 1
     strs <- ngrp
     iters <- 1
     icause <- 1
     maxits <- maxits
     iwork <- rep(0,ndis)
     grad <- matrix(0,nrow=totn,ncol=ndim)
     grlast <- matrix(0,nrow=totn,ncol=ndim)
     out <- .Fortran('monoMDS',
                     nobj=as.integer(totn),
                     nfix=as.integer(nfix),
                     ndim=as.integer(ndim),
                     ndis=as.integer(ndis),
                     ngrp=as.integer(ngrp),
                     diss=as.double(diss),
                     iidx=as.integer(iidx),
                     jidx=as.integer(jidx),
                     xinit=as.double(xinit),
                     istart=as.integer(istart),
                     isform=as.integer(isform),
                     ities=as.integer(ities),
                     iregn=as.integer(iregn),
                     iscal=as.integer(iscal),
                     maxits=as.integer(maxits),
                     sratmx=as.double(sratmx),
                     strmin=as.double(strmin),
                     sfgrmn=as.double(sfgrmn),
                     dist=as.double(dist),
                     dhat=as.double(dhat),
                     points=as.double(x),
                     stress=as.double(stress),
                     strs=as.double(strs),
                     iters=as.integer(iters),
                     cause=as.integer(icause),
                     ##################
                     PACKAGE='vegan'
                     ##################
     )
     res <- list(points=matrix(out$points,ncol=ndim),
                 newpoints=matrix(out$points,
                                  ncol=ndim)[(oldn+1):totn,],
                 stress=out$stress,
                 iters=out$iters,
                 cause=out$cause)
     dimnames(res$points)[[1]] <- attr(dis,'Labels')
     dimnames(res$newpoints)[[1]] <- attr(dis,'Labels')[(oldn+1):totn]
     class(res) <- 'nmds'
     res
}
