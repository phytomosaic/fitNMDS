#' @name smoky
#' @title Whittaker's Great Smoky Mountains Vegetation Data (USA)
#' @aliases smoky
#' @docType data
#' @description
#' A list of two data.frames \sQuote{spe} and \sQuote{env} derived
#'     from Whittaker's (1956) Table 3, showing vegetation changes
#'     along a moisture gradient in Great Smoky Mountains National
#'     Park, USA.
#'
#' @format A list of 2 data.frames:
#' \itemize{
#' \item \code{spe} 12 observations of 41 woody plant species
#' \item \code{env} 12 observations of 7 nvironmental variables
#' }
#' @details Species matrix values are abundance of tree species in 12
#'     stations, (percentages of total woody plant stems per station
#'     > or = 1-inch diameter). Each station is an aggregate of 1-7
#'     plots of variable size and sampling intensity. Presences < 0.5%
#'     are coded as 0.1 in the present matrix.\cr
#'
#'     Environmental matrix values are:\cr
#'     \itemize{
#'     \item \code{mesic}    mesic indicator value
#'     \item \code{submesic} submesic indicator value
#'     \item \code{subxeric} subxeric indicator value
#'     \item \code{xeric}    xeric indicator value
#'     \item \code{treect} count of trees per station
#'     \item \code{nsites} count of sites per station
#'     \item \code{moisture} position on a putative moisture gradient,
#'         ranging from 1 = mesic to 12 = xeric
#'     }
#'     where the first four variables are species weighted averages as
#'     moisture indicator values.\cr
#'
#'     Whittaker's caption verbatim:\cr
#'     \dQuote{Table 3. Composite transect of moisture gradient
#'     between 3500 and 4500 ft, distribution of trees along gradient.
#'     Transect along the moisture gradient from mesic valley sites
#'     (Sta. 1) to xeric southwest slope sites (Sta. 12), based on 46
#'     site counts including 4906 stems from elevations between 3500
#'     ft and 4500 ft. All figures are percentages of total stems in
#'     station from 1-in. diameter class up.
#'     }\cr
#'
#'     Location: Great Smoky Mountains of Tennessee and North
#'     Carolina, USA.
#'
#' @source Table 3 in Whittaker (1956).
#'
#' @references
#' Whittaker, R. H. 1956. Vegetation of the Great Smoky Mountains.
#'     Ecological Monographs 26:2–80.
#'
#' @examples
#' # split into two data.frames
#' data(smoky)
#' spe <- smoky$spe
#' env <- smoky$env
#'
#' # visualize the species matrix
#' op <- par(mfrow=c(1,1), mar=c(0,2,10,0), oma=c(0,0,0,0), cex=0.7)
#' r  <- ecole::colvec(1:99, alpha=1)
#' cc <- log1p(spe[12:1,])
#' image((1:ncol(cc)), (1:nrow(cc)), t(cc), col=r, axes=FALSE)
#' axis(3, at = 1:ncol(cc), labels=colnames(cc), las=3, tick=FALSE)
#' axis(2, at = 1:nrow(cc), labels=rownames(cc), las=1, tick=FALSE)
#' par(op)
#'
#' # following Whittaker's Fig. 4 (top):
#' e <- sapply(env[,1:4], FUN=function(x)(x-min(x))/(max(x)-min(x)))
#' e <- cbind(e, moisture=env$moisture)
#' plot(1:12, ylim=c(0,1), type='n', las=1, bty='l', ylab='')
#' for(i in 1:4){
#'      points(e[,i], pch=16, col=i)
#'      lines(loess(e[,i] ~ e[,'moisture']), col=i)
#' }
#'
#' # following Whittaker's Fig. 4 (middle):
#' nm <- c('tilia_heterophylla','halesia_monticola',
#'          'tsuga_canadensis','quercus_alba','pinus_pungens')
#' e <- sapply(spe[,names(spe)%in%nm],
#'            FUN=function(x)(x-min(x))/(max(x)-min(x)))
#' e <- cbind(e, moisture=env$moisture)
#' plot(1:12, ylim=c(0,1), type='n', las=1, bty='l', ylab='')
#' for(i in 1:5){
#'      points(e[,i], pch=16, col=i)
#'      lines(loess(e[,i] ~ e[,'moisture']), col=i)
#' }
#'
#' @keywords datasets
"smoky"