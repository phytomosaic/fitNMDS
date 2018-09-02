# fitNMDS
Fitting resampled and reciprocal NMDS models.


## Motivation

Combining two different datasets into one nonmetric multidimensional scaling (NMDS) model can be risky if they each cover different attribute spaces (e.g., different species pools in ecology). Therefore, comparing two datasets requires estimating internal sampling variability (using resampled NMDS) relative to external exchangeability (using reciprocal NMDS).


## Installation

Install the package from github as follows:
```r
install.packages('devtools')
devtools::install_github('phytomosaic/fitNMDS')
require(fitNMDS)
```


## Usage

### Prepare two candidate datasets (here we just modify one)
```r
set.seed(231)
data(smoky)
spe1 <- smoky$spe
env1 <- env2 <- smoky$env
spe2 <- spe1 + abs(rnorm(prod(dim(spe1)), 0, 2))  # add noise
tw   <- twin(spe1, spe2, env1, env2)
```


### Resampled NMDS of one dataset
```r
x   <- spe1
res <- resamp_nmds(x, k=2)
summary(res)
plot(res, col='#00000050')
```


### Reciprocal NMDS of both datasets
```r
res <- recip_nmds(tw)
summary(res)
plot(res, noaxes=FALSE)
```


### Correct unequal sample sizes prior to reciprocal NMDS
```r
# two candidate datasets, full and partial
spe2 <- spe1[1:11,]
env2 <- env1[1:11,]

# subset the full matrix, based on compositional nearest neighbors...
(i   <- nearestspecies(spe1, spe2, ties=FALSE))
spe1 <- spe1[i,,]
env1 <- env1[i,,]

# ...then proceed to reciprocal NMDS
tw   <- twin(spe1, spe2, env1, env2)
res  <- recip_nmds(tw)
summary(res)
plot(res)     # expect near-perfect match
```
