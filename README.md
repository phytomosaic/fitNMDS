# fitNMDS
Measures of agreement for bootstrapped and reciprocal NMDS.


## Motivation

Combining two different datasets into one nonmetric multidimensional scaling (NMDS) model can be risky if they each cover different
attribute spaces (e.g., different species pools in ecology). Therefore, comparing two datasets requires estimating internal agreement
(sampling variability) relative to external agreement (exchangeability). Bootstrapped NMDS estimates internal agreement of a candidate dataset, while 
reciprocal NMDS estimates external agreement among two candidate datasets.


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

### Bootstrap NMDS of one dataset
```r
x   <- list(spe=spe1, id=env1)
res <- boot_nmds(x, B=29, k=2, rot=TRUE)
summary(res)
plot(res, col='#00000040')
```


### Reciprocal NMDS of both datasets
```r
res <- recip_nmds(tw)
summary(res)
plot(res)
```


### Correct unequal sample sizes among two datasets
```r
# two candidate datasets, full and partial
spe2 <- spe1[1:11,]
env2 <- env1[1:11,]

<<<<<<< HEAD
# subset the full matrix, based on compositional nearest neighbors...
=======
# subset the full matrix, based on compositional nearest neighbors
>>>>>>> 60fe352f2f7f7f1d367228b86662ad4f281f9dc3
(i   <- nearestspecies(spe1, spe2, ties=FALSE))
spe1 <- spe1[i,,]
env1 <- env1[i,,]

<<<<<<< HEAD
# ...then proceed to reciprocal NMDS
=======
# then proceed to reciprocal NMDS
>>>>>>> 60fe352f2f7f7f1d367228b86662ad4f281f9dc3
tw   <- twin(spe1, spe2, env1, env2)
res  <- recip_nmds(tw)
summary(res)
plot(res)
<<<<<<< HEAD
=======

>>>>>>> 60fe352f2f7f7f1d367228b86662ad4f281f9dc3
```
