# fitNMDS
Measures of fit for bootstrapped and reciprocal NMDS.


## Motivation

Combining two different datasets into one nonmetric multidimensional scaling (NMDS) model can be risky if they each cover different
attribute spaces (e.g., different species pools in ecology). Therefore, comparing two datasets requires estimating internal fit
(stability) relative to external fit (exchangeability). Bootstrapped NMDS estimates internal fit of a candidate dataset, while 
reciprocal NMDS estimates external fit of two candidate datasets.


## Installation

Install the package from github as follows:
```r
install.packages('devtools')
devtools::install_github('phytomosaic/fitNMDS')
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
plot_boot_spider(res, col='#00000040')
```


### Reciprocal NMDS of both datasets
```r
r <- recip_nmds(tw)
summary(r)
plot(r)
```


