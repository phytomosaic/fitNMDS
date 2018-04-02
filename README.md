Using the R package fitNMDS
================
Robert J. Smith
31 March 2018

Motivation
==========

Combining two different datasets into one nonmetric multidimensional scaling (NMDS) model can be risky if they each cover different attribute spaces (e.g., different species pools in ecology). Therefore, comparing two datasets requires estimating internal sampling variability (using bootstrapped NMDS) relative to external exchangeability (using reciprocal NMDS).

Installation
============

Install the package from github as follows:

``` r
## install.packages('devtools')
## devtools::install_github('phytomosaic/fitNMDS')
require(fitNMDS)
```

    Loading required package: fitNMDS

Usage
=====

Prepare two candidate datasets
------------------------------

``` r
# here we just modify one dataset by adding noise to create a 'second'
set.seed(231)
data(smoky)
spe1 <- smoky$spe
env1 <- env2 <- smoky$env
spe2 <- spe1 + abs(rnorm(prod(dim(spe1)), 0, 2))  # add noise
tw   <- twin(spe1, spe2, env1, env2)
```

Bootstrapped NMDS of one dataset
--------------------------------

``` r
x   <- list(spe=spe1, id=env1)
res <- boot_nmds(x, B=29, k=2, rot=TRUE)
```


    time elapsed for bootstrap: 0.08252718 minutes

``` r
summary(res)
```

          SRV rP_median 
        0.981     0.967 

``` r
plot(res, col='#00000040')
```

![Bootstrapped NMDS.](Appendix_S3_tutorial_fitNMDS_files/figure-markdown_github/bootstrap-1.png)

Reciprocal NMDS of both datasets
--------------------------------

``` r
res <- recip_nmds(tw)
summary(res)
```

       d1_rP    d2_rP M1vM2_rP  stress1  stress2  varexp1  varexp2 
       0.957    0.942    0.949    0.172    0.166    0.784    0.759 

``` r
plot(res, noaxes=FALSE)
```

![Reciprocal NMDS.](Appendix_S3_tutorial_fitNMDS_files/figure-markdown_github/reciprocal-1.png)

### Correct unequal sample sizes among two datasets

``` r
# two candidate datasets, full and partial
spe2 <- spe1[1:11,]
env2 <- env1[1:11,]

# subset the full matrix, based on compositional nearest neighbors...
(i   <- nearestspecies(spe1, spe2, ties=FALSE))
```

    _13 _14 _15 _16 _17 _18 _19 _20 _21 _22 _23 
      1   2   3   4   5   6   7   8   9  10  11 

``` r
spe1 <- spe1[i,,]
env1 <- env1[i,,]

# ...then proceed to reciprocal NMDS
tw   <- twin(spe1, spe2, env1, env2)
res  <- recip_nmds(tw)
summary(res)
```

       d1_rP    d2_rP M1vM2_rP  stress1  stress2  varexp1  varexp2 
       1.000    1.000    1.000    0.048    0.048    0.896    0.896 

``` r
plot(res, noaxes=FALSE)
```

![Corrected reciprocal NMDS.](Appendix_S3_tutorial_fitNMDS_files/figure-markdown_github/corr_uneq-1.png)
