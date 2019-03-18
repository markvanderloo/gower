
[![Build Status](https://travis-ci.org/markvanderloo/gower.svg?branch=master)](https://travis-ci.org/markvanderloo/gower)
[![Coverage Status](https://coveralls.io/repos/markvanderloo/gower/badge.svg)](https://coveralls.io/r/markvanderloo/gower) 
[![CRAN](http://www.r-pkg.org/badges/version/gower)](https://CRAN.R-project.org/package=gower)
[![status](https://tinyverse.netlify.com/badge/gower)](https://CRAN.R-project.org/package=gower)
[![Downloads](http://cranlogs.r-pkg.org/badges/gower)](http://cran.r-project.org/package=gower) 


### gower

Gower's distance for R. Based in C, using openMP for parallelization.

### Usage

```
library(gower)
reviris <- iris[rev(seq_len(nrow(iris))),,drop=FALSE]
# compute distances
d <- gower_dist(iris, reviris)

# data.frame with less records is recycled
d <- gower_dist(iris[1:3,,drop=FALSE], reviris)

# compute top-n matches
mat <- gower_topn(iris, reviris, n=5)

## mat$index   : Array of indices (sorted from better to worse match)
## mat$distance: Array of distances (sorted from small to large)

```

More info in the [vignette](https://github.com/markvanderloo/gower/blob/master/pkg/vignettes/intro.Rmd)


### Installation

#### From CRAN

```
install.packages("gower")
```


#### Beta versions

Made available through my drat repo.

First, install the [drat](https://cran.r-project.org/package=drat) package. Users of the OS who's name shall not be mentioned need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first.
```
if(!require(drat)) install.packages('drat')
drat::addRepo('markvanderloo')
install.packages('gower',type='source')

```


### Reference

Gower (1971) A general coefficient of similarity and some of its properties. _Biometrics_ **27** 857-874 [pdf](http://venus.unive.it/romanaz/modstat_ba/gowdis.pdf)
