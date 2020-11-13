<!-- PRIME -->

# Projective Resampling Imputation Estimation Method for Missing Covariates Problems

In this work, motivated by the spirit of projection resampling/random projection, we propose a Projection
Resampling Imputation Mean Estimation (PRIME) to project the covariates along randomly sampled
directions to obtain samples of scalar-valued predictors (dimension reduction). And then a simple geometric
average is taken on the scalar predictors to impute the missing parts (using all-sided information).
In general, several advantages of our method can be concluded as follows. First, the PRIME is capable
of dealing with missing data containing no complete observations, while most existed methods require
a fraction of subjects to have fully completed observations. Second, we can integrate the averaging of
the imputation estimation of multiple projection directions and fully utilize the available information to
reach a more reliable and insightful result.


##  Reproduce code for case1 in our paper

- `help_func.R` includes some helpful functions about PRIME 

- `example.R` is the reproducing code for case1 in our paper

```r
#you should install some depedency packages, such as MASS,misaem,ggplot2,gridExtra
source("./help_func.R")
source("./example.R")
```

# Qustions

Please contact [13023988263@163.com](mailto:13023988263@163.com) or [ele717@163.com](mailto:ele717@163.com) with any questions about the PRIME method.
