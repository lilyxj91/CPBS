
# Cured Proportion and Biased Sampling

<!-- badges: start -->
<!-- badges: end -->

Evaluating and understanding the risk and safety of using medications
for autoimmune disease in a woman during her pregnancy will help both
clinicians and pregnant women to make better treatment decisions.
However, utilizing spontaneous abortion (SAB) data collected in
observational studies of pregnancy to derive valid inference poses two
major challenges. First, the data from the observational cohort are not
random samples of the target population due to the sampling mechanism.
Pregnant women with early SAB are more likely to be excluded from the
cohort, and there may be substantial differences between the observed
SAB time and those in the target population. Second, the observed data
are heterogeneous and contain a “cured” proportion. In this article, we
consider semiparametric models to simultaneously estimate the
probability of being cured and the distribution of time to SAB for the
uncured subgroup. To derive the maximum likelihood estimators, we
appropriately adjust the sampling bias in the likelihood function and
develop an expectation-maximization algorithm to overcome the
computational challenge. We apply the empirical process theory to prove
the consistency and asymptotic normality of the estimators. We examine
the finite sample performance of the proposed estimators in simulation
studies and illustrate the proposed method through an application to SAB
data from pregnant women. \[1\]

## Installation

You can install the development version of cpbs like so:

``` r
library(devtools)
install_github("lilyxj91/cpbs",force = TRUE)
library(cpbs)
```

## Example

R file contians functions used to implement the proposed method and
naive method.

``` r
library(cpbs)
## basic example code
data(example_data)
data=data.matrix(example_data)
data = data[order(data[,1]),]
output=cpbs_em(tol=1e-4,maxit=1000,data)
#> 
#> Successfully Converged
output$theta
#> [1] 1.017538 2.598353
output$alpha
#> [1] 1.535884 1.011702 1.206555
output$beta
#>       cov3       cov4 
#> -0.5696793  1.4060352

output.naive = naive.method(T= data[,1],IND = data[,2],
                            X = data[,3:4],A = data[,5],Y_SAB = data[,6])
output.naive
#> $alpha
#>          (Intercept) X[!is.na(Y_SAB), ]z1 X[!is.na(Y_SAB), ]z2 
#>            2.3445080            0.8438043            0.7968939 
#> 
#> $alpha.se
#>          (Intercept) X[!is.na(Y_SAB), ]z1 X[!is.na(Y_SAB), ]z2 
#>            0.2532021            0.5046373            0.7859039 
#> 
#> $beta
#> X[Y_SAB == 1, ]z1 X[Y_SAB == 1, ]z2 
#>        -0.4501259         0.9851500 
#> 
#> $beta.se
#> X[Y_SAB == 1, ]z1 X[Y_SAB == 1, ]z2 
#>         0.1491193         0.2553199 
#> 
#> $theta
#> [1] 0.9670423 2.1970326
#> 
#> $theta.se
#> [1] 0.04316005 0.13828019
```

## Reference

Piao, J., Ning, J., Chambers, C. D., & Xu, R. (2018). Semiparametric
model and inference for spontaneous abortion data with a cured
proportion and biased sampling. Biostatistics (Oxford, England), 19(1),
54–70. <https://doi.org/10.1093/biostatistics/kxx024>
