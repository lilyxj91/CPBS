
# cpbs

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

This is a basic example which shows the results of cpbs and a naive
method.

``` r
library(cpbs)
## basic example code
```

## Reference

Piao, J., Ning, J., Chambers, C. D., & Xu, R. (2018). Semiparametric
model and inference for spontaneous abortion data with a cured
proportion and biased sampling. Biostatistics (Oxford, England), 19(1),
54–70. <https://doi.org/10.1093/biostatistics/kxx024>
