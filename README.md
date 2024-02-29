# TestsSymmetry

Tests for symmetry when the center of symmetry is unknown

Provides functionality of implementation of statistical testing whether a dataset comes from a symmetric distribution when the center of symmetry is unknown, including Wilcoxon test and sign test procedure. In addition, sample size determination for both tests is provided. 

To install this package, use the statement below in R:

 ```r
 devtools::install_github("jzhou54/TestsSymmetry")
 # or install.packages("TestsSymmetry")
 
 library(TestsSymmetry)
 ```
 
 Here are some examples for the functions
 
 1. Symmetry tests
 
 ```r
  x <- rchisq(50, df = 5)  # asymmtric case
  y <- rnorm(n=50)         # symmetric case
  mod.symm.test(x, alternative="two.sided", method="wilcox")
  mod.symm.test(y, alternative="two.sided", method="wilcox")
 ```
 
 2. Sample size determination
 
 ```r
  x <- rchisq(30, df = 5)  # asymmtric case
  n.symm.test(x, power = 0.6, method = "wilcox")
  n.symm.test(x, power = 0.6, method = "sign")

 ```
 
 