# TestsSymmetry

Tests for symmetry when the center of symmetry is unknown

Provides functionality of implementation of statistical testing whether a dataset comes from a symmetric distribution when the center of symmetry is unknown, including wilcoxon test and sign test procedure. In addition, sample size determination for both tests is provided. 

To install this package, use the statement below in R:

 ```r
 devtools::install_github("jzhou54/TestsSymmetry")
 
 library(TestsSymmetry)
 ```
 
 The function is 
 
 ```r
  mod.symm.test(x, y=NULL, alternative="two.sided", method="wilcox")
  n.symm.test(diff, power = 0.6, method = "wilcox")
 ```
 
 In this function, two methods are incorporated, one is modified wilcoxon rank test and the other is modified sign test. 
