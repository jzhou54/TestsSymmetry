# modsymmtest

Modified symmetry test

Provide functionality to test the symmetry characteristic of a dataset, or paired datasets.

To install this package, use the statement below in R:

 ```r
 devtools::install_github("jzhou54/modsymmtest")
 
 library(modsymmtest)
 ```
 
 The function is 
 
 ```r
  mod.symm.test(x, y=NULL, alternative="two.sided", method="wilcox")
 ```
 
 In this function, two methods are incorporated, one is modified wilcoxon rank test and the other is modified sign test. 
