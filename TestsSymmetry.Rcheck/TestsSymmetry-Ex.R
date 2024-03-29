pkgname <- "TestsSymmetry"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "TestsSymmetry-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('TestsSymmetry')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("mod.symm.test")
### * mod.symm.test

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mod.symm.test
### Title: Wilcoxon and Sign tests for symmetry about an unknown center
### Aliases: mod.symm.test

### ** Examples

## No test: 
  # A study measures the plasma silicon levels before and after silicone implants surgery in 30
  # women to evaluate the effect of the surgery. Informally speaking, we can be interested
  # in that there is an unknown constant shift such that the the plasma silicon level of 
  # post-surgery can be explained completely based on that of pre-surgery. This can be stated 
  # as the null hypothesis `H_0` The difference of plasma silicon level between post-surgery and 
  # pre-surgery has a symmetric distribution around a shift that is unknown.  
  data("plasma.silicon")
  post <- plasma.silicon$postoperative 
  pre <- plasma.silicon$preoperative
  # post <- c(0.21,0.24,0.1,0.12,0.28,0.25,0.22,0.21,0.22,0.23,0.22,0.24,0.45,0.38,
  #           0.23,0.22,0.18,0.15,0.04,0.14,0.24,0.2,0.24,0.18,0.19,0.15,0.26,0.3,0.22,0.24)
  # pre <- c(0.15,0.13,0.39,0.2,0.39,0.42,0.24,0.18,0.26,0.12,0.1,0.11,0.19, 0.15,0.27,
  #          0.28,0.11,0.11,0.18,0.18,0.24,0.48,0.27,0.22,0.18,0.19,0.32,0.31,0.19,0.21)
  mod.symm.test(x=post, y=pre, alternative ="two.sided", method = "wilcox")
  
  # Result:
  # Modified Wilcoxon signed-rank test
  # data:  post and pre
  # W = 238, p-value = 0.767
  # alternative hypothesis: two.sided
  
  # Interpretation:
  # Test statistic `W` is the number of walsh average higher than sample mean, see more details 
  # in paper authored by Vexler, etc. 
  # p-value is 0.767, which implies there is no clue to reject the null hypothesis that
  # the distribution of the difference of plasma silicon levels before and after 
  # silicone implants surgery is symmetric. 
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mod.symm.test", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("n.symm.test")
### * n.symm.test

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: n.symm.test
### Title: Sample size determination for nonparametric tests of symmetry
###   when the center is unknown
### Aliases: n.symm.test

### ** Examples

## No test: 
 data("plasma.silicon")
 post <- plasma.silicon$postoperative 
 pre <- plasma.silicon$preoperative
 diff <- post - pre
 n.symm.test(diff, sig.level = 0.05, power = 0.5, method = "wilcox", alternative ="two.sided" )

# Result:
# Sample size calculation under wilcox procedure 

#           N = 83
#   sig.level = 0.05
#       power = 0.5
#        type = wilcox
# alternative = two.sided

# Interpretation: 
# Given the pilot sample `diff` and the significance level 0.05. The sample size of 
# the data that is expected toprovide the target power 0.5 of the Wilcoxon test procedure 
# is computed as 83.
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("n.symm.test", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
