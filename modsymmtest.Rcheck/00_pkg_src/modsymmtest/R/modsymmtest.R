#' @title
#' Wilcoxon and Sign tests for symmetry about an unknown center
#'
#' @description
#' R built-in function `wilcox.test()` is designed to perform both the one- and two-sample Wilcoxon test for symmetry 
#' under the assumption that the center of symmetry is specified. The procedure `mod.symm.test()` extends the capabilities of 
#' `wilcox.test()` for situations where the center of symmetry is unknown. Such cases can be found in e.g., regression residuals evaluations, 
#' as well as in the book `Nonparametric statistical methods using R` by Kloke and McKean, and in the scholarly work of Gastwirth.
#'
#'
#' @param x numeric vector of data values. Non-finite (e.g., infinite or missing) values will be omitted.
#' @param y an optional numeric vector of data values: as with x non-finite values will be omitted.
#' @param alternative
#'   a character string specifying the alternative hypothesis, must be one of "two sided" (default), "right.skewed", or "left.skewed".
#'   You can specify just the initial letter.
#'   "right.skewed": test whether positively skewed, 
#'   "left.skewed" : test whether negatively skewed.
#' @param method a character string specifying which symmetry test to be used, "wilcox" refers to Wilcoxon signed-rank test,
#'   and "sign" is sign test.
#'   
#'
#' @details 
#' When "wilcox", the default method, is used, the test statistic W has the form 
#' of the Wilcoxon test statistic with unknown center of symmetry replaced by its sample estimator.
#' For details, see Vexler et al. (2023)
#' 
#' When method="sign", the test statistic S is the total number of 
#' the observations that smaller than sample mean.For more details, see Gastwitrh (1971).    
#'   
#' @returns A list of class "htest" containing the following components:
#'   \itemize{
#'   \item statistic - the value of the test statistic.
#'   \item var - the asymptotic variance.
#'   \item alternative - a character string describing the alternative hypothesis.
#'   \item p.value - the p-value for the test.
#'   \item method - the type of test applied.
#'
#'   }
#'  
#'   
#' @author Jiaojiao Zhou, Xinyu Gao, Albert Vexler
#' 
#' @references
#' Vexler, A., Gao, X., & Zhou, J. (2023). How to implement signed-rank `wilcox.test()` type procedures when a center of symmetry is unknown. Computational Statistics & Data Analysis, 107746. 
#' 
#' Gastwirth, J. L. (1971). On the Sign Test for Symmetry. Journal of the American Statistical Association, 66(336), 821-823.
#' 
#'  
#' @importFrom stats pnorm wilcox.test var setNames complete.cases
#' @importFrom Rcpp evalCpp 
#' @useDynLib modsymmtest
#' @examples 
#' \dontrun{
#'   # A study measures the plasma silicon levels before and after silicone implants surgery
#'   # in 30 women to evaluate the effect of the surgery. Informally speaking, we can be interested
#'   # in that there is an unknown constant shift such that the the plasma silicon level of post-surgery
#'   # can be explained completely based on that of pre-surgery. This can be stated as the null hypothesis
#'   # `H_0` The difference of plasma silicon level between post-surgery pre-surgery has a symmetric distribution
#'   # around a shift that is unknown.  
#'   data("plasma.silicon")
#'   post <- plasma.silicon$postoperative 
#'   pre <- plasma.silicon$preoperative
#'   # post <- c(0.21,0.24,0.1,0.12,0.28,0.25,0.22,0.21,0.22,0.23,0.22,0.24,0.45,0.38,
#'   #           0.23,0.22,0.18,0.15,0.04,0.14,0.24,0.2,0.24,0.18,0.19,0.15,0.26,0.3,0.22,0.24)
#'   # pre <- c(0.15,0.13,0.39,0.2,0.39,0.42,0.24,0.18,0.26,0.12,0.1,0.11,0.19, 0.15,0.27,
#'   #          0.28,0.11,0.11,0.18,0.18,0.24,0.48,0.27,0.22,0.18,0.19,0.32,0.31,0.19,0.21)
#'   mod.symm.test(x=post, y=pre, alternative ="two.sided", method = "wilcox")
#'   
#'   Result:
#'   Modified Wilcoxon signed-rank test
#'   data:  post and pre
#'   W = 238, p-value = 0.767
#'   alternative hypothesis: two.sided
#'   
#'   Interpretation:
#'   Test statistic `W` is the number of walsh average higher than sample mean, see more details 
#'   in paper authored by Vexler, etc. 
#'   p-value is 0.767, which implies there is no clue to reject the null hypothesis that
#'   the distribution of the difference of plasma silicon levels before and after 
#'   silicone implants surgery is symmetric. 
#' }
#' @export
mod.symm.test <- function(x, y=NULL,
                          alternative = c("two.sided", "left.skewed", "right.skewed"),
                          method = "wilcox")
{
  alternative <- match.arg(alternative)
  # alt.text <- switch (alternative,
  #                     "two.sided" = "not symmetric",
  #                     "less" = "left-skewed",
  #                     "greater" = "right-skewed" )
  if(!is.numeric(x)) stop("'x' must be numeric")
  if(!is.null(y)) {
    if(!is.numeric(y)) stop("'y' must be numeric")
    DNAME <- paste(deparse(substitute(x)), "and",
                   deparse(substitute(y)))
    
    if(length(x) != length(y))
      stop("'x' and 'y' must have the same length")
    OK <- complete.cases(x, y)
    x <- x[OK] - y[OK]
    y <- NULL
  }else {
    DNAME <- deparse(substitute(x))
    x <- x[is.finite(x)]
  }
  
  if(length(x) < 1L)
    stop("not enough (finite) 'x' observations")
  n <- as.double(length(x))
  m <- mean(x)
  sigma2 <- var(x)
  xc <- x-m
  
  
  if (method == "wilcox"){
    METHOD <- "Modified Wilcoxon signed-rank test"
    
    ## Test statistic
    # STATISTIC <- setNames(as.numeric(wilcox.test(x,mu=m)$statistic), "W")
    
    xc_abs <- abs(xc)
    xc_abs_rank <- rank(xc_abs)
    STATISTIC <- sum(xc_abs_rank*(xc>0))
    STATISTIC <- setNames(STATISTIC, "W")
    
    # Asymptotic mean and variance
    E <- n*(n+1)/4
    
    ## call get_V function
    V = getV(x)
    
    ## The resulting p-value
    pval <- switch (alternative,
                    "two.sided" =  2 * ( 1 - pnorm(abs(E-STATISTIC)/sqrt(V))),
                    "right.skewed" = 1 - pnorm( (E-STATISTIC)/sqrt(V) ),
                    "left.skewed" = pnorm( (E-STATISTIC)/sqrt(V))
    )
    
  } else if (method == "sign") {
    METHOD <- "Modified sign test"
    
    ## Test statistic
    STATISTIC <- setNames(sum(1*(x<m)), "S")
    
    ## Estimation of w
    a <- 1 * (x>m-n^(-1/5)) * (x < m+n^(-1/5))
    D <- sum(a)
    A <- max(c(1,D))
    VHW <- (n^(3/10)/A)^2
    hat_w <- sqrt(1/(4*n*VHW))
    E <- n/2
    
    ## probability weight moment for x
    CE <- mean(x*(x<m)) - m/2
    
    ## Asymptotic variance
    V <- 1/4 + var(x)*(hat_w)^2 + 2 *hat_w*CE
    
    pval <- switch (alternative,
                    "two.sided" =  2 * ( 1 - pnorm(abs(STATISTIC-E)/sqrt(n*V))),
                    "right.skewed" = 1 - pnorm((STATISTIC-E)/sqrt(n*V)),
                    "left.skewed" = 1 - pnorm( -(STATISTIC-E)/sqrt(n*V) )
    )
    
  }
  
  
  pval <- setNames(pval, "p.value")
  RVAL <- list(statistic = STATISTIC,
               var = V,
               p.value = as.numeric(pval),
               method = METHOD,
               alternative = alternative,
               data.name = DNAME )
  class(RVAL) <- "htest"
  RVAL
  
  
}