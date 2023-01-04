#' @title
#' Modified symmetry test when the center of symmetry is unknown
#'
#' @description
#'   Provides one unified symmetry test that incorporates two methods, "wilcox" and "sign". Also, paired data can be tested.
#'
#'
#' @param x data set to be tested
#' @param y another data set, default is NULL
#' @param alternative
#'   a character string specifying the alternative hypothesis, must be one of "two sided" (default), "greater", or "less".
#'   greater: test whether positively skewed (right-skewed),
#'   less : test whether negatively skewed (left-skewed)
#'
#'
#' @param method a character string specifying which symmetry test to be used, "wilcox" is modified wilcoxon sign rank test,
#'   and "sign" is modified sign test.
#' @returns A list of class "htest" containing the following components:
#'   \itemize{
#'   \item statistic - the value of the test statistic.
#'   \item var - the asymptotic variance.
#'   \item alternative - a character string describing the alternative hypothesis.
#'   \item p.value - the p-value for the test.
#'   \item method - the type of test applied.
#'
#'   }
#' @importFrom stats pnorm wilcox.test var setNames complete.cases
#' @importFrom Rcpp evalCpp
#' @useDynLib modsymmtest
#' @examples 
#' \dontrun{
#'   x <- rnorm(100, 1, 1)
#'   y <- rnorm(100)
#'   mod.symm.test(x)
#'   mod.symm.test(x,y)
#' 
#' }
#' @export
mod.symm.test <- function(x, y=NULL,
                          alternative = c("two.sided", "less", "greater"),
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
    METHOD <- "Modified wilcoxon signed rank test"
    
    ## Test statistic
    STATISTIC <- setNames(as.numeric(wilcox.test(x,mu=m)$statistic), "W")
    
    # ## Tn selection
    # r <- quantile(x, c(0.25, 0.75))
    # h <- (r[2] - r[1])/1.34
    # Tn<-log(n)/( 3 * 1.06 * min(sqrt(var(x)), h))
    #
    # ## Estimation of theta
    # S <- function(u) sum( sin(2*pi*(xc[xc!=u]-u)*Tn)/(2*pi*(xc[xc!=u]-u)))+
    #   sum(sin(2*pi*(xc+u)*Tn)/(2*pi*(xc+u)))
    # SV <- Vectorize(S)
    # hat_theta <- 2*sum(SV(xc))/n^2 + 2*Tn/n
    #
    # ## Estimation of tau
    # xs <- sort(xc)  # V(i)
    # S1 <- seq(from=1, to=n, by=1) # i
    # hat_tau <- sum(xs*S1)/n^2
    #
    # ## Asymptotic mean and variance
    E <- n*(n+1)/4
    # V <- n*(n+1)*(2*n+1)/24 - n*(n-1)*(n-3) * hat_theta * hat_tau +
    #   (n-1)*(n-2)*(n-3)*(n-4)*sigma2/(4*n)*(hat_theta)^2
    
    ## call get_V function
    V = getV(x)
    
    ## The resulting p-value
    pval <- switch (alternative,
                    "two.sided" =  2 * ( 1 - pnorm(abs(E-STATISTIC)/sqrt(V))),
                    "greater" = 1 - pnorm( (STATISTIC-E)/sqrt(V) ),
                    "less" = pnorm( (STATISTIC-E)/sqrt(V))
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
    CE <- mean((x-m)*(x<m))
    
    ## Asymptotic variance
    V <- 1/4 + var(x)*(hat_w)^2 + 2 *hat_w*CE
    
    pval <- switch (alternative,
                    "two.sided" =  2 * ( 1 - pnorm(abs(STATISTIC-E)/sqrt(n*V))),
                    "greater" = 1 - pnorm( -(STATISTIC-E)/sqrt(n*V) ),
                    "less" = 1-pnorm((STATISTIC-E)/sqrt(n*V))
    )
    
  }
  
  
  pval <- setNames(pval, "p.value")
  # pval <- 2*(1-pnorm(abs(STATISTIC-n/2)/sqrt(n*V)))
  RVAL <- list(statistic = STATISTIC,
               var = V,
               p.value = as.numeric(pval),
               method = METHOD,
               alternative = alternative,
               data.name = DNAME )
  class(RVAL) <- "htest"
  RVAL
  
  
}