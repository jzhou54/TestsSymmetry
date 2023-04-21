#' Title
#' @title Sample size calculation for Wilcoxon test of symmetry with estimated location 
#' @description 
#' Determine the sample size required for one-sample wilcoxon signed-rank test or sign test 
#' with unknown center of symmetry estimated by sample mean.
#' 
#' @details 
#' 
#' 
#' @param x pilot data, numeric vector of data values. Non-finite (e.g., infinite or missing) values will be omitted.
#' @param sig.level Significance level (Type I error probability)
#' @param power Power of test (1 minus Type II error probability)
#' @param method a character string specifying which symmetry test to be used, "wilcox" refers to Wilcoxon signed-rank test,
#'   and "sign" is sign test.
#' @param alternative
#'   a character string specifying the alternative hypothesis, must be one of "two sided" (default), "right.skewed", or "left.skewed".
#'   You can specify just the initial letter.
#'   "right.skewed": test whether positively skewed, 
#'   "left.skewed" : test whether negatively skewed.
#'   
#' @returns A list of class "power.htest" containing the following components:
#'    \itemize{
#'   \item n - sample size 
#'   \item sig.level - Significance level (Type I error probability)
#'   \item power - Power of test (1 minus Type II error probability)
#'   \item method - the type of test applied.
#'   }
#' @export
#' @importFrom MASS boxcox
#' @importFrom stats qnorm lm var optimize 
#' 
#' @references{
#' Vexler, A., Gao, X., & Zhou, J. (2023). How to implement signed-rank wilcox. test () type procedures when a center of symmetry is unknown. Computational Statistics & Data Analysis, 107746. }
#'
#' @examples
#' x <- rgamma(30, shape=2, scale=2)
#' pwr.symm.test(x, sig.level = 0.05, power = 0.8, method="wilcox")
#' 
pwr.symm.test<- function(x, sig.level = 0.05, power = 0.8, method="wilcox",
                         alternative = c("two.sided", "right.skewed", "left.skewed"),
                         plot.extraplation = FALSE){
  
  if(!is.numeric(x)) stop("'x' must be numeric")
  if(length(x) < 1L)
    stop("not enough (finite) 'x' observations")
  
  DNAME <- deparse(substitute(x))
  alternative <- match.arg(alternative)
  x <- x[is.finite(x)]
  m <- mean(x)
  n <- length(x)
  sigma2.x <- var(x)
  z_alpha <- qnorm(1-sig.level/2)
  beta = 1 - power
  
  if (method == "wilcox"){
    METHOD <- "Modified Wilcoxon signed-rank test"
   
    ########## estimates under H0 ##########
    
    theta0 <- get_quant_H0(x)$theta0
    tau0 <- get_quant_H0(x)$tau0
    mu0 <- function(N) N*(N+1)/4
    sigma0 <- function(N) sqrt(N*(N+1)*(2*N+1)/24 - N*(N-1)*(N-3) * theta0 * tau0 +
                                 (N-1)*(N-2)*(N-3)*(N-4)*sigma2.x/(4*N)*(theta0)^2)
    
    ##### estimates under H1 #####
    
    est_quant_H1 <- get_quant_H1(x)
    qx <- est_quant_H1$qx
    qz <- est_quant_H1$qz
    L1 <- est_quant_H1$L1
    L2 <- est_quant_H1$L2
    theta1 <- est_quant_H1$theta1
    tau1 <- est_quant_H1$tau1
    # print(data.frame("theta1"=theta1, "tau1"=tau1, "sigma2.x"=sigma2.x))
    
    mu1 <- function(N) N*(N-1)/2*qz + N*qx
    sigma1 <- function(N) {
      Kn <- (1 - 2*qz - qz^2 - 2*L1) * N^3 # +
        # (5/2*qz^2 +13/2*qz - 2*qz*qx + 2*(1 - qx) + 6*L1 - 2*L2 - 3) * N^2 +
        # (2*qx*qz -3/2*qz^2 - qx^2 - 9/2*qz + 3*qx - 4*L1 + 2*L2) * N
      
      sig1.square <- Kn - N*(N-1)*(N-3) * theta1 * tau1 +
        (N-1)*(N-2)*(N-3)*(N-4)*sigma2.x/(4*N)*(theta1)^2
      # sigma1.value <- ifelse(sig1.square < 0, sigma0(N), sqrt(sig1.square))
      
      if (sig1.square > 0) {
        sigma1.value <- sqrt(sig1.square)
      }else{
        # print("sig1.square<=0")
        sigma1.value <- sqrt(N*(N+1)*(2*N+1)/24 - N*(N-1)*(N-3) * theta0 * tau0 +
                               (N-1)*(N-2)*(N-3)*(N-4)*sigma2.x/(4*N)*(theta0)^2)
      }
      
      return(sigma1.value)
    }
    
    sigma1 <- Vectorize(sigma1, "N")
    
    
    # Normal approximation
    z1 <- function(N) (mu0(N)-mu1(N)+sigma0(N)*z_alpha)/sigma1(N)
    z2 <- function(N) (mu0(N)-mu1(N)-sigma0(N)*z_alpha)/sigma1(N)  
    z3 <- function(N) (mu0(N)-mu1(N)+sigma0(N)*qnorm(1-sig.level))/sigma1(N) # one-sided
    z4 <- function(N) (mu0(N)-mu1(N)-sigma0(N)*qnorm(1-sig.level))/sigma1(N) # one-sided
    
    if (alternative == "two.sided"){
      # f <- function(N)  (1/2*( sqrt(1-exp(-2/pi*z1(N)^2)) - sqrt(1-exp(-2/pi*z2(N)^2)) ) - beta)^2 # Polya approximation to standard normal
      f <- function(N) (pnorm(z1(N)) - pnorm(z2(N))-beta)^2
      power.fn <-  function(N) pnorm(z1(N)) - pnorm(z2(N))-beta
    }else if(alternative == "right.skewed"){
      f <- function(N)  (z3(N) - qnorm(beta) )^2 # one-sided
      power.fn <-  function(N) z3(N) - qnorm(beta)
    }else if (alternative == "left.skewed"){
      f <- function(N) (z4(N) - qnorm(1 - beta))^2
      power.fn <-  function(N) z4(N) - qnorm(1 - beta)
    }
    
    f <- Vectorize(f, "N")
    power.fn <- Vectorize(power.fn, "N")
    
    getN <-  optimize(f, c(n, 600), maximum = F)$minimum
    getN.int <- ceiling(getN)
    
    # getN.initial <-  getN.int
    # DeltaN <- 10
    # findN <- FALSE
    # 
    # while (getN.initial - DeltaN > 0 && getN.initial + DeltaN < 1000 && findN == FALSE) {  # n = n.pilot
    #   power.fn.left <- power.fn(getN.initial - DeltaN)
    #   power.fn.right  <- power.fn(getN.initial + DeltaN)
    #   
    #   if (power.fn.left * power.fn.right < 0) {
    #     getN.int.res <- uniroot(power.fn, lower = getN.initial - DeltaN, upper=getN.initial + DeltaN )
    #     getN.int <- ceiling(getN.int.res$root)
    #     findN = TRUE
    #   }else{
    #     DeltaN = DeltaN + 10
    #   }
    # }
    
    
  }else if (method == "sign"){
    METHOD <- "Modified sign test"
    
    }
  
  ## return list of 'power.htest' class
  RVAL <- list(N = getN.int,
               sig.level = sig.level,
               power = power,
               method = METHOD,
               alternative = alternative,
               data.name = DNAME
  )
  ### extraplation plot ###
  if (plot.extraplation == TRUE){
    MC <- 1000
    n.seq <- seq(5, n, by=5)
    p.val <- matrix(NA_real_, nrow = MC, ncol = length(n.seq))
    colnames(p.val) <- paste0("n",n.seq)
    
    ind.col <- 1
    for (n in n.seq){
      for (mc in 1:MC){
        set.seed(mc+n+20230410)
        x.sub <- sample(x, size=n, replace = T)
        tryCatch({
          p.val[mc, ind.col] <- suppressWarnings(mod.symm.test(x.sub)$p.value)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        
      }
      ind.col <- ind.col + 1
    }
    
    p.avg <- colMeans(p.val, na.rm = T) 
    
    ## fit model ##
    mod <- lm(log(p.avg)~n.seq)
    extrap.n <- 100
    n.seq.extrap <- seq(from=max(n.seq), to=extrap.n)
    y.predict <- predict(mod, data.frame(x=n.seq))
    y.predict.extrap <- exp(mod$coefficients[[1]] + mod$coefficients[[2]]*n.seq.extrap)
    
    plot(n.seq, p.avg, pch=16, ylim=c(0,1), xlim=c(5, extrap.n))
    lines(n.seq, exp(y.predict), col='blue')
    lines(n.seq.extrap, y.predict.extrap, col="blue", lty=2)
    abline(h = 0.05, col = "red", lty=2)
    
    N.extraplation <- ceiling((log(0.05) - mod$coefficients[[1]] )/mod$coefficients[[2]])
    RVAL <- c(RVAL, list("N.extraplation" = N.extraplation))
  }
  
  ######################################
 
  class(RVAL) <- "power.htest"
  RVAL
}
