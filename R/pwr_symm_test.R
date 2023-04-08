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
#' 
pwr.symm.test<- function(x, sig.level = 0.05, power = 0.8, method="wilcox"){
  
  if(!is.numeric(x)) stop("'x' must be numeric")
  if(length(x) < 1L)
    stop("not enough (finite) 'x' observations")
  
  DNAME <- deparse(substitute(x))
  x <- x[is.finite(x)]
  m <- mean(x)
  n <- length(x)
  sigma2.x <- var(x)
  z_alpha <- qnorm(1-sig.level/2)
  beta = 1 - power
  
  if (method == "wilcox"){
    METHOD <- "Modified Wilcoxon signed-rank test"
   
    ########## estimates under H0 ##########
    # Box-cox transformation
    b <- boxcox(lm(x ~ 1), plotit = F)
    lambda <- b$x[which.max(b$y)]
    x.sym = (x^lambda - 1)/lambda

    
    # a simple method
    # x.sym <- rnorm(n, mean = m, sd=sqrt(sigma2.x))
    # std = sqrt(sigma2.x)
    # theta0 <- dnorm(0, 0, sd = std/sqrt(2))
    # tau0 <- std / (2*sqrt(pi))
    
    sigma2.x.sym <- var(x.sym)
    
    theta0 <- get_quant_H0(x.sym)$theta0
    tau0 <- get_quant_H0(x.sym)$tau0
    mu0 <- function(N) N*(N+1)/4
    # sigma0 <- function(N) sqrt(N*(N+1)*(2*N+1)/24 - N*(N-1)*(N-3) * theta0 * tau0 + 
    #                              (N-1)*(N-2)*(N-3)*(N-4)*std^2/(4*N)*(theta0)^2)
    sigma0 <- function(N) sqrt(N*(N+1)*(2*N+1)/24 - N*(N-1)*(N-3) * theta0 * tau0 +
                                 (N-1)*(N-2)*(N-3)*(N-4)*sigma2.x.sym/(4*N)*(theta0)^2)
    

    # print(data.frame("theta0"=theta0, "tau0"=tau0, "sigma2.x"=std^2))
    
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
                               (N-1)*(N-2)*(N-3)*(N-4)*sigma2.x.sym/(4*N)*(theta0)^2)
      }
      
      return(sigma1.value)
    }
    
    sigma1 <- Vectorize(sigma1, "N")
    
    
    # Normal approximation
    z1 <- function(N) (mu0(N)-mu1(N)+sigma0(N)*z_alpha)/sigma1(N)
    z2 <- function(N) (mu0(N)-mu1(N)-sigma0(N)*z_alpha)/sigma1(N)  
    
    # solve N for function Phi(z1)-Phi(z2)-beta = 0
    # f <- function(N)  (1/2*( sqrt(1-exp(-2/pi*z1(N)^2)) - sqrt(1-exp(-2/pi*z2(N)^2)) ) - beta)^2 # Polya approximation to standard normal
    f <- function(N) (pnorm(z1(N)) - pnorm(z2(N))-beta)^2
    f <- Vectorize(f, "N")
    
    getN <-  optimize(f, c(n, 10*n), maximum = F)$minimum
    getN.int <- ceiling(getN)
    
  }else if (method == "sign"){
    METHOD <- "Modified sign test"
    
    }
  
  
  ######################################
  ## return list of 'power.htest' class
  RVAL <- list(N = getN.int,
               sig.level = sig.level,
               power = power,
               method = METHOD,
               data.name = DNAME
               )
  class(RVAL) <- "power.htest"
  RVAL
}
