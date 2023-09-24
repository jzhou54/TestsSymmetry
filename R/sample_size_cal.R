#' Title
#' @title Sample size determination for a nonparametric tests of symmetry
#' @description Determine the sample size required for one-sample wilcoxon signed-rank test or sign test 
#' with unknown center of symmetry estimated by sample mean.
#' @param x pilot data, numeric vector of data values.
#' @param sig.level significance level (Type I error probability)
#' @param power power of test (1 minus Type II error probability)
#' @param method a character string specifying which symmetry test to be used, "wilcox" refers to Wilcoxon signed-rank test,
#'   and "sign" is sign test.
#' @param alternative 
#'  a character string specifying the alternative hypothesis, must be one of "two sided" (default), "right.skewed", or "left.skewed".
#'   You can specify just the initial letter.
#'   "right.skewed": test whether positively skewed, 
#'   "left.skewed" : test whether negatively skewed.
#' 
#' @importFrom e1071 skewness kurtosis 
#'
#' @returns A list of class "power.htest" containing the following components:
#'    \itemize{
#'   \item N - sample size estimated
#'   \item sig.level - significance level (Type I error probability)
#'   \item power - power of test (1 minus Type II error probability)
#'   \item method - the type of test applied
#'   \item alternative - one- or two-sided test. Can be abbreviated
#'   }
#' @export
#'
#' @examples
n.symm.test <- function(x, sig.level = 0.05, power = 0.8, method="wilcox",
                        alternative = c("two.sided", "right.skewed", "left.skewed")
                        ){
 
  logit.rev <- function(u) exp(u) / (1+exp(u))
  B <- 1000
  upper.cutoff <- 10000
  sigma2.x <- var(x)
  n <- length(x)
  z_alpha <- qnorm(1-sig.level/2)
  alternative <- match.arg(alternative)
  
  if (method == "wilcox"){
    METHOD <- "Sample size calculation under wilcox procedure"
    # quantities under H0
    est_quant_H0 <- get_quant_H0(x)
    theta0 <- est_quant_H0$theta0
    tau0 <- est_quant_H0$tau0
    mu0 <- function(N) N*(N+1)/4
    sigma0 <- function(N) sqrt(N*(N+1)*(2*N+1)/24 - N*(N-1)*(N-3) * theta0 * tau0 +
                                 (N-1)*(N-2)*(N-3)*(N-4)*sigma2.x/(4*N)*(theta0)^2)
    # quantities under H1  
    est_quant_H1 <- get_quant_H1(x)
    p1 <- est_quant_H1$p1
    p2hat <- est_quant_H1$p2
    p3 <- est_quant_H1$p3
    p4 <- est_quant_H1$p4
    theta1 <- est_quant_H1$theta1
    tau1 <- est_quant_H1$tau1
    
    
    ## Goal: Determine whether a value of qz is close to 0.5 or not ##
    ## Method: Construct CI for qz using bootstrap
    
    Wb <- array()
    p1b <- array()
    for(b in 1:B){
      xb <- sample(x, n ,replace =TRUE)
      mb <- mean(xb)
      p1b[b] <- mean(xb> mb)
      Wb[b]<- wilcox.test(xb, mu=mb)$statistic # Wb[b]<-wilcox.test(xb,mu=m)$statistic
    }
    
    p2B<-(Wb-n*p1b)*2/(n*(n-1)) #qzB<-(Wb-n*qx)*2/(n*(n-1))
    
    p2BS<-sort(p2B)
    new_data <- data.frame(n.pilot = n, 
                           distance_p2_half = abs(p2hat-0.5), 
                           abs_obs_skewness = abs(skewness(x)), 
                           obs_sd=sd(x), 
                           obs_kurtosis = kurtosis(x)
    )
    
    pred <- as.numeric(predict(model, newdata = new_data, method = "anova"))
    alp <- logit.rev(pred)
    
    
    aa<-p2BS[round((1-alp/2)*B)]
    bb<-p2BS[round((1/2)*B)]
    cc<-p2BS[round((alp/2)*B)]
    
    Com1<-c(aa,bb,cc)
    Com2<-c(abs(aa-0.5),abs(bb-0.5),abs(cc-0.5))
    Index<-which.max(Com2)
    
    p2.farest <- Com1[Index]
    p2 <- p2.farest
    
    # mu1() and sigma1() function
    mu1 <- function(N) N*(N-1)/2*p2 + N*p1
    sigma1 <- function(N) {
      Kn <- N*p1*(1-p1) + N*(N-1)/2 * (p2*(1-p2) + 4*(p3-p1*p2)) +
        N*(N-1)*(N-2)*(p4-p2^2)
      
      sig1.square <- Kn - N*(N-1)*(N-3) * theta1 * tau1 +
        (N-1)*(N-2)*(N-3)*(N-4)*sigma2.x/(4*N)*(theta1)^2
      
      if (sig1.square > 0) {
        sigma1.value <- sqrt(sig1.square)
      }else{
        sigma1.value <- sqrt(N*(N+1)*(2*N+1)/24 - N*(N-1)*(N-3) * theta0 * tau0 +
                               (N-1)*(N-2)*(N-3)*(N-4)*sigma2.x/(4*N)*(theta0)^2)
      }
      return(sigma1.value)
    }
    sigma1 <- Vectorize(sigma1, "N")
    
    
    z1 <- function(N) (mu0(N)-mu1(N)+sigma0(N)*z_alpha)/sigma1(N)
    z2 <- function(N) (mu0(N)-mu1(N)-sigma0(N)*z_alpha)/sigma1(N)
    power.fn <-  function(N) pnorm(z2(N)) - pnorm(z1(N)) + 1
    power.fn <- Vectorize(power.fn, "N")
    
    pow.piliot <- power.fn(n)
    
    # if (pow.piliot >= power) {
    #   print("Power evaluated at the pilot sample size exceeds the target power.")
    #   n.out = n
    #   warning("The power of the given pilot sample exceeds the target power!")
    #   
    # }else{
    #   n.out <- 10
    #   while (n.out < upper.cutoff & power.fn(n.out) < power){
    #     n.out <- n.out + 1
    #     if (n.out >= upper.cutoff) {
    #       n.out = upper.cutoff
    #       break
    #     }
    #   }
    # }
    
    if (pow.piliot >= power){
      warning("The power of the given pilot sample exceeds the target power!") }
    
    n.out <- 10
    while (n.out < upper.cutoff & power.fn(n.out) < power){
      n.out <- n.out + 1
      if (n.out >= upper.cutoff) {
        n.out = upper.cutoff
        break
      }
    }
    
  }
  
  ## return list of 'power.htest' class
  RVAL <- list(method = METHOD,
               N = n.out,
               sig.level = sig.level,
               power = power,
               type = method,
               alternative = alternative
  )
  
  
  class(RVAL) <- "power.htest"
  RVAL
  
}
