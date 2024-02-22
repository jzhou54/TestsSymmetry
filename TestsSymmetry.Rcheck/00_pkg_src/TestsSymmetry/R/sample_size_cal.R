#' @title Sample size determination for nonparametric tests of symmetry when the center is unknown
#' @description Determine the sample size required for one-sample Wilcoxon signed-rank test or Sign test 
#' with an unknown symmetry center estimated by sample mean given a target power.The function uses learning sample 
#' in order to predict (calculate) sample size needed to reach a preassumed power based on the underlying data that
#' is exemplified by the learning sample.
#' 
#' @param x learning sample data, numeric vector of data values
#' @param sig.level significance level (Type I error probability), the default value is 0.05
#' @param power power of test (1 minus Type II error probability), the default value is 0.8
#' @param method a character string specifying which symmetry test to be used, "wilcox" refers to Wilcoxon signed-rank test,
#'   and "sign" sign test
#' @param alternative 
#'  a character string specifying the alternative hypothesis
#'   "two.sided" : test whether skewed
#'  
#' 
#' @importFrom e1071 skewness kurtosis
#' @importFrom rpart rpart
#' @details A Normal approximation to the power requires specification of some unknown quantities in the nonparametric context. 
#' In this regard, empirical smoothed CDF and Bootstrap methods were leveraged to estimate these quantities using learning sample `x`.  
#' 
#' Remark: If the test provides a power, say, P, based on the learning data and P is higher than the target power, a warning
#' message will be shown. However, a needed sample size N to reach the target power will be conducted. 
#' 
#' @returns A list of class "power.htest" containing the following components:
#'    \itemize{
#'   \item N - sample size estimated
#'   \item sig.level - significance level (Type I error probability)
#'   \item power - power of test (1 minus Type II error probability)
#'   \item method - the test method applied
#'   \item alternative - two-sided test. Can be abbreviated
#'   }
#' @references
#' Chakraborti, S., Hong, B., & van de Wiel, M. A. (2006). A Note on Sample Size Determination for a Nonparametric Test of Location. Technometrics, 48(1), 88-94.
#' 
#' Vexler, A., Gao, X., & Zhou, J. (2023). How to implement signed-rank `wilcox.test()` type procedures when a center of symmetry is unknown. Computational Statistics & Data Analysis, 107746. 
#' 
#' Gastwirth, J. L. (1971). On the Sign Test for Symmetry. Journal of the American Statistical Association, 66(336), 821-823.
#'
#'
#' @examples
#' \donttest{
#'  data("plasma.silicon")
#'  post <- plasma.silicon$postoperative 
#'  pre <- plasma.silicon$preoperative
#'  diff <- post - pre
#'  n.symm.test(diff, sig.level = 0.05, power = 0.5, method = "wilcox", alternative ="two.sided" )
#' 
#' # Result:
#' # Sample size calculation under wilcox procedure 
#' 
#' #           N = 83
#' #   sig.level = 0.05
#' #       power = 0.5
#' #        type = wilcox
#' # alternative = two.sided
#' 
#' # Interpretation: 
#' # Given the pilot sample `diff` and the significance level 0.05. The sample size of 
#' # the data that is expected toprovide the target power 0.5 of the Wilcoxon test procedure 
#' # is computed as 83.
#' }
#' @export
n.symm.test <- function(x, sig.level = 0.05, power = 0.8, method="wilcox",
                        alternative = c("two.sided")
                        ){
  
  logit.rev <- function(u) exp(u) / (1+exp(u))
  B <- 3000
  upper.cutoff <- 3000
  sigma2.x <- var(x)
  n <- length(x)
  m <- mean(x)
  z_alpha <- qnorm(1-sig.level/2)
  alternative <- match.arg(alternative)
  
  if (method == "wilcox"){
    METHOD <- "Sample size calculation under Wilcoxon procedure"
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
      p1b[b] <- mean(xb > mb)
      
      xbc_abs <- abs(xb - mb)
      xbc_abs_rank <- rank(xbc_abs)
      # walsh.avg.mat <- outer(xb, xb, FUN = function(x, y) (x+y)/2)
      # Wb[b] <- sum(walsh.avg.mat > mb)
      Wb[b] <- sum(xbc_abs_rank*(xb>mb))
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
    
  }else if (method == "sign") {
    METHOD <- "Sample size calculation under sign test procedure"
    
    Sb <- array()
    for (b in 1:B) {
      xb <- sample(x, size = n, replace = T)
      mb <- mean(xb)
      Sb[b] <- sum(xb < mb) #  test statistic of sign test
    }
    
    pxB <- Sb / n # estimation of px, E(S*) = n*px
    pxBS <- sort(pxB)
    
    
    ## quantity under H0 ##
    quant_H0 <- get_quant_H0(x)
    w0 <- quant_H0$w
    CE0 <- quant_H0$CE
    
    mu0 <-  function(N) N / 2
    sigma0 <- function(N) sqrt((1 / 4 + var(x) * w0 ^ 2 + 2 * w0 * CE0) * N)
    
    ## quantities under H1 ##
    quant_H1 <- get_quant_H1(x)
    w1 <- quant_H1$w
    CE1 <- quant_H1$CE
  
    new_data <- data.frame(n.pilot = n, 
                           sd=sd(x),
                           abs_sk = abs(skewness(x)), 
                           kur = kurtosis(x))
    
    pred <- as.numeric(predict(tree_model_sign, newdata = new_data, method = "anova"))
    alp <- logit.rev(pred)
    
    px <- pxBS[alp * B]
    mu1 <-  function(N) N * px
    sigma1 <- function(N) sqrt(((1 - px) * px + var(x) * w1 ^ 2 + 2 * w1 * CE1) * N)
    if (1 / 4 + var(x) * w0 ^ 2 + 2 * w0 * CE0 < 0)  {
      sigma0 <- function(N)  sqrt(((1 - px) * px + var(x)*w1^2 + 2*w1*CE1) * N)}
    
    z1 <- function(N) (mu0(N) - mu1(N) + sigma0(N) * z_alpha) / sigma1(N)
    z2 <- function(N) (mu0(N) - mu1(N) - sigma0(N) * z_alpha) / sigma1(N)
    power.fn <-  function(N) pnorm(z2(N)) - pnorm(z1(N)) + 1
    power.fn <- Vectorize(power.fn, "N")
    pow.piliot <- power.fn(n)
    
    if (pow.piliot >= power) {
      warning("The power of the given pilot sample exceeds the target power!") }
    
    n.out <- 10
    while (n.out < upper.cutoff & power.fn(n.out) < power) {
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
               method = method,
               alternative = alternative
  )
  
  
  class(RVAL) <- "power.htest"
  RVAL
  
}
