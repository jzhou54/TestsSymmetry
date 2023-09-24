#' Title
#' @title Sample size calculation for Wilcoxon test of symmetry with estimated location 
#' @description 
#' Determine the sample size required for one-sample wilcoxon signed-rank test or sign test 
#' with unknown center of symmetry estimated by sample mean.
#' 
#' @param x pilot data, numeric vector of data values. Non-finite (e.g., infinite or missing) values will be omitted.
#' @param sig.level significance level (Type I error probability)
#' @param power power of test (1 minus Type II error probability)
#' @param method a character string specifying which symmetry test to be used, "wilcox" refers to Wilcoxon signed-rank test,
#'   and "sign" is sign test.
#' @param alternative
#'   a character string specifying the alternative hypothesis, must be one of "two sided" (default), "right.skewed", or "left.skewed".
#'   You can specify just the initial letter.
#'   "right.skewed": test whether positively skewed, 
#'   "left.skewed" : test whether negatively skewed.
#' @param plot.extrapolation a logical variable indicating whether to show the data-driven extrapolation plot and the corresponding sample size estimated. 
#'   
#' @returns A list of class "power.htest" containing the following components:
#'    \itemize{
#'   \item N - sample size estimated
#'   \item sig.level - significance level (Type I error probability)
#'   \item power - power of test (1 minus Type II error probability)
#'   \item method - the type of test applied.
#'   }
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
#' y <- rexp(35, rate=3)
#' pwr.symm.test(y, plot.extrapolation = T)
#' @export
pwr.symm.test<- function(x, sig.level = 0.05, power = 0.8, method="wilcox",
                         alternative = c("two.sided", "right.skewed", "left.skewed"),
                         plot.extrapolation = FALSE){
  
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
    METHOD <- "Modified Wilcoxon signed-rank test power calculation"
   
    ########## estimates under H0 ##########
    
    theta0 <- get_quant_H0(x)$theta0
    tau0 <- get_quant_H0(x)$tau0
    mu0 <- function(N) N*(N+1)/4
    sigma0 <- function(N) sqrt(N*(N+1)*(2*N+1)/24 - N*(N-1)*(N-3) * theta0 * tau0 +
                                 (N-1)*(N-2)*(N-3)*(N-4)*sigma2.x/(4*N)*(theta0)^2)
    
    ##### estimates under H1 #####
    est_quant_H1 <- get_quant_H1(x)
    p1 <- est_quant_H1$p1
    p2 <- est_quant_H1$p2
    p3 <- est_quant_H1$p3
    p4 <- est_quant_H1$p4
    theta1 <- est_quant_H1$theta1
    tau1 <- est_quant_H1$tau1
    
    
    ### exact estimation of qz ######
    sig <- sd(x)
    Y <- x - m
    CO <- function(u) mean(cos(Y*u))
    COV <- Vectorize(CO)
    SI <- function(u) mean(sin(Y*u))
    SIV <- Vectorize(SI)
    Intg <- function(t) COV(t*(2-n)/n)*SIV(t*(2-n)/n)*exp(-2*t^2*sig^2*(n-2)/n^2)/t
    
    r = quantile(Y, c(0.25, 0.75))
    h = (r[2] - r[1])/ 1.34
    Tn <- log(n) / (3*1.06 * min(sqrt(var(Y)),h ))
    qz <- 1/2 - 2/pi*integrate(Intg, 1/Tn, Tn, stop.on.error = FALSE)$value
    
    
    ## mu1 using p2 ##
    #mu1 <- function(N) N*(N-1)/2*p2 + N*p1
    
    ## mu1 using qz ##
    mu1 <- function(N) N*(N-1)/2*qz + N*p1
    
    sigma1 <- function(N) {
      Kn <- N*p1*(1-p1) + N*(N-1)/2*(p2*(1-p2) + 4*(p3-p1*p2)) +
        N*(N-1)*(N-2)*(p4-p2^2)
      # Kn <- N*(N+1)*(2*N+1)/24
      
      sig1.square <- Kn - N*(N-1)*(N-3) * theta1 * tau1 +
        (N-1)*(N-2)*(N-3)*(N-4)*sigma2.x/(4*N)*(theta1)^2
      # sig1.square <- Kn - N*(N-1)*(N-3) * theta0 * tau0 +
      #   (N-1)*(N-2)*(N-3)*(N-4)*sigma2.x/(4*N)*(theta0)^2
      
      if (sig1.square > 0) {
        sigma1.value <- sqrt(sig1.square)
      }else{
        sigma1.value <- sqrt(N*(N+1)*(2*N+1)/24 - N*(N-1)*(N-3) * theta0 * tau0 +
                               (N-1)*(N-2)*(N-3)*(N-4)*sigma2.x/(4*N)*(theta0)^2)
      }
      return(sigma1.value)
    }
    sigma1 <- Vectorize(sigma1, "N")
    
    
    # Normal approximation
    z_fns <- function(N) {
      z1 <- (mu0(N)-mu1(N)+sigma0(N)*z_alpha)/sigma1(N)
      z2 <- (mu0(N)-mu1(N)-sigma0(N)*z_alpha)/sigma1(N)  
      z3 <- (mu0(N)-mu1(N)+sigma0(N)*qnorm(1-sig.level))/sigma1(N) # right-skewed
      z4 <- (mu0(N)-mu1(N)-sigma0(N)*qnorm(1-sig.level))/sigma1(N) # one-sided
      return(list("z1"=z1, "z2"=z2, "z3"=z3, "z4"=z4))
    }
    # z1 <- function(N) (mu0(N)-mu1(N)+sigma0(N)*z_alpha)/sigma1(N)
    # z2 <- function(N) (mu0(N)-mu1(N)-sigma0(N)*z_alpha)/sigma1(N)  
    # z3 <- function(N) (mu0(N)-mu1(N)+sigma0(N)*qnorm(1-sig.level))/sigma1(N) # right-skewed
    # z4 <- function(N) (mu0(N)-mu1(N)-sigma0(N)*qnorm(1-sig.level))/sigma1(N) # one-sided
    
    if (alternative == "two.sided"){
      f <- function(N) (pnorm(z_fns(N)$z1) - pnorm(z_fns(N)$z2)-beta)^2
      power.fn <-  function(N) pnorm(z_fns(N)$z2) +1 - pnorm(z_fns(N)$z1)
    }else if(alternative == "right.skewed"){
      f <- function(N)  (z_fns(N)$z3 - qnorm(beta) )^2 
      power.fn <-  function(N) 1-pnorm(z_fns(N)$z3) 
    }else if (alternative == "left.skewed"){
      f <- function(N) (z_fns(N)$z4 - qnorm(1 - beta))^2
      power.fn <-  function(N) pnorm(z_fns(N)$z4)
    }
    
    f <- Vectorize(f, "N")
    power.fn <- Vectorize(power.fn, "N")
    
    if (power.fn(n) >= power ){
      getN <- n
    }else{
      getN <-  optimize(f, c(10, 1000), maximum = F)$minimum
    }
    
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
    METHOD <- "Modified sign test power calculation"
    ## Estimation of w
    a <- 1 * (x>m-n^(-1/5)) * (x < m+n^(-1/5))
    D <- sum(a)
    A <- max(c(1,D))
    VHW <- (n^(3/10)/A)^2
    hat_w <- sqrt(1/(4*n*VHW))
    E <- n/2
    CE <- mean((x-m)*(x<m))
    V0 <- function(n) (1/4 + var(x)*(hat_w)^2 + 2 *hat_w*CE)*n
    mu0 <- function(n) n/2
    
    px <- 1-qx
    f.x <- function(u) {density(x, from = u, to = u)$y[1] }
    f.x.mu <- f.x(m)
    mu1 <- function(n) n*px
    V1 <- function(n) (px*qx + var(x)*f.x.mu^2 + 2*f.x.mu*CE)*n
    
    
    
    z1 <- function(N) (mu0(N)-mu1(N)+sigma0(N)*z_alpha)/sigma1(N)
    z2 <- function(N) (mu0(N)-mu1(N)-sigma0(N)*z_alpha)/sigma1(N)  
    z3 <- function(N) (mu0(N)-mu1(N)+sigma0(N)*qnorm(1-sig.level))/sigma1(N) # right-skewed
    z4 <- function(N) (mu0(N)-mu1(N)-sigma0(N)*qnorm(1-sig.level))/sigma1(N) # one-sided
      
      
      
      
    if (alternative == "two.sided"){
      f <- function(N) (pnorm(z1(N)) - pnorm(z2(N))-beta)^2
      power.fn <-  function(N) pnorm(z1(N)) - pnorm(z2(N))-beta
    }else if(alternative == "right.skewed"){
      f <- function(N)  (z3(N) - qnorm(beta) )^2 
      power.fn <-  function(N) z3(N) - qnorm(beta)
    }else if (alternative == "left.skewed"){
      f <- function(N) (z4(N) - qnorm(1 - beta))^2
      power.fn <-  function(N) z4(N) - qnorm(1 - beta)
    }
    
    f <- Vectorize(f, "N")
    power.fn <- Vectorize(power.fn, "N")
    
    
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
  if (plot.extrapolation == TRUE){
    MC <- 1000
    n.seq <- seq(10, n, by=5)
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
    extrap.n <- 200
    n.seq.extrap <- seq(from=max(n.seq), to=extrap.n)
    y.predict <- predict(mod, data.frame(x=n.seq))
    y.predict.extrap <- exp(mod$coefficients[[1]] + mod$coefficients[[2]]*n.seq.extrap)
    
    plot(n.seq, p.avg, pch=16, ylim=c(0,1), xlim=c(5, extrap.n), xlab="n")
    lines(n.seq, exp(y.predict), col='blue')
    lines(n.seq.extrap, y.predict.extrap, col="blue", lty=2)
    abline(h = 0.05, col = "red", lty=2)
    
    N.extrapolation <- ceiling((log(0.05) - mod$coefficients[[1]] )/mod$coefficients[[2]])

    RVAL <- c(RVAL, list("N.extrapolation" = N.extrapolation))
  }
  
  ######################################
 
  class(RVAL) <- "power.htest"
  RVAL
}
