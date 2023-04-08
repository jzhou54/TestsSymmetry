############### function for getting estimates under H0 ########
get_quant_H0 <- function(x)  {
  n <- as.double(length(x))
  m <- mean(x)
  sigma2 <- var(x)
  xc <- x-m
  
  ## Tn selection
  r <- quantile(x, c(0.25, 0.75))
  h <- (r[2] - r[1])/1.34
  Tn<-log(n)/( 3 * 1.06 * min(sqrt(var(x)), h))
  
  ## Estimation of theta
  S <- function(u) sum( sin(2*pi*(xc[xc!=u]-u)*Tn)/(2*pi*(xc[xc!=u]-u)))+
    sum(sin(2*pi*(xc+u)*Tn)/(2*pi*(xc+u)))
  SV <- Vectorize(S)
  hat_theta <- 2*sum(SV(xc))/n^2 + 2*Tn/n
  
  ## Estimation of tau
  xs <- sort(xc)  # V(i)
  S1 <- seq(from=1, to=n, by=1) # i
  hat_tau <- sum(xs*S1)/n^2
  return(list(tau0 = hat_tau,
              theta0 = hat_theta) )
}

############### function for getting estimates under H1 ########
get_quant_H1 <- function(x) {
  mean.X <- mean(x)
  z = Rfit::walsh(x)
  n = length(x)
  
  # 1. nonparametric estimate qX = Pr(X> E(X))
  qx <- mean(x > mean.X)
  
  # 2. nonparametric estimate qZ = Pr(Z12> E(Z12))
  qz <- mean(z > mean.X)
  # 3. nonparametric estimation of L1
  
  Y = x - mean.X
  Fy = ecdf(Y)
  L1 <- mean(Fy(Y)*Fy(-Y) - (1 - Fy(-Y))^2)
  
  # 4. nonparametric estimate L2
  L2 <- mean(Fy(-Y) * (Y>0))
  
  # 5. nonparametric estimate theta = f_Z(mu)
  f.Z <- function(u) {density(z, from = u, to = u)$y[1] }
  theta1 <- f.Z(mean.X)
  
  # 6. nonparametric estimation of tau
  P_z = ecdf(z)
  Zs <- sort(z[z>=mean.X])
  P_Zs <- P_z(Zs)
  tau1 <- sum(diff(Zs) * 0.5*( 2 - P_Zs[2:length(P_Zs)] - P_Zs[1:(length(P_Zs)-1)]))
  return(list("qx"=qx, "qz"=qz, "L1"=L1, "L2"=L2, "theta1"=theta1, "tau1"=tau1))
}
