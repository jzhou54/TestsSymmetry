############### function for getting estimates under H0 ########
get_quant_H0 <- function(x)  {
  n <- as.double(length(x))
  m <- mean(x)
  sigma2 <- var(x)
  xc <- x-m
  
  ## Tn selection
  r <- quantile(x, c(0.25, 0.75))
  h <- (r[2] - r[1])/1.34
  Tn <- log(n)/( 3 * 1.06 * min(sqrt(var(x)), h))
  
  ## Estimation of theta
  S <- function(u) sum( sin(2*pi*(xc[xc!=u]-u)*Tn)/(2*pi*(xc[xc!=u]-u)))+
    sum(sin(2*pi*(xc+u)*Tn)/(2*pi*(xc+u)))
  SV <- Vectorize(S)
  hat_theta <- 2*sum(SV(xc))/n^2 + 2*Tn/n
  
  ## Estimation of tau
  xs <- sort(xc)  # V(i)
  S1 <- seq(from=1, to=n, by=1) # i
  hat_tau <- sum(xs*S1)/n^2
  
  ## quantity for sign test ##
  a <- 1 * (x > m - n ^ (-1 / 5)) * (x < m + n ^ (-1 / 5))
  D <- sum(a)
  A <- max(c(1, D))
  VHW <- (n ^ (3 / 10) / A) ^ 2
  hat_w <- sqrt(1 / (4 * n * VHW))
  CE <- mean(x * (x < m)) - m / 2
  
  return(list("tau0" = hat_tau,
              "theta0" = hat_theta,
              "w" = hat_w, 
              "CE" = CE) )
}

############### function for getting estimates under H1 ########
get_quant_H1 <- function(x) {
  #------------------------------------------------------#
  m <- mean(x)
  n <- length(x)
  y <- x - m
  ys <- sort(y)
  y.star <- c(2*ys[1]-ys[2], ys, 2*ys[n] - ys[n-1])
  # Gy
  G.y.star <- function(u) {
    if (u <= y.star[1]) {
      return(0)
    } else if (u >= y.star[length(y.star)]) {
      return(1)
    } else {
      interval_index <- findInterval(u, y.star) # find which interval that u belongs to
      Yi <- y.star[interval_index]; Yii <- y.star[interval_index+1]
      val <- (interval_index-1) / (n+1) + (u - Yi) / ((n+1)*(Yii - Yi))
      return(val)
    }
  }
  G.y.star.V <- Vectorize(G.y.star)
  #------------------------------------------------------#
  # estimators under H1
  Gy_y <- G.y.star.V(-y)
  # p1=Pr(X1>mu)
  p1 <- 1-G.y.star.V(0)
  # p2=Pr(X1+X2>2mu)
  p2 <- mean(1 - Gy_y)
  # p3=Pr(X1+X2>2mu, X1>mu)
  p3 <- mean( (1 - Gy_y) * (y > 0))
  # p4=Pr(X1+X2>2mu, X1+X3>2mu)
  p4 <- mean( (1 - Gy_y)^2 )
  #------------------------------------------------------#
  
  # theta
  f.x <- function(u) density(x, from = u, to=u)$y[1]
  f.x.vec <- Vectorize(f.x)
  f.Z <- function (u) 2*mean(f.x.vec(2*u-x))
  theta1 <- f.Z(m)
  
  # tau
  FX <- ecdf(x)
  tau1 <- -mean(x*FX(2*m-x)) + m*mean(FX(2*m-x)) # based on ecdf
  
  
  ## quantity for sign test ##
  hat_w <- f.x.vec(m)
  CE <- mean((x - m) * (x < m))
  
  return(list("p1"=p1, "p2"=p2, "p3"=p3, "p4"=p4, "theta1"=theta1, "tau1"=tau1, "w"= hat_w, "CE" = CE))
}

