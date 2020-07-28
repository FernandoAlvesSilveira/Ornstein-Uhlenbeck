"""
References:

ALATON, Peter; DJEHICHE, Boualem; STILLBERGER, David. 
On modelling and pricing weather derivatives. Applied mathematical finance, v. 9, n. 1, p. 1-20, 2002.

KORDI, KONSTANTINA. Pricing weather derivatives (Unpublished doctoral dissertation). University of Piraeus, 2012.
"""
library(tidyverse)

### Ornstein–Uhlenbeck model to estimate temperature and pricing weather derivative on Monte Carlo simulations

#Import a time series object, is suggest use a na.omit() or replace_na() to avoid "NA" effects

#tetha function to estimate A, B, C and phi params
tetha <- function(tmedia) {
  time <- rep((seq(1:365)), length(tmedia)/365)
  serie <- data.frame(tmedia[1:length(time)], time)
  
  tetha_model <- lm(serie$tmedia ~ serie$time + sin(2*pi*(serie$time)/365) +
                      cos(2*pi*(serie$time/365)), data = serie)
  tpredict <- predict.lm(tetha_model)
  Parameters = coef(tetha_model)
  A = Parameters[1]
  B = Parameters[2]
  C = sqrt(Parameters[3]^2+Parameters[4]^2)
  phi = (atan(Parameters[4]/Parameters[3]) - pi)
  
  print(A)
  print(B)
  print(C)
  print(phi)
  
  return(tpredict)
}

#Speed of mean reversion
sigma <- function(tmedia) {
  counter <- 1
  sum_temp <- NA
  while (counter < length(tmedia)) {
    sum_temp[counter] <- (tmedia[counter + 1] - tmedia[counter])^2
    counter <- counter + 1
  } 
  return( (1/length(tmedia))*(sum(sum_temp, na.rm = TRUE)) )
}

#Volatility estimation
alpha <- function(tmedia, tpredict, sigma) {
  counter <- 2
  y <- NA
  Y <- NA
  while (counter < length(tmedia)) {
    y[counter] <- (tpredict[counter-1] - tmedia[counter -1] / sigma) * (tmedia[counter] - tpredict[counter])
    Y[counter] <- (tpredict[counter-1] - tmedia[counter -1] / sigma) * (tmedia[counter-1] - tpredict[counter-1])
    counter <- counter + 1
  }
  return(-log(sum(y, na.rm = TRUE)/sum(Y, na.rm = TRUE)))
  
}

#Ornstein–Uhlenbeck simulation process

#1 steep: simple simulation of OU process
ou_simulation <- function(To, A, B, C, alpha, omega, sigma, path, lambda, phi) {
  paths <- seq(1:path)
  i <- 2
  paths[1] <- To
  while (i <= path) {
    paths[i] <- A + B*(paths[i]+1) + C*sin((omega*(paths[i] + 1) + phi)) + 
      (1 - alpha) * paths - (A + B* paths[i] + C*sin(omega+paths[i]+phi)) +
      sigma*rnorm(1) - lambda*sigma
    i <- i + 1
  }
  return(paths)
}

#2 steep: repeat and storage results of OU simulation
ou_process <- function(To, A, B, C, alpha, omega, sigma, path, lambda, phi, repets) {
  j <- 1
  results <- data.frame(seq(1:365))
  for (j in 1:repets) {
    results[,j] <- ou_simulation(To, A, B, C, alpha, omega, sigma, path, lambda, phi)
    j <- j + 1
  }
  p <- 1
  z <- NA
  for (p in 1:365) {
    z[p] <-(sum(results[p,1:repets])/repets)
    p <- p + 1
  }
  return(z)
}

#Pricing with Monte Carlo simulation
monte_carlo <- function(To, A, B, C, alpha, omega, sigma, path, lambda, phi, repets, K, r, tmax, tmin, repet, tick) {
  PayoffCall <- data.frame(seq(1:repet))
  PayoffPut <- data.frame(seq(1:repet))
  i <- 1
  H <- NA
  SimPayoffCall <- NA
  SimPayoffPut <- NA
  while (i <= repet) {
    values <- ou_process(To, A, B, C, alpha, omega, sigma, path, lambda, phi, repets)
    H[i] <- sum(max(tmax - values[1:path]))
    SimPayoffCall[i] <- exp(-r*(repet))*max(H[i] - K)
    SimPayoffPut[i] <- exp(-r*(repet))*max(K - H[i])
    i <- i + 1
  }
  print(Callprice <- (tick*mean(PayoffCall)))
  print(Putprice <- (tick*mean(SimPayoffPut)))
}
