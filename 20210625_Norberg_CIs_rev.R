# Define roll function
roll <- function( x , n ){
  if( n == 0 )
    return( x )
  c( tail(x,n) , head(x,-n) )
}

# load libraries
library(REBayes)

# Simulation parameters
NN <- 72
alpha <- 1
beta <- 1
signif <- 3*0.05/72

storage <- matrix(NA,72,3)

for (ii in 1:72) {
  # Simulate data
  THETA <- rgamma(NN,alpha,beta)
  EE <- Norberg[,2]/344
  XX <- Norberg[,3]
  ROLLED <- ii - 1
  EE <- roll(EE,ROLLED)
  XX <- roll(XX,ROLLED)
  
  
  # Integrated likelihood functional
  int_log_like <- function(para) {
    hold <- 0
    for (ii in 1:(NN-1)) {
      hold <- hold + log(choose(XX[ii]+exp(para[1])-1,XX[ii])) +
        exp(para[1])*log(exp(para[2])) - exp(para[1])*log(EE[ii]+exp(para[2])) +
        XX[ii]*log(EE[ii]) - XX[ii]*log(EE[ii]+exp(para[2]))
    }
    return(-hold)
  }
  
  # MLE
  mle_optim <- optim(log(c(alpha,beta)),int_log_like)
  alpha_hat <- exp(mle_optim$par[1])
  beta_hat <- exp(mle_optim$par[2])
  
  # THETA_HAT <- (XX+alpha_hat)/(EE+beta_hat)
  
  # plot(THETA,THETA_HAT)
  # abline(a=0,b=1)
  
  ## Ln
  Ln <- choose(XX[NN] + alpha_hat  - 1, XX[NN]) *
    (beta_hat/(EE[NN]+beta_hat))^alpha_hat *
    (EE[NN]/(EE[NN]+beta_hat))^XX[NN]
  elln <- function(theta) {
    dpois(XX[NN],EE[NN]*theta)
  }
  Ratio <- function(theta) {
    Ln/elln(theta)
  }
  
  tt <- seq(0,500,length.out = 100000)
  rr <- Ratio(tt)
  
  # 
  # plot(tt,log(rr),type='l')
  # abline(h = log(1/signif),col='red')
  # abline(v = THETA_HAT[NN],col='blue')
  
  lower_b <- tt[min(which(log(rr)<log(1/signif)))]
  upper_b <- tt[max(which(log(rr)<log(1/signif)))]
  
  # abline(v = lower_b,col='green')
  # abline(v = upper_b,col='green')
  
  # (THETA[NN]>=lower_b)*(THETA[NN]<=upper_b)
  storage[72-ii+1,1] <- lower_b
  storage[72-ii+1,3] <- upper_b
}

# Simulate data
THETA <- rgamma(NN,alpha,beta)
EE <- Norberg[,2]/344
XX <- Norberg[,3]

# Integrated likelihood functional
int_log_like <- function(para) {
  hold <- 0
  for (ii in 1:(NN)) {
    hold <- hold + log(choose(XX[ii]+exp(para[1])-1,XX[ii])) +
      exp(para[1])*log(exp(para[2])) - exp(para[1])*log(EE[ii]+exp(para[2])) +
      XX[ii]*log(EE[ii]) - XX[ii]*log(EE[ii]+exp(para[2]))
  }
  return(-hold)
}
# MLE
mle_optim <- optim(log(c(alpha,beta)),int_log_like)
alpha_hat <- exp(mle_optim$par[1])
beta_hat <- exp(mle_optim$par[2])

THETA_HAT <- (XX+alpha_hat)/(EE+beta_hat)

storage[,2] <- THETA_HAT
par(mfrow=c(2,1),mar=c(4,4,1,2)+0.1)
plot(1:72,c(max(10),rep(0,71)),type='n',xlab='Occupation group',ylab='Risk factor')
points(storage[,2],pch=4,col='red',cex=1)
for (ii in 1:72) {
  lines(c(ii,ii),storage[ii,c(1,3)],col='blue')
}
grid()

