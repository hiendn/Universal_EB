reject <- c()
for (kk in 1:1000) {
  NN_vec <- c(10,100,1000)
  alpha_vec <- c(2,2,5)
  beta_vec <- c(2,5,2)
  signif_vec <- c(0.05,0.005,0.0005)
  diff_vec <- c(0,1,5,10)
  
  # Simulation parameters
  NN <- NN_vec[2]
  alpha <- alpha_vec[2]
  beta <- beta_vec[2]
  signif <- signif_vec[1]
  
  # Simulate data
  THETA <- rgamma(NN,alpha,beta)
  THETA[NN] <- THETA[NN-1] + diff_vec[3] 
  
  EE <- runif(NN)*10
  XX <- rpois(NN,THETA*EE)
  
  # Integrated likelihood functional
  int_log_like <- function(para) {
    hold <- 0
    for (ii in 1:(NN-2)) {
      hold <- hold + log(choose(XX[ii]+exp(para[1])-1,XX[ii])) +
        exp(para[1])*log(exp(para[2])) - exp(para[1])*log(EE[ii]+exp(para[2])) +
        XX[ii]*log(EE[ii]) - XX[ii]*log(EE[ii]+exp(para[2]))
    }
    return(-hold)
  }
  
  # MLE
  mle_optim <- optim(log(c(alpha,beta)),int_log_like,method='BFGS')
  alpha_hat <- exp(mle_optim$par[1])
  beta_hat <- exp(mle_optim$par[2])
  
  THETA_HAT <- (XX+alpha_hat)/(EE+beta_hat)
  
  
  ## Ln
  Ln <- exp(log(choose(XX[NN] + alpha_hat  - 1, XX[NN])) +
    log((beta_hat/(EE[NN]+beta_hat))^alpha_hat) +
    log((EE[NN]/(EE[NN]+beta_hat))^XX[NN]) +
    log(choose(XX[NN-1] + alpha_hat  - 1, XX[NN-1])) +
    log((beta_hat/(EE[NN-1]+beta_hat))^alpha_hat) +
    log((EE[NN-1]/(EE[NN-1]+beta_hat))^XX[NN-1]))
  
  elln <- function(theta) {
    exp(dpois(XX[NN],EE[NN]*theta,log = TRUE) + dpois(XX[NN-1],EE[NN-1]*theta,log = TRUE))
  }
  Neg_denom <- function(para) {
    -elln(para)
  }
  
  denom_optim <- optimize(Neg_denom,c(0,100))
  theta_denom <- denom_optim$minimum
  
  reject[kk] <- ( (elln(theta_denom)/Ln)<=signif )
  print(c(kk,reject[kk],elln(theta_denom),Ln,alpha_hat,beta_hat))
}

mean(reject,na.rm=TRUE)
sum(is.na(reject))