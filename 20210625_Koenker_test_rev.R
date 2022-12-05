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
  
  Ln_vec <- c()
  
  for (rr in 1:5) {
    
    # Sub sampling
    sub_obs <- NN
    sam_index <- sample(1:(NN-2),size=sub_obs,replace = TRUE)
    XX_sam <- XX[sam_index]
    EE_sam <- EE[sam_index]
    
    # Integrated likelihood functional
    int_log_like <- function(para) {
      hold <- 0
      for (ii in 1:sub_obs) {
        hold <- hold + log(choose(XX_sam[ii]+exp(para[1])-1,XX_sam[ii])) +
          exp(para[1])*log(exp(para[2])) - exp(para[1])*log(EE_sam[ii]+exp(para[2])) +
          XX_sam[ii]*log(EE_sam[ii]) - XX_sam[ii]*log(EE_sam[ii]+exp(para[2]))
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
    
    Ln_vec[rr] <- Ln
  }


  
  elln <- function(theta) {
    exp(dpois(XX[NN],EE[NN]*theta,log = TRUE) + dpois(XX[NN-1],EE[NN-1]*theta,log = TRUE))
  }
  Neg_denom <- function(para) {
    -elln(para)
  }
  
  denom_optim <- optimize(Neg_denom,c(0,100))
  theta_denom <- denom_optim$minimum
  
  reject[kk] <- ( (elln(theta_denom)/mean(Ln))<=signif )
  print(c(kk,reject[kk],elln(theta_denom),mean(Ln),alpha_hat,beta_hat))
}

mean(reject,na.rm=TRUE)
sum(is.na(reject))