N_vec <- c(10,100,1000)
sig_vec <- c(1,5,10)
alpha_vec <- c(0.05,0.005,0.0005)

Results <- matrix(NA,3^3,4)

count <- 0
for (NN in 1:3) {
  for (SS in 1:3) {
    for (AA in 1:3) {
      ### Sim
      count <- count + 1
      efr_cov <- c()
      was_cov <- c()
      efr_len <- c()
      was_len <- c()
      sigma <- sig_vec[SS]
        alpha <- alpha_vec[AA]
        for (nn in 1:1000) {
          N <- N_vec[NN]
          THETA <- rnorm(N,0,sigma)
          X <- rnorm(N,THETA,1)
          EST <- max(var(X[-1])-1,0)
          
          Bhat <- 1 - (N-2)/sum(X*X)
          efrL <- Bhat*X[1]-qnorm(1-alpha/2)*sqrt(Bhat+2/(N-2)*X[1]^2*(1-Bhat)^2)
          efrU <- Bhat*X[1]+qnorm(1-alpha/2)*sqrt(Bhat+2/(N-2)*X[1]^2*(1-Bhat)^2)
          
          
          wasL <- X[1] - sqrt(2*log(1/alpha) + 2*log(1+EST) + X[1]^2/(1+EST))
          wasU <- X[1] + sqrt(2*log(1/alpha) + 2*log(1+EST) + X[1]^2/(1+EST))
          
          efr_cov[nn] <- (efrL<THETA[1])*(THETA[1]<efrU)
          was_cov[nn] <- (wasL<THETA[1])*(THETA[1]<wasU)
          efr_len[nn] <- efrU-efrL
          was_len[nn] <- wasU-wasL
          
        }
      Results[count,1] <- mean(efr_cov,na.rm=T)
      Results[count,2] <- mean(was_cov,na.rm=T)
      Results[count,3] <- sum(is.na(efr_cov))
      Results[count,4] <- mean(was_len/efr_len,na.rm=T)
    }
  }
}


