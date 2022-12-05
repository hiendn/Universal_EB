### Sim
efr_cov <- c()
was_cov <- c()
efr_len <- c()
was_len <- c()
sigma <- 5
alpha <- 0.01
for (nn in 1:100) {
  N <- 10
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

crit2 <- function(x) {log(dnorm(X[1],0,sqrt(1+EST)))-log(dnorm(X[1],x,1))}
curve(crit2,from = -30,to = 30,n=1000)
abline(h=log(1/alpha),col='blue')
abline(v=efrL)
abline(v=efrU)
abline(v=wasL,col='green')
abline(v=wasU,col='green')
abline(v=THETA[1],col='red')

mean(efr_cov)
mean(was_cov)
mean(efr_len)
mean(was_len)
mean(was_len/efr_len,na.rm=T)