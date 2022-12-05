N <- 3
THETA <- rnorm(N,0,1)
X <- rnorm(N,THETA,1)

EST <- max(var(X[-1])-1,0.00001)

crit2 <- function(x) {log(dnorm(X[1],0,sqrt(EST+1)))-log(dnorm(X[1],x,1))}

Bhat <- 1 - (N-2)/sum(X*X)
efrL <- Bhat*X[1]-1.96*sqrt(Bhat+2/(N-2)*X[1]^2*(1-Bhat)^2)
efrU <- Bhat*X[1]+1.96*sqrt(Bhat+2/(N-2)*X[1]^2*(1-Bhat)^2)


wasL <- X[1] - sqrt(2*log(1/0.05) + 2*log(EST+1) + X[1]^2/(1+EST))
wasU <- X[1] + sqrt(2*log(1/0.05) + 2*log(EST+1) + X[1]^2/(1+EST))

curve(crit2,from = -10,to = 10,n=1000)
abline(h=log(1/0.05),col='blue')
abline(v=efrL)
abline(v=efrU)
abline(v=wasL,col='green')
abline(v=wasU,col='green')
abline(v=THETA[1],col='red')

### Sim
efr_cov <- c()
was_cov <- c()
efr_len <- c()
was_len <- c()
sigma <- 1
alpha <- 0.05
for (nn in 1:1000) {
  N <- 3
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

efr_cov
was_cov
efr_len
was_len