UfuncSing <- function(xvec,mvec,abvec){
  ## function for beta-binomial, returns log U 

  xm <- (mvec-xvec)
  xsum <- sum(xvec)
  msum <-  sum(xm)
  thetatilde <- xsum/sum(mvec)
  
  logR <- lbeta(xvec[1]+abvec[1], xm[1]+abvec[2]) + 
    lbeta(xvec[2]+abvec[1], xm[2]+abvec[2]) -
    2*lbeta(abvec[1],abvec[2])    - 
    xsum*log(thetatilde)- 
    msum*log(1-thetatilde)
  
  return(logR)
}


nn <- 20  

## simulation 1: theta1=theta2 for all i

NN <- 100  ## number of replications
Umat <- matrix(NA, NN, nn)

for(k in 1: NN){
  
  mm1 <- sample(seq(15,40, by=1),nn,  replace=T) 
  mm2 <- sample(seq(15,40, by=1),nn,  replace=T) 
  theta.true <- seq(0.1,0.9,length.out=nn) 
  xx <- rbinom(nn, mm1,theta.true) 
  yy <- rbinom(nn, mm2,theta.true)
 
  for(j in 1:nn){ 
    xx2 <- xx[j]; mx2 <- mm1[j]
    yy2 <- yy[j]; my2 <- mm2[j]
    compset <- setdiff(1:nn,j)
    xx1 <- xx[compset]; mx1 <- mm1[compset]
    yy1 <- yy[compset]; my1 <- mm2[compset]
    
    n <- length(xx1)
    x.st <- xx1/mx1; y.st <- yy1/my1
    allxy <- c(x.st,y.st)
    allm <- c(mx1,my1)
    
    mubar <- mean(allxy)
    xymvar <- (2*n-1)*var(allxy)/(2*n) 
    mbar <- mean(allm)
    mu.eb <- mubar
    phi.eb <- ((mbar*xymvar)/(mu.eb*(1-mu.eb))  - 1)/(mbar - 1)
    
    aa.hat <- (1/phi.eb -1)*mu.eb
    bb.hat <- (1/phi.eb -1)*(1 - mu.eb)
    abvec <- c(aa.hat, bb.hat)
    
    if(abvec[1]>0 & abvec[2]>0){
      xvec2 <- c(xx2,yy2) 
      mvec2 <- c(mx2,my2)  
      Umat[k,j] <- UfuncSing(xvec2,mvec2,abvec)
    }
  } 
} 

testsizemax <- 0.01
ylims <- range(c(na.omit(Umat),log(1/testsizemax)))

plot(theta.true,Umat[1,],type="n", ylim=ylims, xlab=expression(theta), ylab=expression(log ~ T(D[n])), 
     main="",cex.lab=1.5,cex.axis=1.1)
abline(h=log(1/0.05), col=3, lty=2,lwd=2)
abline(h=log(1/0.02), col=4, lty=2,lwd=2)
abline(h=log(1/0.01), col=5, lty=2,lwd=2)

for(k in 1:NN){
  points(theta.true,Umat[k,],pch=19, col=2)
 }

