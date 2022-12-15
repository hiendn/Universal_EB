RfuncBBsing <- function(theta,xv,mv,abvec){
  ## function to calculate R for single \theta_i
  xm <- (mv-xv)
  logR <- lbeta(xv+abvec[1], xm+abvec[2]) -
    lbeta(abvec[1],abvec[2]) - xv*log(theta) - xm*log(1-theta)
  return(logR)
}

nn<-10
mm <- sample(seq(15,40, by=1),nn,  replace=T) 

theta.true <- seq(0.1,0.9,length.out=nn) 
thetavec <- seq(0.0001,0.9999,by=0.0005); ntheta <- length(thetavec)
Rcurvemat <- matrix(NA,ntheta, nn)
xx <- rbinom(nn, mm,theta.true) 

for(j in 1:nn){ 
  xx2 <- xx[j]; mm2 <- mm[j]
  compset <- setdiff(1:nn,j)
  xx1 <- xx[compset]; mm1 <- mm[compset]
  
  n <- length(xx1)
  xm <- xx1/mm1 
  xmbar <- mean(xm)
  xmvar <- (n-1)*var(xm)/n 
  mbar <- mean(mm)
  mu.eb <- xmbar
  phi.eb <- ((mbar*xmvar)/(mu.eb*(1-mu.eb))  - 1)/(mbar - 1)
  
  aa.hat <- (1/phi.eb -1)*mu.eb
  bb.hat <- (1/phi.eb -1)*(1 - mu.eb)
  abvec <- c(aa.hat, bb.hat)
  
  if(abvec[1]>0 & abvec[2]>0){
    Rvec <- rep(0, ntheta)
    for(i in 1:ntheta){
      Rvec[i] <- RfuncBBsing(thetavec[i],xx2,mm2,abvec)
    }
    Rcurvemat[,j] <- Rvec
  }
}


## plotting CIs

testsize <- 0.05 
finalR <- na.omit(t(Rcurvemat))
na.indices <- na.action(na.omit(t(Rcurvemat)))[1]
ylims <- range(finalR) 
ylims <- c(-5,20)

plot(thetavec,Rcurvemat[,1],type="n", ylim=ylims, xlab=expression(theta), 
     ylab=expression(log ~ R[n]), main="", cex.lab=1.5, cex.axis=1.2)

theta.truef <- theta.true[setdiff(1:nn,na.indices)]
nn.f <- nrow(finalR)
abline(h=log(1/testsize), col=1, lty=1,lwd=2)

Allcol <- rainbow(nn.f)
for(j in 1:nn.f){
  lines(thetavec,finalR[j,],lty=2, col=Allcol[j])
  abline(v=theta.truef[j], col=Allcol[j], lwd=1.5)
}





