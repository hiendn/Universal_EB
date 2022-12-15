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
Powermat <- matrix(0,10,3)
effdiffvec <- seq(0,0.9,by=0.1)

for(eff.ind in 1:10){
  
  nn <- 50  
  effdiff <- effdiffvec[eff.ind]  
  NN <- 20  
  Umat <- matrix(NA, nn, NN)
  theta.true <- seq(0.05,0.95,length.out=nn) 
  NumAbove0.01 <- 0
  NumAbove0.02 <- 0
  NumAbove0.05 <- 0
  
  for(k in 1: nn){
    
    mm1 <- sample(seq(15,40, by=1),NN,  replace=T) 
    mm2 <- sample(seq(15,40, by=1),NN,  replace=T) 
    theta1 <- theta.true[k]
    theta2 <- theta1 + effdiff 
    xx <- rbinom(NN, mm1,theta1) 
    summary(xx/mm1-theta1); 
    
    if(theta2<1){
      yy <- rbinom(NN, mm2,theta2)
      summary(yy/mm2-theta2)
      
      for(j in 1:NN){ 
        xx2 <- xx[j]; mx2 <- mm1[j]
        yy2 <- yy[j]; my2 <- mm2[j]
        compset <- setdiff(1:NN,j)
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
  } 
  
  NumAbove0.01 <- length(which(Umat > log(1/0.01)))
  NumAbove0.02 <- length(which(Umat > log(1/0.02)))
  NumAbove0.05 <- length(which(Umat > log(1/0.05)))
  totalN <- length(na.omit(Umat))
  Powermat[eff.ind,] <- c(NumAbove0.01,NumAbove0.02,NumAbove0.05)/totalN
} 

plot(effdiffvec,Powermat[,1], type="n",xlab= expression(theta[2]-theta[1]),
     ylab="Power", main="",cex.lab=1.5,cex.axis=1.2)

lines(effdiffvec, Powermat[,1], col="red", lty=2) 
lines(effdiffvec, Powermat[,2], col="green", lty=3)
lines(effdiffvec, Powermat[,3], col="blue", lty=4) 

legend(0.6,0.6, 
       c(expression(paste(alpha, " = ", 0.01)),
         expression(paste(alpha, " = ", 0.02)),
         expression(paste(alpha, " = ", 0.05))),
       col=c("red","green","blue"), lty=2:4, lwd=1.5,cex=1.2)
