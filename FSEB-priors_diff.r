UfuncSing <- function(xvec,mvec,abvec){
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

effdiffvec <- c(0, 0.2, 0.5, 0.9)
Neff <- length(effdiffvec)
NN.all <- c(10,50,100)  
nn <- 100 
Powermat <- matrix(0,Neff,3) 
setbeta <- matrix(c(2,2,2,5,5,2),3,2,byrow=T) 
powermat.all <- NULL

for(NN.chosen in 1:3){
    NN <- NN.all[NN.chosen]
    AllBetaMat <- NULL
 
for(betaind in 1:3){
  selbeta <- setbeta[betaind,]
  
  for(eff.ind in 1:Neff){
    effdiff <- effdiffvec[eff.ind] 
    theta.true <- rbeta(nn, selbeta[1], selbeta[2])
    NumAbove0.0005 <- 0
    NumAbove0.005 <- 0
    NumAbove0.05 <- 0
    Umat <- matrix(NA, nn, NN)
  
    for(k in 1: nn){
      mm1 <- sample(seq(15,40, by=1),NN,  replace=T) 
      mm2 <- mm1 
      theta1 <- theta.true[k]
      theta2 <- theta1 + effdiff 
      if(theta2 >= 1){
        if(theta1 - effdiff > 0)
          theta2 <- theta1 - effdiff
      }
      xx <- rbinom(NN, mm1,theta1) 

      if(theta2<1){
        yy <-  rbinom(NN, mm2,theta2) 
         
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
    
    NumAbove0.0005 <- length(which(Umat > log(1/0.0005)))
    NumAbove0.005 <- length(which(Umat > log(1/0.005)))
    NumAbove0.05 <- length(which(Umat > log(1/0.05)))
    totalN <- length(na.omit(Umat))
    Powermat[eff.ind,1:3] <- c(NumAbove0.0005,NumAbove0.005,NumAbove0.05)/totalN
  } 
  AllBetaMat <- rbind(AllBetaMat, Powermat)
 } 
  powermat.all <- rbind(powermat.all, round(AllBetaMat,3))
}
