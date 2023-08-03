###################################################
# Let (\lambda_1, \lambda_2, \cdots, \lambda_n\right) be the beta-Jacobi ensemble J_n(p_1, p_2).
# Set X_i = \frac{\left(p_1+p_2\right) \lambda_i-p_1}{\sqrt{np_1 }}.
# In this R code, we construct  an asymptotically efficient  importance Sampling estimator for P(x_{(n)} > x).
# Let X_{(1)}< \cdots <X_{(n)} be the order statistics of X_1, \cdots, X_n.
###################################################



###########################################
Ratefun <- function(beta, n, p1, p2, x){
  if( n/p1 < 0.01 ){
    gam <- 0
    sig <- p1/p2
    r <- beta*(n-1)*(1+sig)*sqrt(x^2-4/(1+sig))
    return( (n-1)*r/n  )
  } else if(  p1/p2 < 0.01){
    gam <- n/p1
    sig <- 0
    r <- beta*(n-1)*sqrt((gam-x)^2-4)/(1+sqrt(gam)*x)
    return((n-1)*r/n )
  }else{
    gam <- n/p1
    sig <- p1/p2
    u1 <- ((1-sig)*sqrt(gam) - 2*sqrt(1+sig-sig*gam))/(1+sig)
    u2 <- ((1-sig)*sqrt(gam) + 2*sqrt(1+sig-sig*gam))/(1+sig)
    r <-2*beta*(n-1)*gam*sig*sqrt((u2-x)*(u1-x))/(1+sqrt(gam)*x)/(1-sqrt(gam)*sig*x)
    return((n-1)*r/n )
  }
}

Ffun <- function(beta, n, p1 ,p2, x, xn,xn1, eigenvalues){
  if(xn < max(xn1,x)){
    return(0)
  }else{
    s1n <- sqrt(n*p1)/p1
    s2n <- sqrt(n*p1)/p2
    r1n <- beta*(p1-n+1)/2
    r2n <- beta*(p2-n+1)/2
    logbn <- log(n) + (beta*(n-1)+r2n)*log(s1n/(s1n+s2n)) + (beta*(n-1)+r1n)*log(s2n/(s1n+s2n)) + lgamma(1+beta/2) + lgamma(beta*(p1+p2)/2) + lgamma(beta*(p1+p2-1)/2) - lgamma(1+beta*n/2) - lgamma(beta*p1/2) - lgamma(beta*p2/2)-  - lgamma(beta*(p1+p2-n)/2)
    logunx <- (r1n-1)*log(1+s1n*x) + (r2n-1)*log(1 - s2n*x)
    r <- Ratefun(beta, n, p1, p2, x)
    loghnx <- log(r) - r*xn + r* max(xn1,x)
    prob <- logbn + sum(beta*log(xn-eigenvalues)) + logunx - loghnx
    #    return(exp(-prob))
    return(prob/2)
  }
}

# \widetilde{\mathcal{J}}_{n}^{p_1, p_2} :=\frac{(p_1+p_2)\mathcal{J}_{n}^{p_1, p_2}-p_1  I_{n-1}}{\sqrt{n p_1}}
RescaleJ <- function(beta, n, p1, p2){
  J <- array(0, c(n,n))
  c <- rep(0,n)
  s <- rep(0,n)
  for (i in 1:n){
    c[i] <- rbeta(1,beta*(p1-i+1)/2 , beta*(p2-i+1)/2 )
    s[i] <- rbeta(1,beta*(n-i)/2 , beta*(p1+p2 -n -i+1)/2 )
  }
  
  for (i in 2:n){
    J[i,i] <-((p1+p2+2)*(s[i-1]*(1-c[i-1]) + c[i]*(1-s[i-1]))-(p1+1))/sqrt((n+1)*(p1+1)) 
  }
  J[1,1] <-  ((p1+p2+2)*c[1]-p1)/sqrt((n+1)*(p1+1)) 
  for (i in 3:n){
    J[i-1,i] <- sqrt(c[i-1]*(1- c[i-1])*s[i-1]*(1-s[i-2]))*(p1+p2+2)/sqrt((n+1)*(p1+1)) 
    J[i,i-1] <- J[i-1,i] 
  }
  J[1,2] <- sqrt(c[1]*(1- c[1])*s[1])  *(p1+p2+2)/sqrt((n+1)*(p1+1)) 
  J[2,1] <- J[1,2] 
  return(J)
}


EffSim <- function(beta,n,p1, p2,x, iter_num){
  LL <- rep(0, iter_num)
  
  for (iter in 1: iter_num){
    J <- RescaleJ(beta, n-1,p1-1,p2-1)
    eigenvalues <- eigen(J)$values
    xn1 <-eigenvalues[1]
    
    r <- Ratefun(beta, n, p1, p2, x)
    
    xn <- rexp(1,r) + max(xn1,x)
    
    LL[iter] <- Ffun(1, n, p1, p2 ,x,xn, xn1, eigenvalues)
  }
  return(c(mean(LL), 
           sqrt(mean(LL^2)-mean(LL)^2),
           sqrt(mean(LL^2)-mean(LL)^2)/mean(LL)))
}



####################################################
#Simulation
beta <- 2 
x <- 3
iter_num <- 10000
nseq <- c(seq(50,400,50))

y1 <- array(0, c(length(nseq),4 ))
for ( i in 1: length(nseq)  ){
  n <- nseq[i]
  p1 <- 3*n
  p2 <- 3*p1
  
  y1[i,1] <- n
  y1[i,2:4] <- EffSim(beta,n,p1, p2,x, iter_num)
  
}
#############
y2 <- array(0, c(length(nseq),4 ))
for ( i in 1: length(nseq)  ){
  n <- nseq[i]
  p1 <- 50*n
  p2 <- 3*p1
  
  y2[i,1] <- n
  y2[i,2:4] <- EffSim(beta,n,p1, p2,x, iter_num)
}
#############
y3 <- array(0, c(length(nseq),4 ))
for ( i in 1: length(nseq)  ){
  n <- nseq[i]
  p1 <- 3*n
  p2 <- 50*p1
  
  y3[i,1] <- n
  y3[i,2:4] <- EffSim(beta,n,p1, p2,x, iter_num)
}
#############
y4 <- array(0, c(length(nseq),4 ))
for ( i in 1: length(nseq)  ){
  n <- nseq[i]
  p1 <- 50*n
  p2 <- 50*p1
  
  y4[i,1] <- n
  y4[i,2:4] <- EffSim(beta,n,p1, p2,x, iter_num)
}

