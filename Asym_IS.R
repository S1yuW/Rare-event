###################################################
# Let (\lambda_1, \lambda_2, \cdots, \lambda_n\right) be the beta-Jacobi ensemble J_n(p_1, p_2).
# Set X_i = \frac{\left(p_1+p_2\right) \lambda_i-p_1}{\sqrt{np_1 }}.
# In this R code, we construct  an asymptotically efficient  importance Sampling estimator for P(x_{(n)} > x).
# Let X_{(1)}< \cdots <X_{(n)} be the order statistics of X_1, \cdots, X_n.
###################################################


# Generate the tridiagonal matrix J. (The eigenvalues of J) ~ J_n(p_1, p_2).
FunJ <- function(beta, n, p1, p2) {
  # Initialize matrix J
  J <- matrix(0, n, n)
  
  # Initialize vectors c and s
  c <- numeric(n)
  s <- numeric(n)
  
  # Generate random values for c and s
  for (i in 1:n) {
    c[i] <- rbeta(1, beta * (p1 - i + 1) / 2, beta * (p2 - i + 1) / 2)
    s[i] <- rbeta(1, beta * (n - i) / 2, beta * (p1 + p2 - n - i + 1) / 2)
  }
  
  # Fill in the diagonal and off-diagonal elements of J
  for (i in 2:n) {
    J[i, i] <- s[i - 1] * (1 - c[i - 1]) + c[i] * (1 - s[i - 1])
  }
  J[1, 1] <- c[1]
  
  for (i in 3:n) {
    J[i - 1, i] <- sqrt(c[i - 1] * (1 - c[i - 1]) * s[i - 1] * (1 - s[i - 2]))
    J[i, i - 1] <- J[i - 1, i]
  }
  
  J[1, 2] <- sqrt(c[1] * (1 - c[1]) * s[1])
  J[2, 1] <- J[1, 2]
  
  return(J)
}

###########################################
Ratefun <- function(beta, n, p1, p2, x) {
  if (n/p1 < 0.01) {
    gam <- 0
    sig <- p1/p2
    J <- (1 + sig) * sqrt(x^2 - 4/(1 + sig)) / 2
    return(beta * (n - 1) * J)
  } else if (p1/p2 < 0.01) {
    gam <- n/p1
    sig <- 0
    J <- sqrt((sqrt(gam) - x)^2 - 4) / (2 * (1 + sqrt(gam) * x))
    return(beta * (n - 1) * J)
  } else {
    gam <- n/p1
    sig <- p1/p2
    u1 <- ((1 - sig) * sqrt(gam) - 2 * sqrt(1 + sig - sig * gam)) / (1 + sig)
    u2 <- ((1 - sig) * sqrt(gam) + 2 * sqrt(1 + sig - sig * gam)) / (1 + sig)
    J <- gam * sig * sqrt((u2 - x) * (u1 - x)) / (1 + sqrt(gam) * x) / (1 - sqrt(gam) * sig * x)
    return(beta * (n - 1) * J)
  }
}




#Calculate Fn
Ffun <- function(beta, n, p1, p2, x) {
  J <- FunJ(beta, n - 1, p1 - 1, p2 - 1)
  eigenvalues <- eigen(J)$values
  eigenvalues <- ((p1 + p2) * eigenvalues - p1) / sqrt(n * p1)
  xn1 <- eigenvalues[1] 
  r <- Ratefun(beta, n, p1, p2, x)
  
  while (TRUE) {
    xn <- rexp(1, rate = r) + max(xn1, x)
    if (xn < p2 / sqrt(n * p1)) {
      break
    }
  }
  
  s1n <- sqrt(n * p1) / p1
  s2n <- sqrt(n * p1) / p2
  r1n <- beta * (p1 - n + 1) / 2
  r2n <- beta * (p2 - n + 1) / 2
  logbn <- log(n) + (beta * (n - 1) + r2n) * log(s1n / (s1n + s2n)) + (beta * (n - 1) + r1n) * log(s2n / (s1n + s2n)) + lgamma(1 + beta / 2) + lgamma(beta * (p1 + p2) / 2) + lgamma(beta * (p1 + p2 - 1) / 2) - lgamma(1 + beta * n / 2) - lgamma(beta * p1 / 2) - lgamma(beta * p2 / 2) - lgamma(beta * (p1 + p2 - n) / 2)
  logunx <- (r1n - 1) * log(1 + s1n * xn) + (r2n - 1) * log(1 - s2n * xn)
  r <- Ratefun(beta, n, p1, p2, x)
  loghnx <- log(r) - r * (xn - max(xn1, x)) - log(1 - exp(-r * (1 / s2n - max(xn1, x))))
  logf <- logbn + sum(beta * log(xn - eigenvalues)) + logunx - loghnx
  return(exp(logf))
}


# Simulate iter_num i.i.d. copies using Ffun -> LL[iter]
CopiesIS <- function(beta, n, p1, p2, x, iter_num) {
  LL <- numeric(iter_num)
  
  for (iter in 1:iter_num) {
    if (is.na(x)) {
      # Handle case where x is NA
      break
    }
    
    LL[iter] <- Ffun(beta, n, p1, p2, x)
  }
  
  return(LL)
}

####### ggplot
# Given beta, n, p1, p2, and different x, the transformation curve of C.O.V. with respect to N
# N <- c(seq(N1, N2, gap))

COVNSim <- function(beta, n, p1, p2, x, N1, N2, gap) {
  Nseq <- seq(N1, N2, gap)
  CN <- matrix(0, nrow = length(Nseq), ncol = 2)
  
  for (i in 1:length(Nseq)) {
    iter_num <- Nseq[i]
    CN[i, 1] <- iter_num
    LL <- CopiesIS(beta, n, p1, p2, x, iter_num)
    CN[i, 2] <- sqrt(mean(LL^2) - mean(LL)^2) / sqrt(iter_num) / mean(LL)
  }
  
  return(CN)
}


####################################################
#Simulation

library(ggplot2)

library(ggpubr)
#Given beta,n,p1,p2 and different x, the transformation curve of C.O.V. with respect to N
#  N <- c(seq(N1,N2,gap))

beta <- 2
n <- 10
p1 <- 20
p2 <- 40
N1 <- 500
N2 <- 10000
gap <- 500

x <- c(seq(0.75,0.85,0.02))
y <- array(0, c(length(x)*length(c(seq(N1,N2,gap))), 3))


for (i in 1:length(x)) {
  i1 <- (i-1)*length(c(seq(N1,N2,gap)))+1
  i2 <- i*length(c(seq(N1,N2,gap)))
  y[i1:i2,1] <- c(seq(N1,N2,gap))
  a <- ((p1 + p2) * x[i] - p1) / sqrt(n * p1)
  
  CNSim <- COVNSim(beta,n,p1, p2,a,N1, N2,gap)
  y[i1:i2,2]  <- CNSim[,2]
  y[i1:i2,3] <- x[i]
}

data <- data.frame(y)
colnames(data) <- c('N','C.O.V.','x')


cols <- c(rgb(168,3,38, maxColorValue = 255),rgb(251,132,2, maxColorValue = 255),rgb(255,195,0, maxColorValue = 255),rgb(95,198,201, maxColorValue = 255),rgb(170,89,246, maxColorValue = 255),rgb(0,0,3, maxColorValue = 255) )          
data$x <- as.character(data$x)

p <- ggplot(data, aes(N, C.O.V.)) +
  geom_point(aes(color = x), size = 2, shape = 16) +
  theme_bw() +
  scale_color_manual(values = cols)

p

ggsave('COVN.png', plot = p, width = 15, height = 5, dpi = 600 )



##################################

#################
beta <- 2
n <- 10
p1 <- 20
p2 <- 40
x <- c(seq(0.75,0.85,0.001))
iter_num <- 5000
y <- array(0, c(length(x),2 ) )
for (i in 1: length(x) ){

  y[i,1] <- x[i]
  a <- ((p1 + p2) * x[i] - p1) / sqrt(n * p1)
  
  LL <- CopiesIS(beta,n,p1, p2,a,iter_num)
  y[i,2]  <- mean(LL)
}
data <- data.frame(y)
colnames(data) <- c('x','logP')


p<-ggplot(data, aes(x, logP)) +
  geom_point(size = 0.6) +
  theme_bw() +
  scale_color_manual(values = cols) +
  scale_y_log10()


ggsave('logP.png', plot = p, width = 5, height = 5, dpi = 600 )
##############################
