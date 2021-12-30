library(jagsUI)

path ='C:/Users/ProDesk/Dropbox/2020_2/simulation/final/data_new'
setwd(path)

analysis <- function(filename, iter, burnin){
  it <- iter 
  bn <- burnin
  df1 <- read.table(file=filename, header=T)
  df <- as.data.frame(df1)
  
  Y <- as.matrix(df[,-c(1:2)])
  n <- nrow(Y)
  p <- ncol(Y)
  K <- apply(Y, 2, max, na.rm=TRUE)
  group <- as.numeric(df$group)
  
  data <- list("Y","n","p","K","group")
  
  
  modfile <- tempfile()
  
  writeLines("
model{
  for(i in 1:n){
    for(j in 1:p){
      Y[i, j]~dcat(prob[i,j,1:K[j]])
    }
    theta[i]~dnorm(mu[i],1)
    mu[i] <- b10 + b11 * group[i]
    
    res[i] <- theta[i] - mu[i]
    the.new[i]~dnorm(mu[i], 1)
    res.new[i] <- the.new[i] - mu[i]
  }
  b10 <- 0.00
  b11 ~ dnorm(0.0, 100)
  
  for(i in 1:n){
    for(j in 1:p){
      for(k in 1:K[j]) {
        eta[i, j, k] <- alpha[j] * (theta[i] - (b[j]+tau[j, k]))
        psum[i, j, k] <- sum(eta[i, j, 1:k])
        exp.psum[i, j, k] <- exp(psum[i,  j, k])
        prob[i, j, k] <- exp.psum[i,j,k] / sum(exp.psum[i,j,1:K[j]])
      }
    }
  }
  
  for (j in 1:p){
    alpha[j] ~ dlnorm(0, pow((1/sqrt(2)),2))
    b[j] ~ dnorm(0, pow((1/sqrt(2)),2))
    tau[j, 1] <- 0.0
    for(k in 2:4) {
      tau[j, k] ~ dnorm(0, pow((1/sqrt(10)),2))
    }
    tau[j, 5] <- -(tau[j, 2] + tau[j, 3] + tau[j, 4])
  }
  
  fit <- sum(res[])
  fit.new <- sum(res.new[])
  
}
", con=modfile)
  
  
  param <- c('alpha','b','b10','b11','tau','mu','theta','fit','fit.new')
  inits <- function(){list(alpha=rep(1, p), b=rep(0, p), tau=matrix(data=(c(rep(NA, p),rep(0, 3*p),rep(NA, p))),nrow=p,ncol=5))}
  
  out <- jags(data = data,
              inits = inits,
              parameters.to.save = param,
              model.file = modfile,
              n.chains = 2,
              n.iter = it,
              n.burnin = bn,
              n.thin = 10)
  return(out)
}

mu1 = c(0, 0.3)
mu2 = c(0, 0)
momega1 = c(0.3, 1, 0)
momega2 = c(-0.3, -1, 0)
momname1 = c(1.34, 2.71, 1)
momname2 = c(0.74, 0.36, 1)

rep = c(1:10)
b = 2

for(i in rep){
  for(a in 2:2){
    filename = paste0("mu_",mu1[a],"_",mu2[a],"_momega_",momname1[b],"_",momname2[b],"_data",i,".txt")
    out <- analysis(filename=filename, iter=15000, burnin=5000)
    path = paste0("mu_",mu1[a],"_",mu2[a],"_momega_",momname1[b],"_",momname2[b],"_out_gpcm",i,".RData")
    save.image(path)}
}
