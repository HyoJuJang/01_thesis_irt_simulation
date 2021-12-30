library(mirtCAT)

###simul function

simul <- function(seed1, seed2, theta1_m, theta2_m, theta_sd, omega1_m, omega2_m, psi, nitems, n){
  
  n <- n
  n2 <- n/2
  sigma <- 1/psi
  seed1 <- seed1
  seed2 <- seed2
  
  set.seed(seed1)
  theta_g1 <- matrix(rnorm(n2, theta1_m, theta_sd))
  theta_g2 <- matrix(rnorm(n2, theta2_m, theta_sd))
  theta <- rbind(theta_g1, theta_g2)
  omega_g1 <- matrix(rlnorm(n2, omega1_m, sigma))
  omega_g2 <- matrix(rlnorm(n2, omega2_m, sigma))
  omega <- rbind(omega_g1, omega_g2)
  
  set.seed(seed2)
  nitems <- nitems
  alpha <- rlnorm(nitems, 0, 0.3) 
  b <- runif(nitems, -2, 2)
  tau=matrix(data=(c(rep(-0.9, nitems),rep(-0.3, nitems),rep(0.3, nitems),rep(0.9, nitems))),nrow=nitems,ncol=4)
  
  q <- matrix(ncol=5, nrow=nitems)
  p <- matrix(ncol=5, nrow=nitems)
  denom <- c()
  pp1 <- list()
  prob.list <- list()
  
  i = 1
  for(j in 1:nitems){
    q[j, 1] <- c(1)
    q[j, 2] <- q[j, 1]*exp(alpha[j] * (theta[i] - (b[j] + omega[i] * tau[j, 1])))
    q[j, 3] <- q[j, 2]*exp(alpha[j] * (theta[i] - (b[j] + omega[i] * tau[j, 2])))
    q[j, 4] <- q[j, 3]*exp(alpha[j] * (theta[i] - (b[j] + omega[i] * tau[j, 3])))
    q[j, 5] <- q[j, 4]*exp(alpha[j] * (theta[i] - (b[j] + omega[i] * tau[j, 4])))
    
    denom[j] <- q[j, 1] + q[j, 2] + q[j, 3] + q[j, 4] + q[j, 5]
    
    p[j, 1] <- q[j, 1]/denom[j]
    p[j, 2] <- q[j, 2]/denom[j]
    p[j, 3] <- q[j, 3]/denom[j]
    p[j, 4] <- q[j, 4]/denom[j]
    p[j, 5] <- q[j, 5]/denom[j]
    
    pp1[[j]] <- p[j, ]
  }
  
  prob.list <- list()
  for(s in 1:nitems){
    prob.list[[s]] <- pp1[[s]]
  }
  
  for(i in 2:n){
    for(j in 1:nitems){
      q[j, 1] <- c(1)
      q[j, 2] <- q[j, 1]*exp(alpha[j] * (theta[i] - (b[j] + omega[i] * tau[j, 1])))
      q[j, 3] <- q[j, 2]*exp(alpha[j] * (theta[i] - (b[j] + omega[i] * tau[j, 2])))
      q[j, 4] <- q[j, 3]*exp(alpha[j] * (theta[i] - (b[j] + omega[i] * tau[j, 3])))
      q[j, 5] <- q[j, 4]*exp(alpha[j] * (theta[i] - (b[j] + omega[i] * tau[j, 4])))
      
      denom[j] <- q[j, 1] + q[j, 2] + q[j, 3] + q[j, 4] + q[j, 5]
      
      p[j, 1] <- q[j, 1]/denom[j]
      p[j, 2] <- q[j, 2]/denom[j]
      p[j, 3] <- q[j, 3]/denom[j]
      p[j, 4] <- q[j, 4]/denom[j]
      p[j, 5] <- q[j, 5]/denom[j]
      
      pp1[[j]] <- p[j, ]
    }
    
    for(s in 1:nitems){
      prob.list[[s]] <- rbind(prob.list[[s]],pp1[[s]])
    }
  }
  
  data0 <- simdata(prob.list = prob.list)
  data1 <- data0+1
  data <- cbind(c(rep(0,n2),rep(1,n2)),data1)
  
  a <- as.matrix(alpha, nrow=1)
  aname <- paste0('a',1:nitems)
  df <- cbind(aname, a)
  b <- as.matrix(b, nrow=1)
  bname <- paste0('b',1:nitems)
  bb <- cbind(bname, b)
  mu <- as.matrix(c(theta1_m,theta2_m), nrow=1)
  muname <- c('mu_g1','mu_g2')
  mm <- cbind(muname, mu)
  m_omega <- as.matrix(c(omega1_m,omega2_m), nrow=1)
  mname <- c('m_omega_g1','m_omega_g2')
  oo <- cbind(mname, m_omega)
  psi <- as.matrix(psi, nrow=1)
  psiname <- 'psi'
  pp <- cbind(psiname, psi)
  tauuu <- c(rep(c(-0.9,-0.3,0.3,0.9),nitems))
  tauname <- paste0('tau',rep(1:nitems, each=4))
  tname <- paste0(tauname,'_',1:4)
  tt <- cbind(tname, tauuu)
  thname <- paste0('theta',1:n)
  the <- as.matrix(theta, nrow=1)
  thth <- cbind(thname, the)
  omename <- paste0('omega',1:n)
  ome <- as.matrix(omega, nrow=1)
  omom <- cbind(omename, ome)
  df <- rbind(df, bb, mm, oo, pp, tt, thth, omom)
  param<- as.data.frame(df)
  
  return(list(data, param))
}

###dir
dir ='C:/Users/commend/Dropbox/2020_2/simulation/final/data_new'
setwd(dir)


###conditions
itemseed <- 8520
seed2 <- itemseed
n <- 1000
nitem <- 10

mu1 = c(0.3, 0)
mu2 = c(0, 0)
momega1 = c(0.3, 1, 0)
momega2 = c(-0.3, -1, 0)
momname1 = c(1.34, 2.71, 1)
momname2 = c(0.74, 0.36, 1)
psi = c(1/0.4)
#seed = c(seq(from=27, to=2700, by=27)) #for rep 1:20
seed = c(seq(from=28, to=2800, by=28))  #for rep 31:50

filename = c()
filename2 = c()


###generate data

for(b in 2:2){
  for(a in 1:1){
    for(i in 31:50){
      filename = paste0("mu_",mu1[a],"_",mu2[a],"_momega_",momname1[b],"_",momname2[b],"_data",i,".txt")
      filename2 = paste0("mu_",mu1[a],"_",mu2[a],"_momega_",momname1[b],"_",momname2[b],"_param",i,".txt")
      mydata <- simul(seed1=seed[i], seed2=itemseed, theta1_m=mu1[a], theta2_m=mu2[a], theta_sd=1, omega1_m=momega1[b], omega2_m=momega2[b], psi=psi, nitems=nitem, n=n)
      dt0 <- mydata[[1]]
      param <- mydata[[2]] 
      id <- c(1:nrow(dt0))
      dt <- cbind(id,dt0)
      colnames(dt)[2] = c('group')
      write.table(x=dt,file=filename, row.names=F)
      write.table(x=param,file=filename2, row.names=F)
    }
  }}

