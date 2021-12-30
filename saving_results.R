library(jagsUI)
path1 = file.path('C:/Users/commend/Dropbox/2020_2/simulation/final/data')
path2 = file.path('C:/Users/commend/Dropbox/2020_2/simulation/final/data_rep50')
path3 = file.path('C:/Users/commend/Dropbox/2020_2/simulation/final/data_new')

rslt <- function(nitems, n){
  rs <- round(out$summary,4)
  nitems <- nitems
  n <- n
  n2 <- n/2
  
  
  r1 <- rs[c(1:(nitems*2+4)),]
  mn <- c(paste0("mu[",c(1,(n2+1)),']'))
  r2 <- rs[mn,]
  tn1 <- c(paste0("theta[",1:n2,']'))
  theta1 <- colMeans(rs[tn1,])
  tn2 <- c(paste0("theta[",(n2+1):n,']'))
  theta2 <- colMeans(rs[tn2,])
  mon <- c(paste0("m[",c(1,(n2+1)),']'))
  r3 <- rs[mon,]
  on1 <- c(paste0("omega[",1:n2,']'))
  o1_1 <- colMeans(log(rs[on1,c(1,3:7)]))
  o1_2 <- colMeans(rs[on1,-c(1,3:7)])
  omega1 <- c(o1_1[1],o1_2[1],o1_1[2:6],o1_2[2:3])
  on2 <- c(paste0("omega[",(n2+1):n,']'))
  o2_1 <- colMeans(log(rs[on2,c(1,3:7)]))
  o2_2 <- colMeans(rs[on2,-c(1,3:7)])
  omega2 <- c(o2_1[1],o2_2[1],o2_1[2:6],o2_2[2:3])
  pn <- c("psi")
  psi <- rs[pn,]
  sigma <- sqrt(psi)
  sn <- c("sigma")
  names(sigma) <- sn
  tn <- c(paste0("tau[",rep(1:nitems,5),',',c(rep(1,nitems),rep(2,nitems),rep(3,nitems),rep(4,nitems),rep(5,nitems)),']'))
  r5 <- rs[tn,]
  dn <- c("deviance")
  deviance <- rs[dn,]
  
  
  rslt <- round(rbind(r1,r2,theta1,theta2,r3,omega1,omega2,sigma,r5,deviance),4)
  return(rslt)
}
parm <- function(nitems, path_pr){
  pr <- read.table(path_pr, header=T)
  prs <- matrix(round(pr[,2],4), ncol=1)
  rownames(prs) = pr[,1]
  nitems <- nitems
  
  pr1 <- prs[1:(nitems*2),]
  prb <- c(rep(NA,4))
  pmn <- c(paste0("mu_g",c(1,2)))
  pr2 <- prs[pmn,]
  pr3 <- c(NA, NA)
  pon <- c(paste0("m_omega_g",c(1,2)))
  pomega <- prs[pon,]
  pr4 <- c(NA, NA)
  ppn <- c("psi")
  ppsi <- prs[ppn,]
  psig <- 1/ppsi
  names(psig) <- c('sigma')
  pt <- c(rep(NA,nitems),rep(-0.9,nitems),rep(-0.3,nitems),rep(0.3,nitems),rep(0.9,nitems))
  ptn <- c(paste0("tau[",rep(1:nitems,5),',',c(rep(1,nitems),rep(2,nitems),rep(3,nitems),rep(4,nitems),rep(5,nitems)),']'))
  names(pt) <- ptn
  pr5 <- c(NA)
  
  parm <- round(c(pr1,prb,pr2,pr3,pomega,pr4,psig,pt,pr5),4)
  return(parm)
}
rslt2 <- function(nitems, n){
  ws <- round(out$summary,4)
  nitems <- nitems
  n <- n
  n2 <- n/2
  
  w1 <- ws[c(1:(nitems*2+2)),]
  w1_1 <- c(rep(NA,2))
  mn <- c(paste0("mu[",c(1,(n2+1)),']'))
  w2 <- ws[mn,]
  tn1 <- c(paste0("theta[",1:n2,']'))
  theta1 <- colMeans(ws[tn1,])
  tn2 <- c(paste0("theta[",(n2+1):n,']'))
  theta2 <- colMeans(ws[tn2,])
  w3 <- c(rep(NA,2))
  omega1 <- c(rep(NA,2))
  omega2 <- c(rep(NA,2))
  sigma <- c(NA)
  sn <- c("sigma")
  names(sigma) <- sn
  tn <- c(paste0("tau[",rep(1:nitems,5),',',c(rep(1,nitems),rep(2,nitems),rep(3,nitems),rep(4,nitems),rep(5,nitems)),']'))
  w5 <- ws[tn,]
  dn <- c("deviance")
  deviance <- ws[dn,]
  
  
  wslt <- round(rbind(w1,w1_1,w1_1,w2,theta1,theta2,w3,omega1,omega2,sigma,w1_1,w5,deviance),4)
  return(wslt)
}


n <- 1000
nitems <- 10

mu1 = c(0, 0, 0.3)
mu2 = c(0, 0.3, 0)
momega1 = c(0.3, 1, 0)
momega2 = c(-0.3, -1, 0)
momname1 = c(1.34, 2.71, 1)
momname2 = c(0.74, 0.36, 1)
psi = c(1/0.4)

result <- c()
rst <- c()


### ERS GPCM ###

  setwd(path3)
  repl = c(1:20, 31:50)
    for(b in 2:2){
      for(a in 3:3){
        for(i in repl){
          path <- paste0("mu_",mu1[a],"_",mu2[a],"_momega_",momname1[b],"_",momname2[b],"_out",i,".RData")
          load(path)
          n <- 1000
          nitems <- 10
          rst <- rslt(nitems = nitems, n=n)
          assign(paste0("result_",mu1[a],"_",mu2[a],"_",momname1[b],"_",momname2[b],"_",i),rst)
        }
      }
    }
  setwd(path2)
    repl = c(31:50)
    for(b in 1:length(momega1)){
      for(a in 1:length(mu1)){
        for(i in repl){
          path <- paste0("mu_",mu1[a],"_",mu2[a],"_momega_",momname1[b],"_",momname2[b],"_out",i,".RData")
          load(path)
          n <- 1000
          nitems <- 10
          rst <- rslt(nitems = nitems, n=n)
          assign(paste0("result_",mu1[a],"_",mu2[a],"_",momname1[b],"_",momname2[b],"_",i),rst)
        }
      }
    }
    
  #ls(pattern = "result_0_0_*")
  #ls(pattern = "result_0.3_0_*")


### Param ###

  setwd(path1)
  repl = c(1:20)
  for(y in 2:2){
    for(x in 3:3){
      for(z in repl){
        path_pr <- paste0("mu_",mu1[x],"_",mu2[x],"_momega_",momname1[y],"_",momname2[y],"_param",z,".txt")
        prm <- parm(nitems=10, path_pr=path_pr)
        prm <- as.matrix(prm, ncol=1)
        assign(paste0("param_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y],"_",z),prm)
      }
    }
  }
  
  setwd(path2)
  repl = c(31:50)
  for(y in 1:length(momega1)){
    for(x in 1:length(mu1)){
      for(z in repl){
        path_pr <- paste0("mu_",mu1[x],"_",mu2[x],"_momega_",momname1[y],"_",momname2[y],"_param",z,".txt")
        prm <- parm(nitems=10, path_pr=path_pr)
        prm <- as.matrix(prm, ncol=1)
        assign(paste0("param_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y],"_",z),prm)
      }
    }
  }
  
  
  #ls(pattern = "param*")
  #ls(pattern = "param_0.3_0_*")


### GPCM ###


  setwd(path3)
  repl = c(31:47)
  for(b in 2:2){
    for(a in 3:3){
      for(i in repl){
        path <- paste0("mu_",mu1[a],"_",mu2[a],"_momega_",momname1[b],"_",momname2[b],"_out_gpcm",i,".RData")
        load(path)
        n <- 1000
        nitems <- 10
        rst <- rslt2(nitems = nitems, n=n)
        assign(paste0("result_",mu1[a],"_",mu2[a],"_",momname1[b],"_",momname2[b],"_gpcm_",i),rst)
      }
    }
  }
  setwd(path2)
  repl = c(31:50)
  for(b in 1:length(momega1)){
    for(a in 1:length(mu1)){
      for(i in repl){
        path <- paste0("mu_",mu1[a],"_",mu2[a],"_momega_",momname1[b],"_",momname2[b],"_out_gpcm",i,".RData")
        load(path)
        n <- 1000
        nitems <- 10
        rst <- rslt2(nitems = nitems, n=n)
        assign(paste0("result_",mu1[a],"_",mu2[a],"_",momname1[b],"_",momname2[b],"_gpcm_",i),rst)
      }
    }
  }
  
  #ls(pattern = "_gpcm_*")
  #ls(pattern = "result*")
  
  exsubset = ls()[!grepl("result", ls()) & !grepl("param", ls())]
  rm(list=exsubset)

  
  
###################################
###  Bias & RMSE  #################
###################################

  mu1 = c(0, 0, 0.3)
  mu2 = c(0, 0.3, 0)
  momega1 = c(0.3, 1, 0)
  momega2 = c(-0.3, -1, 0)
  momname1 = c(1.34, 2.71, 1)
  momname2 = c(0.74, 0.36, 1)
  
  reslt <- matrix(nrow=84,ncol=50)
  repli = c(1:20,31:50)
  for(y in 1:length(momega1)){
    for(x in 1:length(mu1)){
      for(z in repli){
        reslt[,z] <- get(paste0("result_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y],"_",z))[,1] - get(paste0("param_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y],"_",z))[,1]
      }
      result <- reslt[-c(21:24,27:28,31:32,34:43,84),c(1:20,31:50)]
      
      bias_a <- sum(colMeans(result[1:10,]))/40
      bias_b <- sum(colMeans(result[11:20,]))/40
      bias_mu_g2 <- sum(result[22,])/40
      bias_momega_g1 <- sum(result[23,])/40
      bias_momega_g2 <- sum(result[24,])/40
      bias_tau <- sum(colMeans(result[26:65,]))/40
      
      bias <- rbind(bias_a, bias_b, bias_mu_g2, bias_momega_g1, bias_momega_g2, bias_tau)
      
      rmse_a <- sqrt(sum(colMeans(result[1:10,]^2))/40)
      rmse_b <- sqrt(sum(colMeans(result[11:20,]^2))/40)
      rmse_mu_g2 <- sqrt(sum(result[22,]^2)/40)
      rmse_momega_g1 <- sqrt(sum(result[23,]^2)/40)
      rmse_momega_g2 <- sqrt(sum(result[24,]^2)/40)
      rmse_tau <- sqrt(sum(colMeans(result[26:65,]^2))/40)
      
      rmse <- rbind(rmse_a, rmse_b, rmse_mu_g2, rmse_momega_g1, rmse_momega_g2, rmse_tau)
      
      r <- rbind(bias,rmse)
      colnames(r) <- paste0("rslt_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y])
      assign(paste0("rt_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y]),r)
    }}
  
  
  
  rslt_ers_gcpm <- cbind(rt_0_0_1_1, rt_0_0_1.34_0.74, rt_0_0_2.71_0.36, rt_0_0.3_1_1, rt_0_0.3_1.34_0.74, rt_0_0.3_2.71_0.36)
  
  #write.csv(rslt_ers_gcpm, file='C:/Users/commend/Dropbox/2020_2/simulation/final/results/ers_gpcm.csv')
  
  
  reslt <- matrix(nrow=84,ncol=50)
  
  for(y in 1:length(momega1)){
    for(x in 1:length(mu1)){
      for(z in repli){
        reslt[,z] <- get(paste0("result_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y],"_gpcm_",z))[,1] - get(paste0("param_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y],"_",z))[,1]
      }
      result <- reslt[-c(21:24,27:28,31:32,34:43,84),c(1:20,31:50)]
      
      bias_a <- sum(colMeans(result[1:10,]))/40
      bias_b <- sum(colMeans(result[11:20,]))/40
      bias_mu_g2 <- sum(result[22,])/40
      bias_tau <- sum(colMeans(result[26:65,]))/40
      
      bias <- rbind(bias_a, bias_b, bias_mu_g2, NA, NA, bias_tau)
      
      rmse_a <- sqrt(sum(colMeans(result[1:10,]^2))/40)
      rmse_b <- sqrt(sum(colMeans(result[11:20,]^2))/40)
      rmse_mu_g2 <- sqrt(sum(result[22,]^2)/40)
      rmse_tau <- sqrt(sum(colMeans(result[26:65,]^2))/40)
      
      rmse <- rbind(rmse_a, rmse_b, rmse_mu_g2, NA, NA, rmse_tau)
      
      r <- rbind(bias,rmse)
      colnames(r) <- paste0("rslt_gpcm_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y])
      assign(paste0("rt_gpcm_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y]),r)
    }}
  
  
  rslt_gpcm <- cbind(rt_gpcm_0_0_1_1, rt_gpcm_0_0_1.34_0.74, rt_gpcm_0_0_2.71_0.36, rt_gpcm_0_0.3_1_1, rt_gpcm_0_0.3_1.34_0.74, rt_gpcm_0_0.3_2.71_0.36)
  
  
  #write.csv(rslt_gpcm, file='C:/Users/commend/Dropbox/2020_2/simulation/final/results/gpcm.csv')
  
  
  
  ### ERS-GPCM & GPCM
  
  rslt <- cbind(rt_0_0_1_1, rt_gpcm_0_0_1_1, rt_0_0_1.34_0.74, rt_gpcm_0_0_1.34_0.74, rt_0_0_2.71_0.36, rt_gpcm_0_0_2.71_0.36, rt_0_0.3_1_1, rt_gpcm_0_0.3_1_1, rt_0_0.3_1.34_0.74, rt_gpcm_0_0.3_1.34_0.74, rt_0_0.3_2.71_0.36, rt_gpcm_0_0.3_2.71_0.36)
  
  rslt <- cbind(rt_0.3_0_2.71_0.36, rt_gpcm_0.3_0_2.71_0.36)
  
  write.csv(rslt, file='C:/Users/commend/Dropbox/2020_2/simulation/final/results_new/bias_rmse_40.csv')
  
  
  
  
  
###################################
###  variance  #################
###################################


library('matrixStats')


###1. ERS-GPCM

reslt <- matrix(nrow=84,ncol=50)

for(y in 1:length(momega1)){
  for(x in 1:length(mu1)){
    for(z in repli){
      reslt[,z] <- get(paste0("result_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y],"_",z))[,1] - get(paste0("param_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y],"_",z))[,1]
    }
    result <- reslt[-c(21:24,27:28,31:32,34:43,84),c(1:20,31:50)]
    
    var_a <- sum(rowVars(result[1:10,]))/10
    var_b <- sum(rowVars(result[11:20,]))/10
    var_mu_g2 <- var(result[22,])
    var_momega_g1 <- var(result[23,])
    var_momega_g2 <- var(result[24,])
    var_tau <- sum(rowVars(result[26:65,]))/40
    
    var <- rbind(var_a, var_b, var_mu_g2, var_momega_g1, var_momega_g2, var_tau)
    
    r <- var
    colnames(r) <- paste0("Var_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y])
    assign(paste0("var_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y]),r)
  }}


var_ers_gcpm <- cbind(var_0_0_1_1, var_0_0_1.34_0.74, var_0_0_2.71_0.36, var_0_0.3_1_1, var_0_0.3_1.34_0.74, var_0_0.3_2.71_0.36)

for(y in 1:length(momega1)){
  for(x in 1:length(mu1)){
    for(z in repli){
      reslt[,z] <- get(paste0("result_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y],"_gpcm_",z))[,1] - get(paste0("param_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y],"_",z))[,1]
    }
    result <- reslt[-c(21:24,27:28,31:32,34:43,84),c(1:20,31:50)]
    
    var_a <- sum(rowVars(result[1:10,]))/10
    var_b <- sum(rowVars(result[11:20,]))/10
    var_mu_g2 <- var(result[22,])
    var_tau <- sum(rowVars(result[26:65,]))/40
    
    var <- rbind(var_a, var_b, var_mu_g2, NA, NA, var_tau)
    
    r <- var
    colnames(r) <- paste0("Var_gpcm_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y])
    assign(paste0("var_gpcm_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y]),r)
  }}


var_gpcm <- cbind(var_gpcm_0_0_1_1, var_gpcm_0_0_1.34_0.74, var_gpcm_0_0_2.71_0.36, var_gpcm_0_0.3_1_1, var_gpcm_0_0.3_1.34_0.74, var_gpcm_0_0.3_2.71_0.36)



### ERS-GPCM & GPCM

var <- cbind(var_0_0_1_1, var_gpcm_0_0_1_1, var_0_0_1.34_0.74, var_gpcm_0_0_1.34_0.74, var_0_0_2.71_0.36, var_gpcm_0_0_2.71_0.36, var_0_0.3_1_1, var_gpcm_0_0.3_1_1, var_0_0.3_1.34_0.74, var_gpcm_0_0.3_1.34_0.74, var_0_0.3_2.71_0.36, var_gpcm_0_0.3_2.71_0.36)

var <- cbind(var_0.3_0_2.71_0.36, var_gpcm_0.3_0_2.71_0.36)

write.csv(var, file='C:/Users/commend/Dropbox/2020_2/simulation/final/results_new/var_40.csv')




###################################
###  False & True Positive  #######
###################################


library('matrixStats')


###1. ERS-GPCM

reslt1 <- matrix(nrow=50,ncol=11)
reslt2 <- matrix(nrow=50,ncol=11)


for(y in 2:2){
  for(x in 3:3){
    for(z in repl){
      reslt1[z,] <- get(paste0("result_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y],"_",z))[22,]
      reslt2[z,] <- get(paste0("result_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y],"_",z))[24,]
    }
    colnames(reslt1) = colnames(result_0.3_0_2.71_0.36_31)
    colnames(reslt2) = colnames(result_0.3_0_2.71_0.36_31)
    
    reslt1 <- reslt1[c(1:20, 31:50),]
    reslt2 <- reslt2[c(1:20, 31:50),]
    
    b11_overlap0 <- sum(reslt1[,10], na.rm = T)
    b21_overlap0 <- sum(reslt2[,10], na.rm = T)
    
    overlap <- rbind(b11_overlap0, b21_overlap0)
    
    r <- overlap
    colnames(r) <- paste0("Overlap_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y])
    assign(paste0("overlap0_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y]),r)
  }}


overlap0_ers_gcpm <- cbind(overlap0_0_0_1_1, overlap0_0_0_1.34_0.74, overlap0_0_0_2.71_0.36, overlap0_0_0.3_1_1, overlap0_0_0.3_1.34_0.74, overlap0_0_0.3_2.71_0.36)


###2. GPCM


reslt3 <- matrix(nrow=50,ncol=11)
repli = c(1:20, 31:50)

for(y in 2:2){
  for(x in 3:3){
    for(z in repl){
      reslt3[z,] <- get(paste0("result_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y],"_gpcm_",z))[22,]
    }
    colnames(reslt3) = colnames(result_0.3_0_2.71_0.36_31)
    
    b11_overlap0 <- sum(reslt3[,10], na.rm = T)
    b21_overlap0 <- NA
    
    overlap <- rbind(b11_overlap0, b21_overlap0)
    
    r <- overlap
    
    colnames(r) <- paste0("Overlap_gpcm_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y])
    assign(paste0("overlap0_gpcm_",mu1[x],"_",mu2[x],"_",momname1[y],"_",momname2[y]),r)
  }}


overlap0_gpcm <- cbind(overlap0_gpcm_0_0_1_1, overlap0_gpcm_0_0_1.34_0.74, overlap0_gpcm_0_0_2.71_0.36, overlap0_gpcm_0_0.3_1_1, overlap0_gpcm_0_0.3_1.34_0.74, overlap0_gpcm_0_0.3_2.71_0.36)



### ERS-GPCM & GPCM

overlap0 <- cbind(overlap0_0_0_1_1, overlap0_gpcm_0_0_1_1, overlap0_0_0_1.34_0.74, overlap0_gpcm_0_0_1.34_0.74, overlap0_0_0_2.71_0.36, overlap0_gpcm_0_0_2.71_0.36, overlap0_0_0.3_1_1, overlap0_gpcm_0_0.3_1_1, overlap0_0_0.3_1.34_0.74, overlap0_gpcm_0_0.3_1.34_0.74, overlap0_0_0.3_2.71_0.36, overlap0_gpcm_0_0.3_2.71_0.36)


overlap0 <- cbind(overlap0_0.3_0_2.71_0.36, overlap0_gpcm_0.3_0_2.71_0.36)
reslt1
reslt2
reslt3

write.csv(overlap0, file='C:/Users/commend/Dropbox/2020_2/simulation/final/results_new/overlap0_40.csv')


