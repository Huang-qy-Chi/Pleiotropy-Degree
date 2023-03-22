rm(list = ls())
library(LaplacesDemon)
library(foreach)
library(doParallel)
library(beepr)
setwd("C:/Users/Administrator/Desktop/Pleio Degree/programme")
source("Gen_data.R")
source("Vec.R")
source("Msq.R")
source("DGP.R")

S1 = 0.2*matrix(rep(1,100),ncol = 10)+0.8*diag(10)
S2 = 0.5*matrix(rep(1,100),ncol = 10)+0.5*diag(10)
S3 = 0.8*matrix(rep(1,100),ncol = 10)+0.2*diag(10)

cl <- makeCluster(6)
registerDoParallel(cl)
RR <- foreach(i = 1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
  dat = Gen_data(n = 500, beta = c(0.25,0.25,0,0,rep(0,6)), Sigma0 = S3)
  y = dat$Y
  x = dat$X
  DEDGP(y,x,alpha = 0.05)$DGP
}
stopImplicitCluster()
stopCluster(cl)
1-max(table(RR))/1000
beep(1)
table(RR)


S1 = 0.2*matrix(rep(1,16),ncol = 4)+0.8*diag(4)
S2 = 0.5*matrix(rep(1,16),ncol = 4)+0.5*diag(4)
S3 = 0.8*matrix(rep(1,16),ncol = 4)+0.2*diag(4)
cl <- makeCluster(6)
registerDoParallel(cl)
RR <- foreach(i = 1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
  dat = Gen_data(n = 1000, beta = c(1,1,1,0), Sigma0 = S1)
  y = dat$Y
  x = dat$X
  DEDGP(y,x,alpha = 0.01)$DGP
}
stopImplicitCluster()
stopCluster(cl)
1-max(table(RR))/1000
beep(1)
mean(RR)




################################################################


cl <- makeCluster(6)
registerDoParallel(cl)
RS = matrix(rep(1,8),ncol = 4)
n = c(500,1000)
corr = c(0.2,0.5,0.8)
alpha = c(0.05,0.01)
Beta = matrix(c(rep(0,4),c(1,rep(0,3)),c(rep(1,2),rep(0,2)),c(rep(1,3),0)),ncol = 4)

#Normal,q=4
for(p in 1:2){ #alpha
  for(r in 1:3){ #COrr
    S = corr[r]*matrix(rep(1,16),ncol = 4)+(1-corr[r])*diag(4)
    for(m in 1:2){ #sample size
      for(u in 1:4){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
          #result = Gen_data3(n = n[m], beta = Beta[,u],mu = c(3,3,3), Sigma1 = S,Sigma2 = 3*diag(3))
          #result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          y = result$Y
          x = result$X
          DEDGP(y,x,alpha = alpha[p])$DGP
        }
        RS[m,u]=1-max(table(RR))/1000
      }
    }
    rownames(RS) = c("n=500","n=1000")
    colnames(RS) = c("K^*=0","K^*=1","K^*=2","K^*=3")  
    write.table(RS,sep=",",file =paste("C:/Users/Administrator/Desktop/Pleio Degree/DGP-4-",100*alpha[p],"/PDEDGP",10*corr[r],".csv",sep=""),row.names = T)
    beep(1)
    Sys.sleep(20)
  }
}

#q=10
RS = matrix(rep(1,8),ncol = 4)
Beta = matrix(c(rep(0,10),c(1,rep(0,9)),c(rep(1,2),rep(0,8)),c(rep(1,3),rep(0,7))),ncol = 4)
#Normal
for(p in 1:2){ #alpha
  for(r in 1:3){ #COrr
    S = corr[r]*matrix(rep(1,100),ncol = 10)+(1-corr[r])*diag(10)
    for(m in 1:2){ #sample size
      for(u in 1:4){ #different hypothesis
        RR <- foreach(i=1:1000,.combine = c,.packages = "LaplacesDemon")%dopar%{
          #result = Gen_data3(n = n[m], beta = Beta[,u],mu = c(3,3,3), Sigma1 = S,Sigma2 = 3*diag(3))
          #result = Gen_data5(n = n[m], beta = Beta[,u], Sigma0 = S)
          result = Gen_data(n = n[m], beta = Beta[,u], Sigma0 = S)
          y = result$Y
          x = result$X
          DEDGP(y,x,alpha = alpha[p])$DGP
        }
        RS[m,u]=1-max(table(RR))/1000
      }
    }
    rownames(RS) = c("n=500","n=1000")
    colnames(RS) = c("K^*=0","K^*=1","K^*=2","K^*=3")  
    write.table(RS,sep=",",file =paste("C:/Users/Administrator/Desktop/Pleio Degree/DGP-10-",100*alpha[p],"/PDEDGP10_",10*corr[r],".csv",sep=""),row.names = T)
    beep(1)
    Sys.sleep(20)
  }
}
stopImplicitCluster()
stopCluster(cl)
beep(3)




















