
DGP = function(y,x,K,alpha=0.05){
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  V = diag(q)
  DG = combn(1:q,K)
  DG1 = combn(1:q,K-1)
  x = x - mean(x)
  sX = solve(t(x)%*%(x))
  beta = sX%*%(t(x)%*%y)
  Sigma = (1/n)*t(y)%*%(diag(n)-x%*%sX%*%t(x))%*%y
  Omegasq = kronecker(Msq(Sigma,-0.5),diag(n))
  X = kronecker(diag(q),x)
  Y = Vec(y,byrow = F)
  #beta1 = solve(t(X)%*%Omegasq%*%(X))%*%(t(X)%*%Omegasq%*%(Y))
  #standalized
  X1 = Omegasq%*%X
  Y1 = Omegasq%*%Y
  V = diag(q)
  beta2 = solve(t(X1)%*%(X1))%*%(t(X1)%*%(Y1))
  HX = t(X1)%*%(X1)
  #Trecord = t(beta2)%*%HX%*%beta2
  Trecord = c()
  r = ncol(DG)
  for(k in 1:r){
    sf = DG[,k]
    Trecord = c(Trecord,t(beta2)%*%t(V[-sf,])%*%solve(V[-sf,]%*%solve(HX)%*%t(V[-sf,]))%*%(V[-sf,])%*%beta2)
  }
  Tmin = min(Trecord)
  location = DG[,which.min(Trecord)] 
  
  Trecord1 = c()
  r1 = ncol(DG1)
  for(k1 in 1:r1){
    sf1 = DG1[,k1]
    Trecord1 = c(Trecord1,t(beta2)%*%t(V[-sf1,])%*%solve(V[-sf1,]%*%solve(HX)%*%t(V[-sf1,]))%*%(V[-sf1,])%*%beta2)
  }
  Tmin1 = min(Trecord1)
  
  #index = 0
  #if(Trecord[1]>qchisq(df = q,1-alpha)){
  #if(Tmin>qchisq(df = q-1,1-alpha)){
  #index = 1
  #}
  #}
  index = 1 - (Tmin<qchisq(df = q-K,1-alpha)&Tmin1>qchisq(df = q-K+1,1-alpha))
  
  
  return(list(index = index, location = location,beta = beta2,test = Trecord,quan = qchisq(df = q-K,1-alpha),test1 = Trecord1,quan1 = qchisq(df = q-K+1,1-alpha)))
}

DEDGP = function(y,x,alpha=0.05){ #detect K^*
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  V = diag(q)
  x = x - mean(x)
  sX = solve(t(x)%*%(x))
  beta = sX%*%(t(x)%*%y)
  Sigma = (1/n)*t(y)%*%(diag(n)-x%*%sX%*%t(x))%*%y
  Omegasq = kronecker(Msq(Sigma,-0.5),diag(n))
  X = kronecker(diag(q),x)
  Y = Vec(y,byrow = F)
  #beta1 = solve(t(X)%*%Omegasq%*%(X))%*%(t(X)%*%Omegasq%*%(Y))
  #standalized
  X1 = Omegasq%*%X
  Y1 = Omegasq%*%Y
  
  beta2 = solve(t(X1)%*%(X1))%*%(t(X1)%*%(Y1))
  HX = t(X1)%*%(X1)
  
  
  T0 = t(beta2)%*%HX%*%beta2
  if(T0 < qchisq(df = q, 1-alpha)){
    return(list(DGP=0,location = "No"))
  }
  
 
  for(k in 1:(q-2)){  #TEST FOR DGP=k
    DG = combn(1:q,k)
    v = ncol(DG)
    Trecord = rep(0,v)
    for(m in 1:v){
      sf = DG[,m]
      V1 = as.matrix(V[-sf,])
      Trecord[m] = t(beta2)%*%t(V1)%*%solve(V1%*%solve(HX)%*%t(V1))%*%(V1)%*%beta2
    }
    Tmin = min(Trecord)
    loc = which.min(Trecord)
    if(Tmin < qchisq(df = q - k, 1-alpha)){
      return(list(DGP=k,location = DG[,loc]))
    }
  }
  
  DG = combn(1:q,q-1)
  v = ncol(DG)
  Trecord = rep(0,v)
  for(m in 1:v){
    sf = DG[,m]
    V1 = t(as.matrix(V[-sf,]))
    Trecord[m] = t(beta2)%*%t(V1)%*%solve(V1%*%solve(HX)%*%t(V1))%*%(V1)%*%beta2
  }
  Tmin = min(Trecord)
  loc = which.min(Trecord)
  if(Tmin < qchisq(df = 1, 1-alpha)){
    return(list(DGP=q-1,location = DG[,loc]))
  }
  
  
  
  
  
  
  #index = 0
  #if(Trecord[1]>qchisq(df = q,1-alpha)){
  #if(Tmin>qchisq(df = q-1,1-alpha)){
  #index = 1
  #}
  #}
   
  return(list(DGP = q, location = 1:q))
}



DEDGP1 = function(y,x){ #detect K^*
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  V = diag(q)
  x = x - mean(x)
  sX = solve(t(x)%*%(x))
  beta = sX%*%(t(x)%*%y)
  Sigma = (1/n)*t(y)%*%(diag(n)-x%*%sX%*%t(x))%*%y
  Omegasq = kronecker(Msq(Sigma,-0.5),diag(n))
  X = kronecker(diag(q),x)
  Y = Vec(y,byrow = F)
  #beta1 = solve(t(X)%*%Omegasq%*%(X))%*%(t(X)%*%Omegasq%*%(Y))
  #standalized
  X1 = Omegasq%*%X
  Y1 = Omegasq%*%Y
  
  beta2 = solve(t(X1)%*%(X1))%*%(t(X1)%*%(Y1))
  HX = t(X1)%*%(X1)
  
  Tre = rep(0, q)
  pvalue = rep(0,q)
  T1 = t(beta2)%*%HX%*%beta2
  Tre[1] = T1
  pvalue = 1 - pchisq(T1,q)
  
  
  for(k in 1:(q-2)){  #TEST FOR DGP=k
    DG = combn(1:q,k)
    v = ncol(DG)
    Trecord = rep(0,v)
    for(m in 1:v){
      sf = DG[,m]
      V1 = as.matrix(V[-sf,])
      Trecord[m] = t(beta2)%*%t(V1)%*%solve(V1%*%solve(HX)%*%t(V1))%*%(V1)%*%beta2
    }
    Tmin = min(Trecord)
    Tre[k+1] = Tmin
    pvalue[k+1] = 1 - pchisq(Tmin,q-k)
    #loc = which.min(Trecord)
  }
 
  DG = combn(1:q,q-1)
  v = ncol(DG)
  Trecord = rep(0,v)
  for(m in 1:v){
    sf = DG[,m]
    V1 = t(as.matrix(V[-sf,]))
    Trecord[m] = t(beta2)%*%t(V1)%*%solve(V1%*%solve(HX)%*%t(V1))%*%(V1)%*%beta2
  }
  Tmin = min(Trecord)
  Tre[q] = Tmin
  pvalue[q] = 1 - pchisq(Tmin, 1)
  #loc = which.min(Trecord)
  
  
  
  
  
  
  
  #index = 0
  #if(Trecord[1]>qchisq(df = q,1-alpha)){
  #if(Tmin>qchisq(df = q-1,1-alpha)){
  #index = 1
  #}
  #}
  
  return(list(Tre = Tre,pvalue=pvalue,beta = beta2))
}


