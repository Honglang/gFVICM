#############################################################################
##### estimation
library(MASS)
library(mvtnorm)
library(Matrix)

########## p-spline functions ################################

spline = function(xx, K, d){
  Jn = d+K+1
  
  S = matrix(0, nrow=length(xx), ncol=Jn)	 
  if (K == 0){ 
    for(i in 1:(d+1)){
      S[,i]=xx^(i-1)
    }
  }  else {
    for(i in 1:(d+1)){
      S[,i]=xx^(i-1)
    }
    
    knots=quantile(xx, probs = (1:K)/(K+1))
    
    for(j in 1:K){
      
      S[,d+1+j]=((xx-knots[j])^d)*(xx>=knots[j])
    } 
  }
  S
}

### first derivative of splines

spline_dev = function(xx, K, d){
  
  Jn = d+K+1
  
  S_dev = matrix(0, nrow=length(xx), ncol=Jn)
  if (K==0){
    for(i in 2:(d+1)){
      S_dev[,i]=(i-1)*(xx^(i-2))
    }
  }  else{ for(i in 2:(d+1)){
    S_dev[,i]=(i-1)*(xx^(i-2))
  }
    
    knots=quantile(xx, probs = (1:K)/(K+1))
    
    for(j in 1:K){
      
      S_dev[,d+1+j]=d*((xx-knots[j])^(d-1))*(xx>=knots[j]) 
    }
  } 
  S_dev
}

##############################################################
#### estimation function ##########
estimation = function(x, y, G, K, d, lambda){
  
  K0=K1=K
  d0=d1=d
  
  D = diag(c(rep(0,d0+1),rep(1,K0),rep(0,d1+1),rep(1,K1))) # K: num of knots, d: order   
  p = ncol(x)
  
  ##step 0: initial values
  beta.old = c(sqrt(1/3),sqrt(1/3),sqrt(1/3), sqrt(5/13),sqrt(4/13),sqrt(4/13))
  beta.old0 = beta.old[1:p]
  beta.old1 = beta.old[(p+1):(2*p)]
  u0 = x%*%beta.old0
  u1 = x%*%beta.old1 
  u = cbind(u0, u1)  
  
  B0 = spline(u0, K0, d0)
  B1 = spline(u1, K1, d1)
  
  G_ext = matrix(rep(G,K1+d1+1), nrow=Nobs, ncol=K1+d1+1)
  B = cbind(B0, B1*G_ext)
  
  my.data = data.frame(y,B)  
  #res = geeglm(y~-1+., family=binomial(link="logit"), data=my.data, id=id, corstr = "ar1", std.err="san.se")   
  res = glm(y~-1+., family=binomial, data=my.data)
  gamma.old = as.vector(res$coefficients)
  
  maxstep = 0   ##loop for gamma and beta
  while(maxstep <= 50){
    maxstep = maxstep+1
    
    ##################
    #step 1: spline approximation
    beta.old0 = beta.old[1:p]
    beta.old1 = beta.old[(p+1):(2*p)]
    u0 = x%*%beta.old0
    u1 = x%*%beta.old1 
    u = cbind(u0, u1)  
    
    B0 = spline(u0, K0, d0)
    B1 = spline(u1, K1, d1)
    
    B = cbind(B0, B1*G_ext)
    
    ###########################
    ## step 2: gamma est
    ### Newton-Raphson for estimating gamma
    
    run = 0  
    while(run <= 50){
      #cat("run==============================================",run)
      run = run+1
      
      mu = B %*% gamma.old    
      mudot_gamma = B
      ngamma = (K0+d0+1+K1+d1+1)
      arsumg = matrix(rep(0,2*ngamma),nrow=2*ngamma)
      arsumc = matrix(rep(0,2*ngamma*2*ngamma),nrow=2*ngamma)
      #gi = matrix(rep(0,2*ngamma),nrow=2*ngamma)
      arsumgfirstdev1 = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
      #firstdev = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
      
      for(i in 1:N){
        seq = c(mvec0[i]:mvec1[i])
        ni = mvec[i]
        yi = y[seq]
        xi = x[seq,]
        mui = mu[seq]
        mudoti = mudot_gamma[seq,]
        
        ui = 1 /(1 + exp(- mui))
        #fui = log(ui) - log(1-ui)
        fui=mui
        fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
        #vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1+exp(mui))))#diag(as.vector(sqrt(1/(1-ui))))
        vui = diag(as.vector(exp(mui/2)+exp(-mui/2))) 
        
        # cat("yi",yi)
        # cat("xi",xi)
        # cat("ui",ui)
        # cat("mui",mui)
        # cat("fui",fui)
        # cat("vui",vui)
        # 
        wi = t(mudoti) %*% fui_dev %*% vui %*% m0[[i]] %*% vui
        zi = t(mudoti) %*% fui_dev %*% vui %*% m1[[i]] %*% vui
        
        gi0 = (1/N)* wi %*% (yi-ui)
        gi1 = (1/N)* zi %*% (yi-ui)
        
        gi = c(gi0, gi1)
        # cat("gi")
        # print(gi)
        arsumc = arsumc + gi %*% t(gi) 
        arsumg = arsumg + gi 
        
        di0 = - (1/N)* wi %*% fui_dev %*% mudoti #first dev of gi0
        di1 = - (1/N)* zi %*% fui_dev %*% mudoti #first dev of gi1
        
        firstdev = rbind(di0,di1)
        arsumgfirstdev1 = arsumgfirstdev1 + firstdev
        
      }
      
      # cat("arsumc")
      # print(arsumc)
      # svdresult=svd(arsumc)
      # arcinv10=svdresult$v%*%diag(1/svdresult$d)%*%t(svdresult$u)
      # cat("arcinv10")
      # print(arcinv10)
      arcinv1 = ginv(arsumc)
      # cat("arcinv1")
      # print(arcinv1)
      
      Q1 = t(arsumg) %*% arcinv1 %*% arsumg 
      
      arqif1dev1 = (2*t(arsumgfirstdev1) %*% arcinv1 %*% arsumg)/N + 2*lambda*D%*%gamma.old
      arqif2dev1 = (2*t(arsumgfirstdev1) %*% arcinv1 %*% arsumgfirstdev1)/N + 2*lambda*D 	
      
      invarqif2dev1 = ginv(arqif2dev1)
      
      gamma.new = gamma.old - invarqif2dev1 %*% arqif1dev1
      gammadiff = max(abs(gamma.new - gamma.old))
      if(gammadiff<1e-6){break}
      
      gamma.old = gamma.new
      #print(gamma.old)
    } #loop for gamma
    
    
    #########################               
    # step 3: est for beta
    ### Newton-Raphson for estimating beta
    
    NRbeta.old=beta.old
    run2 = 0  
    while(run2 <= 50){
      #cat("run2==============================================",run2)
      run2 = run2+1
      
      NRbeta.old0 = NRbeta.old[1:p]
      NRbeta.old1 = NRbeta.old[(p+1):(2*p)]
      u0 = x%*%NRbeta.old0
      u1 = x%*%NRbeta.old1 
      u = cbind(u0, u1)  
      
      B0 = spline(u0, K0, d0)
      B1 = spline(u1, K1, d1)
      B = cbind(B0, B1*G_ext)
      mu = B %*% gamma.old
      
      B0d = spline_dev(u0, K0, d0)
      B1d = spline_dev(u1, K1, d1)
      
      gamma0.old = gamma.old[1:(K0+d0+1)]
      gamma1.old = gamma.old[(K0+d0+1+1):(K0+d0+1+K1+d1+1)]
      
      # first derivative of mu 
      mudot_beta = cbind(matrix(rep(B0d%*%gamma0.old,p),nrow=Nobs, ncol=p)*x, matrix(rep(B1d%*%gamma1.old*G,p),nrow=Nobs, ncol=p)*x)        
      
      nbeta = 2*p
      arsumg = matrix(rep(0,2*nbeta),nrow=2*nbeta)
      arsumc = matrix(rep(0,2*nbeta*2*nbeta),nrow=2*nbeta)
      #gi = matrix(rep(0,2*nbeta),nrow=2*nbeta)
      arsumgfirstdev = matrix(rep(0,2*nbeta*nbeta),nrow=2*nbeta)
      #firstdev = matrix(rep(0,2*nbeta*nbeta),nrow=2*nbeta)
      
      for(i in 1:N){
        seq = c(mvec0[i]:mvec1[i])
        ni = mvec[i]
        yi = y[seq]
        xi = x[seq,]
        mui = mu[seq]
        mudoti = mudot_beta[seq,]
        
        ui = 1 /(1 + exp(- mui))
        #fui = log(ui) - log(1-ui)
        fui=mui
        fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
        #vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1+exp(mui))))#diag(as.vector(sqrt(1/(1-ui))))
        vui = diag(as.vector(exp(mui/2)+exp(-mui/2))) 
        
        wi = t(mudoti) %*% fui_dev %*% vui %*% m0[[i]] %*% vui
        zi = t(mudoti) %*% fui_dev %*% vui %*% m1[[i]] %*% vui
        
        gi0 = (1/N)* wi %*% (yi-ui)
        gi1 = (1/N)* zi %*% (yi-ui)
        # if(i==13){
        #   cat("yi",yi)
        #   cat("xi",xi)
        #   cat("ui",ui)
        #   cat("mui",mui)
        #   cat("fui",fui)
        #   cat("vui",diag(vui))
        # }
        gi = c(gi0, gi1)
        # cat("gi")
        # print(gi)
        arsumc = arsumc + gi %*% t(gi) 
        arsumg = arsumg + gi 
        
        di0 = - (1/N)* wi %*% fui_dev %*% mudoti #first dev of gi0
        di1 = - (1/N)* zi %*% fui_dev %*% mudoti #first dev of gi1
        
        firstdev = rbind(di0,di1)
        arsumgfirstdev = arsumgfirstdev + firstdev
        
      }
      # cat("arsumc")
      # print(arsumc)
      # svdresult=svd(arsumc)
      # arcinv10=svdresult$v%*%diag(1/svdresult$d)%*%t(svdresult$u)
      # cat("arcinv10")
      # print(arcinv10)
      arcinv=ginv(arsumc)
      # cat("arcinv")
      # print(arcinv)
      
      
      
      Q2 = t(arsumg) %*% arcinv %*% arsumg
      
      arqif1dev2 = 2*t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev2 = 2*t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev   	
      
      invarqif2dev2 = ginv(arqif2dev2)
      
      NRbeta.new = NRbeta.old - invarqif2dev2 %*% arqif1dev2
      
      NRbeta.new1 = sign(NRbeta.new[1])*NRbeta.new[1:p]/sqrt(sum(NRbeta.new[1:p]^2))
      NRbeta.new2 = sign(NRbeta.new[1+p])*NRbeta.new[(p+1):(2*p)]/sqrt(sum(NRbeta.new[(p+1):(2*p)]^2))
      NRbeta.new_norm = c(NRbeta.new1, NRbeta.new2)
      
      betadiff = max(abs(NRbeta.new_norm - NRbeta.old))
      if(betadiff<1e-6){break}
      
      NRbeta.old = NRbeta.new_norm
      #print(NRbeta.old)
    } #loop for beta
    
    beta.new = NRbeta.old
    
    dif = max(abs(beta.new-beta.old))
    #dif = sqrt(sum((beta.new-beta.old)^2))
    if(dif<1e-6){break}
    #print(dif)
    #print(gamma.new)
    #print(beta.new)
    beta.old = beta.new
    
  }#outer iteration for beta and gamma
  
  
  
  # based on the estimation, calculate Q, covariance
  QNdotdot = 2*t(arsumgfirstdev1) %*% arcinv1 %*% arsumgfirstdev1
  QN = Q1
  df = sum(diag(ginv(QNdotdot + 2*N*lambda*D)%*%QNdotdot))
  GCV = (QN/N)/(1-df/N)^2
  
  lam_gam = lambda*D%*%gamma.old
  arsums = matrix(rep(0,ngamma*ngamma),nrow=ngamma)
  
  for (i in 1:N){
    seq = c(mvec0[i]:mvec1[i])
    ni = mvec[i]
    yi = y[seq]
    xi = x[seq,]
    mui = mu[seq]
    mudoti = mudot_gamma[seq,]
    
    # ui = 1 /(1 + exp(- mui))
    # fui = log(ui) - log(1-ui)
    # fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
    # vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1-ui))))
    ui = 1 /(1 + exp(- mui))
    #fui = log(ui) - log(1-ui)
    fui=mui
    fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
    #vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1+exp(mui))))#diag(as.vector(sqrt(1/(1-ui))))
    vui = diag(as.vector(exp(mui/2)+exp(-mui/2))) 
    
    wi = t(mudoti) %*% fui_dev %*% vui %*% m0[[i]] %*% vui
    zi = t(mudoti) %*% fui_dev %*% vui %*% m1[[i]] %*% vui
    
    gi0 = (1/N)* wi %*% (yi-ui)
    gi1 = (1/N)* zi %*% (yi-ui)
    
    gi = c(gi0, gi1)
    si = t(arsumgfirstdev1) %*% (arcinv1/N) %*% gi + lam_gam
    
    arsums = arsums + si%*%t(si)
    
  }
  
  cov_part = ginv(arqif2dev1*N/2)
  cov_gamma = cov_part%*%arsums%*%cov_part
  cov_gamma1 = cov_gamma[(K0+d0+1+1):(K0+d0+1+K1+d1+1),(K0+d0+1+1):(K0+d0+1+K1+d1+1)]
  cov_gamma0 = cov_gamma[1:(K0+d0+1),1:(K0+d0+1)]
  
  
  sp_matrix0 = spline(w,K0,d0)
  sp_matrix1 = spline(w,K1,d1)
  
  m1SE = sqrt(diag(sp_matrix1%*%cov_gamma1%*%t(sp_matrix1)))
  m0SE = sqrt(diag(sp_matrix0%*%cov_gamma0%*%t(sp_matrix0)))
  
  asy.cov = ginv(arqif2dev2/2)
  SE_beta = sqrt(diag(asy.cov))
  
  
  
  list(GCV, beta.old, gamma.old, SE_beta, m0SE, m1SE, QN)
  
  
} #estimation function

### golden search method

golden = function(f, lower, upper, K, d) { 
  
  golden.ratio = 2/(sqrt(5) + 1)
  x1 = upper - golden.ratio*(upper - lower)
  x2 = lower + golden.ratio*(upper - lower)
  
  f1 = f(x1)
  f2 = f(x2)
  
  iteration = 0
  
  while (abs(upper - lower) > 0.1)
  {
    iteration = iteration + 1
    
    if (f2 > f1)
    {
      ### Set the new upper bound
      upper = x2
      ### Set the new upper test point
      ### Use the special result of the golden ratio
      x2 = x1
      f2 = f1
      
      ### Set the new lower test point
      x1 = upper - golden.ratio*(upper - lower)
      f1 = f(x1)
    } 
    else 
    {
      # the minimum is to the right of x1
      # let x1 be the new lower bound
      # let x2 be the new lower test point
      
      ### Set the new lower bound
      lower = x1
      
      ### Set the new lower test point
      x1 = x2
      
      f1 = f2
      
      ### Set the new upper test point
      x2 = lower + golden.ratio*(upper - lower)
      f2 = f(x2)
    }
  }
  
  ### Use the mid-point of the final interval as the estimate of the optimzer
  estimated.minimizer = (lower + upper)/2
  estimated.minimizer
}


########################################################################
########################################################################
### testing the linearity of m1



#### estimation function ##########
estimation1 = function(x, y, G, K0, d0, K1, d1, lambda){
  
  D = diag(c(rep(0,d0+1),rep(1,K0),rep(0,d1+1),rep(1,K1))) # K: num of knots, d: order
  #Jn = d+K+1
  p = ncol(x)
  
  ##step 0: initial values
  beta.old = c(sqrt(1/3),sqrt(1/3),sqrt(1/3), sqrt(5/13),sqrt(4/13),sqrt(4/13))
  beta.old0 = beta.old[1:p]
  beta.old1 = beta.old[(p+1):(2*p)]
  u0 = x%*%beta.old0
  u1 = x%*%beta.old1 
  u = cbind(u0, u1)  
  
  B0 = spline(u0, K0, d0)
  B1 = spline(u1, K1, d1)
  
  G_ext = matrix(rep(G,K1+d1+1), nrow=Nobs, ncol=K1+d1+1)
  B = cbind(B0, B1*G_ext)
  
  my.data = data.frame(y,B)  
  #res = geeglm(y~-1+., family=binomial(link="logit"), data=my.data, id=id, corstr = "ar1", std.err="san.se")   
  res = glm(y~-1+., family=binomial, data=my.data)
  gamma.old = as.vector(res$coefficients)
  
  maxstep = 0   ##loop for gamma and beta
  while(maxstep <= 50){
    maxstep = maxstep+1
    
    ##################
    #step 1: b-spline approximation
    beta.old0 = beta.old[1:p]
    beta.old1 = beta.old[(p+1):(2*p)]
    u0 = x%*%beta.old0
    u1 = x%*%beta.old1 
    u = cbind(u0, u1)  
    
    B0 = spline(u0, K0, d0)
    B1 = spline(u1, K1, d1)
    
    B = cbind(B0, B1*G_ext)
    
    ###########################
    ## step 2: gamma est
    ### Newton-Raphson for estimating gamma
    
    run = 0  
    while(run <= 50){
      
      run = run+1
      
      mu = B %*% gamma.old    
      mudot_gamma = B
      ngamma = K0+d0+1+K1+d1+1
      arsumg = matrix(rep(0,2*ngamma),nrow=2*ngamma)
      arsumc = matrix(rep(0,2*ngamma*2*ngamma),nrow=2*ngamma)
      #gi = matrix(rep(0,2*ngamma),nrow=2*ngamma)
      arsumgfirstdev1 = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
      #firstdev = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
      
      for(i in 1:N){
        seq = c(mvec0[i]:mvec1[i])
        ni = mvec[i]
        yi = y[seq]
        xi = x[seq,]
        mui = mu[seq]
        mudoti = mudot_gamma[seq,]
        
        # ui = 1 /(1 + exp(- mui))
        # fui = log(ui) - log(1-ui)
        # fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
        # vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1-ui))))
        
        ui = 1 /(1 + exp(- mui))
        #fui = log(ui) - log(1-ui)
        fui=mui
        fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
        #vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1+exp(mui))))#diag(as.vector(sqrt(1/(1-ui))))
        vui = diag(as.vector(exp(mui/2)+exp(-mui/2))) 
        
        wi = t(mudoti) %*% fui_dev %*% vui %*% m0[[i]] %*% vui
        zi = t(mudoti) %*% fui_dev %*% vui %*% m1[[i]] %*% vui
        
        gi0 = (1/N)* wi %*% (yi-ui)
        gi1 = (1/N)* zi %*% (yi-ui)
        
        gi = c(gi0, gi1)
        
        arsumc = arsumc + gi %*% t(gi) 
        arsumg = arsumg + gi 
        
        di0 = - (1/N)* wi %*% fui_dev %*% mudoti #first dev of gi0
        di1 = - (1/N)* zi %*% fui_dev %*% mudoti #first dev of gi1
        
        firstdev = rbind(di0,di1)
        arsumgfirstdev1 = arsumgfirstdev1 + firstdev
        
      }
      
      
      arcinv1 = ginv(arsumc)
      
      Q1 = t(arsumg) %*% arcinv1 %*% arsumg 
      
      arqif1dev1 = (2*t(arsumgfirstdev1) %*% arcinv1 %*% arsumg)/N + 2*lambda*D%*%gamma.old
      arqif2dev1 = (2*t(arsumgfirstdev1) %*% arcinv1 %*% arsumgfirstdev1)/N + 2*lambda*D 	
      
      invarqif2dev1 = ginv(arqif2dev1)
      
      gamma.new = gamma.old - invarqif2dev1 %*% arqif1dev1
      gammadiff = max(abs(gamma.new - gamma.old))
      if(gammadiff<1e-6){break}
      
      gamma.old = gamma.new
      #print(gamma.old)
    } #loop for gamma
    
    
    #########################               
    # step 3: est for beta
    ### Newton-Raphson for estimating beta
    
    NRbeta.old=beta.old
    run2 = 0  
    while(run2 <= 50){
      
      run2 = run2+1
      
      NRbeta.old0 = NRbeta.old[1:p]
      NRbeta.old1 = NRbeta.old[(p+1):(2*p)]
      u0 = x%*%NRbeta.old0
      u1 = x%*%NRbeta.old1 
      u = cbind(u0, u1)  
      
      B0 = spline(u0, K0, d0)
      B1 = spline(u1, K1, d1)
      B = cbind(B0, B1*G_ext)
      mu = B %*% gamma.old
      
      B0d = spline_dev(u0, K0, d0)
      B1d = spline_dev(u1, K1, d1)
      
      gamma0.old = gamma.old[1:(K0+d0+1)]
      gamma1.old = gamma.old[(K0+d0+1+1):(K0+d0+1+K1+d1+1)]
      
      # first derivative of mu 
      mudot_beta = cbind(matrix(rep(B0d%*%gamma0.old,p),nrow=Nobs, ncol=p)*x, matrix(rep(B1d%*%gamma1.old*G,p),nrow=Nobs, ncol=p)*x)        
      
      nbeta = 2*p
      arsumg = matrix(rep(0,2*nbeta),nrow=2*nbeta)
      arsumc = matrix(rep(0,2*nbeta*2*nbeta),nrow=2*nbeta)
      #gi = matrix(rep(0,2*nbeta),nrow=2*nbeta)
      arsumgfirstdev = matrix(rep(0,2*nbeta*nbeta),nrow=2*nbeta)
      #firstdev = matrix(rep(0,2*nbeta*nbeta),nrow=2*nbeta)
      
      for(i in 1:N){
        seq = c(mvec0[i]:mvec1[i])
        ni = mvec[i]
        yi = y[seq]
        xi = x[seq,]
        mui = mu[seq]
        mudoti = mudot_beta[seq,]
        
        # ui = 1 /(1 + exp(- mui))
        # fui = log(ui) - log(1-ui)
        # fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
        # vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1-ui))))
        
        ui = 1 /(1 + exp(- mui))
        #fui = log(ui) - log(1-ui)
        fui=mui
        fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
        #vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1+exp(mui))))#diag(as.vector(sqrt(1/(1-ui))))
        vui = diag(as.vector(exp(mui/2)+exp(-mui/2))) 
        
        
        wi = t(mudoti) %*% fui_dev %*% vui %*% m0[[i]] %*% vui
        zi = t(mudoti) %*% fui_dev %*% vui %*% m1[[i]] %*% vui
        
        gi0 = (1/N)* wi %*% (yi-ui)
        gi1 = (1/N)* zi %*% (yi-ui)
        
        gi = c(gi0, gi1)
        
        arsumc = arsumc + gi %*% t(gi) 
        arsumg = arsumg + gi 
        
        di0 = - (1/N)* wi %*% fui_dev %*% mudoti #first dev of gi0
        di1 = - (1/N)* zi %*% fui_dev %*% mudoti #first dev of gi1
        
        firstdev = rbind(di0,di1)
        arsumgfirstdev = arsumgfirstdev + firstdev
        
      }
      
      
      
      arcinv=ginv(arsumc)
      
      Q2 = t(arsumg) %*% arcinv %*% arsumg
      
      arqif1dev2 = 2*t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev2 = 2*t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev   	
      
      invarqif2dev2 = ginv(arqif2dev2)
      
      NRbeta.new = NRbeta.old - invarqif2dev2 %*% arqif1dev2
      
      NRbeta.new1 = sign(NRbeta.new[1])*NRbeta.new[1:p]/sqrt(sum(NRbeta.new[1:p]^2))
      NRbeta.new2 = sign(NRbeta.new[1+p])*NRbeta.new[(p+1):(2*p)]/sqrt(sum(NRbeta.new[(p+1):(2*p)]^2))
      NRbeta.new_norm = c(NRbeta.new1, NRbeta.new2)
      
      betadiff = max(abs(NRbeta.new_norm - NRbeta.old))
      if(betadiff<1e-6){break}
      
      NRbeta.old = NRbeta.new_norm
      #print(NRbeta.old)
    } #loop for beta
    
    beta.new = NRbeta.old
    
    dif = max(abs(beta.new-beta.old))
    #dif = sqrt(sum((beta.new-beta.old)^2))
    if(dif<1e-6){break}
    #print(dif)
    #print(gamma.new)
    #print(beta.new)
    beta.old = beta.new
    
  }#outer iteration for beta and gamma
  
  
  # based on the estimation, calculate Q, covariance
  QNdotdot = 2*t(arsumgfirstdev1) %*% arcinv1 %*% arsumgfirstdev1
  QN = Q1
  df = sum(diag(ginv(QNdotdot + 2*N*lambda*D)%*%QNdotdot))
  GCV = (QN/N)/(1-df/N)^2
  
  
  list(GCV, beta.old, gamma.old, Q1)
  
  
} #estimation function

########################################################################  


########################################################################

Qstat = function(x, y, G, K0, d0, K1, d1, beta, gamma){      
  # based on the estimation, calculate Q
  beta.old0 = beta[1:p]
  beta.old1 = beta[(p+1):(2*p)]
  u0 = x%*%beta.old0
  u1 = x%*%beta.old1 
  u = cbind(u0, u1)  
  
  B0 = spline(u0, K0, d0)
  B1 = spline(u1, K1, d1)
  G_ext = matrix(rep(G,K1+d1+1), nrow=Nobs, ncol=K1+d1+1)
  B = cbind(B0, B1*G_ext)
  mu = B %*% gamma
  
  B0d = spline_dev(u0, K0, d0)
  B1d = spline_dev(u1, K1, d1)
  
  gamma0 = gamma[1:(K0+d0+1)]
  gamma1 = gamma[(K0+d0+1+1):length(gamma)]
  
  J01 = -beta.old0[-1]/beta.old0[1]
  J0 = rbind(J01,diag(p-1))
  
  J11 = -beta.old1[-1]/beta.old1[1]
  J1 = rbind(J11,diag(p-1))
  
  
  # first derivative of mu wrt beta   
  mudot0 = (matrix(rep(B0d%*%gamma0,p),nrow=Nobs, ncol=p)*x)%*%J0
  mudot1 = (matrix(rep(B1d%*%gamma1*G,p),nrow=Nobs, ncol=p)*x)%*%J1     
  
  # first derivative of mu 
  mudot_gamma = B
  mudot_beta = cbind(mudot0, mudot1)        
  
  mudot = cbind(mudot_beta, mudot_gamma)
  
  npar = length(beta)+length(gamma)-2  #number of parameters
  arsumg = matrix(rep(0,2*npar),nrow=2*npar)
  arsumc = matrix(rep(0,2*npar*2*npar),nrow=2*npar)
  #gi = matrix(rep(0,2*npar),nrow=2*npar)
  arsumgfirstdev = matrix(rep(0,2*npar*npar),nrow=2*npar)
  #firstdev = matrix(rep(0,2*npar*npar),nrow=2*npar)
  
  for(i in 1:N){
    seq = c(mvec0[i]:mvec1[i])
    ni = mvec[i]
    yi = y[seq]
    xi = x[seq,]
    mui = mu[seq]
    mudoti = mudot[seq,]
    
    # ui = 1 /(1 + exp(- mui))
    # fui = log(ui) - log(1-ui)
    # fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
    # vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1-ui))))
    
    ui = 1 /(1 + exp(- mui))
    #fui = log(ui) - log(1-ui)
    fui=mui
    fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
    #vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1+exp(mui))))#diag(as.vector(sqrt(1/(1-ui))))
    vui = diag(as.vector(exp(mui/2)+exp(-mui/2))) 
    
    wi = t(mudoti) %*% fui_dev %*% vui %*% m0[[i]] %*% vui
    zi = t(mudoti) %*% fui_dev %*% vui %*% m1[[i]] %*% vui
    
    gi0 = (1/N)* wi %*% (yi-ui)
    gi1 = (1/N)* zi %*% (yi-ui)
    
    gi = c(gi0, gi1)
    
    arsumc = arsumc + gi %*% t(gi) 
    arsumg = arsumg + gi 
    
    di0 = - (1/N)* wi %*% fui_dev %*% mudoti #first dev of gi0
    di1 = - (1/N)* zi %*% fui_dev %*% mudoti #first dev of gi1
    
    firstdev = rbind(di0,di1)
    arsumgfirstdev = arsumgfirstdev + firstdev
    
  }
  
  
  arcinv=ginv(arsumc)
  
  Q = t(arsumg) %*% arcinv %*% arsumg
  Q
  
} #Qstat function

###################################################################
