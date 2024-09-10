library("MASS")
library("mvtnorm")
library("Matrix")
library("bindata")

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

###################################################################
### golden search method

#f0 = function(x){estimation_0(X, Y_0, K, d, x)[[1]]}

golden = function(f, lower, upper) { 
  
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

#########################################################################################

##############################################################
#################### estimation function ###############

#### estimation function ##########
estimation = function(x, y, G, N,mvec0,mvec1,mvec,K, d, lambda){
  beta0_true = c(sqrt(5/13),sqrt(4/13),sqrt(4/13))  
  beta1_true = c(sqrt(1/3),sqrt(1/3),sqrt(1/3))
  beta_true=c(beta0_true, beta1_true)
  
  m0 = list()
  m1 = list()
  for (i in 1:N){
    ni = mvec[i]
    ## corr matrix: AR(1) 
    m0[[i]] = diag(ni)           
    m1[[i]] = matrix(rep(0,ni*ni),ni)
    for (k in 1:ni) {
      for (l in 1:ni) {
        if (abs(k-l)==1) m1[[i]][k,l] <-1
      }
    }
    
  }
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
  Nobs=sum(mvec)
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
        fui = log(ui) - log(1-ui)
        fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
        vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1-ui))))
        
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
        
        ui = 1 /(1 + exp(- mui))
        fui = log(ui) - log(1-ui)
        fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
        vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1-ui))))
        
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
  
  lam_gam = lambda*D%*%gamma.old
  arsums = matrix(rep(0,ngamma*ngamma),nrow=ngamma)
  
  for (i in 1:N){
    seq = c(mvec0[i]:mvec1[i])
    ni = mvec[i]
    yi = y[seq]
    xi = x[seq,]
    mui = mu[seq]
    mudoti = mudot_gamma[seq,]
    
    ui = 1 /(1 + exp(- mui))
    fui = log(ui) - log(1-ui)
    fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
    vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1-ui))))
    
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
  
  
  w = seq(0.25,1.4,0.01)
  sp_matrix0 = spline(w,K0,d0)
  sp_matrix1 = spline(w,K1,d1)
  
  m1SE = sqrt(diag(sp_matrix1%*%cov_gamma1%*%t(sp_matrix1)))
  m0SE = sqrt(diag(sp_matrix0%*%cov_gamma0%*%t(sp_matrix0)))
  
  asy.cov = ginv(arqif2dev2/2)
  SE_beta = sqrt(diag(asy.cov))
  
  CIlower = beta.old - 1.96*SE_beta   #confidence interval
  CIupper = beta.old + 1.96*SE_beta
  
  cp = (beta_true<=CIupper) & (beta_true>=CIlower)
  
  
  list(GCV, beta.old, gamma.old, SE_beta, m0SE, m1SE, cp, QN)
  
  
} #estimation function

########################################################################

Qstat = function(x, y, G, N,mvec0,mvec1,mvec,K0, d0, K1, d1, beta, gamma){      
  # based on the estimation, calculate Q
  m0 = list()
  m1 = list()
  for (i in 1:N){
    ni = mvec[i]
    ## corr matrix: AR(1) 
    m0[[i]] = diag(ni)           
    m1[[i]] = matrix(rep(0,ni*ni),ni)
    for (k in 1:ni) {
      for (l in 1:ni) {
        if (abs(k-l)==1) m1[[i]][k,l] <-1
      }
    }
    
  }
  p = ncol(x)
  
  beta.old0 = beta[1:p]
  beta.old1 = beta[(p+1):(2*p)]
  u0 = x%*%beta.old0
  u1 = x%*%beta.old1 
  u = cbind(u0, u1)  
  
  B0 = spline(u0, K0, d0)
  B1 = spline(u1, K1, d1)
  Nobs = sum(mvec)
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
    
    ui = 1 /(1 + exp(- mui))
    fui = log(ui) - log(1-ui)
    fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
    vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1-ui))))
    
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

##############################################################
#### estimation function ##########
estimation1 = function(x, y, G,N, mvec0,mvec1,mvec,K0, d0, K1, d1, lambda){
  m0 = list()
  m1 = list()
  for (i in 1:N){
    ni = mvec[i]
    ## corr matrix: AR(1) 
    m0[[i]] = diag(ni)           
    m1[[i]] = matrix(rep(0,ni*ni),ni)
    for (k in 1:ni) {
      for (l in 1:ni) {
        if (abs(k-l)==1) m1[[i]][k,l] <-1
      }
    }
    
  }
  
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
  Nobs = sum(mvec)
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
        
        ui = 1 /(1 + exp(- mui))
        fui = log(ui) - log(1-ui)
        fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
        vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1-ui))))
        
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
        
        ui = 1 /(1 + exp(- mui))
        fui = log(ui) - log(1-ui)
        fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
        vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1-ui))))
        
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

#### estimation0 function for m0 only ##########
estimation0 = function(x, y, G,N,mvec0,mvec1,mvec, K0, d0, lambda){
  m0 = list()
  m1 = list()
  for (i in 1:N){
    ni = mvec[i]
    ## corr matrix: AR(1) 
    m0[[i]] = diag(ni)           
    m1[[i]] = matrix(rep(0,ni*ni),ni)
    for (k in 1:ni) {
      for (l in 1:ni) {
        if (abs(k-l)==1) m1[[i]][k,l] <-1
      }
    }
    
  }
  
  D = diag(c(rep(0,d0+1),rep(1,K0))) # K: num of knots, d: order
  #Jn = d+K+1
  p = ncol(x)
  
  ##step 0: initial values
  beta.old = c(sqrt(1/3),sqrt(1/3),sqrt(1/3))
  u0 = x%*%beta.old     
  
  B0 = spline(u0, K0, d0)
  
  B = B0
  
  gamma.old = ginv(t(B)%*%B)%*%t(B)%*%y  ##simple LS, initial value 
  
  maxstep = 0   ##loop for gamma and beta
  while(maxstep <= 50){
    maxstep = maxstep+1
    
    ##################
    #step 1: b-spline approximation
    
    u0 = x%*%beta.old
    
    B0 = spline(u0, K0, d0)
    
    B = B0
    
    ###########################
    ## step 2: gamma est
    ### Newton-Raphson for estimating gamma
    
    run = 0  
    while(run <= 50){
      
      run = run+1
      
      mu = B %*% gamma.old    
      mudot_gamma = B
      ngamma = K0+d0+1
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
        fui = log(ui) - log(1-ui)
        fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
        vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1-ui))))
        
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
      
      u0 = x%*%NRbeta.old
      
      B0 = spline(u0, K0, d0)
      B = B0
      mu = B %*% gamma.old
      
      B0d = spline_dev(u0, K0, d0)   
      
      # first derivative of mu 
      mudot_beta = matrix(rep(B0d%*%gamma.old,p),nrow=Nobs, ncol=p)*x       
      
      nbeta = p
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
        fui = log(ui) - log(1-ui)
        fui_dev = diag(as.vector(ui)) %*% diag(as.vector(1-ui))
        vui = diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1-ui))))
        
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
      
      NRbeta.new_norm = sign(NRbeta.new[1])*NRbeta.new/sqrt(sum(NRbeta.new^2))
      
      betadiff = max(abs(NRbeta.new_norm - NRbeta.old))
      if(betadiff<1e-6){break}
      
      NRbeta.old = NRbeta.new_norm
      #print(NRbeta.old)
    } #loop for beta
    
    beta.new = NRbeta.old
    
    dif = max(abs(beta.new-beta.old))
    #dif = sqrt(sum((beta.new-beta.old)^2))
    if(dif<1e-4){break}
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
  
  
} #estimation0 function, for m0 only

estimation_main=function(N = 200,m = 10,pr = 0.5,seednumber){
  set.seed(seednumber)
  ####true values
  beta0_true = c(sqrt(5/13),sqrt(4/13),sqrt(4/13))  
  beta1_true = c(sqrt(1/3),sqrt(1/3),sqrt(1/3))
  beta_true=c(beta0_true, beta1_true)
  
  #genarate data
  mvec = rep(0,N)
  indx = vector("list", N)
  for (i in 1:N){
    visit = c(1,rbinom((m-1),1,1)) # 30% missing except 1st time
    mvec[i] = sum(visit)  
    indx[[i]] = which(visit==1)
  }
  Nobs = sum(mvec)
  mvec0 = cumsum(c(1, mvec[-N])) #start loc for each id
  mvec1 = cumsum(mvec) #end loc for each id
  
  N0 = N*2
  Nobs0 = N0*m
  x1 = runif(Nobs0)
  x2 = runif(Nobs0)
  x3 = runif(Nobs0)
  x = cbind(x1,x2,x3)
  p = ncol(x)
  
  
  G0 = sample(c(0,1,2), size = N0, replace = T, prob=c((1-pr)^2, 2*pr*(1-pr), pr^2))
  G = rep(G0,rep(m,N0)) 
  
  
  ### AR(1)
  rho = 0.2
  R = matrix(0,m,m)
  for (k in 1:m){
    for (l in 1:m){
      R[k,l]=rho^(abs(k-l))
    }
  }
  
  A1 = sqrt(3)/2-1.645/sqrt(12)
  A2 = sqrt(3)/2+1.645/sqrt(12)
  
  u0 = x%*%beta0_true
  u1 = x%*%beta1_true
  
  mean_f = cos(pi*u0) + sin(pi*(u1-A1)/(A2-A1))*G 
  
  P_1 = exp(mean_f)/(1+exp(mean_f))
  
  y_temp = vector("list", N0)
  for(i in 1:N0){
    seq = c((1+(i-1)*m):(i*m))
    meani = mean_f[seq]
    P_1i = P_1[seq]
    #y.old = try(ep(mu=P_1i, R=R, nRep=1, seed=NULL))
    y.old = try(rmvbin(n=1, margprob=P_1i, bincorr=R)) 
    if (!inherits(y.old,"try-error")){
      y_temp[[i]] = y.old
    } else (y_temp[[i]] = rep(NA,m))
  }
  
  y.new = unlist(y_temp)
  my_indx = which(y.new!='NA')
  y_in = y.new[my_indx]     
  x_in = x[my_indx,]
  G_in = G[my_indx]
  
  #new data
  y = y_in[1:Nobs]
  x = x_in[1:Nobs,]
  G = G_in[1:Nobs]
  id = rep(1:N, mvec)
  
  ### corr matrix: AR(1) #####
  m0 = list()
  m1 = list()
  for (i in 1:N){
    ni = mvec[i]
    ## corr matrix: AR(1) 
    m0[[i]] = diag(ni)           
    m1[[i]] = matrix(rep(0,ni*ni),ni)
    for (k in 1:ni) {
      for (l in 1:ni) {
        if (abs(k-l)==1) m1[[i]][k,l] <-1
      }
    }
    
  }
  ###################################
  ############ determine K, d ######################
  
  d_list = c(3,4)
  K_list = c(1,2,3)
  all_list = expand.grid(d_list,K_list)
  
  #p.value1 = rep(0,nrow(all_list))
  #BIC1 = rep(0,nrow(all_list))
  BIC1_1 = rep(0,nrow(all_list))
  
  for (kk in 1:nrow(all_list)){
    
    d = all_list[kk,1]    
    K = all_list[kk,2]
    est = function(lam){estimation(x, y, G,N,mvec0,mvec1,mvec, K, d, lam)[[1]]}
    lamb_est = golden(est, 0.01, 2)
    fit1 = estimation(x, y, G, N,mvec0,mvec1,mvec, K, d, lamb_est)
    
    beta = fit1[[2]]
    gamma = fit1[[3]]
    Q = Qstat(x, y, G, N,mvec0,mvec1,mvec, K, d, K, d, beta, gamma)
    
    npar = length(gamma) 
    
    BIC1_1[kk] = Q+npar*log(N)
    
  }
  
  min.row = which.min(BIC1_1)
  
  d = all_list[min.row,1]
  K = all_list[min.row,2]
  
  
  
  #######################################################
  ########## after determination of knots and order
  
  
  f_1 = function(z){estimation(x, y, G,N,mvec0,mvec1,mvec,  K, d, z)[[1]]}
  
  lambda_est = golden(f_1, 0.01, 2)
  
  result2 = estimation(x, y, G, N,mvec0,mvec1,mvec, K, d, lambda_est)
  return(list(beta.old = result2[[2]],
              gamma.old = result2[[3]],
              SE_beta = result2[[4]],
              m0SE = result2[[5]],
              m1SE = result2[[6]],
              cp = result2[[7]]))
}



tesing_main=function(N = 200,m = 10,pA = 0.5,abc = 0.2,seednumber){
  set.seed(seednumber)
  N0 = N*2
  Nobs0 = N0*m
  
  mvec = rep(m,N)
  #mvec = sample(c(10:20),size=N,replace=T)
  Nobs = sum(mvec)
  mvec0 = cumsum(c(1, mvec[-N])) #start location for each id
  mvec1 = cumsum(mvec) #end location for each id
  
  N0 = N*2
  Nobs0 = N0*m
  x1 = runif(Nobs0)
  x2 = runif(Nobs0)
  x3 = runif(Nobs0)
  
  x = cbind(x1,x2,x3)
  p = ncol(x)
  
  
  G0 = sample(c(0,1,2), size = N0, replace = T, prob=c((1-pA)^2, 2*pA*(1-pA), pA^2))
  G = rep(G0,rep(m,N0)) 
  
  
  ### AR(1)
  rho = 0.2
  R = matrix(0,m,m)
  for (k in 1:m){
    for (l in 1:m){
      R[k,l]=rho^(abs(k-l))
    }
  }
  
  
  ####true values
  beta0_true = c(sqrt(5/13),sqrt(4/13),sqrt(4/13))  
  beta1_true = c(sqrt(1/3),sqrt(1/3),sqrt(1/3))
  beta_true=c(beta0_true, beta1_true)
  A1 = sqrt(3)/2-1.645/sqrt(12)
  A2 = sqrt(3)/2+1.645/sqrt(12)
  
  u0 = x%*%beta0_true
  u1 = x%*%beta1_true
  
  
  m1_full = sin(pi*(u1-A1)/(A2-A1))
  m1_H0 = (0.3+u1)
  mean_f = cos(pi*u0) + (m1_H0 + abc*(m1_full-m1_H0))*G 
  
  P_1 = exp(mean_f)/(1+exp(mean_f))
  
  y_temp = vector("list", N0)
  for(i in 1:N0){
    seq = c((1+(i-1)*m):(i*m))
    meani = mean_f[seq]
    P_1i = P_1[seq]
    #y.old = try(ep(mu=P_1i, R=R, nRep=1, seed=NULL))
    y.old = try(rmvbin(n=1, margprob=P_1i, bincorr=R)) 
    if (!inherits(y.old,"try-error")){
      y_temp[[i]] = y.old
    } else (y_temp[[i]] = rep(NA,m))
  }
  
  y.new = unlist(y_temp)
  my_indx = which(y.new!='NA')
  y_in = y.new[my_indx]     
  x_in = x[my_indx,]
  G_in = G[my_indx]
  
  #new data
  y = y_in[1:Nobs]
  x = x_in[1:Nobs,]
  G = G_in[1:Nobs]
  id = rep(1:N, mvec)
  
  ### corr matrix: AR(1) #######################################
  m0 = list()
  m1 = list()
  for (i in 1:N){
    ni = mvec[i]
    m0[[i]] = diag(ni)
    m1[[i]] = matrix(rep(0,ni*ni),ni)
    for (k in 1:ni) {
      for (l in 1:ni) {
        if (abs(k-l)==1) m1[[i]][k,l] =1
      }
    }
  }
  ##############################################################
  
  ##### determine K0, d0 under y=m0+eps #################
  
  d0=4
  K0=1
  #K0_list = c(1,2,3)
  #p.value0 = rep(0,length(K0_list))
  #BIC0 = rep(0,length(K0_list))
  #BIC0_1 = rep(0,length(K0_list))
  
  #for (jj in 1:length(K0_list)){
  
  #    K0 = K0_list[jj]
  #    est_H0 = function(lam){estimation0(x, y_m0, G, K0, d0, lam)[[1]]}
  #    lamb_est = golden(est_H0, 0.01, 2)
  #    fit0 = estimation0(x, y_m0, G, K0, d0, lamb_est)
  
  #    beta = fit0[[2]]
  #    gamma = fit0[[3]]
  #    Q = fit0[[4]]
  
  #    npar = length(gamma) 
  #    #p.value0[jj] = pchisq(Q,npar,lower.tail=F)
  
  #    #BIC0[jj] = Q+N*lamb_est*t(gamma)%*%DD%*%gamma + npar*log(N) 
  #    BIC0_1[jj] = Q+npar*log(N) 
  
  #}
  
  
  #K0 = K0_list[which.min(BIC0_1)]
  
  
  ############ determine K1, d1 ######################
  
  d1_list = c(2,3)
  K1_list = c(1,2,3)
  all_list = expand.grid(d1_list,K1_list)
  
  #p.value1 = rep(0,nrow(all_list))
  #BIC1 = rep(0,nrow(all_list))
  BIC1_1 = rep(0,nrow(all_list))
  
  for (kk in 1:nrow(all_list)){
    
    d1 = all_list[kk,1]    
    K1 = all_list[kk,2]
    est_H1 = function(lam){estimation1(x, y, G, N,mvec0,mvec1,mvec, K0, d0, K1, d1, lam)[[1]]}
    lamb_est1 = golden(est_H1, 0.01, 2)
    fit1 = estimation1(x, y, G, N,mvec0,mvec1,mvec, K0, d0, K1, d1, lamb_est1)
    
    beta = fit1[[2]]
    gamma = fit1[[3]]
    Q = Qstat(x, y, G, N,mvec0,mvec1,mvec, K0, d0, K1, d1, beta, gamma)
    
    npar = length(gamma) 
    
    BIC1_1[kk] = Q+npar*log(N)
    
  }
  
  min.row = which.min(BIC1_1)
  
  d1 = all_list[min.row,1]
  K1 = all_list[min.row,2]
  
  
  
  ########################################################################################
  
  
  result1 = estimation1(x,y,G,N,mvec0,mvec1,mvec, K0,d0,K1,d1,0)
  beta1 = result1[[2]]
  gamma1 = result1[[3]]
  Q1 = Qstat(x,y,G,N,mvec0,mvec1,mvec, K0,d0,K1,d1,beta1,gamma1)
  
  result0 = estimation1(x,y,G,N,mvec0,mvec1,mvec, K0,d0,0,1,0)
  beta0 = result0[[2]]
  gamma0 = c(result0[[3]], rep(0,K1+d1-1))
  Q0 = Qstat(x,y,G,N,mvec0,mvec1,mvec, K0,d0,K1,d1,beta0,gamma0)
  
  test = Q0 - Q1
  df_test = K1+d1-1
  cutoff = qchisq(0.95, df_test)
  
  return(rej = (test > cutoff)*1)
}


