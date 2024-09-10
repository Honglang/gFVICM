source(".\\functions.R")


gFVICM_simulation=function(N,m,pA,tau,S){
  #N is the sample size
  #m is the number of repeated measurements for each individual
  #pA is the minor alle frequency for the genetic data 
  #tau is the parameter for the alternative hypothesis in the testing of function forms
  #S is the number of simulation repetition 
  beta0_true = c(sqrt(5/13),sqrt(4/13),sqrt(4/13))  
  beta1_true = c(sqrt(1/3),sqrt(1/3),sqrt(1/3))
  beta_true=c(beta0_true, beta1_true)
  
  beta_est=matrix(0,ncol=6,nrow=S)
  se_beta=matrix(0,ncol=6,nrow=S)
  coverage=matrix(0,ncol=6,nrow=S)
  test_res=rep(0,S)
  for(i in 1:S){
    print(i)
    result=estimation_main(N,m,pA,i)
    beta_est[i,]=result$beta.old
    se_beta[i,]=result$SE_beta
    coverage[i,]=result$cp
    test_res[i]=tesing_main(N,m,pA,tau,i)
  }
  

  beta_bias=apply(beta_est,2,mean)-beta_true
  beta_sd=apply(beta_est,2,sd)
  beta_se=apply(se_beta,2,mean)
  beta_coverage_prob=apply(coverage,2,mean)
  #############################################
  power=mean(test_res)
  
  
  return(list(beta_bias=beta_bias,beta_sd=beta_sd,beta_se=beta_se,
              beta_coverage_prob=beta_coverage_prob,power=power))
}


#############################################
### Table 1 and 2                        ####
## Table 1 from the following  res[,,,,1] ###
## Table 1 from the following  res[,,,,2] ###
#############################################
S=500
Ns=c(200,500)
ms=c(10,20)
pAs=c(0.1,0.3,0.5)
res=array(0,c(6,4,3,2,2)) # first two dim for 6 betas and their c(bias,sd,se,cp), the last three dims for pA,N,m
for(i in 1:2){#for m
  for(j in 1:2){# for N
    for(l in 1:3){# for pA
      result_temp=gFVICM_simulation(Ns[j],ms[i],pAs[l],0,S=S)
      res[,,l,j,i]=matrix(c(result_temp$beta_bias,result_temp$beta_sd,result_temp$beta_se,
                            result_temp$beta_coverage_prob),ncol=4)
    }
  }
}

##############################################
### Figure 5                              ####
##############################################
tau_seq=(0:5)/10
power_res=array(0,c(length(tau_seq),length(Ns),length(ms)))
for(i in 1:length(tau_seq)){
  for(j in 1:length(Ns)){
    for(l in 1:length(ms)){
      res=gFVICM_simulation(Ns[j],ms[l],pA=0.3,tau_seq[i],S)
      power_res[i,j,l]=res$power
    }
  }
}

plot(tau_seq,power_res[,1,1],xlab=expression(tau),ylab="Power",type='l',col='red')
lines(tau_seq,power_res[,2,1],col='blue')

plot(tau_seq,power_res[,1,2],xlab=expression(tau),ylab="Power",type='l',col='red')
lines(tau_seq,power_res[,2,2],col='blue')

#############################################
##############################################
### Figure 6                              ####
##############################################
tau_seq=(0:5)/10
power_res=matrix(0,c(length(tau_seq),length(pAs)))
for(i in 1:length(tau_seq)){
  for(j in 1:length(pAs)){
    res=gFVICM_simulation(500,10,pAs[j],tau_seq[i],S)
    power_res[i,j]=res$power
  }
}

plot(tau_seq,power_res[,1],xlab=expression(tau),ylab="Power",type='l',col='red')
lines(tau_seq,power_res[,2],col='blue')
lines(tau_seq,power_res[,3],col='green')


