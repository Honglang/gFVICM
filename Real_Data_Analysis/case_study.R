source(".\\functions_real_data.R")
load(".\\codon492.RData")
########################################################################################
### This is for codon 492. For other codons, need to change two places for this code:
### 1) the above data loading
### 2) the G=data$G492 in the following last line of the data prepration part
########################################################################################
#### prepare for analysis
dose = c(0,5,10,20,30,40)
#rescale of dosage
dose1 = (dose-min(dose))/(max(dose)-min(dose))

age = data$AGE
bmi = data$BMI
id = data$ID
uniq_id = unique(id)
N = length(uniq_id)
m = 6
Nobs = nrow(data)

## rescale of X 
x2 = (age - min(age))/(max(age)-min(age))
x3 = (bmi - min(bmi))/(max(bmi)-min(bmi))
dose2 = rep(dose1,N)
x1 = dose2
x = cbind(x1,x2,x3)
p = ncol(x)

mvec = rep(m,N)
mvec0 = cumsum(c(1, mvec[-N])) #start loc for each id
mvec1 = cumsum(mvec) #end loc for each id

y = data$BP
G = data$G492

##############################################################

### corr matrix: AR(1) #####
m0 = list()
m1 = list()
for (i in 1:N){
  ni = mvec[i]
  m0[[i]] = diag(ni)           
  m1[[i]] = matrix(rep(0,ni*ni),ni)
  for (k in 1:ni) {
    for (l in 1:ni) {
      if (abs(k-l)==1) m1[[i]][k,l] <-1
    }
  }
  
}
###################################

set.seed(12345)
w = sort(rnorm(100,0.47,0.213))

f1 = function(z){estimation(x, y, G, K, d, z)[[1]]}

####### choose order and number of knots ######################################

knot = c(1,2,3)
order = 3
comb = as.matrix(expand.grid(knot,order))

p_list = rep(0,nrow(comb))
BIC = rep(0,nrow(comb))
BIC0 = rep(0,nrow(comb))

for(k in 1:nrow(comb)){
  K = comb[k,1]
  d = comb[k,2]
  
  lambda_est = golden(f1, 0.05, 2, K, d)
  
  fit = estimation(x,y,G,K,d,lambda_est)
  
  beta = fit[[2]]
  gamma = fit[[3]]
  Q = fit[[7]]
  
  lamb = lambda_est
  D = diag(rep(c(rep(0,d+1),rep(1,K)),2))
  
  npar = length(gamma)
  p_list[k] = pchisq(Q+N*lamb*t(gamma)%*%D%*%gamma,npar,lower.tail=F)
  
  BIC[k] = Q+N*lamb*t(gamma)%*%D%*%gamma + npar*log(N) 
  BIC0[k] = Q + npar*log(N)
}


knots = comb[which.max(p_list),]
#knots = comb[which.min(BIC0),]
K = knots[1]
d = knots[2]


#################################Estimation##############################
lambda_est = golden(f1, 0.01, 2, K, d)

result = estimation(x, y, G, K, d, lambda_est)

beta.old = result[[2]]
gamma.old = result[[3]]
SE_beta = result[[4]]
m0SE = result[[5]]*(Nobs/6)
m1SE = result[[6]]*(Nobs/6)

u0=x%*%beta.old[1:p]
u1=x%*%beta.old[(p+1):(2*p)]
B0=spline(u0,K,d)
B1=spline(u1,K,d)
m0hat=B0%*%gamma.old[1:(K+d+1)]
m1hat=B1%*%gamma.old[(K+d+1+1):(2*(K+d+1))]
yhat=m0hat+m1hat*G

err=y-yhat

plot(u0,m0hat)
plot(u1,m1hat)

MSE = (1/Nobs)*sum(err^2)

### test significance of beta1
z.test = beta.old/SE_beta
p.beta11 = 2*pnorm(-abs(z.test[4]))
p.beta12 = 2*pnorm(-abs(z.test[5]))
p.beta13 = 2*pnorm(-abs(z.test[6]))

m1hat_new = spline(w,K,d)%*%gamma.old[(K+d+1+1):(2*(K+d+1))]

CB_up = m1hat_new + 1.96*m1SE
CB_low = m1hat_new - 1.96*m1SE
###################################################################
##############     Figure 7 for one codon          ################
###################################################################
plot(w[5:90],m1hat_new[5:90],type='l',col='red',xlim=c(0.1,0.85),ylim=c(-1,5.5), main=expression(codon492),xlab=expression(u[1]), ylab=expression(m[1](u[1])) )
lines(w[5:90],CB_up[5:90],lty=2,col='blue')
lines(w[5:90],CB_low[5:90],lty=2,col='blue')
###################################################################

##########################Testing#########################################
result1 = estimation1(x,y,G,K,d,K,d,lambda_est)
beta1 = result1[[2]]
gamma1 = result1[[3]]
Q1 = Qstat(x,y,G,K,d,K,d,beta1,gamma1)

result0 = estimation1(x,y,G,K,d,0,1,lambda_est)
beta0 = result0[[2]]
gamma0 = c(result0[[3]], rep(0,K+d-1))
Q0 = Qstat(x,y,G,K,d,K,d,beta0,gamma0)

test = Q0 - Q1
df_test = K+d-1
cutoff = qchisq(0.95, df_test)

rej = (test > cutoff)*1
p.value = pchisq(test,df_test,lower.tail = F)


testresults=list(MAF=min(2*table(G)[1]+table(G)[2],2*table(G)[3]+table(G)[2])/(2*sum(table(G))),
                       p.beta11=p.beta11,p.beta12=p.beta12,p.beta13=p.beta13,pm=1-pchisq(test,df_test),Q=Q1)

odds_calcution=function(dosage,G){
  u0=c(dosage,mean(x[,2]),mean(x[,3]))%*%beta.old[1:p]
  u1=c(dosage,mean(x[,2]),mean(x[,3]))%*%beta.old[(p+1):(2*p)]
  B0=spline(u0,K,d)
  B1=spline(u1,K,d)
  m0hat=B0%*%gamma.old[1:(K+d+1)]
  m1hat=B1%*%gamma.old[(K+d+1+1):(2*(K+d+1))]
  odds=exp(m0hat+m1hat*G)
}

odds_results=matrix(0,nrow=3,ncol=6)
dosage_levels=c(0,5,10,20,30,40)/40
g_levels=c(0,1,2)
for(do in 1:6){
  for(g in 1:3){
    odds_results[g,do]=odds_calcution(dosage_levels[do],g_levels[g])
  }
}

###################################################################
#########     Table 3 and Table 4 for one codon          ##########
###################################################################
testresults
odds_results
###################################################################

