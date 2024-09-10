####################################################
## Example figures in Figure 1,2,3,4 in the paper ##
####################################################

#Using intermediate data
beta=read.table(".\\beta.txt")
gamma=read.table(".\\gamma.txt")
cp=read.table(".\\cp.txt")
SE=read.table(".\\SE.txt")
m0SE=read.table(".\\m0SE.txt")
m1SE=read.table(".\\m1SE.txt")

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


beta=beta[1:500,]
gamma=gamma[1:500,]
cp=cp[1:500,]
SE=SE[1:500,]
m0SE=m0SE[1:500,]
m1SE=m1SE[1:500,]

beta0_true = c(sqrt(5/13),sqrt(4/13),sqrt(4/13))  
beta1_true = c(sqrt(1/3),sqrt(1/3),sqrt(1/3))
beta_true=c(beta0_true, beta1_true)

beta_est = apply(beta,2,mean)
CP = apply(cp,2,mean)

bias = beta_est - beta_true    #bias
SD = apply(beta, 2, sd)    #standard deviation

SE = apply(SE, 2, mean)

MSE = bias^2 + SE^2

m1SE = apply(m1SE,2,mean)
m0SE = apply(m0SE,2,mean)

m1SE=m1SE*200
m0SE=m0SE*200

w = seq(0.25,1.4,0.01)
d=4
K=1
Jn=d+K+1
m1hat_rep = spline(w,K,d)%*%t(gamma[,(Jn+1):(2*Jn)])
m1hat = apply(m1hat_rep,1,mean)
CBupper = m1hat + 1.96*m1SE
CBlower = m1hat - 1.96*m1SE


A1 = sqrt(3)/2-1.645/sqrt(12)
A2 = sqrt(3)/2+1.645/sqrt(12)


m1true = sin(pi*(w-A1)/(A2-A1))

plot(w, m1true, type="l", col="red", main=expression(paste('n=200,',~'m=10,',~p[A],'=0.5')),xlab=expression(u[1]), ylab=expression(m[1](u[1])), ylim=c(-1.5,1.5), cex.lab=1.4,cex.axis=1.2,cex.main=1.8)
lines(w, m1hat, lty=2, col="blue")
lines(w, CBupper, lty=6, col="blue")
lines(w, CBlower, lty=6, col="blue")



m0true = cos(pi*w)
m0hat_rep = spline(w,K,d)%*%t(gamma[,1:Jn])
m0hat = apply(m0hat_rep,1,mean)
CBupper0 = m0hat + 1.96*m0SE
CBlower0 = m0hat - 1.96*m0SE

plot(w, m0true, type="l", col="red", main=expression(paste('n=200,',~'m=10,',~p[A],'=0.5')),xlab=expression(u[0]), ylab=expression(m[0](u[0])), ylim=c(-1.5,2), cex.lab=1.4,cex.axis=1.2,cex.main=1.8)
lines(w, m0hat, lty=2, col="blue")
lines(w, CBupper0, lty=6, col="blue")
lines(w, CBlower0, lty=6, col="blue")









