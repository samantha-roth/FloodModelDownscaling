rm(list=ls())

nCh<- 200
nLoc<- 10
nPars=2

library(mvtnorm)
#library(foreach); library(doParallel)

if(dir.exists("/storage/work/svr5482/Reification/Philly/data/Norristown/nCh/calibration")==F){
  dir.create("/storage/work/svr5482/Reification/Philly/data/Norristown/nCh/calibration")}

if(dir.exists("/storage/work/svr5482/Reification/Philly/data/Norristown/nCh/calibration/10m")==F){
  dir.create("/storage/work/svr5482/Reification/Philly/data/Norristown/nCh/calibration/10m")}

if(dir.exists("/storage/work/svr5482/Reification/Philly/data/Norristown/nCh/calibration/10m/ST10")==F){
  dir.create("/storage/work/svr5482/Reification/Philly/data/Norristown/nCh/calibration/10m/ST10")}


#take another random sample
set.seed(56)
RS10inds<- sample(1:30,nLoc)

################################################################################
#load the HWMs
load("/storage/work/svr5482/Reification/Philly/data/Norristown/nCh/obsWE_RS.RData")
Z<- obsWE_RS<- obsWE_RS[RS10inds]
s_ob<- .03
ig_b=.5
ig_a= (ig_b/(s_ob^2))-1

#load parameter settings
load("/storage/work/svr5482/FloodingModelCalibrationProject/parameterSamples/Philly/nCh_brkunif.RData")

#load the projections closest to the HWMs
load("/storage/work/svr5482/Reification/Philly/data/Norristown/nCh/closest10m_RS.RData")
closest10m_RS<- closest10m_RS[RS10inds,]

#load emulator parameters from the emulator fit with R
load("/storage/work/svr5482/Reification/Philly/data/Norristown/nCh/ST/parsFN1e-5ExpGPbyLoc.10m.RData")
kappa.vec10<- kappa.vec10[RS10inds]
phi.vec10<- phi.vec10[RS10inds]
zeta.vec10<- rep(1e-5,length(RS10inds))

#load parameters
load("/storage/work/svr5482/FloodingModelCalibrationProject/parameterSamples/Philly/nCh_brkunif.RData")

################################################################################

x.std<- (samp10$ch-.01)/(.1-.01)

y.c<- apply(closest10m_RS,1,function(y) y- mean(y))
mean.y<- rowMeans(closest10m_RS)

HWM.c<- obsWE_RS-mean.y

y.c.Diffs<- apply(y.c,1,function(y) HWM.c-y)
y.c.meanAbsDiffs<- colMeans(abs(y.c.Diffs))
ascMAB<- order(y.c.meanAbsDiffs)
best3<- ascMAB[1:3]
best3y.c.Diffs<- y.c.Diffs[,best3]
delta<- rowMeans(best3y.c.Diffs)
#delta<- rep(0,nLoc)

################################################################################

Sig2inv <- array(dim = c(nCh, nCh, nLoc))

sigma_sq= kappa.vec10
omega= 1/(phi.vec10^2)
nugget=  zeta.vec10

for(loc in 1:nLoc){
  
  Sig2= matrix(NA, nrow= nCh, ncol= nCh)
  for(i in 1:nCh){
    for(j in 1:nCh){
      Sig2[i,j]<- sigma_sq[loc]*exp(-omega[loc]*(abs(x.std[i] - x.std[j])))
    }
  }
  
  Sig2= Sig2 + diag(nugget[loc],nrow=nCh)
  Sig2inv[,,loc]<- solve(Sig2)
}

################################################################################

predMeanCov<- function(prop_ch,loc){
  #SS=c(prop,x.std)
  #AA = matrix(SS,length(SS),length(SS))
  #BB = t(AA)
  #C1 = abs(AA-BB)
  
  zeta=zeta.vec10[loc]
  phi=phi.vec10[loc] #retrieve the estimated parameter from the matrix
  kappa=kappa.vec10[loc]
  
  Sigma11<- zeta+kappa
  Sigma12<- matrix(kappa*exp(-(abs(x.std-prop_ch)/phi^2)),1,nCh)
  #cov.red=zeta*diag(1,1+nCh)+kappa*exp(-(C1/phi^2)) #compute the covariance function
  #Sigma22=cov.red[-1,-1] #Sigma_22 matrix
  #Sigma12=matrix(cov.red[1,-1],1,nCh) #Sigma matrix
  #Sigma11=cov.red[1,1]
  #Sigma22inv=solve(Sigma22)
  pred10=Sigma12%*%Sig2inv[,,loc]%*%y.c[,loc] #predicted values from the emulator
  var10=diag(Sigma11-Sigma12%*%Sig2inv[,,loc]%*%t(Sigma12))
  list(pred=pred10,var=var10)
}

loglik<- function(prop_ch,loc,prop_s2){
  -.5*(log(predMeanCov(prop_ch,loc)$var+prop_s2)+
         (HWM.c[loc]-predMeanCov(prop_ch,loc)$pred-delta[loc])^2/(predMeanCov(prop_ch,loc)$var+prop_s2))
}

sum_ll_func<- function(prop_ch,prop_s2){
  sum_ll<- 0
  for(loc in 1:nLoc){
    sum_ll<- sum_ll + loglik(prop_ch,loc,prop_s2)
  }
  return(as.numeric(sum_ll))
}
################################################################################
logDetJ<- function(p,a=0,b=1){
  -sum(log(p-a)+log(b-p))
}

UtoR<- function(p,a=0,b=1){
  log((p-a)/(b-a))-log(1-((p-a)/(b-a)))
}

RtoU<- function(x,a=0,b=1){
  (b-a)*(1/(1+exp(-x)))+a
}


################################################################################
#when fixing observation variance
#z_s2= s_ob^2
################################################################################
#Code the Gibbs sampler.
mh.alg<- function(init, n.sample) {
  x.t <- init[1] #ch,rwe
  s2.t <- init[2] #sigma^2
  
  x.out <- rep(NA,n.sample+1)
  s2.out <- rep(NA,n.sample+1)
  
  x.out[1]= x.t
  s2.out[1]= s2.t
  
  x.test<- x.out[1]
  s2.test<- s2.out[1]
  #predMeanCov<- pred_mean_cov(x.test)
  
  #predMean<- predMeanCov[[1]]
  #predCov<- predMeanCov[[2]]
  sum_ll_test<- sum_ll_func(x.test,s2.test)
  
  proposals<- matrix(NA, nrow= n.sample, ncol= 2*nPars)
  
  for (i in 1 : n.sample) {
    #if(i%%1000==0){print(paste("Iteration",i))}
    
    ############################################################################
    
    #n_ch
    u<- runif(1)
    
    propDistMean<- UtoR(as.numeric(x.out[i]),a=0,b=1)
    
    propX<- rnorm(1, mean = propDistMean, sd = .1) #was .5
    
    prop<- RtoU(propX,a=0,b=1)
    
    #record proposals
    proposals[i,1]<- prop
    
    logMetRatio1<- sum_ll_func(prop,s2.out[i])
    logMetRatio2<- sum_ll_test
    
    
    logMetRatio<- logMetRatio1-logMetRatio2+
      logDetJ(as.numeric(x.out[i]),a= 0,b= 1)-
      logDetJ(prop,a= 0,b= 1)
    
    #metropolis ratio in the metropolis hastings algorithm 
    if(log(u)<=logMetRatio){
      x.out[i+1]<- prop
      
      x.test<- x.out[i+1]
      sum_ll_test<- sum_ll_func(x.test,s2.out[i])
      #record acceptance
      proposals[i,2]<- 1
    }
    if(log(u)>logMetRatio){
      x.out[i+1]<- x.out[i]
      #record rejection
      proposals[i,2]<- 0
    } 
    
    
    ############################################################################
    #MH algorithm for sigma^2
    u<- runif(1)
    propX<- rnorm(1,log(s2.out[i]),sd= 1) #was 3
    prop<- exp(propX); 
    
    ##record proposal
    proposals[i,3]<- prop
    
    
    logMetRatio1<- sum_ll_func(x.out[i],prop) - 
      (ig_a+1)*log(prop) - (ig_b/prop)
    logMetRatio2<- sum_ll_func(x.out[i],s2.out[i]) - 
      (ig_a+1)*log(s2.out[i]) - (ig_b/s2.out[i])
    
    logMetRatio<- logMetRatio1-logMetRatio2+propX-log(s2.out[i])
    
    
    if(log(u)<=logMetRatio){
      s2.out[i+1]<- prop
      s2.test<- prop
      #record acceptance
      proposals[i,4]=1
    } 
    
    if(log(u)>logMetRatio){
      s2.out[i+1]<- s2.out[i]
      
      #record rejection
      proposals[i,4]=0
    } 
    
  }
  out <- cbind(x.out, s2.out)
  #out<- x.out
  list(out=out,proposals=proposals)
}



#set.seed(73)
#init_vals<- c(runif(1),rinvgamma(1,shape=ig_a,rate=ig_b))
init_vals<- c(0.4423369169,0.0009579518)
#init_vals<- (.03-.01)/(.1-.01)

#load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m10m/outputData/nCh_RWE/100e400c/just10m/cal/mcmc_res_homGP10.RData")

#init_vals<- c(res[nrow(res),]) #use last step of first Markov chain
#rm(res)

set.seed(51) #for last 200k steps
pt<-proc.time() # Start Time
#output <- mh.alg(init = init_vals, n.sample = 200000)
#output <- mh.alg(init = init_vals, n.sample = 100000)
output <- mh.alg(init = init_vals, n.sample = 10000)
#output <- mh.alg(init = init_vals, n.sample = 1000)

proposals<- output$proposals
res<- output$out

ptFinal<-proc.time()-pt ; ptFinal<-ptFinal[3] # End Time to be used in Effective Samples per Second Calculation

plot(1:nrow(res),res[,1]*(.1-.01)+.01)

plot(1:nrow(res),res[,2])
#res<- res[-1,]

#setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m10m/outputData/nCh_RWE/100e400c/just10m/cal/.25/3cm")
setwd("/storage/work/svr5482/Reification/Philly/data/Norristown/nCh/calibration/10m/ST10") 
save(res,file="chain_MH_ch_2ds2.RData")
save(proposals,file="props_MH_ch_2ds2.RData")
save(ptFinal,file="time_MH_ch_2ds2.RData")

