library(parallel)
library(boot)


pval<-function(x){
  a=mean(x)
  if(a<0){p=length(which(x<0))
  }else{p=length(which(x>0))}
  return(p/length(x))
}

###Function to estimate Rhat (Gelman and Rubin convergence diagnostic)
Rhatfun<-function(data1,data2,data3,m){
  var2=apply(data2,2,var)
  var3=apply(data3,2,var)
  long=dim(data1)[2]
  var1=apply(data1,2,var)
  mea1=apply(data1,2,mean)
  mea2=apply(data2,2,mean)
  mea3=apply(data3,2,mean)
  r=matrix(c(t(data3),t(data2),t(data1)), byrow=T,nrow=m*3)
  mea4=apply(r,2,mean)
  
  W=apply(matrix(c(var1,var2,var3), byrow=T, nrow=3),2,mean)
  W1=apply(matrix(c(var1,var2,var3), byrow=T, nrow=3),2,var)
  B=apply(matrix(c(mea1,mea2,mea3), byrow=T, nrow=3),2,var)
  cov1=cov2=rep(0,long)
  for (i in 1:(long)){
    cov1[i]=cov(matrix(c(var1,var2,var3), byrow=T, nrow=3)[,i],y=matrix(c(mea1,mea2,mea3)^2, byrow=T, nrow=3)[,i])
    cov2[i]=cov(matrix(c(var1,var2,var3), byrow=T, nrow=3)[,i],y=matrix(c(mea1,mea2,mea3), byrow=T, nrow=3)[,i])
  }
  sig2=((m-1)/m)*W+B
  V=sqrt(sig2+B/3)^2
  varV=((m-1)/m)^2/3*W1+(4/3)^2*B^2+2*(m-1)*4/(9*m)*(cov1-2*mea4*cov2)
  df=2*V^2/varV
  Rhat=abs((V/W*df)/(df-2))
  return(Rhat)
}


runim_ageC<-function(val,parameters,consts,data,inits,m){
  library(nimble)
  #   ###
  js_ms_age1 <- nimbleCode({
    # -----------------------------------
    # Parameters:
    # imi: immigration probability
    # gamma: removal entry probability 
    # S: survival probability
    # g: probability of starting depradation
    # p: depredating capture probability
    # q: from coast capture probability
    # -----------------------------------
    # States (S):
    # 1 not yet entered
    # 2 alive not depredating
    # 3 alive depredating
    # 4 dead
    # Observations (O):
    # 1 seen from the coast
    # 2 seen depredating
    # 3 not seen
    # -----------------------------------
    # Priors and constraints
    
    beta_effort~dnorm(0,0.001)
    mean.g~dnorm(0,0.001)
    mean.gamma~dunif(0,0.4)
    beta_s~dnorm(0,0.001)
    beta_state~dnorm(0,0.001)
    
    mean.S[2]~dnorm(0,0.001)
    mean.S[1]~dnorm(0,0.001)

    S[2]<-mean.S[1]
    S[3]<-mean.S[2]
    S[1]<--1000

    mean.p~dnorm(0,0.001)
    mean.q~dnorm(0,0.001)
    for (t in 1:(n.occasions-1)){
      g[t]<-ilogit(mean.g+beta_s*time[t])
      p[t]<-ilogit(mean.p+beta_effort*effort[t+1])
    }
    S1[1]<-ilogit(mean.S[1])
    S1[2]<-ilogit(mean.S[2])
    S2[1]<-ilogit(mean.S[1]+beta_state)
    S2[2]<-ilogit(mean.S[2]+beta_state)
    q<-ilogit(mean.q)
      # Define probabilities of state S(t+1) given S(t)
      for (t in 1:(n.occasions-1)){
        for (i in 1:3){ 
          ps[1,i,t,1] <- 1-mean.gamma
          ps[1,i,t,2] <- mean.gamma*(1-ilogit(mean.g+beta_s*time[t]))
          ps[1,i,t,3] <- ilogit(mean.g+beta_s*time[t])*mean.gamma
          ps[1,i,t,4] <- 0
          ps[2,i,t,1] <- 0
          ps[2,i,t,2] <- ilogit(S[i])*(1-ilogit(mean.g+beta_s*time[t]))
          ps[2,i,t,3] <- ilogit(mean.g+beta_s*time[t])*ilogit(S[i])
          ps[2,i,t,4] <- 1-ilogit(S[i])
          ps[3,i,t,1] <- 0
          ps[3,i,t,2] <- 0
          ps[3,i,t,3] <- ilogit(S[i]+beta_state)
          ps[3,i,t,4] <- 1-ilogit(S[i]+beta_state)
          ps[4,i,t,1] <- 0
          ps[4,i,t,2] <- 0
          ps[4,i,t,3] <- 0
          ps[4,i,t,4] <- 1
          }#j
        # Define probabilities of O(t) given S(t)
        po[1,t,1] <- 0
        po[1,t,2] <- 0
        po[1,t,3] <- 1
        po[2,t,1] <- ilogit(mean.q)
        po[2,t,2] <- 0
        po[2,t,3] <- 1-ilogit(mean.q)
        po[3,t,1] <- ilogit(mean.q)*(1-ilogit(mean.p+beta_effort*effort[t+1]))
        po[3,t,2] <- ilogit(mean.p+beta_effort*effort[t+1])
        po[3,t,3] <- 1-ilogit(mean.p+beta_effort*effort[t+1])- ilogit(mean.q)*(1-ilogit(mean.p+beta_effort*effort[t+1]))
        po[4,t,1] <- 0
        po[4,t,2] <- 0
        po[4,t,3] <- 1
      } #t
  
    
    
    # Likelihood
    for (i in 1:M){
      # Define latent state at first occasion
      z[i,1] <-1  # Make sure that all M individuals are in state 1 at t=1
      for (t in (2:n.occasions)){
        # State process: draw S(t) given S(t-1)
        z[i,t] ~ dcat(ps[z[i,t-1], stage[i,t-1], t-1,1:4])
        # Observation process: draw O(t) given S(t)
        y[i,t] ~ dcat(po[z[i,t], t-1,1:3])
      } #t
    } #i
    for (i in 1:M){
      for (t in 2:n.occasions){
        al.c[i,t-1] <- equals(z[i,t], 2) # alive non depre
        al.d[i,t-1] <- equals(z[i,t], 3) # alive depre
      } #t
      alive.c[i]<- sum(al.c[i,1:(n.occasions-1)])
      alive.d[i] <- sum(al.d[i,1:(n.occasions-1)])}
    
    for (t in 1:(n.occasions-1)){
      NC[t] <- sum(al.c[1:M,t]) # Actual population size non depre
      ND[t] <- sum(al.d[1:M,t]) # Actual population size depre
    } #t
    for (i in 1:M){
      w[i] <- 1-equals(alive.c[i]+alive.d[i],0)
    } #i
    
    Nsuper <- sum(w[1:M])
  })
  
  Sim_age    <- nimbleModel(code = js_ms_age1, name = 'Sim_age', constants = consts, data = data, inits = inits, check = F)
  CSim_age  <- compileNimble(Sim_age)
  specSim_age <- configureMCMC(Sim_age, monitors=parameters,thin=75, useConjugacy=FALSE)
  specSim_age$addSampler(target = c("mean.S",'beta_state'), type = 'AF_slice')
  specSim_age$addSampler(target = c("mean.g",'beta_s'), type = 'AF_slice')
  specSim_age$addSampler(target = c("mean.p","mean.q","mean.gamma"), type =  'AF_slice')
  specSim_age$addSampler(target = c("mean.gamma"), type =  'RW')
  Sim_ageMCMC <- buildMCMC(specSim_age)
  CSim_ageMCMC<- compileNimble(Sim_ageMCMC, project=Sim_age, resetFunctions = TRUE)
  CSim_ageMCMC$run(3*m*75)
  samplesSim_age=as.matrix(CSim_ageMCMC$mvSamples)[(m+1):(m*3),]
  return(samplesSim_age)}
  
  runim_ageD<-function(val,parameters,consts,data,inits,m){
    library(nimble)
    #   ###
    js_ms_age1 <- nimbleCode({
      # -----------------------------------
      # Parameters:
      # S: survival probability
      # g: entry & probability of starting depradation
      # p: depredating capture probability
      # q: from coast capture probability
      # -----------------------------------
      # States (S):
      # 1 not yet entered
      # 2 alive not depredating
      # 3 alive depredating
      # 4 dead
      # Observations (O):
      # 1 seen from the coast
      # 2 seen depredating
      # 3 not seen
      # -----------------------------------
      # Priors and constraints
      
      beta_effort~dnorm(0,0.001)
      mean.g~dnorm(0,0.001)
      beta_s~dnorm(0,0.001)
      mean.S[2]~dnorm(0,0.001)
      mean.S[1]~dnorm(0,0.001)
      
      S[2]<-mean.S[1]
      S[3]<-mean.S[2]
      S[1]<--1000
      
      mean.p~dnorm(0,0.001)
      for (t in 1:(n.occasions-1)){
        g[t]<-ilogit(mean.g+beta_s*time[t])
        p[t]<-ilogit(mean.p+beta_effort*effort[t+1])
      }
      S1[1]<-ilogit(mean.S[1])
      S1[2]<-ilogit(mean.S[2])

      q<-ilogit(mean.q)
      # Define probabilities of state S(t+1) given S(t)
      for (t in 1:(n.occasions-1)){
        for (i in 1:3){ 
          ps[1,i,t,1] <- 1-ilogit(mean.g+beta_s*time[t])
          ps[1,i,t,2] <- 0
          ps[1,i,t,3] <- ilogit(mean.g+beta_s*time[t])
          ps[1,i,t,4] <- 0
          ps[2,i,t,1] <- 0
          ps[2,i,t,2] <- 0
          ps[2,i,t,3] <- 0
          ps[2,i,t,4] <- 0
          ps[3,i,t,1] <- 0
          ps[3,i,t,2] <- 0
          ps[3,i,t,3] <- ilogit(S[i])
          ps[3,i,t,4] <- 1-ilogit(S[i])
          ps[4,i,t,1] <- 0
          ps[4,i,t,2] <- 0
          ps[4,i,t,3] <- 0
          ps[4,i,t,4] <- 1
        }#j
        # Define probabilities of O(t) given S(t)
        po[1,t,1] <- 0
        po[1,t,2] <- 0
        po[1,t,3] <- 1
        po[2,t,1] <- 0
        po[2,t,2] <- 0
        po[2,t,3] <- 1
        po[3,t,1] <- 0
        po[3,t,2] <- ilogit(mean.p+beta_effort*effort[t+1])
        po[3,t,3] <- 1-ilogit(mean.p+beta_effort*effort[t+1])
        po[4,t,1] <- 0
        po[4,t,2] <- 0
        po[4,t,3] <- 1
      } #t
      
      
      
      # Likelihood
      for (i in 1:M){
        # Define latent state at first occasion
        z[i,1] <-1  # Make sure that all M individuals are in state 1 at t=1
        for (t in (2:n.occasions)){
          # State process: draw S(t) given S(t-1)
          z[i,t] ~ dcat(ps[z[i,t-1], stage[i,t-1], t-1,1:4])
          # Observation process: draw O(t) given S(t)
          y[i,t] ~ dcat(po[z[i,t], t-1,1:3])
        } #t
      } #i
      for (i in 1:M){
        for (t in 2:n.occasions){
          al.d[i,t-1] <- equals(z[i,t], 3) # alive depre
        } #t
        alive.d[i] <- sum(al.d[i,1:(n.occasions-1)])}
      
      for (t in 1:(n.occasions-1)){
        ND[t] <- sum(al.d[1:M,t]) # Actual population size depre
      } #t
      for (i in 1:M){
        w[i] <- 1-equals(alive.d[i],0)
      } #i
      
      Nsuper <- sum(w[1:M])
    })
    
  Sim_age    <- nimbleModel(code = js_ms_age1, name = 'Sim_age', constants = consts, data = data, inits = inits, check = F)
  CSim_age  <- compileNimble(Sim_age)
  specSim_age <- configureMCMC(Sim_age, monitors=parameters,thin=75, useConjugacy=FALSE)
  specSim_age$addSampler(target = c("mean.g",'beta_s'), type = 'AF_slice')
  Sim_ageMCMC <- buildMCMC(specSim_age)
  CSim_ageMCMC<- compileNimble(Sim_ageMCMC, project=Sim_age, resetFunctions = TRUE)
  CSim_ageMCMC$run(3*m*75)
  samplesSim_age=as.matrix(CSim_ageMCMC$mvSamples)[(m+1):(m*3),]
  return(samplesSim_age)}

CH_cro_2705 <- read.csv("CH_2106_m.csv", sep=",")
marking_cro2705 <- read.csv("Mk_2106_m.csv")
statut_cro2705 <- read.csv("St_2106_m.csv")
z_init_cro2705 <- read.csv("z_2106_m.csv")
Effort_stand_boat <- read.csv("Effort_stand_boat.csv", sep="")

CHD <- as.matrix(read.csv("CH_2106_d.csv", sep=","))
markingD <- as.matrix(read.csv("Mk_2106_d.csv"))
ageD <- as.matrix(read.csv("St_2106_d.csv"))
CHD=CHD+1
zD=CHD
for (i in 1: dim(CHD)[1]){
  fi=min(which(CHD[i,]==2))
  la=max(which(CHD[i,]==2))
  if (fi >1){zD[i,1:(fi-1)]=1}
  zD[i,fi:la]=3
  if (la <dim(CHD)[2]){zD[i,(la+1):dim(CHD)[2]]=4}
}




nmax=176
n.occasions=16
CH<-as.matrix(CH_cro_2705)
z<-as.matrix(z_init_cro2705)
marking<-as.matrix(marking_cro2705)
age<-as.matrix(statut_cro2705)
Effort<-as.vector(Effort_stand_boat$x)
CH[101:102,5]=3
nz<-200
age[age==2]<-1
age[age==3]<-2
CH.du=cbind(rep(3,dim(CH)[1]),CH)
z.du=cbind(rep(1,dim(z)[1]),z)
marking.du<-cbind(rep(1,dim(CH)[1]),marking)
effort<-c(0,Effort)
age <- cbind(rep(0,c(dim(CH)[1])), age)
CH.aug<-rbind(CH.du, matrix(3, ncol=dim(CH.du)[2], nrow=nz))
z.aug <- rbind(z.du, matrix(1,ncol=dim(z.du)[2], nrow=nz))
age.aug <- rbind(age, matrix(2, ncol=dim(age)[2], nrow=nz))
marking.aug<-rbind(marking.du, matrix(3, ncol=dim(marking.du)[2], nrow=nz))

nz<-100
ageD[ageD==2]<-1
ageD[ageD==3]<-2
CHD.du=cbind(rep(3,dim(CHD)[1]),CHD)
zD.du=cbind(rep(1,dim(zD)[1]),zD)
markingD.du<-cbind(rep(1,dim(CHD)[1]),markingD)
ageD <- cbind(rep(0,c(dim(CHD)[1])), ageD)
CHD.aug<-rbind(CHD.du, matrix(3, ncol=dim(CHD.du)[2], nrow=nz))
zD.aug <- rbind(zD.du, matrix(1,ncol=dim(zD.du)[2], nrow=nz))
ageD.aug <- rbind(ageD, matrix(2, ncol=dim(ageD)[2], nrow=nz))
markingD.aug<-rbind(markingD.du, matrix(3, ncol=dim(markingD.du)[2], nrow=nz))

zinit=z.aug
zinit[,1]=NA
zinitD=zD.aug
zinitD[,1]=NA



parameters <- c("mean.p","mean.q", "mean.g",
                'beta_effort',
                "mean.gamma", 'beta_s', "beta_state",
                 "Nsuper", "NC","ND", "S1", "S2", "g",'q','p')

parametersD <- c("mean.p", "mean.g",
                'beta_effort',"beta_s",
                "Nsuper", "ND", "S1", "g",'p')


###model cpl
stage=age.aug+1
stageD=ageD.aug+1

consts=list(M = dim(CH.aug)[1],
             n.occasions = dim(CH.aug)[2], stage=stage)
initsage <-list(z=zinit,mean.S=rep(2,2), mean.p=-1, mean.q=-1,
                 mean.g=0,beta_effort=1,
                 mean.gamma=0.05, beta_s=0.25, beta_state=0.5)
dataage <- list(y = as.matrix(CH.aug),
                     effort=effort,time=c(1:17))

constsD=list(M = dim(CHD.aug)[1],
            n.occasions = dim(CHD.aug)[2], stage=stageD)
initsageD <-list(z=zinitD,mean.S=rep(2,2), mean.p=-1,
                mean.g=0,beta_effort=1,
                 beta_s=0.25)
dataageD <- list(y = as.matrix(CHD.aug),
                effort=effort,time=c(1:17))

resage=array(NA,dim=c(77,6))
resageD=array(NA,dim=c(55,6))

m=1500

print(Sys.time())
clu=makeCluster(3)
rbage=parLapply(cl=clu,X=1:3,fun=runim_ageC,parameters,consts,dataage,initsage,m)
rbage=simplify2array(rbage)
Rhatage=Rhatfun(rbage[,,1],rbage[,,2],rbage[,,3],m*2)
###store the results
resage[,1]=   apply(rbage,2, mean)##store the mean of the parameters
resage[,2]= apply(rbage,2, sd)##store the standard deviation of the parameters
resage[,3]= apply(rbage,2, quantile,0.025)
resage[,4]= apply(rbage,2, quantile,0.975)
resage[,5]= Rhatage
resage[,6]=apply(rbage,2, pval)
print(Sys.time())
rownames(resage)=rownames(t(rbage[,,1]))
save(resage,rbage,file=paste("CrozetFin.Rdata",sep=""))
print(Sys.time())

rbageD=parLapply(cl=clu,X=1:3,fun=runim_ageD,parametersD,constsD,dataageD,initsageD,m)
rbageD=simplify2array(rbageD)
RhatageD=Rhatfun(rbageD[,,1],rbageD[,,2],rbageD[,,3],m*2)
###store the results
resageD[,1]=   apply(rbageD,2, mean)##store the mean of the parameters
resageD[,2]= apply(rbageD,2, sd)##store the standard deviation of the parameters
resageD[,3]= apply(rbageD,2, quantile,0.025)
resageD[,4]= apply(rbageD,2, quantile,0.975)
resageD[,5]= RhatageD
resageD[,6]=apply(rbageD,2, pval)
print(Sys.time())
rownames(resageD)=rownames(t(rbageD[,,1]))
save(resageD,rbageD,file=paste("DFin.Rdata",sep=""))

stopCluster(clu)


