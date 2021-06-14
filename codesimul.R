library(parallel)
library(boot)

#####Modifier no ICI
args = commandArgs(TRUE) #ca ca récupère les arguments passer par bash, les mettre a la suite pour en mettre d'autres
no= as.numeric(args[1])
nsimul= as.numeric(args[2])
#no=1



# Define parameter values
phi=c(0.5,1.5,2,0.5)
im<-0.05 ## imigration
recap=c(-1.5,2,1.5,-0.75,1.5)
dep=c(-0.2,0.1)
nsimul=10
m=1000
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


runim_cpl<-function(val,parameters,consts,data,inits,m){
  library(nimble)
  #   ###
  js_ms_marking_age1 <- nimbleCode({
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
    
    beta_marking~dnorm(0,0.001)
    beta_effort~dnorm(0,0.001)
    beta_s~dnorm(0,0.001)
    beta_state~dnorm(0,0.001)
    
    mean.g~dnorm(0,0.001)
    mean.gamma~dunif(0,0.25)
    
    mean.S[2]~dnorm(0,0.001)
    mean.S[1]~dnorm(0,0.001)
    
    S[3]<-mean.S[1]
    S[2]<-mean.S[1]
    S[4]<-mean.S[2]
    S[5]<-mean.S[2]
    S[1]<--1000
    mean.p~dnorm(0,0.001)
    mean.q~dnorm(0,0.001)
    for (i in 1:5){
      # Define probabilities of state S(t+1) given S(t)
      for (t in 1:(n.occasions-1)){
        ps[1,i,t,1] <- 1-mean.gamma
        ps[1,i,t,2] <- mean.gamma*(1-ilogit(-mean.g+beta_s*time[t]))
        ps[1,i,t,3] <- ilogit(-mean.g+beta_s*time[t])*mean.gamma
        ps[1,i,t,4] <- 0
        ps[2,i,t,1] <- 0
        ps[2,i,t,2] <- ilogit(S[i])*(1-ilogit(-mean.g+beta_s*time[t]))
        ps[2,i,t,3] <- ilogit(-mean.g+beta_s*time[t])*ilogit(S[i])
        ps[2,i,t,4] <- 1-ilogit(S[i])
        ps[3,i,t,1] <- 0
        ps[3,i,t,2] <- 0
        ps[3,i,t,3] <- ilogit(S[i]+beta_state)
        ps[3,i,t,4] <- 1-ilogit(S[i]+beta_state)
        ps[4,i,t,1] <- 0
        ps[4,i,t,2] <- 0
        ps[4,i,t,3] <- 0
        ps[4,i,t,4] <- 1

          # Define probabilities of O(t) given S(t)
        po[1,i,t,1] <- 0
        po[1,i,t,2] <- 0
        po[1,i,t,3] <- 1
        po[2,i,t,1] <- ilogit(mean.q+beta_marking*mark[i])
        po[2,i,t,2] <- 0
        po[2,i,t,3] <- 1-ilogit(mean.q+beta_marking*mark[i])
        po[3,i,t,1] <- ilogit(mean.q+beta_marking*mark[i])*(1-ilogit(mean.p+beta_effort*effort[t+1]+beta_marking*mark[i]))
        po[3,i,t,2] <- ilogit(mean.p+beta_effort*effort[t+1]+beta_marking*mark[i])
        po[3,i,t,3] <- 1-ilogit(mean.p+beta_effort*effort[t+1]+beta_marking*mark[i])- ilogit(mean.q+beta_marking*mark[i])*(1-ilogit(mean.p+beta_effort*effort[t+1]+beta_marking*mark[i]))
        po[4,i,t,1] <- 0
        po[4,i,t,2] <- 0
        po[4,i,t,3] <- 1
      } #t
    } #i
    
    
    # Likelihood
    for (i in 1:M){
      # Define latent state at first occasion
      z[i,1] <-1  # Make sure that all M individuals are in state 1 at t=1
      for (t in (2:n.occasions)){
        # State process: draw S(t) given S(t-1)
        z[i,t] ~ dcat(ps[z[i,t-1], stage[i,t-1], t-1,1:4])
        # Observation process: draw O(t) given S(t)
        y[i,t] ~ dcat(po[z[i,t], stage[i,t], t-1,1:3])
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
  print(Sys.time())
   Sim_cpl    <- nimbleModel(code = js_ms_marking_age1, name = 'Sim_cpl', constants = consts, data = data, inits = inits, check = F)
  CSim_cpl  <- compileNimble(Sim_cpl)
  specSim_cpl <- configureMCMC(Sim_cpl, monitors=parameters,thin=75, useConjugacy=FALSE)
  specSim_cpl$addSampler(target = c("mean.S",'beta_state'), type = 'AF_slice')
  specSim_cpl$addSampler(target = c("mean.g",'beta_s'), type = 'AF_slice')
  specSim_cpl$addSampler(target = c("mean.p","mean.q","mean.gamma"), type =  'AF_slice')
  specSim_cpl$addSampler(target = c("mean.gamma"), type =  'RW')
  Sim_cplMCMC <- buildMCMC(specSim_cpl)
  CSim_cplMCMC<- compileNimble(Sim_cplMCMC, project=Sim_cpl, resetFunctions = TRUE)
   CSim_cplMCMC$run(m*3*75)
   samplesSim_cpl=as.matrix(CSim_cplMCMC$mvSamples)[(m+1):(m*3),]
  return(samplesSim_cpl)}

runim_age<-function(val,parameters,consts,data,inits,m){
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
    
    beta_marking~dnorm(0,0.001)
    beta_effort~dnorm(0,0.001)
    beta_s~dnorm(0,0.001)
    beta_state~dnorm(0,0.001)
    
    mean.g~dnorm(0,0.001)
    mean.gamma~dunif(0,0.25)
    
    mean.S[2]~dnorm(0,0.001)
    mean.S[1]~dnorm(0,0.001)
    
    S[3]<-mean.S[1]
    S[2]<-mean.S[1]
    S[4]<-mean.S[2]
    S[5]<-mean.S[2]
    S[1]<--1000
    mean.p~dnorm(0,0.001)
    mean.q~dnorm(0,0.001)
   
      # Define probabilities of state S(t+1) given S(t)
      for (t in 1:(n.occasions-1)){
        for (i in 1:5){ 
          ps[1,i,t,1] <- 1-mean.gamma
        ps[1,i,t,2] <- mean.gamma*(1-ilogit(-mean.g+beta_s*time[t]))
        ps[1,i,t,3] <- ilogit(-mean.g+beta_s*time[t])*mean.gamma
        ps[1,i,t,4] <- 0
        ps[2,i,t,1] <- 0
        ps[2,i,t,2] <- ilogit(S[i])*(1-ilogit(-mean.g+beta_s*time[t]))
        ps[2,i,t,3] <- ilogit(-mean.g+beta_s*time[t])*ilogit(S[i])
        ps[2,i,t,4] <- 1-ilogit(S[i])
        ps[3,i,t,1] <- 0
        ps[3,i,t,2] <- 0
        ps[3,i,t,3] <- ilogit(S[i]+beta_state)
        ps[3,i,t,4] <- 1-ilogit(S[i]+beta_state)
        ps[4,i,t,1] <- 0
        ps[4,i,t,2] <- 0
        ps[4,i,t,3] <- 0
        ps[4,i,t,4] <- 1
        } #i  
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
# Function to simulate capture-recapture data under the JS model
simul.js <- function(PHI,x,PBC,B, Bgroup,b ,PBgroup,PCgroup, group){
     # B=c(BJi,BAi);group=c(group.j,group.a)
    N=dim(PBC)[1]
  n.occasions <- dim(PBC)[2]
  CH.p<-CH.sur<- matrix(0, ncol = n.occasions, nrow = N)
  # Define a vector with the occasion of entering the population
  ent.occ <-B
  ent.occgr<-numeric()
  nbgroup.soc<-20
  depred.group<-p.group<-matrix(0,nrow=nbgroup.soc,ncol=n.occasions)
 
  for (t in 1:(n.occasions-1)){
    ent.occgr <- c(ent.occgr, rep(t, Bgroup[t]))
  }
  # Group loop Simulate Depredation and capture
  for (i in 1: nbgroup.soc){
    depred.group[i, ent.occgr[i]] <- 1 # Write 1 when ind. enters the pop.
    if(depred.group[i,ent.occgr[i]]==1){p.group[i,t] <- rbinom(1, 1,PCgroup[ent.occgr[i]])}  #observed from the coast?
    if (ent.occgr[i] == n.occasions) next
    depred=1 ###1: non depredateur, 2: depredateur, donc on commence par defaut a 1 pour tous les individus
    for (t in (ent.occgr[i]+1):n.occasions){
      ###########probabilite de devenir depredateur
      if(depred==1){depred<- rbinom(1, 1, b[t])+1}
      depred.group[i,t] <- depred
      if(depred.group[i,t]==1){p.group[i,t] <- rbinom(1, 1,PCgroup[t])}  #observed from the coast?
      if(depred.group[i,t]==2){p.group[i,t] <- rbinom(1, 1,PBgroup[t])*2 #observed from the boat?
      if (p.group[i,t]==0){p.group[i,t] <- rbinom(1, 1,PCgroup[t])}} #if not observed from the boat,
    } #t
  }
  # Simulate survival
  for (i in 1:N){
    CH.sur[i, ent.occ[i]] <-  depred.group[group[i],ent.occ[i]]
    if (ent.occ[i] == n.occasions) next
    for (t in (ent.occ[i]+1):n.occasions){
      # Bernoulli trial: has individual survived occasion?
      sur <- rbinom(1, 1, PHI[x[i,t-1],depred.group[group[i],t-1]])
      ifelse (sur==1, CH.sur[i,t] <- depred.group[group[i],t], break)
    } #t
  }
  
  ## test if social group present if yes p increase for all individuals of the group
  for (i in 1:N){
    for (t in ent.occ[i]:n.occasions){
      if(CH.sur[i,t]>0){#s'il est vivant et present
        if(p.group[group[i],t]>0){CH.p[i,t] <- rbinom(1, 1,PBC[i,t])*p.group[group[i],t]}  #observed from marking
        ### PROBA DUN INDIVIDU DETRE CAPTURE EST PBgroup*PB
      }}}
  # Full capture-recapture matrix
  CH <- CH.p
  # Remove individuals never captured
  Ntot<-colSums(CH.sur) # store all individuals including non captured
  cap.sum <- rowSums(CH)
  never <- which(cap.sum == 0)
  CH <- CH[-never,]
   x <- x[-never,]
   CH.sur2=CH.sur[-never,]
  return(list(CH=CH,CH.sur=CH.sur,CH.sur2=CH.sur2,Ntot=Ntot,never=never,p.group=p.group,depred.group=depred.group,x=x))
}



n.occasions <- 20 # Number of capture occasions
Ngroup<-20
NNJ<-40
NNA<-220



# La probabilite que les individus deviennent depredateur 
time=c(rep(-10,8),c(1:12))
b<-inv.logit(dep[1]+dep[2]*time)

#survie
phi.juv <- c(inv.logit(phi[1]),inv.logit(phi[1]+phi[4])) # Juvenile annual survival
phi.aM <- c(inv.logit(phi[2]),inv.logit(phi[2]+phi[4])) # Adult annual survival
phi.aF <- c(inv.logit(phi[3]),inv.logit(phi[3]+phi[4])) # Adult annual survival
PHI=rbind(phi.juv,phi.aM,phi.aF)

# Simulate immigration in the population and group arrival
imi=rep(im,n.occasions-1)
ptot=imi
cptot=rep(0,(n.occasions-1))
cptot[1]=ptot[1]
for (t in 2:(n.occasions-1)){cptot[t]=ptot[t]*prod(1-ptot[1:(t-1)])}
bb<-cptot/sum(cptot)


resage=rescpl=array(NA,dim=c(51,6,nsimul))

 NDepred=NCrozet=array(NA,dim=c(20,nsimul))


for (si in 1:nsimul){

 print(paste(si,Sys.time()))
Bgroup<-as.vector(rmultinom(1,Ngroup,bb))# migration per gorup
Nb.group<-cumsum(Bgroup)
fgroup=c()
for (t in 1:(n.occasions-1)){
  fgroup=c(fgroup,rep(t,Bgroup[t]))
}

# Associate social group to each individuals
group.j=sort(round(runif(NNJ,1,Ngroup)))
group.a=sort(round(runif(NNA,1,Ngroup)))
BJ=fgroup[group.j]
BA=fgroup[group.a]
BJ=tabulate(BJ,20)
BA=tabulate(BA,20)

# Simulate sex
sex.j<-rbinom(NNJ,1,0.5)+1
sex.a<-rbinom(NNA,1,0.5)+1

# Simulate marking
### marquage doit correspondre a l'annee de naissance de l'individu
level_marking<-function(marked,n.occasions,mark1){
  marking<-matrix(mark1, nrow=sum(marked), ncol=n.occasions)
  cmarked=c(0,cumsum(marked))
  maxocc=min(which(cmarked==sum(marked)))
  minocc=max(which(cmarked==0))
  # maxocc=n.occasions
for(t in minocc:(maxocc-1)){
    for (i in (cmarked[t]+1):cmarked[t+1]){
      marking[i,t]=mark1
      for (j in t:(n.occasions-1)){
        if(marking[i,j]==1){
          marking[i,j+1]<-round(rbinom(1,1,0.1), digits=0)+1}
        if(marking[i,j]==2){
          marking[i,j+1]<-round(rbinom(1,1,0.3), digits=0)+2}
        if(marking[i,j]==3){
          marking[i,j+1]<-3}
      }}}
  return(marking)
}
marking_J<-level_marking(BJ, n.occasions,1)
marking_A<-level_marking(BA, n.occasions,2)
marking=rbind(marking_J,marking_A)


# Simulate photo effort
Effort<-NULL
#### Simulation de l'effort et des proba de captures crozet
# et bateau au cours du temps 
for ( t in 1:n.occasions){
  Effort[t]<-round(runif(1,500,5000), digits=0)
}
Effort<-Effort/max(Effort)
PBC<-inv.logit(recap[1]+recap[2]*marking)
# Simulate capture probability per group 
PBgroup<-inv.logit(recap[3]+recap[5]*Effort)
PCgroup<-inv.logit(recap[4]+recap[5]*Effort)



# Create matrices X indicating age classes
x.j <- matrix(0, ncol = n.occasions, nrow = NNJ)
x.a <- matrix(0, ncol = n.occasions, nrow = NNA)

#juvenile
BJi=rep(1:n.occasions, BJ)
  for (i in 1:NNJ){
    if(BJi[i]>1){x.j[i,BJi[i]-1] <- 1}
   if(BJi[i]<n.occasions){ x.j[i,BJi[i]:(n.occasions-1)] <- 1}
    if(BJi[i]+10<=n.occasions){
      x.j[i,(BJi[i]+10):(n.occasions)] <- 1+sex.j[i]
    }
  }
#Adulte 
BAi=rep(1:n.occasions, BA)
  for (i in 1:NNA){
    if(BAi[i]>1){x.a[i,(1:BAi[i]-1)] <- 1}
    if(BAi[i]>10){x.a[i,BAi[i]-11] <- 1}
    
    if(BAi[i]<=n.occasions){
      x.a[i,BAi[i]:(n.occasions)] <- 1+sex.a[i]
    }}

x=rbind(x.j,x.a)
# Execute the Capture History simulation 
simNB <- simul.js(PHI,x, PBC ,c(BJi,BAi),Bgroup, b, PBgroup, PCgroup,c(group.j,group.a))
CH <- simNB$CH
CH_sur <- simNB$CH.sur
z<-simNB$CH.sur2
x<-simNB$x
marking=marking[-simNB$never,] 
never=simNB$never

for(t in 1:n.occasions){
NDepred[t,si]<-sum(CH_sur[,t]==2)
NCrozet[t,si]<-sum(CH_sur[,t]==1)
}

#### changement des indicateurs d'état absent avant d'avoir était observe pour la premiere fois 1
CH[CH==0]<-3

z[z==2]<-3
z[z==1]<-2
# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
zf<-apply(z,1,get.first)
for (i in which(zf>1)){z[i,1:(zf[i]-1)]<-1}
z[z==0]<-4



###Regrouper l'ensemble
 CH.du=cbind(rep(3,dim(CH)[1]),CH)
# CH.du=CH
z.du=cbind(rep(1,dim(z)[1]),z)
age <- cbind(rep(0,c(length(c(sex.j,sex.a))))[-never], x)
marking=cbind(rep(1,c(length(c(sex.j,sex.a))))[-never],marking) 
Effort<-c(0,Effort)
time=c(-1000,time)


# Data augmentation 
nz<-100
CH.aug<-rbind(CH.du, matrix(3, ncol=dim(CH.du)[2], nrow=nz))
z.aug <- rbind(z.du, matrix(1,ncol=dim(z.du)[2], nrow=nz))
# f.aug <- c(f, rep(20,nz))
age.aug <- rbind(age, matrix(2, ncol=dim(age)[2], nrow=nz))
marking.aug<-rbind(marking, matrix(sample(c(2,3),1), ncol=dim(marking)[2], nrow=nz))
# state.aug<-rbind(state,matrix(1,ncol=dim(state)[2], nrow=nz))

age.aug[age.aug>2]=2
zinit=z.aug
zinit[,1]=NA
m1=1
m2=dim(CH.aug)[1]


parameters <- c("mean.p","mean.q", "mean.S", "mean.g",
                'beta_marking','beta_effort',
                "mean.gamma", 'beta_s', "beta_state",
                 "Nsuper", "NC","ND")


###model cpl
stage=age.aug+1
stage[age.aug==2]=4
stage[age.aug==2&marking.aug==3]=5
stage[age.aug==1&marking.aug>1]=3
###model cpl
consts=list(M = dim(CH.aug[m1:m2,])[1],
             n.occasions = dim(CH.aug[m1:m2,])[2], stage=stage)
datacpl <- list(y = as.matrix(CH.aug[m1:m2,]), mark=c(0,1,2,2,3),effort=Effort,time=time)
initsage <-list(z=zinit[m1:m2,],mean.S=rep(2,2), mean.p=-1, mean.q=-1
                 , mean.g=0,beta_marking=0.05,beta_effort=1,
                 mean.gamma=0.05, beta_s=0.25, beta_state=0.5)
dataage <- list(y = as.matrix(CH.aug[m1:m2,]),
                     effort=Effort,time=time)



clu=makeCluster(3)
rbcpl=parLapply(cl=clu,X=1:3,fun=runim_cpl,parameters,consts,datacpl,initsage,m)
rbcpl=simplify2array(rbcpl)
Rhatcpl=Rhatfun(rbcpl[,,1],rbcpl[,,2],rbcpl[,,3],m*2)
###store the results
rescpl[,1,si]= apply(rbcpl,2, mean)##store the mean of the parameters
rescpl[,2,si]= apply(rbcpl,2, sd)##store the standard deviation of the parameters
rescpl[,3,si]= apply(rbcpl,2, quantile,0.025)
rescpl[,4,si]= apply(rbcpl,2, quantile,0.975)
rescpl[,5,si]= Rhatcpl
rescpl[,6,si]=apply(rbcpl,2, pval)

rbage=parLapply(cl=clu,X=1:3,fun=runim_age,parameters,consts,dataage,initsage,m)
rbage=simplify2array(rbage)
Rhatage=Rhatfun(rbage[,,1],rbage[,,2],rbage[,,3],m*2)
###store the results
resage[,1,si]=   apply(rbage,2, mean)##store the mean of the parameters
resage[,2,si]= apply(rbage,2, sd)##store the standard deviation of the parameters
resage[,3,si]= apply(rbage,2, quantile,0.025)
resage[,4,si]= apply(rbage,2, quantile,0.975)
resage[,5,si]= Rhatage
resage[,6,si]=apply(rbage,2, pval)
rownames(resage)=rownames(t(rbage[,,1]))
rownames(rescpl)=rownames(t(rbcpl[,,1]))
save(rescpl,resage,phi,im,recap,dep,NDepred,NCrozet, file=paste("Rsim_bstate",phi[4],"bs",dep[2],"_",no,".rdata",sep=""))
   stopCluster(clu)
   }

