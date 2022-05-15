library(dplyr)
library(sampling)
library(MASS)
N<-100000
rep<-1000

############### SAMPLING  FUNCTION ###############

generate.sample <- function(
  #Randomly generates PS with equal inclusion
  #Systemic random samples NPS
  population, size_ps = size_ps, size_nps = size_nps,k=k){
  PS<-NULL
  NPS<-NULL #starts with empty 
  for (i in (1:k)) {
    domain<-population[population$dom==i,] #select cases that belong to i domain
    PS <- rbind(PS,domain[sample.int(nrow(domain), size_ps),]) #add domains
    lam <- domain$lambda
    lam<-inclusionprobabilities(lam,size_nps)
    NPS <- rbind(NPS,domain[as.logical(UPrandomsystematic(lam)), ])
  }
  PS$sample <- 1
  NPS$sample <- 0
  combined_sample <- full_join(PS, NPS)
  return(combined_sample)
}

#Elements that will be saved in a list by each function:

#1=estimates
#2=mean squared error of probability sample
#3=mean squared error of non-probability sample
#4=mean squared error of combined estimator
#5=root mean squared error of probability sample
#6=root mean squared error of non-probability sample
#7=root mean squared error of combined estimator
#8=average root mean squared error of probability sample 
#9= average root mean squared error of non-probability sample 
#10= average root mean squared error of combined estimator
#11= bias of non-probability sample 
#12=bias of combined estimator

######## FUNCTION FOR CATEGORIES OF EQUAL SIZE - MODEL A ########

equal.size.1<-function(k,c,n.pk,n.npk,corr){
  seed<-as.integer(paste0(k,c,10))
  set.seed(seed)
  pop<-mvrnorm(N,c(0,0,0),matrix(c(1,0,corr,0,1,corr,corr,corr,1),nrow=3,ncol=3)) 
  d<-sqrt(80/3)
  pop<-as.data.frame(pop)
  names(pop)<-c("X1","X2","W")
  pop$yc<-10+d*pop$X1+d*pop$X2+rnorm(N,0,d) 
  #two rnorm= per domain and observation (N in the pop). Still rnorm per obs., 
  #Total contribution of two rnorm to be the same as is it now. Total variance to be= d-power-2
  #choose= which percentage in domain specific and which in the individual.
  br<-quantile(pop$yc,probs = (0:(c))/c)      
  pop$y<-cut(pop$yc,breaks=br,labels=1:c,include.lowest = T) 
  pop$dom<-as.factor(sample(1:k,N,replace=T)) #             
  pop$lambda<-1/(1+exp(-(2*(pop$W-0.75)))) 
  table(pop$dom,pop$y) 
  
  estim<-matrix(0,rep,k*c*4)
  mse.ps<-matrix(0,1,k*c)
  mse.nps<-matrix(0,1,k*c)
  mse.comb<-matrix(0,1,k*c)
  bias.nps<-matrix(0,1,k*c)
  bias.comb<-matrix(0,1,k*c)
  
  for(i in 1:rep){
    samp<-generate.sample(pop,n.pk,n.npk,k)
    prob<-samp[which(samp$sample=="1"),]
    np<-samp[which(samp$sample=="0"),]
    probtab<-prop.table(table(prob$dom,prob$y),margin=1)
    nprobtab<-prop.table(table(np$dom,np$y),margin=1)
    true<-prop.table(table(pop$dom,pop$y),margin=1)
    #sample for each domain separately 
    
    v.kc<-nprobtab*(1-nprobtab)
    sigma2.k<-(n.pk/(n.pk-1))*((nprobtab-probtab)^2-v.kc*(1+n.npk/n.pk)/(n.npk-1))
    sigma2.k<-apply(sigma2.k,MARGIN = 1,FUN = mean) 
    new.sigma<-matrix(sigma2.k,nrow=k,ncol=c)
    #new.sigma[new.sigma<0]<-0
    wkc<-((n.npk-1)*new.sigma+v.kc)/((n.npk-1)*(1-1/n.pk)*new.sigma+(1+n.npk/n.pk)*v.kc)
    wkc[wkc<0]<-0
    wkc[wkc>1]<-1
    ccc<-wkc*probtab+(1-wkc)*nprobtab
    
    estim[i,1:(k*c)]<-c(t(probtab))            #PROBABILITY SAMPLE
    estim[i,(k*c+1):(2*k*c)]<-c(t(nprobtab))   #NON-PROBABILITY SAMPLE
    estim[i,(2*k*c+1):(3*k*c)]<-c(t(ccc))      #COMBINED ESTIMATOR
    estim[i,(3*k*c+1):(4*k*c)]<-c(t(true))     #TRUE VALUES
    
    mse.ps[1,1:(k*c)]<-mse.ps[1,1:(k*c)]+((estim[i,1:(k*c)]-estim[i,(3*k*c+1):(4*k*c)])^2)/rep
    mse.nps[1,1:(k*c)]<-mse.nps[1,1:(k*c)]+((estim[i,(k*c+1):(2*k*c)]-estim[i,(3*k*c+1):(4*k*c)])^2)/rep
    mse.comb[1,1:(k*c)]<-mse.comb[1,1:(k*c)]+((estim[i,(2*k*c+1):(3*k*c)]-estim[i,(3*k*c+1):(4*k*c)])^2)/rep
    
    bias.nps[1,1:(k*c)]<-estim[i,(k*c+1):(2*k*c)]-estim[i,(3*k*c+1):(4*k*c)]
    bias.comb[1,1:(k*c)]<-estim[i,(2*k*c+1):(3*k*c)]-estim[i,(3*k*c+1):(4*k*c)]
  }
  
  rmse.ps<-sqrt(mse.ps)
  rmse.nps<-sqrt(mse.nps)
  rmse.comb<-sqrt(mse.comb)
  armse.ps<-matrix(0,1,k)
  armse.nps<-matrix(0,1,k)
  armse.comb<-matrix(0,1,k)
  
  for(n in 1:(k)){armse.ps[1,n]<-mean(rmse.ps[1,(n*c-c+1):(c*n)])}
  
  for(n in 1:(k)){armse.nps[1,n]<-mean(rmse.nps[1,(n*c-c+1):(c*n)])}
  
  for(n in 1:(k)){armse.comb[1,n]<-mean(rmse.comb[1,(n*c-c+1):(c*n)])}

  return(list(estim,mse.ps,mse.nps,mse.comb,
              rmse.ps,rmse.nps,rmse.comb,
              armse.ps,armse.nps,armse.comb,
              bias.nps,bias.comb))
  
}


######## FUNCTION FOR CATEGORIES OF EQUAL SIZE - MODEL B ########

equal.size.2<-function(k,c,n.pk,n.npk,corr){
  seed<-as.integer(paste0(n.pk*2,c+3,k))
  set.seed(seed)
  pop<-mvrnorm(N,c(0,0,0),matrix(c(1,0,corr,0,1,corr,corr,corr,1),nrow=3,ncol=3)) #create population of 3 variables
  d<-sqrt(80/3)
  pop<-as.data.frame(pop)
  names(pop)<-c("X1","X2","W")
  pop$yc<-10+d*pop$X1+d*pop$X2+rnorm(N,0,d)   
  br<-quantile(pop$yc,probs = (0:(c))/c)      
  pop$y<-cut(pop$yc,breaks=br,labels=1:c,include.lowest = T) 
  pop$dom<-as.factor(sample(1:k,N,replace=T))              
  pop$lambda<-1/(1+exp(-(2*(pop$W-0.75)))) 
  table(pop$dom,pop$y) 
  
  estim<-matrix(0,rep,k*c*4)
  
  mse.ps<-matrix(0,1,k*c)
  mse.nps<-matrix(0,1,k*c)
  mse.comb<-matrix(0,1,k*c)
  
  bias.nps<-matrix(0,1,k*c)
  bias.comb<-matrix(0,1,k*c)
  
  for(i in 1:rep){
    samp<-generate.sample(pop,n.pk,n.npk,k)
    prob<-samp[which(samp$sample=="1"),]
    np<-samp[which(samp$sample=="0"),]
    probtab<-prop.table(table(prob$dom,prob$y),margin=1)
    nprobtab<-prop.table(table(np$dom,np$y),margin=1)
    true<-prop.table(table(pop$dom,pop$y),margin=1)
    #sample for each domain separately 
    
    v.kc<-nprobtab*(1-nprobtab)
    bc<-matrix(colMeans(nprobtab-probtab),nrow=k,ncol=c,byrow=T) #repeat for every row
    sigma2.k<-(n.pk/(n.pk-1))*((nprobtab-probtab)^2-(1-1/n.pk)*(bc^2)-(2/n.pk)*bc*nprobtab)-v.kc*(1+n.npk/n.pk)/(n.npk-1)
    sigma2.k<-apply(sigma2.k,MARGIN = 1,FUN = mean)
    new.sigma<-matrix(sigma2.k,nrow=k,ncol=c) #repeat for every column
    #new.sigma[new.sigma<0]<-0
    wkc<-((n.npk-1)*(bc^2+new.sigma)+v.kc)/((n.npk-1)*(1-1/n.pk)*(bc^2+new.sigma)+(1+n.npk/n.pk)*v.kc+((n.npk-1)/n.pk)*bc*(2*nprobtab-1))
    wkc[wkc<0]<-0
    wkc[wkc>1]<-1
    ccc<-wkc*probtab+(1-wkc)*nprobtab
    
    estim[i,1:(k*c)]<-c(t(probtab))            #PROBABILITY SAMPLE
    estim[i,(k*c+1):(2*k*c)]<-c(t(nprobtab))   #NON-PROBABILITY SAMPLE
    estim[i,(2*k*c+1):(3*k*c)]<-c(t(ccc))      #COMBINED ESTIMATOR
    estim[i,(3*k*c+1):(4*k*c)]<-c(t(true))     #TRUE VALUES
    
    mse.ps[1,1:(k*c)]<-mse.ps[1,1:(k*c)]+((estim[i,1:(k*c)]-estim[i,(3*k*c+1):(4*k*c)])^2)/rep
    mse.nps[1,1:(k*c)]<-mse.nps[1,1:(k*c)]+((estim[i,(k*c+1):(2*k*c)]-estim[i,(3*k*c+1):(4*k*c)])^2)/rep
    mse.comb[1,1:(k*c)]<-mse.comb[1,1:(k*c)]+((estim[i,(2*k*c+1):(3*k*c)]-estim[i,(3*k*c+1):(4*k*c)])^2)/rep
    
    bias.nps[1,1:(k*c)]<-estim[i,(k*c+1):(2*k*c)]-estim[i,(3*k*c+1):(4*k*c)]
    bias.comb[1,1:(k*c)]<-estim[i,(2*k*c+1):(3*k*c)]-estim[i,(3*k*c+1):(4*k*c)]
    
  }
  
  rmse.ps<-sqrt(mse.ps)
  rmse.nps<-sqrt(mse.nps)
  rmse.comb<-sqrt(mse.comb)
  
  armse.ps<-matrix(0,1,k)
  armse.nps<-matrix(0,1,k)
  armse.comb<-matrix(0,1,k)
  
  for(i in 1:(k)){armse.ps[1,i]<-mean(rmse.ps[1,(i*c-c+1):(c*i)])}
  
  for(i in 1:(k)){armse.nps[1,i]<-mean(rmse.nps[1,(i*c-c+1):(c*i)])}
  
  for(i in 1:(k)){armse.comb[1,i]<-mean(rmse.comb[1,(i*c-c+1):(c*i)])}
  
  return(list(estim,mse.ps,mse.nps,mse.comb,
              rmse.ps,rmse.nps,rmse.comb,
              armse.ps,armse.nps,armse.comb,
              bias.nps,bias.comb))
  
}



######## FUNCTION FOR CATEGORIES OF UNEQUAL SIZE - MODEL A ########

unequal.size.1<-function(k,c,n.pk,n.npk,corr){
  seed<-as.integer(paste0(k,3,c,c*2))
  set.seed(seed)
  pop<-mvrnorm(N,c(0,0,0),matrix(c(1,0,corr,0,1,corr,corr,corr,1),nrow=3,ncol=3)) 
  d<-sqrt(80/3)
  pop<-as.data.frame(pop)
  names(pop)<-c("X1","X2","W")
  pop$yc<-10+d*pop$X1+d*pop$X2+rnorm(N,0,d)   
  
  br<-NA     
  if (c=="3"){
    br=c(0,.29,.65,1)
  } else if (c=="5") {
    br=c(0,.09,.29,.65,.87,1) 
  } else if (c=="8") {
    br=c(0,.09,.20,.28,.42,.55,.65,.87,1) 
  } else if (c=="15") {
    br=c(0,.03,.06,.11,.18,.27,.37,.49,.58,.68,.77,.85,.93,.95,.97,1) 
  }
  
  br<-quantile(pop$yc,probs=br)
  pop$y<-cut(pop$yc,breaks=br,labels=1:c,include.lowest = T) 
  pop$dom<-as.factor(sample(1:k,N,replace=T)) #             
  pop$lambda<-1/(1+exp(-(2*(pop$W-0.75)))) 
  table(pop$dom,pop$y) 
  
  estim<-matrix(0,rep,k*c*4)
  mse.ps<-matrix(0,1,k*c)
  mse.nps<-matrix(0,1,k*c)
  mse.comb<-matrix(0,1,k*c)
  bias.nps<-matrix(0,1,k*c)
  bias.comb<-matrix(0,1,k*c)
  
  for(i in 1:rep){
    samp<-generate.sample(pop,n.pk,n.npk,k)
    prob<-samp[which(samp$sample=="1"),]
    np<-samp[which(samp$sample=="0"),]
    probtab<-prop.table(table(prob$dom,prob$y),margin=1)
    nprobtab<-prop.table(table(np$dom,np$y),margin=1)
    true<-prop.table(table(pop$dom,pop$y),margin=1)
    #sample for each domain separately 
    
    v.kc<-nprobtab*(1-nprobtab)
    sigma2.k<-(n.pk/(n.pk-1))*((nprobtab-probtab)^2-v.kc*(1+n.npk/n.pk)/(n.npk-1))
    sigma2.k<-apply(sigma2.k,MARGIN = 1,FUN = mean) 
    new.sigma<-matrix(sigma2.k,nrow=k,ncol=c)
    #new.sigma[new.sigma<0]<-0
    wkc<-((n.npk-1)*new.sigma+v.kc)/((n.npk-1)*(1-1/n.pk)*new.sigma+(1+n.npk/n.pk)*v.kc)
    wkc[wkc<0]<-0
    wkc[wkc>1]<-1
    ccc<-wkc*probtab+(1-wkc)*nprobtab
    
    estim[i,1:(k*c)]<-c(t(probtab))            #PROBABILITY SAMPLE
    estim[i,(k*c+1):(2*k*c)]<-c(t(nprobtab))   #NON-PROBABILITY SAMPLE
    estim[i,(2*k*c+1):(3*k*c)]<-c(t(ccc))      #COMBINED ESTIMATOR
    estim[i,(3*k*c+1):(4*k*c)]<-c(t(true))     #TRUE VALUES
    
    mse.ps[1,1:(k*c)]<-mse.ps[1,1:(k*c)]+((estim[i,1:(k*c)]-estim[i,(3*k*c+1):(4*k*c)])^2)/rep
    mse.nps[1,1:(k*c)]<-mse.nps[1,1:(k*c)]+((estim[i,(k*c+1):(2*k*c)]-estim[i,(3*k*c+1):(4*k*c)])^2)/rep
    mse.comb[1,1:(k*c)]<-mse.comb[1,1:(k*c)]+((estim[i,(2*k*c+1):(3*k*c)]-estim[i,(3*k*c+1):(4*k*c)])^2)/rep
    
    bias.nps[1,1:(k*c)]<-estim[i,(k*c+1):(2*k*c)]-estim[i,(3*k*c+1):(4*k*c)]
    bias.comb[1,1:(k*c)]<-estim[i,(2*k*c+1):(3*k*c)]-estim[i,(3*k*c+1):(4*k*c)]
  }
  
  rmse.ps<-sqrt(mse.ps)
  rmse.nps<-sqrt(mse.nps)
  rmse.comb<-sqrt(mse.comb)
  armse.ps<-matrix(0,1,k)
  armse.nps<-matrix(0,1,k)
  armse.comb<-matrix(0,1,k)
  
  for(n in 1:(k)){armse.ps[1,n]<-mean(rmse.ps[1,(n*c-c+1):(c*n)])}
  
  for(n in 1:(k)){armse.nps[1,n]<-mean(rmse.nps[1,(n*c-c+1):(c*n)])}
  
  for(n in 1:(k)){armse.comb[1,n]<-mean(rmse.comb[1,(n*c-c+1):(c*n)])}
  
  return(list(estim,mse.ps,mse.nps,mse.comb,
              rmse.ps,rmse.nps,rmse.comb,
              armse.ps,armse.nps,armse.comb,
              bias.nps,bias.comb))
  
}





######## FUNCTION FOR CATEGORIES OF UNEQUAL SIZE - MODEL B ########

unequal.size.2<-function(k,c,n.pk,n.npk,corr){
  seed<-as.integer(paste0(k*3,n.pk,c*4))
  set.seed(seed)
  pop<-mvrnorm(N,c(0,0,0),matrix(c(1,0,corr,0,1,corr,corr,corr,1),nrow=3,ncol=3)) 
  d<-sqrt(80/3)
  pop<-as.data.frame(pop)
  names(pop)<-c("X1","X2","W")
  pop$yc<-10+d*pop$X1+d*pop$X2+rnorm(N,0,d)   
  
  br<-NA     
  if (c=="3"){
    br=c(0,.29,.65,1)
  } else if (c=="5") {
    br=c(0,.09,.29,.65,.87,1) 
  } else if (c=="8") {
    br=c(0,.09,.20,.28,.42,.55,.65,.87,1) 
  } else if (c=="15") {
    br=c(0,.03,.06,.11,.18,.27,.37,.49,.58,.68,.77,.85,.93,.95,.97,1) 
  }
  
  br<-quantile(pop$yc,probs=br)
  pop$y<-cut(pop$yc,breaks=br,labels=1:c,include.lowest = T) 
  pop$dom<-as.factor(sample(1:k,N,replace=T)) #             
  pop$lambda<-1/(1+exp(-(2*(pop$W-0.75)))) 
  table(pop$dom,pop$y) 
  
  estim<-matrix(0,rep,k*c*4)
  mse.ps<-matrix(0,1,k*c)
  mse.nps<-matrix(0,1,k*c)
  mse.comb<-matrix(0,1,k*c)
  bias.nps<-matrix(0,1,k*c)
  bias.comb<-matrix(0,1,k*c)
  
  for(i in 1:rep){
    samp<-generate.sample(pop,n.pk,n.npk,k)
    prob<-samp[which(samp$sample=="1"),]
    np<-samp[which(samp$sample=="0"),]
    probtab<-prop.table(table(prob$dom,prob$y),margin=1)
    nprobtab<-prop.table(table(np$dom,np$y),margin=1)
    true<-prop.table(table(pop$dom,pop$y),margin=1)
    #sample for each domain separately 
    
    v.kc<-nprobtab*(1-nprobtab)
    bc<-matrix(colMeans(nprobtab-probtab),nrow=k,ncol=c,byrow=T) #repeat for every row
    sigma2.k<-(n.pk/(n.pk-1))*((nprobtab-probtab)^2-(1-1/n.pk)*(bc^2)-(2/n.pk)*bc*nprobtab)-v.kc*(1+n.npk/n.pk)/(n.npk-1)
    sigma2.k<-apply(sigma2.k,MARGIN = 1,FUN = mean)
    new.sigma<-matrix(sigma2.k,nrow=k,ncol=c) #repeat for every column
    #new.sigma[new.sigma<0]<-0
    wkc<-((n.npk-1)*(bc^2+new.sigma)+v.kc)/((n.npk-1)*(1-1/n.pk)*(bc^2+new.sigma)+(1+n.npk/n.pk)*v.kc+((n.npk-1)/n.pk)*bc*(2*nprobtab-1))
    wkc[wkc<0]<-0
    wkc[wkc>1]<-1
    ccc<-wkc*probtab+(1-wkc)*nprobtab
    
    estim[i,1:(k*c)]<-c(t(probtab))            #PROBABILITY SAMPLE
    estim[i,(k*c+1):(2*k*c)]<-c(t(nprobtab))   #NON-PROBABILITY SAMPLE
    estim[i,(2*k*c+1):(3*k*c)]<-c(t(ccc))      #COMBINED ESTIMATOR
    estim[i,(3*k*c+1):(4*k*c)]<-c(t(true))     #TRUE VALUES
    
    mse.ps[1,1:(k*c)]<-mse.ps[1,1:(k*c)]+((estim[i,1:(k*c)]-estim[i,(3*k*c+1):(4*k*c)])^2)/rep
    mse.nps[1,1:(k*c)]<-mse.nps[1,1:(k*c)]+((estim[i,(k*c+1):(2*k*c)]-estim[i,(3*k*c+1):(4*k*c)])^2)/rep
    mse.comb[1,1:(k*c)]<-mse.comb[1,1:(k*c)]+((estim[i,(2*k*c+1):(3*k*c)]-estim[i,(3*k*c+1):(4*k*c)])^2)/rep
    
    bias.nps[1,1:(k*c)]<-estim[i,(k*c+1):(2*k*c)]-estim[i,(3*k*c+1):(4*k*c)]
    bias.comb[1,1:(k*c)]<-estim[i,(2*k*c+1):(3*k*c)]-estim[i,(3*k*c+1):(4*k*c)]
  }
  
  rmse.ps<-sqrt(mse.ps)
  rmse.nps<-sqrt(mse.nps)
  rmse.comb<-sqrt(mse.comb)
  armse.ps<-matrix(0,1,k)
  armse.nps<-matrix(0,1,k)
  armse.comb<-matrix(0,1,k)
  
  for(n in 1:(k)){armse.ps[1,n]<-mean(rmse.ps[1,(n*c-c+1):(c*n)])}
  
  for(n in 1:(k)){armse.nps[1,n]<-mean(rmse.nps[1,(n*c-c+1):(c*n)])}
  
  for(n in 1:(k)){armse.comb[1,n]<-mean(rmse.comb[1,(n*c-c+1):(c*n)])}
  
  return(list(estim,mse.ps,mse.nps,mse.comb,
              rmse.ps,rmse.nps,rmse.comb,
              armse.ps,armse.nps,armse.comb,
              bias.nps,bias.comb))
}
