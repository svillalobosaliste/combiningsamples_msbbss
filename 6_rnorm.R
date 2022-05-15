library(dplyr)
library(sampling)
library(MASS)
N<-100000
rep<-1000

######################################################################
######################## CREATE FUNCTIONS ############################
######################################################################

############### SAMPLING  FUNCTION ###############

generate.sample <- function(
  #Randomly generates PS with equal inclusion
  #Systemic random samples NPS
  population, size_ps = size_ps, size_dnpsdcxe = size_nps,k=k){
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


######## Function for categories of equal sizes in the approach of bias=0 ########

equal.size.1<-function(k,c,n.pk,n.npk,corr){
  seed<-as.integer(paste0(k,c,10),3,2)
  set.seed(seed)
  pop<-mvrnorm(N,c(0,0,0),matrix(c(1,0,corr,0,1,corr,corr,corr,1),nrow=3,ncol=3)) 
  d<-sqrt(80/3)
  fraction<-0.5 #assign in so 10% comes from domain specific and the rest still random (add to function)
  pop<-as.data.frame(pop)
  names(pop)<-c("X1","X2","W")
  pop$dom<-as.factor(sample(1:k,N,replace=T)) #   
  pop$yc<-10+d*pop$X1+d*pop$X2+rnorm(N,0,d*(sqrt(1-fraction))) #we get 90% of the variance
  
  for (x in 1:k){
    pop$yc[pop$dom==x]<-pop$yc[pop$dom==x]+rnorm(1,0,d*sqrt(fraction)) #every unit inside domain x gets same contribution. 
    #Units will b more similar
  }
  
  br<-quantile(pop$yc,probs = (0:(c))/c)      
  pop$y<-cut(pop$yc,breaks=br,labels=1:c,include.lowest = T) 
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
    wkc[is.nan(wkc)]<-(n.npk)/((n.npk-1)*(1-1/n.pk)+(1+n.npk/n.pk)) #this or 1
    
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
              bias.nps,bias.comb,true,ccc))
  
}


######## Function for categories of equal sizes in the approach of bias=a number ########

equal.size.2<-function(k,c,n.pk,n.npk,corr){
  seed<-as.integer(paste0(5,n.pk*2,c+3))
  set.seed(seed)
  pop<-mvrnorm(N,c(0,0,0),matrix(c(1,0,corr,0,1,corr,corr,corr,1),nrow=3,ncol=3)) #create population of 3 variables
  d<-sqrt(80/3)
  pop<-mvrnorm(N,c(0,0,0),matrix(c(1,0,corr,0,1,corr,corr,corr,1),nrow=3,ncol=3)) 
  d<-sqrt(80/3)
  fraction<-0.5 #assign in so 10% comes from domain specific and the rest still random (add to function)
  pop<-as.data.frame(pop)
  names(pop)<-c("X1","X2","W")
  pop$dom<-as.factor(sample(1:k,N,replace=T)) #   
  pop$yc<-10+d*pop$X1+d*pop$X2+rnorm(N,0,d*(sqrt(1-fraction))) #we get 90% of the variance
  
  for (x in 1:k){
    pop$yc[pop$dom==x]<-pop$yc[pop$dom==x]+rnorm(1,0,d*sqrt(fraction)) #every unit inside domain x gets same contribution. 
    #Units will b more similar
  }
  
  br<-quantile(pop$yc,probs = (0:(c))/c)      
  pop$y<-cut(pop$yc,breaks=br,labels=1:c,include.lowest = T) 
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
    wkc[is.nan(wkc)]<-(n.npk)/((n.npk-1)*(1-1/n.pk)+(1+n.npk/n.pk)) #this or 1
    
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
              bias.nps,bias.comb,true,ccc))
  
}



######## Function for categories of unequal sizes in the approach of bias=0 ########

unequal.size.1<-function(k,c,n.pk,n.npk,corr){
  seed<-as.integer(paste0(k,7,1,3,c,c*2))
  set.seed(seed)
  pop<-mvrnorm(N,c(0,0,0),matrix(c(1,0,corr,0,1,corr,corr,corr,1),nrow=3,ncol=3)) 
  d<-sqrt(80/3)
  fraction<-0.5 #assign in so 10% comes from domain specific and the rest still random (add to function)
  pop<-as.data.frame(pop)
  names(pop)<-c("X1","X2","W")
  pop$dom<-as.factor(sample(1:k,N,replace=T)) #   
  pop$yc<-10+d*pop$X1+d*pop$X2+rnorm(N,0,d*(sqrt(1-fraction))) #we get 90% of the variance
  
  for (x in 1:k){
    pop$yc[pop$dom==x]<-pop$yc[pop$dom==x]+rnorm(1,0,d*sqrt(fraction)) #every unit inside domain x gets same contribution. 
    #Units will b more similar
  }
  
  
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
    wkc[is.nan(wkc)]<-(n.npk)/((n.npk-1)*(1-1/n.pk)+(1+n.npk/n.pk)) #this or 1
    
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
              bias.nps,bias.comb,true,ccc))
  
}





######## Function for categories of unequal sizes in the approach of bias=a number ########

unequal.size.2<-function(k,c,n.pk,n.npk,corr){
  seed<-as.integer(paste0(k*3,n.pk,0,9,c*4))
  set.seed(seed)
  pop<-mvrnorm(N,c(0,0,0),matrix(c(1,0,corr,0,1,corr,corr,corr,1),nrow=3,ncol=3)) 
  d<-sqrt(80/3)
  pop<-mvrnorm(N,c(0,0,0),matrix(c(1,0,corr,0,1,corr,corr,corr,1),nrow=3,ncol=3)) 
  d<-sqrt(80/3)
  fraction<-0.5 #assign in so 10% comes from domain specific and the rest still random (add to function)
  pop<-as.data.frame(pop)
  names(pop)<-c("X1","X2","W")
  pop$dom<-as.factor(sample(1:k,N,replace=T)) #   
  pop$yc<-10+d*pop$X1+d*pop$X2+rnorm(N,0,d*(sqrt(1-fraction))) #we get 90% of the variance
  
  for (x in 1:k){
    pop$yc[pop$dom==x]<-pop$yc[pop$dom==x]+rnorm(1,0,d*sqrt(fraction)) #every unit inside domain x gets same contribution. 
    #Units will b more similar
  }
  
  
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
    wkc[is.nan(wkc)]<-(n.npk)/((n.npk-1)*(1-1/n.pk)+(1+n.npk/n.pk)) #this or 1
    
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
              bias.nps,bias.comb,true,ccc))
}


#######################################################################
############## RUN SIMULATIONS AND STORAGE RESULTS ####################
#######################################################################

#Different grids to add columns with armse later on
cond.eq.rnorm<-expand.grid(k=c(1,4,10,15),
                           c=c(3,5,8,15),
                           n.pk=c(10,100,400,900),
                           n.npk=c(100,1000,2000,6000),
                           corr=c(0.228,0.632))

cond.uneq.rnorm<-expand.grid(k=c(1,4,10,15),
                             c=c(3,5,8,15),
                             n.pk=c(10,100,400,900),
                             n.npk=c(100,1000,2000,6000),
                             corr=c(0.228,0.632))


####### For loop function equal size - MODEL A #######

for (t in 1:nrow(cond.eq)){
  x<-equal.size.1(cond.eq$k[t],cond.eq$c[t],cond.eq$n.pk[t],cond.eq$n.npk[t],cond.eq$corr[t])
  saveRDS(x[[1]],file=paste0("est_",t,"_eq1_rnorm.rds"))
  saveRDS(x[[2]],file=paste0("mse_ps_",t,"_eq1_rnorm.rds"))
  saveRDS(x[[3]],file=paste0("mse_nps_",t,"_eq1_rnorm.rds"))
  saveRDS(x[[4]],file=paste0("mse_comb_",t,"_eq1_rnorm.rds"))
  saveRDS(x[[8]],file=paste0("armse_ps_",t,"_eq1_rnorm.rds"))
  saveRDS(x[[9]],file=paste0("armse_nps_",t,"_eq1_rnorm.rds"))
  saveRDS(x[[10]],file=paste0("armse_comb_",t,"_eq1_rnorm.rds"))
  saveRDS(x[[11]],file=paste0("bias_nps_",t,"_eq1_rnorm.rds"))
  saveRDS(x[[12]],file=paste0("bias_comb_",t,"_eq1_rnorm.rds"))  
}

####### For loop function equal size - MODEL B #######

for (t in 1:nrow(cond.eq)){
  x<-equal.size.2(cond.eq$k[t],cond.eq$c[t],cond.eq$n.pk[t],cond.eq$n.npk[t],cond.eq$corr[t])
  saveRDS(x[[1]],file=paste0("est_",t,"_eq2_rnorm.rds"))
  saveRDS(x[[2]],file=paste0("mse_ps_",t,"_eq2_rnorm.rds"))
  saveRDS(x[[3]],file=paste0("mse_nps_",t,"_eq2_rnorm.rds"))
  saveRDS(x[[4]],file=paste0("mse_comb_",t,"_eq2_rnorm.rds"))
  saveRDS(x[[8]],file=paste0("armse_ps_",t,"_eq2_rnorm.rds"))
  saveRDS(x[[9]],file=paste0("armse_nps_",t,"_eq2_rnorm.rds"))
  saveRDS(x[[10]],file=paste0("armse_comb_",t,"_eq2_rnorm.rds"))
  saveRDS(x[[11]],file=paste0("bias_nps_",t,"_eq2_rnorm.rds"))
  saveRDS(x[[12]],file=paste0("bias_comb_",t,"_eq2_rnorm.rds"))  
}

####### For loop function unequal size - MODEL A #######

for (t in 1:nrow(cond.uneq)){
  x<-unequal.size.1(cond.uneq$k[t],cond.uneq$c[t],cond.uneq$n.pk[t],cond.uneq$n.npk[t],cond.uneq$corr[t])
  saveRDS(x[[1]],file=paste0("est_",t,"_uneq1_rnorm.rds"))
  saveRDS(x[[2]],file=paste0("mse_ps_",t,"_uneq1_rnorm.rds"))
  saveRDS(x[[3]],file=paste0("mse_nps_",t,"_uneq1_rnorm.rds"))
  saveRDS(x[[4]],file=paste0("mse_comb_",t,"_uneq1_rnorm.rds"))
  saveRDS(x[[8]],file=paste0("armse_ps_",t,"_uneq1_rnorm.rds"))
  saveRDS(x[[9]],file=paste0("armse_nps_",t,"_uneq1_rnorm.rds"))
  saveRDS(x[[10]],file=paste0("armse_comb_",t,"_uneq1_rnorm.rds"))
  saveRDS(x[[11]],file=paste0("bias_nps_",t,"_uneq1_rnorm.rds"))
  saveRDS(x[[12]],file=paste0("bias_comb_",t,"_uneq1_rnorm.rds"))  
}

####### For loop function unequal size - MODEL B #######

for (t in 1:nrow(cond.uneq)){
  x<-unequal.size.2(cond.uneq$k[t],cond.uneq$c[t],cond.uneq$n.pk[t],cond.uneq$n.npk[t],cond.uneq$corr[t])
  saveRDS(x[[1]],file=paste0("est_",t,"_uneq2_rnorm.rds"))
  saveRDS(x[[2]],file=paste0("mse_ps_",t,"_uneq2_rnorm.rds"))
  saveRDS(x[[3]],file=paste0("mse_nps_",t,"_uneq2_rnorm.rds"))
  saveRDS(x[[4]],file=paste0("mse_comb_",t,"_uneq2_rnorm.rds"))
  saveRDS(x[[8]],file=paste0("armse_ps_",t,"_uneq2_rnorm.rds"))
  saveRDS(x[[9]],file=paste0("armse_nps_",t,"_uneq2_rnorm.rds"))
  saveRDS(x[[10]],file=paste0("armse_comb_",t,"_uneq2_rnorm.rds"))
  saveRDS(x[[11]],file=paste0("bias_nps_",t,"_uneq2_rnorm.rds"))
  saveRDS(x[[12]],file=paste0("bias_comb_",t,"_uneq_rnorm.rds"))  
}


#######################################################################
########################## COMPUTING MARMSE ###########################
#######################################################################

#Different grids to add columns with armse later on
cond.eq.rnorm<-expand.grid(k=c(1,4,10,15),
                           c=c(3,5,8,15),
                           n.pk=c(10,100,400,900),
                           n.npk=c(100,1000,2000,6000),
                           corr=c(0.228,0.632))

cond.uneq.rnorm<-expand.grid(k=c(1,4,10,15),
                             c=c(3,5,8,15),
                             n.pk=c(10,100,400,900),
                             n.npk=c(100,1000,2000,6000),
                             corr=c(0.228,0.632))
####### Storage the average armse in table of 512 conditions #########
for (i in 1:512){
  X1<-readRDS(file=paste0("armse_ps_",i,"_eq1_rnorm.rds"))
  X2<-readRDS(file=paste0("armse_nps_",i,"_eq1_rnorm.rds"))
  X3<-readRDS(file=paste0("armse_comb_",i,"_eq1_rnorm.rds"))
  cond.eq.rnorm$armse.ps1[i]<-mean(X1)
  cond.eq.rnorm$armse.nps1[i]<-mean(X2)
  cond.eq.rnorm$armse.comb1[i]<-mean(X3)
  
  X4<-readRDS(file=paste0("armse_ps_",i,"_eq2_rnorm.rds"))
  X5<-readRDS(file=paste0("armse_nps_",i,"_eq2_rnorm.rds"))
  X6<-readRDS(file=paste0("armse_comb_",i,"_eq2_rnorm.rds"))
  cond.eq.rnorm$armse.ps2[i]<-mean(X4)
  cond.eq.rnorm$armse.nps2[i]<-mean(X5)
  cond.eq.rnorm$armse.comb2[i]<-mean(X6)
}

for (i in 1:512){
  X1<-readRDS(file=paste0("armse_ps_",i,"_uneq1_rnorm.rds"))
  X2<-readRDS(file=paste0("armse_nps_",i,"_uneq1_rnorm.rds"))
  X3<-readRDS(file=paste0("armse_comb_",i,"_uneq1_rnorm.rds"))
  cond.uneq.rnorm$armse.ps1[i]<-mean(X1)
  cond.uneq.rnorm$armse.nps1[i]<-mean(X2)
  cond.uneq.rnorm$armse.comb1[i]<-mean(X3)
  
  X4<-readRDS(file=paste0("armse_ps_",i,"_uneq2_rnorm.rds"))
  X5<-readRDS(file=paste0("armse_nps_",i,"_uneq2_rnorm.rds"))
  X6<-readRDS(file=paste0("armse_comb_",i,"_uneq2_rnorm.rds"))
  cond.uneq.rnorm$armse.ps2[i]<-mean(X4)
  cond.uneq.rnorm$armse.nps2[i]<-mean(X5)
  cond.uneq.rnorm$armse.comb2[i]<-mean(X6)
}

saveRDS(cond.eq.rnorm,file="cond.eq.norm.rds")
saveRDS(cond.uneq.rnorm,file="cond.uneq.norm.rds")

################################################
### SPLIT COND_EQ AND COND_UNEQ PER APPROACH ###
################################################

# Open cond.eq.norm.rds and cond.uneq.norm.rds

cond_eq_1<-as.data.frame(cbind(cond.eq.rnorm$k,cond.eq.rnorm$c,cond.eq.rnorm$n.pk,cond.eq.rnorm$n.npk,cond.eq.rnorm$corr,
                               cond.eq.rnorm$armse.ps1,cond.eq.rnorm$armse.nps1,cond.eq.rnorm$armse.comb1))


cond_eq_2<-as.data.frame(cbind(cond.eq.rnorm$k,cond.eq.rnorm$c,cond.eq.rnorm$n.pk,cond.eq.rnorm$n.npk,cond.eq.rnorm$corr,
                               cond.eq.rnorm$armse.ps2,cond.eq.rnorm$armse.nps2,cond.eq.rnorm$armse.comb2))


cond_uneq_1<-as.data.frame(cbind(cond.uneq.rnorm$k,cond.uneq.rnorm$c,cond.uneq.rnorm$n.pk,cond.uneq.rnorm$n.npk,cond.uneq.rnorm$corr,
                                 cond.uneq.rnorm$armse.ps1,cond.uneq.rnorm$armse.nps1,cond.uneq.rnorm$armse.comb1))


cond_uneq_2<-as.data.frame(cbind(cond.uneq.rnorm$k,cond.uneq.rnorm$c,cond.uneq.rnorm$n.pk,cond.uneq.rnorm$n.npk,cond.uneq.rnorm$corr,
                                 cond.uneq.rnorm$armse.ps2,cond.uneq.rnorm$armse.nps2,cond.uneq.rnorm$armse.comb2))

colnames(cond_eq_1)<-c("k","c","n.pk","n.npk","corr","armse.ps","armse.nps","armse.comb")
colnames(cond_eq_2)<-c("k","c","n.pk","n.npk","corr","armse.ps","armse.nps","armse.comb")
colnames(cond_uneq_1)<-c("k","c","n.pk","n.npk","corr","armse.ps","armse.nps","armse.comb")
colnames(cond_uneq_2)<-c("k","c","n.pk","n.npk","corr","armse.ps","armse.nps","armse.comb")



eq_1_weak<-cond_eq_1[which(cond_eq_1$corr==0.228),]
eq_1_sev<-cond_eq_1[which(cond_eq_1$corr==0.632),]
eq_2_weak<-cond_eq_2[which(cond_eq_2$corr==0.228),]
eq_2_sev<-cond_eq_2[which(cond_eq_2$corr==0.632),]
uneq_1_weak<-cond_uneq_1[which(cond_uneq_1$corr==0.228),]
uneq_1_sev<-cond_uneq_1[which(cond_uneq_1$corr==0.632),]
uneq_2_weak<-cond_uneq_2[which(cond_uneq_2$corr==0.228),]
uneq_2_sev<-cond_uneq_2[which(cond_uneq_2$corr==0.632),]


#############################################################
### CREATING FUNCTION THAT STORAGES IF COMBINED ESTIMATOR ###
#######   IS LOWER AS A DUMMY VARIABLE YES=1 NO=0     #######
#############################################################

less.than.function<-function(data) {
  
  data$less_ps<-NA
  data$less_nps<-NA
  data$less_both<-NA
  
  for (i in 1:nrow(data)) {
    
    if (data[i,8] < data[i,6]){
      data[i,9]<-1
    } else {data[i,9]<-0}
    
    if (data[i,8] < data[i,7]){
      data[i,10]<-1
    } else {data[i,10]<-0}
    
    if (data[i,8] < data[i,6] & data[i,8] < data[i,7] ){
      data[i,11]<-1
    } else {data[i,11]<-0}}
  
  return(data)
  
}

eq_1_weak<-less.than.function(eq_1_weak)
eq_1_sev<-less.than.function(eq_1_sev)
eq_2_weak<-less.than.function(eq_2_weak)
eq_2_sev<-less.than.function(eq_2_sev)
uneq_1_weak<-less.than.function(uneq_1_weak)
uneq_1_sev<-less.than.function(uneq_1_sev)
uneq_2_weak<-less.than.function(uneq_2_weak)
uneq_2_sev<-less.than.function(uneq_2_weak)


##########################################################################
#########  COMPUTE PROPORTION OF TIMES COMBINED HAS LOWER MSE    ######### 
##########################################################################
sum.less<-matrix(0,8,3)
sum.less[1,1]<-sum(eq_1_weak$less_ps)
sum.less[1,2]<-sum(eq_1_weak$less_nps)
sum.less[1,3]<-sum(eq_1_weak$less_both)

sum.less[2,1]<-sum(eq_2_weak$less_ps)
sum.less[2,2]<-sum(eq_2_weak$less_nps)
sum.less[2,3]<-sum(eq_2_weak$less_both)

sum.less[3,1]<-sum(uneq_1_weak$less_ps)
sum.less[3,2]<-sum(uneq_1_weak$less_nps)
sum.less[3,3]<-sum(uneq_1_weak$less_both)

sum.less[4,1]<-sum(uneq_2_weak$less_ps)
sum.less[4,2]<-sum(uneq_2_weak$less_nps)
sum.less[4,3]<-sum(uneq_2_weak$less_both)

sum.less[5,1]<-sum(eq_1_sev$less_ps)
sum.less[5,2]<-sum(eq_1_sev$less_nps)
sum.less[5,3]<-sum(eq_1_sev$less_both)

sum.less[6,1]<-sum(eq_2_sev$less_ps)
sum.less[6,2]<-sum(eq_2_sev$less_nps)
sum.less[6,3]<-sum(eq_2_sev$less_both)

sum.less[7,1]<-sum(uneq_1_sev$less_ps)
sum.less[7,2]<-sum(uneq_1_sev$less_nps)
sum.less[7,3]<-sum(uneq_1_sev$less_both)

sum.less[8,1]<-sum(uneq_2_sev$less_ps)
sum.less[8,2]<-sum(uneq_2_sev$less_nps)
sum.less[8,3]<-sum(uneq_2_sev$less_both)


colnames(sum.less)<-c("less than ps","less than nps","less than both")
rownames(sum.less)<-c("eq_1_weak","eq_2_weak","uneq_1_weak","uneq_2_weak",
                      "eq_1_sev","eq_2_sev","uneq_1_sev","uneq_2_sev")

library("xlsx")
sum.less.1<-(sum.less*100/256)/100
sum.less.1<-round(sum.less.1,digits=2)
write.xlsx(sum.less.1, file = "sum_rnorm_table.xlsx")

