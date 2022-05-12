
################ Create grids to run the functions ################
#When size of categories is equal
cond.eq<-expand.grid(k=c(1,4,10,15),
                  c=c(3,5,8,15),
                  n.pk=c(10,100,400,900),
                  n.npk=c(100,1000,2000,6000),
                  corr=c(0.228,0.632))

#When size of categories is unequal
cond.uneq<-expand.grid(k=c(1,4,10,15),
                     c=c(3,5,8,15),
                     n.pk=c(10,100,400,900),
                     n.npk=c(100,1000,2000,6000),
                     corr=c(0.228,0.632))

####### Run simulation and storage = FUNCTION FOR CATEGORIES OF EQUAL SIZE - MODEL A #######

for (t in 1:nrow(cond.eq)){
  x<-equal.size.1(cond.eq$k[t],cond.eq$c[t],cond.eq$n.pk[t],cond.eq$n.npk[t],cond.eq$corr[t])
  saveRDS(x[[1]],file=paste0("est_",t,"_eq1.rds"))
  saveRDS(x[[2]],file=paste0("mse_ps_",t,"_eq1.rds"))
  saveRDS(x[[3]],file=paste0("mse_nps_",t,"_eq1.rds"))
  saveRDS(x[[4]],file=paste0("mse_comb_",t,"_eq1.rds"))
  saveRDS(x[[8]],file=paste0("armse_ps_",t,"_eq1.rds"))
  saveRDS(x[[9]],file=paste0("armse_nps_",t,"_eq1.rds"))
  saveRDS(x[[10]],file=paste0("armse_comb_",t,"_eq1.rds"))
  saveRDS(x[[11]],file=paste0("bias_nps_",t,"_eq1.rds"))
  saveRDS(x[[12]],file=paste0("bias_comb_",t,"_eq1.rds"))  
}

#######  Run simulation  and storage = FUNCTION FOR CATEGORIES OF EQUAL SIZE - MODEL B #######

for (t in 1:nrow(cond.eq)){
  x<-equal.size.2(cond.eq$k[t],cond.eq$c[t],cond.eq$n.pk[t],cond.eq$n.npk[t],cond.eq$corr[t])
  saveRDS(x[[1]],file=paste0("est_",t,"_eq2.rds"))
  saveRDS(x[[2]],file=paste0("mse_ps_",t,"_eq2.rds"))
  saveRDS(x[[3]],file=paste0("mse_nps_",t,"_eq2.rds"))
  saveRDS(x[[4]],file=paste0("mse_comb_",t,"_eq2.rds"))
  saveRDS(x[[8]],file=paste0("armse_ps_",t,"_eq2.rds"))
  saveRDS(x[[9]],file=paste0("armse_nps_",t,"_eq2.rds"))
  saveRDS(x[[10]],file=paste0("armse_comb_",t,"_eq2.rds"))
  saveRDS(x[[11]],file=paste0("bias_nps_",t,"_eq2.rds"))
  saveRDS(x[[12]],file=paste0("bias_comb_",t,"_eq2.rds"))  
}

####### Run simulation  and storage = FUNCTION FOR CATEGORIES OF UNEQUAL SIZE - MODEL A #######

for (t in 1:nrow(cond.uneq)){
  x<-unequal.size.1(cond.uneq$k[t],cond.uneq$c[t],cond.uneq$n.pk[t],cond.uneq$n.npk[t],cond.uneq$corr[t])
  saveRDS(x[[1]],file=paste0("est_",t,"_uneq1.rds"))
  saveRDS(x[[2]],file=paste0("mse_ps_",t,"_uneq1.rds"))
  saveRDS(x[[3]],file=paste0("mse_nps_",t,"_uneq1.rds"))
  saveRDS(x[[4]],file=paste0("mse_comb_",t,"_uneq1.rds"))
  saveRDS(x[[8]],file=paste0("armse_ps_",t,"_uneq1.rds"))
  saveRDS(x[[9]],file=paste0("armse_nps_",t,"_uneq1.rds"))
  saveRDS(x[[10]],file=paste0("armse_comb_",t,"_uneq1.rds"))
  saveRDS(x[[11]],file=paste0("bias_nps_",t,"_uneq1.rds"))
  saveRDS(x[[12]],file=paste0("bias_comb_",t,"_uneq1.rds"))  
}

####### Run simulation  and storage = FUNCTION FOR CATEGORIES OF UNEQUAL SIZE - MODEL B #######

for (t in 1:nrow(cond.uneq)){
  x<-unequal.size.2(cond.uneq$k[t],cond.uneq$c[t],cond.uneq$n.pk[t],cond.uneq$n.npk[t],cond.uneq$corr[t])
  saveRDS(x[[1]],file=paste0("est_",t,"_uneq2.rds"))
  saveRDS(x[[2]],file=paste0("mse_ps_",t,"_uneq2.rds"))
  saveRDS(x[[3]],file=paste0("mse_nps_",t,"_uneq2.rds"))
  saveRDS(x[[4]],file=paste0("mse_comb_",t,"_uneq2.rds"))
  saveRDS(x[[8]],file=paste0("armse_ps_",t,"_uneq2.rds"))
  saveRDS(x[[9]],file=paste0("armse_nps_",t,"_uneq2.rds"))
  saveRDS(x[[10]],file=paste0("armse_comb_",t,"_uneq2.rds"))
  saveRDS(x[[11]],file=paste0("bias_nps_",t,"_uneq2.rds"))
  saveRDS(x[[12]],file=paste0("bias_comb_",t,"_uneq.rds"))  
}


####### Read ARMSE of each simulation and compute the mean per domain ###################
# For equal size categories
for (i in 1:512){
  X1<-readRDS(file=paste0("armse_ps_",i,"_eq1.rds"))
  X2<-readRDS(file=paste0("armse_nps_",i,"_eq1.rds"))
  X3<-readRDS(file=paste0("armse_comb_",i,"_eq1.rds"))
  cond.eq$armse.ps1[i]<-mean(X1)
  cond.eq$armse.nps1[i]<-mean(X2)
  cond.eq$armse.comb1[i]<-mean(X3)
  
  X4<-readRDS(file=paste0("armse_ps_",i,"_eq2.rds"))
  X5<-readRDS(file=paste0("armse_nps_",i,"_eq2.rds"))
  X6<-readRDS(file=paste0("armse_comb_",i,"_eq2.rds"))
  cond.eq$armse.ps2[i]<-mean(X4)
  cond.eq$armse.nps2[i]<-mean(X5)
  cond.eq$armse.comb2[i]<-mean(X6)
}

#For unequal size categories
for (i in 1:512){
  X1<-readRDS(file=paste0("armse_ps_",i,"_uneq1.rds"))
  X2<-readRDS(file=paste0("armse_nps_",i,"_uneq1.rds"))
  X3<-readRDS(file=paste0("armse_comb_",i,"_uneq1.rds"))
  cond.uneq$armse.ps1[i]<-mean(X1)
  cond.uneq$armse.nps1[i]<-mean(X2)
  cond.uneq$armse.comb1[i]<-mean(X3)
  
  X4<-readRDS(file=paste0("armse_ps_",i,"_uneq2.rds"))
  X5<-readRDS(file=paste0("armse_nps_",i,"_uneq2.rds"))
  X6<-readRDS(file=paste0("armse_comb_",i,"_uneq2.rds"))
  cond.uneq$armse.ps2[i]<-mean(X4)
  cond.uneq$armse.nps2[i]<-mean(X5)
  cond.uneq$armse.comb2[i]<-mean(X6)
}

