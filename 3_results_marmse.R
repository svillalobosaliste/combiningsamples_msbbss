library(xlsx)

#open cond_eq.rds and cond_uneq.rds

################################################
### SPLIT COND_EQ AND COND_UNEQ PER APPROACH ###
################################################
cond_eq_1<-as.data.frame(cbind(cond_equal$k,cond_equal$c,cond_equal$n.pk,cond_equal$n.npk,cond_equal$corr,
                               cond_equal$armse.ps1,cond_equal$armse.nps1,cond_equal$armse.comb1))


cond_eq_2<-as.data.frame(cbind(cond_equal$k,cond_equal$c,cond_equal$n.pk,cond_equal$n.npk,cond_equal$corr,
                               cond_equal$armse.ps2,cond_equal$armse.nps2,cond_equal$armse.comb2))


cond_uneq_1<-as.data.frame(cbind(cond_uneq$k,cond_uneq$c,cond_uneq$n.pk,cond_uneq$n.npk,cond_uneq$corr,
                                 cond_uneq$armse.ps1,cond_uneq$armse.nps1,cond_uneq$armse.comb1))


cond_uneq_2<-as.data.frame(cbind(cond_uneq$k,cond_uneq$c,cond_uneq$n.pk,cond_uneq$n.npk,cond_uneq$corr,
                                 cond_uneq$armse.ps2,cond_uneq$armse.nps2,cond_uneq$armse.comb2))

colnames(cond_eq_1)<-c("k","c","n.pk","n.npk","corr","armse.ps","armse.nps","armse.comb")
colnames(cond_eq_2)<-c("k","c","n.pk","n.npk","corr","armse.ps","armse.nps","armse.comb")
colnames(cond_uneq_1)<-c("k","c","n.pk","n.npk","corr","armse.ps","armse.nps","armse.comb")
colnames(cond_uneq_2)<-c("k","c","n.pk","n.npk","corr","armse.ps","armse.nps","armse.comb")

##################################################
### COUNTING HOW MANY TIMES COMBINES ESTIMATOR ###
#####   IS LOWER THAN:PS, NPS, AND BOTH      #####
##################################################

less_ps_eq1<-cond_eq_1[which(cond_eq_1$armse.comb < cond_eq_1$armse.ps),]       #282
less_ps_eq2<-cond_eq_1[which(cond_eq_2$armse.comb < cond_eq_2$armse.ps),]       #351
less_ps_uneq1<-cond_eq_1[which(cond_uneq_1$armse.comb < cond_uneq_1$armse.ps),] #302
less_ps_uneq2<-cond_eq_1[which(cond_uneq_2$armse.comb < cond_uneq_2$armse.ps),] #336

less_nps_eq1<-cond_eq_1[which(cond_eq_1$armse.com < cond_eq_1$armse.nps ),]        #384
less_nps_eq2<-cond_eq_1[which(cond_eq_2$armse.com < cond_eq_2$armse.nps ),]        #391
less_nps_uneq1<-cond_eq_1[which(cond_uneq_1$armse.com < cond_uneq_1$armse.nps ),]  #387  
less_nps_uneq2<-cond_eq_1[which(cond_uneq_2$armse.com < cond_uneq_2$armse.nps ),]  #390  

less_both_eq1<-cond_eq_1[which(cond_eq_1$armse.comb < cond_eq_1$armse.ps 
                               & cond_eq_1$armse.com < cond_eq_1$armse.nps ),]        #155
less_both_eq2<-cond_eq_1[which(cond_eq_2$armse.comb < cond_eq_2$armse.ps 
                               & cond_eq_2$armse.com < cond_eq_2$armse.nps ),]        #230
less_both_uneq1<-cond_eq_1[which(cond_uneq_1$armse.comb < cond_uneq_1$armse.ps 
                                 & cond_uneq_1$armse.com < cond_uneq_1$armse.nps ),]  #177
less_both_uneq2<-cond_eq_1[which(cond_uneq_2$armse.comb < cond_uneq_2$armse.ps 
                                 & cond_uneq_2$armse.com < cond_uneq_2$armse.nps ),]  #214


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

armse.eq1<-less.than.function(cond_eq_1)
armse.eq2<-less.than.function(cond_eq_2)
armse.uneq1<-less.than.function(cond_uneq_1)
armse.uneq2<-less.than.function(cond_uneq_2)

################################################
#########  SPLITING PER SELECTIVITY    ######### 
################################################

eq_1_weak<-round(armse.eq1[which(cond_eq_1$corr=="0.228"),],digits=2)
eq_1_sev<-round(armse.eq1[which(cond_eq_1$corr=="0.632"),],digits=2)
eq_2_weak<-round(armse.eq2[which(cond_eq_1$corr=="0.228"),],digits=2)
eq_2_sev<-round(armse.eq2[which(cond_eq_1$corr=="0.632"),],digits=2)
uneq_1_weak<-round(armse.uneq1[which(cond_uneq_1$corr=="0.228"),],digits=2)
uneq_1_sev<-round(armse.uneq1[which(cond_uneq_1$corr=="0.632"),],digits=2)
uneq_2_weak<-round(armse.uneq2[which(cond_uneq_1$corr=="0.228"),],digits=2)
uneq_2_sev<-round(armse.uneq2[which(cond_uneq_1$corr=="0.632"),],digits=2)

saveRDS(eq_1_weak,file="eq_1_weak.rds")
saveRDS(eq_1_sev,file="eq_1_sev.rds")
saveRDS(eq_2_weak,file="eq_2_weak.rds")
saveRDS(eq_2_sev,file="eq_2_sev.rds")
saveRDS(uneq_1_weak,file="uneq_1_weak.rds")
saveRDS(uneq_1_sev,file="uneq_1_sev.rds")
saveRDS(uneq_2_weak,file="uneq_2_weak.rds")
saveRDS(uneq_2_sev,file="uneq_2_sev.rds")

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


sum.less.1<-(sum.less*100/256)
sum.less.1<-round(sum.less.1,digits=2)
write.xlsx(sum.less.1, file = "summary_marmse.xlsx")

#################################################
##########  CREATE TABLES MARMSE (8)   ##########
#################################################
## This is the position number of vector armse in table marmse
row1<-c(1,17,33,49,2,18,34,50,3,19,35,51,4,20,36,52)
row2<-c(65,81,97,113,66,82,98,114,67,83,99,115,68,84,100,116)
row3<-c(129,145,161,177,130,146,162,178,131,147,163,179,132,148,164,180)
row4<-c(193,209,225,241,194,210,226,242,195,211,227,243,196,212,228,244)
row5<-c(5,21,37,53,6,22,38,54,7,23,39,55,8,24,40,56)
row6<-c(69,85,101,117,70,86,102,118,71,87,103,119,72,88,104,120)
row7<-c(133,149,165,181,134,150,166,182,135,151,167,183,136,152,168,184)
row8<-c(197,213,229,245,198,214,230,246,199,215,231,247,200,216,232,248)
row9<-c(9,25,41,57,10,26,42,58,11,27,43,59,12,28,44,60)
row10<-c(73,89,105,121,74,90,106,122,75,91,107,123,76,92,108,124)
row11<-c(137,153,169,185,138,154,170,186,139,155,171,187,140,156,172,188)
row12<-c(201,217,233,249,202,218,234,250,203,219,235,251,204,220,236,252)
row13<-c(13,29,45,61,14,30,46,62,15,31,47,63,16,32,48,64)
row14<-c(77,93,109,125,78,94,110,126,79,95,111,127,80,96,112,128)
row15<-c(141,157,173,189,142,158,174,190,143,159,175,191,144,160,176,192)
row16<-c(205,221,237,253,206,222,238,254,207,223,239,255,208,224,240,256)

row_it<-cbind(row1,row2,row3,row4,row5,row6,row7,row8,row9,row10,row11,row12,row13,row14,row15,row16)

marmse.table<-function(data,rows){
  
  marmse<-matrix(0,16,16)
  
  for (i in 1:nrow(rows)) {
    
    marmse[i,1]<-data[(rows[i,1]),8]
    marmse[i,2]<-data[(rows[i,2]),8]
    marmse[i,3]<-data[(rows[i,3]),8]
    marmse[i,4]<-data[(rows[i,4]),8]
    marmse[i,5]<-data[(rows[i,5]),8]
    marmse[i,6]<-data[(rows[i,6]),8]
    marmse[i,7]<-data[(rows[i,7]),8]
    marmse[i,8]<-data[(rows[i,8]),8]
    marmse[i,9]<-data[(rows[i,9]),8]
    marmse[i,10]<-data[(rows[i,10]),8]
    marmse[i,11]<-data[(rows[i,11]),8]
    marmse[i,12]<-data[(rows[i,12]),8]
    marmse[i,13]<-data[(rows[i,13]),8]
    marmse[i,14]<-data[(rows[i,14]),8]
    marmse[i,15]<-data[(rows[i,15]),8]
    marmse[i,16]<-data[(rows[i,16]),8]
  }
  
  return(marmse)
}

marmse_eq_1_weak<-marmse.table(eq_1_weak,row_it)
marmse_eq_2_weak<-marmse.table(eq_2_weak,row_it)
marmse_uneq_1_weak<-marmse.table(uneq_1_weak,row_it)
marmse_uneq_2_weak<-marmse.table(uneq_2_weak,row_it)

marmse_eq_1_sev<-marmse.table(eq_1_sev,row_it)
marmse_eq_2_sev<-marmse.table(eq_2_sev,row_it)
marmse_uneq_1_sev<-marmse.table(uneq_1_sev,row_it)
marmse_uneq_2_sev<-marmse.table(uneq_2_sev,row_it)


write.xlsx(marmse_eq_1_weak, file = "marmse_eq_1_weak.xlsx")
write.xlsx(marmse_eq_2_weak, file = "marmse_eq_2_weak.xlsx")
write.xlsx(marmse_uneq_1_weak, file = "marmse_uneq_1_weak.xlsx")
write.xlsx(marmse_uneq_2_weak, file = "marmse_uneq_2_weak.xlsx")
write.xlsx(marmse_eq_1_sev, file = "marmse_eq_1_sev.xlsx")
write.xlsx(marmse_eq_2_sev, file = "marmse_eq_2_sev.xlsx")
write.xlsx(marmse_uneq_1_sev, file = "marmse_uneq_1_sev.xlsx")
write.xlsx(marmse_uneq_2_sev, file = "marmse_uneq_2_sev.xlsx")

#################################################
########  CREATE GROUPS TABLE MARMSE (8) ####### 
#################################################

# This calculation indicates if the value of the MARMSE of the combined estimator is lower than the probability
#sample, lower than the non-probability sample, or lower than both.
# 2 = means that it is lower than the non-probability sample
# 3 = means that it is lower than the probability sample  
# 6 = means that it is lower than both 

eq_1_weak$color<-eq_1_weak$less_ps*2+eq_1_weak$less_nps*3+eq_1_weak$less_both*1 
eq_2_weak$color<-eq_2_weak$less_ps*2+eq_2_weak$less_nps*3+eq_2_weak$less_both*1
eq_1_sev$color<-eq_1_sev$less_ps*2+eq_1_sev$less_nps*3+eq_1_sev$less_both*1 
eq_2_sev$color<-eq_2_sev$less_ps*2+eq_2_sev$less_nps*3+eq_2_sev$less_both*1 
uneq_1_weak$color<-uneq_1_weak$less_ps*2+uneq_1_weak$less_nps*3+uneq_1_weak$less_both*1 
uneq_2_weak$color<-uneq_2_weak$less_ps*2+uneq_2_weak$less_nps*3+uneq_2_weak$less_both*1 
uneq_1_sev$color<-uneq_1_sev$less_ps*2+uneq_1_sev$less_nps*3+uneq_1_sev$less_both*1 
uneq_2_sev$color<-uneq_2_sev$less_ps*2+uneq_2_sev$less_nps*3+uneq_2_sev$less_both*1 


color.marmse.table<-function(data,rows){
  
  marmse<-matrix(0,16,16)
  
  for (i in 1:nrow(rows)) {
    
    marmse[i,1]<-data[(rows[i,1]),12]
    marmse[i,2]<-data[(rows[i,2]),12]
    marmse[i,3]<-data[(rows[i,3]),12]
    marmse[i,4]<-data[(rows[i,4]),12]
    marmse[i,5]<-data[(rows[i,5]),12]
    marmse[i,6]<-data[(rows[i,6]),12]
    marmse[i,7]<-data[(rows[i,7]),12]
    marmse[i,8]<-data[(rows[i,8]),12]
    marmse[i,9]<-data[(rows[i,9]),12]
    marmse[i,10]<-data[(rows[i,10]),12]
    marmse[i,11]<-data[(rows[i,11]),12]
    marmse[i,12]<-data[(rows[i,12]),12]
    marmse[i,13]<-data[(rows[i,13]),12]
    marmse[i,14]<-data[(rows[i,14]),12]
    marmse[i,15]<-data[(rows[i,15]),12]
    marmse[i,16]<-data[(rows[i,16]),12]
  }
  return(marmse)
}


e1w<-color.marmse.table(eq_1_weak,row_it)
e2w<-color.marmse.table(eq_2_weak,row_it)
un1w<-color.marmse.table(uneq_1_weak,row_it)
un2w<-color.marmse.table(uneq_2_weak,row_it)

e1s<-color.marmse.table(eq_1_sev,row_it)
e2s<-color.marmse.table(eq_2_sev,row_it)
un1s<-color.marmse.table(uneq_1_sev,row_it)
un2s<-color.marmse.table(uneq_2_sev,row_it)

write.xlsx(e1w, file = "eq1weak.xlsx")
write.xlsx(e2w, file = "eq2weak.xlsx")
write.xlsx(un1w, file = "uneq1weak.xlsx")
write.xlsx(un2w, file = "uneq2weak.xlsx")
write.xlsx(e1s, file = "eq1sev.xlsx")
write.xlsx(e2s, file = "eq2sev.xlsx")
write.xlsx(un1s, file = "uneq1sev.xlsx")
write.xlsx(un2s, file = "uneq2sev.xlsx")

#These tables are used in excel along with the previous tables to coloured the tables depending
#on the armse being lower. This resulting tables are the ones shown in the report.

###################
p.marmse.table<-function(data,rows){
  
  marmse<-matrix(0,16,16)
  
  for (i in 1:nrow(rows)) {
    
    marmse[i,1]<-data[(rows[i,1]),6]
    marmse[i,2]<-data[(rows[i,2]),6]
    marmse[i,3]<-data[(rows[i,3]),6]
    marmse[i,4]<-data[(rows[i,4]),6]
    marmse[i,5]<-data[(rows[i,5]),6]
    marmse[i,6]<-data[(rows[i,6]),6]
    marmse[i,7]<-data[(rows[i,7]),6]
    marmse[i,8]<-data[(rows[i,8]),6]
    marmse[i,9]<-data[(rows[i,9]),6]
    marmse[i,10]<-data[(rows[i,10]),6]
    marmse[i,11]<-data[(rows[i,11]),6]
    marmse[i,12]<-data[(rows[i,12]),6]
    marmse[i,13]<-data[(rows[i,13]),6]
    marmse[i,14]<-data[(rows[i,14]),6]
    marmse[i,15]<-data[(rows[i,15]),6]
    marmse[i,16]<-data[(rows[i,16]),6]
  }
  
  return(marmse)
}

np.marmse.table<-function(data,rows){
  
  marmse<-matrix(0,16,16)
  
  for (i in 1:nrow(rows)) {
    
    marmse[i,1]<-data[(rows[i,1]),7]
    marmse[i,2]<-data[(rows[i,2]),7]
    marmse[i,3]<-data[(rows[i,3]),7]
    marmse[i,4]<-data[(rows[i,4]),7]
    marmse[i,5]<-data[(rows[i,5]),7]
    marmse[i,6]<-data[(rows[i,6]),7]
    marmse[i,7]<-data[(rows[i,7]),7]
    marmse[i,8]<-data[(rows[i,8]),7]
    marmse[i,9]<-data[(rows[i,9]),7]
    marmse[i,10]<-data[(rows[i,10]),7]
    marmse[i,11]<-data[(rows[i,11]),7]
    marmse[i,12]<-data[(rows[i,12]),7]
    marmse[i,13]<-data[(rows[i,13]),7]
    marmse[i,14]<-data[(rows[i,14]),7]
    marmse[i,15]<-data[(rows[i,15]),7]
    marmse[i,16]<-data[(rows[i,16]),7]
  }
  
  return(marmse)
}

ps_eq_1_weak<-p.marmse.table(eq_1_weak,row_it)
ps_eq_2_weak<-p.marmse.table(eq_2_weak,row_it)
ps_uneq_1_weak<-p.marmse.table(uneq_1_weak,row_it)
ps_uneq_2_weak<-p.marmse.table(uneq_2_weak,row_it)

ps_eq_1_sev<-p.marmse.table(eq_1_sev,row_it)
ps_eq_2_sev<-p.marmse.table(eq_2_sev,row_it)
ps_uneq_1_sev<-p.marmse.table(uneq_1_sev,row_it)
ps_uneq_2_sev<-p.marmse.table(uneq_2_sev,row_it)

write.xlsx(ps_eq_1_weak, file = "ps_eq_1_weak.xlsx")
write.xlsx(ps_eq_2_weak, file = "ps_eq_2_weak.xlsx")
write.xlsx(ps_uneq_1_weak, file = "ps_uneq_1_weak.xlsx")
write.xlsx(ps_uneq_2_weak, file = "ps_uneq_2_weak.xlsx")
write.xlsx(ps_eq_1_sev, file = "ps_eq_1_sev.xlsx")
write.xlsx(ps_eq_2_sev, file = "ps_eq_2_sev.xlsx")
write.xlsx(ps_uneq_1_sev, file = "ps_uneq_1_sev.xlsx")
write.xlsx(ps_uneq_2_sev, file = "ps_uneq_2_sev.xlsx")

nps_eq_1_weak<-np.marmse.table(eq_1_weak,row_it)
nps_eq_2_weak<-np.marmse.table(eq_2_weak,row_it)
nps_uneq_1_weak<-np.marmse.table(uneq_1_weak,row_it)
nps_uneq_2_weak<-np.marmse.table(uneq_2_weak,row_it)

nps_eq_1_sev<-np.marmse.table(eq_1_sev,row_it)
nps_eq_2_sev<-np.marmse.table(eq_2_sev,row_it)
nps_uneq_1_sev<-np.marmse.table(uneq_1_sev,row_it)
nps_uneq_2_sev<-np.marmse.table(uneq_2_sev,row_it)

write.xlsx(nps_eq_1_weak, file = "nps_eq_1_weak.xlsx")
write.xlsx(nps_eq_2_weak, file = "nps_eq_2_weak.xlsx")
write.xlsx(nps_uneq_1_weak, file = "nps_uneq_1_weak.xlsx")
write.xlsx(nps_uneq_2_weak, file = "nps_uneq_2_weak.xlsx")
write.xlsx(nps_eq_1_sev, file = "nps_eq_1_sev.xlsx")
write.xlsx(nps_eq_2_sev, file = "nps_eq_2_sev.xlsx")
write.xlsx(nps_uneq_1_sev, file = "nps_uneq_1_sev.xlsx")
write.xlsx(nps_uneq_2_sev, file = "nps_uneq_2_sev.xlsx")
