library(xlsx)

#open bias_eq_1.rds, bias_eq_2.rds, bias_uneq_1,rds, bias_uneq_2.rds

saveRDS(bias.eq.1,file=paste0("bias_eq_1.rds"))  
saveRDS(bias.eq.2,file=paste0("bias_eq_2.rds")) 
saveRDS(bias.uneq.1,file=paste0("bias_uneq_1.rds")) 
saveRDS(bias.uneq.2,file=paste0("bias_uneq_2.rds"))

##################################################
############## SPLIT PER SELECTIVITY #############
##################################################
bias_eq_1_weak<-bias_eq_1[which(bias_eq_1$corr==0.228),]
bias_eq_1_sev<-bias_eq_1[which(bias_eq_1$corr==0.632),]
bias_eq_2_weak<-bias_eq_2[which(bias_eq_2$corr==0.228),]
bias_eq_2_sev<-bias_eq_2[which(bias_eq_2$corr==0.632),]
bias_uneq_1_weak<-bias_uneq_1[which(bias_uneq_1$corr==0.228),]
bias_uneq_1_sev<-bias_uneq_1[which(bias_uneq_1$corr==0.632),]
bias_uneq_2_weak<-bias_uneq_2[which(bias_uneq_2$corr==0.228),]
bias_uneq_2_sev<-bias_uneq_2[which(bias_uneq_2$corr==0.632),]

##########################################################
#########COMPUTE DIFFERENCE OF BIAS NPS-COMB #############
##########################################################

bias_eq_1_weak$diff<-(bias_eq_1_weak$bias.np)-(bias_eq_1_weak$bias.c)
bias_eq_1_sev$diff<-(bias_eq_1_sev$bias.np)-(bias_eq_1_sev$bias.c)
bias_eq_2_weak$diff<-(bias_eq_2_weak$bias.np)-(bias_eq_2_weak$bias.c)
bias_eq_2_sev$diff<-(bias_eq_2_sev$bias.np)-(bias_eq_2_sev$bias.c)
bias_uneq_1_weak$diff<-(bias_uneq_1_weak$bias.np)-(bias_uneq_1_weak$bias.c)
bias_uneq_1_sev$diff<-(bias_uneq_1_sev$bias.np)-(bias_uneq_1_sev$bias.c)
bias_uneq_2_weak$diff<-(bias_uneq_2_weak$bias.np)-(bias_uneq_2_weak$bias.c)
bias_uneq_2_sev$diff<-(bias_uneq_2_sev$bias.np)-(bias_uneq_2_sev$bias.c)

#negative= np is better
##################################################################### 
################## FINAL TABLE DIFFERENCE OF BIAS ###################  
##################################################################### 

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

bias.diff.table<-function(data,rows){
  
  bias<-matrix(0,16,16)
  for (i in 1:nrow(rows)) {
    
    bias[i,1]<-data[(rows[i,1]),8]
    bias[i,2]<-data[(rows[i,2]),8]
    bias[i,3]<-data[(rows[i,3]),8]
    bias[i,4]<-data[(rows[i,4]),8]
    bias[i,5]<-data[(rows[i,5]),8]
    bias[i,6]<-data[(rows[i,6]),8]
    bias[i,7]<-data[(rows[i,7]),8]
    bias[i,8]<-data[(rows[i,8]),8]
    bias[i,9]<-data[(rows[i,9]),8]
    bias[i,10]<-data[(rows[i,10]),8]
    bias[i,11]<-data[(rows[i,11]),8]
    bias[i,12]<-data[(rows[i,12]),8]
    bias[i,13]<-data[(rows[i,13]),8]
    bias[i,14]<-data[(rows[i,14]),8]
    bias[i,15]<-data[(rows[i,15]),8]
    bias[i,16]<-data[(rows[i,16]),8]
  }
  
  return(bias)
  
}

b_diff_weak_eq1<-round(bias.diff.table(bias_eq_1_weak,row_it),digits=2)
b_diff_sev_eq1<-round(bias.diff.table(bias_eq_1_sev,row_it),digits=2)
b_diff_weak_eq2<-round(bias.diff.table(bias_eq_2_weak,row_it),digits=2)
b_diff_sev_eq2<-round(bias.diff.table(bias_eq_2_sev,row_it),digits=2)
b_diff_weak_un1<-round(bias.diff.table(bias_uneq_1_weak,row_it),digits=2)
b_diff_sev_un1<-round(bias.diff.table(bias_uneq_1_sev,row_it),digits=2)
b_diff_weak_un2<-round(bias.diff.table(bias_uneq_2_weak,row_it),digits=2)
b_diff_sev_un2<-round(bias.diff.table(bias_uneq_2_sev,row_it),digits=2)

write.xlsx(b_diff_weak_eq1, file = "b_diff_weak_eq1.xlsx")
write.xlsx(b_diff_sev_eq1, file = "b_diff_sev_eq1.xlsx")
write.xlsx(b_diff_weak_eq2, file = "b_diff_weak_eq2.xlsx")
write.xlsx(b_diff_sev_eq2, file = "b_diff_sev_eq2.xlsx")
write.xlsx(b_diff_weak_un1, file = "b_diff_weak_un1.xlsx")
write.xlsx(b_diff_sev_un1, file = "b_diff_sev_un1.xlsx")
write.xlsx(b_diff_weak_un2, file = "b_diff_weak_un2.xlsx")
write.xlsx(b_diff_sev_un2, file = "b_diff_sev_un2.xlsx")


###################################################################
########COMPUTE PROPORTION OF COMBINED WITH LOWER BIAS#############  
###################################################################

prop.bias<-function (data){
  
  data$prop<-NA
  
  for (i in 1:256){
    
if (data[i,8]>=0){
  data[i,9]<-1
} else {data[i,9]<-0}
  }
  
  return(data)
}

bias_eq_1_weak<-prop.bias(bias_eq_1_weak)
bias_eq_1_sev<-prop.bias(bias_eq_1_sev)
bias_eq_2_weak<-prop.bias(bias_eq_2_weak)
bias_eq_2_sev<-prop.bias(bias_eq_2_sev)
bias_uneq_1_weak<-prop.bias(bias_uneq_1_weak)
bias_uneq_1_sev<-prop.bias(bias_uneq_1_sev)
bias_uneq_2_weak<-prop.bias(bias_uneq_2_weak)
bias_uneq_2_sev<-prop.bias(bias_uneq_2_sev)

x1<-table(bias_eq_1_weak$prop)
x2<-table(bias_eq_1_sev$prop)
x3<-table(bias_eq_2_weak$prop)
x4<-table(bias_eq_2_sev$prop)
x5<-table(bias_uneq_1_weak$prop)
x6<-table(bias_uneq_1_sev$prop)
x7<-table(bias_uneq_2_weak$prop)
x8<-table(bias_uneq_2_sev$prop)

proportion_bias<-rbind(x1,x2,x3,x4,x5,x6,x7,x8)
write.xlsx(proportion_bias,file="proportion_bias.xlsx")

##################################################################### 
##################### BIAS OF COMBINED ESTIMATOR ####################  
##################################################################### 

bias.comb.table<-function(data,rows){
  
  bias<-matrix(0,16,16)
  for (i in 1:nrow(rows)) {
    
    bias[i,1]<-data[(rows[i,1]),7]
    bias[i,2]<-data[(rows[i,2]),7]
    bias[i,3]<-data[(rows[i,3]),7]
    bias[i,4]<-data[(rows[i,4]),7]
    bias[i,5]<-data[(rows[i,5]),7]
    bias[i,6]<-data[(rows[i,6]),7]
    bias[i,7]<-data[(rows[i,7]),7]
    bias[i,8]<-data[(rows[i,8]),7]
    bias[i,9]<-data[(rows[i,9]),7]
    bias[i,10]<-data[(rows[i,10]),7]
    bias[i,11]<-data[(rows[i,11]),7]
    bias[i,12]<-data[(rows[i,12]),7]
    bias[i,13]<-data[(rows[i,13]),7]
    bias[i,14]<-data[(rows[i,14]),7]
    bias[i,15]<-data[(rows[i,15]),7]
    bias[i,16]<-data[(rows[i,16]),7]
  }
  
  return(bias)
  
}

b_comb_weak_eq1<-round(bias.comb.table(bias_eq_1_weak,row_it),digits=2)
b_comb_sev_eq1<-round(bias.comb.table(bias_eq_1_sev,row_it),digits=2)
b_comb_weak_eq2<-round(bias.comb.table(bias_eq_2_weak,row_it),digits=2)
b_comb_sev_eq2<-round(bias.comb.table(bias_eq_2_sev,row_it),digits=2)
b_comb_weak_un1<-round(bias.comb.table(bias_uneq_1_weak,row_it),digits=2)
b_comb_sev_un1<-round(bias.comb.table(bias_uneq_1_sev,row_it),digits=2)
b_comb_weak_un2<-round(bias.comb.table(bias_uneq_2_weak,row_it),digits=2)
b_comb_sev_un2<-round(bias.comb.table(bias_uneq_2_weak,row_it),digits=2)

write.xlsx(b_comb_weak_eq1, file = "b_comb_weak_eq1.xlsx")
write.xlsx(b_comb_sev_eq1, file = "b_comb_sev_eq1.xlsx")
write.xlsx(b_comb_weak_eq2, file = "b_comb_weak_eq2.xlsx")
write.xlsx(b_comb_sev_eq2, file = "b_comb_sev_eq2.xlsx")
write.xlsx(b_comb_weak_un1, file = "b_comb_weak_un1.xlsx")
write.xlsx(b_comb_sev_un1, file = "b_comb_sev_un1.xlsx")
write.xlsx(b_comb_weak_un2, file = "b_comb_weak_un2.xlsx")
write.xlsx(b_comb_sev_un2, file = "b_comb_sev_un2.xlsx")