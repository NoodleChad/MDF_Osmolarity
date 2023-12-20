#### Figure IV a)  ####
setwd("./Experimental/")
dat<-read.csv("Max_growth_output.csv",sep=",",dec=".",h=T)
dat1<-read.csv("./Results/analysis_II_fitted_kcl.csv",sep=",",dec=".",h=T)
dat2<-read.csv("./Results/analysis_II_fitted_nacl.csv",sep=",",dec=".",h=T)
dat1$k_val = dat1$k_val*1000
dat2$k_val = dat2$k_val*1000
pdf(file = "./Figures/Figure4a.pdf",
    width = 10, height = 6)
par(mar = c(6, 4.1, 4.1, 2.1))
plot(dat1$k_val, dat1$max_growth,col="#1b98e0",type='l',lwd = 4,
     ylim=range(c(0,dat$avg[dat$salt=="NaCl"]+dat$stdev[dat$salt=="NaCl"])), xlab="Osmolarity (mOsm)", ylab="Max Growth Rate (1/h)",
     main="",cex = 1.25,cex.axis=1.5, cex.lab = 1.5,xlim=range(0,1000))
points(dat2$k_val, dat2$max_growth,col="green",type='l',lwd = 4,
       cex = 1.25)
points(dat$osm[dat$salt=="NaCl"], dat$avg[dat$salt=="NaCl"],col="black", pch = 19,
       cex = 1.25)
arrows(dat$osm[dat$salt=="NaCl"], dat$avg[dat$salt=="NaCl"]-dat$stdev[dat$salt=="NaCl"],
       dat$osm[dat$salt=="NaCl"], dat$avg[dat$salt=="NaCl"]+dat$stdev[dat$salt=="NaCl"],
       length=0.05, angle=90, code=3,lwd=2)
points(dat$osm[dat$salt=="KCl"], dat$avg[dat$salt=="KCl"],col="red", pch = 18,
       cex = 1.25)
arrows(dat$osm[dat$salt=="KCl"], dat$avg[dat$salt=="KCl"]-dat$stdev[dat$salt=="KCl"],
       dat$osm[dat$salt=="KCl"], dat$avg[dat$salt=="KCl"]+dat$stdev[dat$salt=="KCl"],
       length=0.05, angle=90, code=3,lwd=2,col="red")
legend(x = "topleft",inset = c(0.25, -0.2), legend = c("NaCl", "KCl"),
       cex = 1.20, pch = c(19, 18), title = "Experimental",
       col = c("black", "red"),xpd = TRUE, horiz = TRUE,bty = 'n')
legend(x = "topright",inset = c(0.25, -0.2), legend = c("NaCl", "KCl"),
       cex = 1.20,lty = c(1,1), lwd = c(3.5,3.5), title = "Simulated",
       col = c("green","#1b98e0"),xpd = TRUE, horiz = TRUE,bty = 'n')
dev.off()

# Stats
  library(energy)
  # NaCl
    dat_nacl<-read.csv("./Experimental/NaCl_max_growth.csv",sep=",",dec=".",h=T)
    # linear
      linear_salt_effect.lm <- lm(dat_nacl$mu_norm ~ dat_nacl$NaCl)
      summary(linear_salt_effect.lm)
    # quadratic
      quad_salt_effect.lm <- lm(dat_nacl$mu ~ poly(dat_nacl$NaCl, 2, raw=TRUE))
      summary(quad_salt_effect.lm)
    # log all
      log_salt_effect.lm <- lm(dat_nacl$mu ~ log(dat_nacl$NaCl))
      summary(log_salt_effect.lm)
    # Distance correlation
      dcor(dat_nacl$NaCl, dat_nacl$mu_norm,index=1)
      
  # KCl
    dat_kcl<-read.csv("./Experimental/KCl_max_growth.csv",sep=",",dec=".",h=T)
    # Quadratic
      quad_salt_effect.lm <- lm(dat_kcl$mu_norm ~ poly(dat_kcl$KCl, 2, raw=TRUE))
      summary(quad_salt_effect.lm)
    # Log
      log_salt_effect.lm <- lm(dat_kcl$mu ~ log(dat_kcl$KCl))
      summary(log_salt_effect.lm)
    
    # Distance correlation
      dcor(dat_kcl$KCl, dat_kcl$mu_norm,index=1)

#### Figure IV b)  ####
dat <- read.csv("./Results/hyper_osm_sim.csv",sep=",",dec=".",h=T) # given by Figure_2_b_c.py
pdf(file = "./Figures/Figure4b.pdf", 
    width = 12, height = 8)
par(mar = c(5.1, 5.1, 4.1, 2.1))  
plot(dat$k_val*1000, dat$max_growth, pch = 19, col="#1b98e0", xlab = "K (mM)", xlim=c(0,50),
     ylab = expression(paste("Max growth rate (h"^"-1",")")),cex.axis=2, cex.lab = 2,cex=2) 
dev.off()

#### Figure IV c) ####
dat<-read.csv("./Experimental/KO_dat_min_max.csv",sep=",",dec=".",h=T)
pdf(file = "./Figures/Figure4c.pdf",
    width = 12, height = 8)
# Plot EDD
par(mar = c(5.1, 5.1, 4.1, 2.1))  
plot(dat$NaCl, dat$EDD_mean,col="green", pch = 15,
     ylim=range(c(0,dat$EDD_mean+dat$EDD_std)), 
     xlab="NaCl concentration (mM)", ylab=expression(paste("Scaled max growth rate (h"^"-1",")")),
     main="",cex = 1.25,cex.axis=2, cex.lab = 2,xlim=range(0,200))
arrows(dat$NaCl, dat$EDD_mean-dat$EDD_std,
       dat$NaCl, dat$EDD_mean+dat$EDD_std,
       length=0.05, angle=90, code=3,lwd=2,col = 'green')

# Plot GND
points(dat$NaCl, dat$GND_mean,col="red", pch = 17,
        cex = 1.25,cex.axis=2, cex.lab = 2)
arrows(dat$NaCl, dat$GND_mean-dat$GND_std,
       dat$NaCl, dat$GND_mean+dat$GND_std,
       length=0.05, angle=90, code=3,lwd=2,col='red')

# Plot ctrl
points(dat$NaCl, dat$CTRL_mean,col="black", pch = 19,
       cex = 1.25,cex.axis=2, cex.lab = 2)
arrows(dat$NaCl, dat$CTRL_mean-dat$CTRL_std,
       dat$NaCl, dat$CTRL_mean+dat$CTRL_std,
       length=0.05, angle=90, code=3,lwd=2,col='black')

# Legend
legend(x = "bottomright", legend = c('Control', expression(paste(Delta, "EDD")),
                                     expression(paste(Delta, "GND"))),
       cex = 2, pch = c(19, 15, 17), title = "Strain", 
       col = c("black","green",'red'),xpd = TRUE,horiz = TRUE)
dev.off()

#### Figure IV d) ####

dat1<-read.csv("./Results/analysis_II_kos.csv",sep=",",dec=".",h=T) # Given by Figure_2b_c with closed reactions
pdf(file = "./Figures/Figure4d.pdf",
    width = 12, height = 8)
par(mar = c(5.1, 5.1, 4.1, 2.1))  
plot(dat1$K, dat1$mu_EDD,col="green", pch = 15,
     ylim=range(c(0.6,0.9)), 
     xlab="K (M)", ylab=expression(paste("Max growth rate (h"^"-1",")")),
     main="",cex = 1.25,cex.axis=2, cex.lab = 2)

points(dat1$K, dat1$mu_GND,col="red", pch = 17,
       cex = 1.25,cex.axis=2, cex.lab = 2)

points(dat1$K, dat1$mu_OG,col="black", pch = 19,
       cex = 1.25,cex.axis=2, cex.lab = 2)

# Stats #
dat<-read.csv("ko_stats.csv",sep=",",dec=".",h=T)
low_0 = c(dat$r1[dat$NaCl==0],dat$r2[dat$NaCl==0],dat$r3[dat$NaCl==0],
          dat$r4[dat$NaCl==0],dat$r5[dat$NaCl==0],dat$r6[dat$NaCl==0])
low_name = c('WT','EDD','GND','WT','EDD','GND','WT','EDD','GND',
               'WT','EDD','GND','WT','EDD','GND','WT','EDD','GND')
kruskal.test(low_0~low_name)

low_1 = c(dat$r1[dat$NaCl==20],dat$r2[dat$NaCl==20],dat$r3[dat$NaCl==20],
          dat$r4[dat$NaCl==20],dat$r5[dat$NaCl==20],dat$r6[dat$NaCl==20])
kruskal.test(low_1~low_name)

low_2 = c(dat$r1[dat$NaCl==40],dat$r2[dat$NaCl==40],dat$r3[dat$NaCl==40],
          dat$r4[dat$NaCl==40],dat$r5[dat$NaCl==40],dat$r6[dat$NaCl==40])
kruskal.test(low_2~low_name)

low_3 = c(dat$r1[dat$NaCl==60],dat$r2[dat$NaCl==60],dat$r3[dat$NaCl==60],
          dat$r4[dat$NaCl==60],dat$r5[dat$NaCl==60],dat$r6[dat$NaCl==60])
kruskal.test(low_3~low_name)

low_4 = c(dat$r1[dat$NaCl==80],dat$r2[dat$NaCl==80],dat$r3[dat$NaCl==80],
          dat$r4[dat$NaCl==80],dat$r5[dat$NaCl==80],dat$r6[dat$NaCl==80])
kruskal.test(low_4~low_name)

low_5 = c(dat$r1[dat$NaCl==100],dat$r2[dat$NaCl==100],dat$r3[dat$NaCl==100],
          dat$r4[dat$NaCl==100],dat$r5[dat$NaCl==100],dat$r6[dat$NaCl==100])
kruskal.test(low_5~low_name)

low_6 = c(dat$r1[dat$NaCl==120],dat$r2[dat$NaCl==120],dat$r3[dat$NaCl==120],
          dat$r4[dat$NaCl==120],dat$r5[dat$NaCl==120],dat$r6[dat$NaCl==120])
kruskal.test(low_6~low_name)

low_7 = c(dat$r1[dat$NaCl==140],dat$r2[dat$NaCl==140],dat$r3[dat$NaCl==140],
          dat$r4[dat$NaCl==140],dat$r5[dat$NaCl==140],dat$r6[dat$NaCl==140])
kruskal.test(low_7~low_name)

low_8 = c(dat$r1[dat$NaCl==180],dat$r2[dat$NaCl==180],dat$r3[dat$NaCl==180],
          dat$r4[dat$NaCl==180],dat$r5[dat$NaCl==180],dat$r6[dat$NaCl==180])
kruskal.test(low_8~low_name)

low_9 = c(dat$r1[dat$NaCl==160],dat$r2[dat$NaCl==160],dat$r3[dat$NaCl==160],
           dat$r4[dat$NaCl==160],dat$r5[dat$NaCl==160],dat$r6[dat$NaCl==160])
kruskal.test(low_9~low_name)

### Figure 4 e) ###

dat<-read.csv("./Experimental/Metabolome_C18.csv",sep=",",dec=".",h=T)
## filter what is too low ##
too_low = 1
for(i in 1:length(dat$Name)){
  temp <- c(dat$Ratio_NaCl_1[i],dat$Ratio_NaCl_25[i],dat$Ratio_NaCl_50[i],dat$Ratio_NaCl_100[i],dat$Ratio_NaCl_400[i])
  if(max(temp)<5){
    too_low[length(too_low)+1] = i
  }
}
too_low_dat = too_low[-1]

Name = dat$Name[-too_low_dat]
Formula = dat$Formula[-too_low_dat]
Area_NaCl_1_1 = dat$Area_1_1[-too_low_dat]
Area_NaCl_1_2 = dat$Area_1_2[-too_low_dat]
Area_NaCl_25_1 = dat$Area_25_1[-too_low_dat]
Area_NaCl_25_2 = dat$Area_25_2[-too_low_dat]
Area_NaCl_50_1 = dat$Area_50_1[-too_low_dat]
Area_NaCl_50_2 = dat$Area_50_2[-too_low_dat]
Area_NaCl_100_1 = dat$Area_100_1[-too_low_dat]
Area_NaCl_100_2 = dat$Area_100_2[-too_low_dat]
Area_NaCl_400_1 = dat$Area_873_1[-too_low_dat]
Area_NaCl_400_2 = dat$Area_873_2[-too_low_dat]

highest_max = 1
second_max = 1
third_max = 1
last_max = 1
highest_min = 1
second_min = 1
third_min = 1
last_min = 1
for(i in 1:length(Name)){
  max_1 = max(c(Area_NaCl_1_1[i],Area_NaCl_1_2[i]))
  max_25 = max(c(Area_NaCl_25_1[i],Area_NaCl_25_2[i]))
  max_50 = max(c(Area_NaCl_50_1[i],Area_NaCl_50_2[i]))
  max_100 = max(c(Area_NaCl_100_1[i],Area_NaCl_100_2[i]))
  max_400 = max(c(Area_NaCl_400_1[i],Area_NaCl_400_2[i]))
  temp <- c(max_1,max_25,max_50,max_100,max_400)
  highest_max[i] = which.max(temp)
  last_max[i] = which.min(temp)
  temp[which.max(temp)] = 0
  second_max[i] = which.max(temp)
  temp[which.max(temp)] = 0
  third_max[i] =  which.max(temp)
  
  min_1 = min(c(Area_NaCl_1_1[i],Area_NaCl_1_2[i]))
  min_25 = min(c(Area_NaCl_25_1[i],Area_NaCl_25_2[i]))
  min_50 = min(c(Area_NaCl_50_1[i],Area_NaCl_50_2[i]))
  min_100 = min(c(Area_NaCl_100_1[i],Area_NaCl_100_2[i]))
  min_400 = min(c(Area_NaCl_400_1[i],Area_NaCl_400_2[i]))
  temp <- c(min_1,min_25,min_50,min_100,min_400)
  highest_min[i] = which.max(temp)
  last_min[i] = which.min(temp)
  temp[which.max(temp)] = 0
  second_min[i] = which.max(temp)
  temp[which.min(temp)] = 0
  third_min[i] =  which.max(temp)
}
NaCl_1_sum_max = (length(which(highest_max[1:length(dat$Name[-too_low_dat])]==1))*5 + length(which(second_max[1:length(dat$Name[-too_low_dat])]==1))*3 + length(which(third_max[1:length(dat$Name[-too_low_dat])]==1))*1)/(length(dat$Name[-too_low_dat])*5)
NaCl_25_sum_max = (length(which(highest_max[1:length(dat$Name[-too_low_dat])]==2))*5 + length(which(second_max[1:length(dat$Name[-too_low_dat])]==2))*3 + length(which(third_max[1:length(dat$Name[-too_low_dat])]==2))*1)/(length(dat$Name[-too_low_dat])*5)
NaCl_50_sum_max = (length(which(highest_max[1:length(dat$Name[-too_low_dat])]==3))*5 + length(which(second_max[1:length(dat$Name[-too_low_dat])]==3))*3 + length(which(third_max[1:length(dat$Name[-too_low_dat])]==3))*1)/(length(dat$Name[-too_low_dat])*5)
NaCl_100_sum_max = (length(which(highest_max[1:length(dat$Name[-too_low_dat])]==4))*5 + length(which(second_max[1:length(dat$Name[-too_low_dat])]==4))*3 + length(which(third_max[1:length(dat$Name[-too_low_dat])]==4))*1)/(length(dat$Name[-too_low_dat])*5)
NaCl_400_sum_max = (length(which(highest_max[1:length(dat$Name[-too_low_dat])]==5))*5 + length(which(second_max[1:length(dat$Name[-too_low_dat])]==5))*3 + length(which(third_max[1:length(dat$Name[-too_low_dat])]==5))*1)/(length(dat$Name[-too_low_dat])*5)

NaCl_1_sum_min = (length(which(highest_min[1:length(dat$Name[-too_low_dat])]==1))*5 + length(which(second_min[1:length(dat$Name[-too_low_dat])]==1))*3 + length(which(third_min[1:length(dat$Name[-too_low_dat])]==1))*1)/(length(dat$Name[-too_low_dat])*5)
NaCl_25_sum_min = (length(which(highest_min[1:length(dat$Name[-too_low_dat])]==2))*5 + length(which(second_min[1:length(dat$Name[-too_low_dat])]==2))*3 + length(which(third_min[1:length(dat$Name[-too_low_dat])]==2))*1)/(length(dat$Name[-too_low_dat])*5)
NaCl_50_sum_min = (length(which(highest_min[1:length(dat$Name[-too_low_dat])]==3))*5 + length(which(second_min[1:length(dat$Name[-too_low_dat])]==3))*3 + length(which(third_min[1:length(dat$Name[-too_low_dat])]==3))*1)/(length(dat$Name[-too_low_dat])*5)
NaCl_100_sum_min = (length(which(highest_min[1:length(dat$Name[-too_low_dat])]==4))*5 + length(which(second_min[1:length(dat$Name[-too_low_dat])]==4))*3 + length(which(third_min[1:length(dat$Name[-too_low_dat])]==4))*1)/(length(dat$Name[-too_low_dat])*5)
NaCl_400_sum_min = (length(which(highest_min[1:length(dat$Name[-too_low_dat])]==5))*5 + length(which(second_min[1:length(dat$Name[-too_low_dat])]==5))*3 + length(which(third_min[1:length(dat$Name[-too_low_dat])]==5))*1)/(length(dat$Name[-too_low_dat])*5)

score_max <- c(NaCl_1_sum_max,NaCl_25_sum_max,NaCl_50_sum_max,NaCl_100_sum_max,NaCl_400_sum_max)
score_min <- c(NaCl_1_sum_min,NaCl_25_sum_min,NaCl_50_sum_min,NaCl_100_sum_min,NaCl_400_sum_min)
score <-c(mean(c(NaCl_1_sum_max,NaCl_1_sum_min)),mean(c(NaCl_25_sum_max,NaCl_25_sum_min)),
          mean(c(NaCl_50_sum_max,NaCl_50_sum_min)),mean(c(NaCl_100_sum_max,NaCl_100_sum_min)),
          mean(c(NaCl_400_sum_max,NaCl_400_sum_min)))
mx_max <- as.matrix(score_max)
mx_min <- as.matrix(score_min)
mx <- as.matrix(score)
pdf(file = "./Figures/Figure4e.pdf",
    width = 8, height = 6)
par(mar = c(5.1, 5.1, 4.1, 2.1))  
barplot_2<-barplot(score,xlab = "NaCl concentration (mM)",ylab = "Abundance index",
                   col = c("red","forestgreen","#1b98e0","purple","grey"),cex.names=1.75,
                   names.arg =c("1","25","50","100","400"),cex.lab=1.75,cex.axis=1.75,
                   ylim=c(0,0.6))
arrows(x0 = barplot_2,
       y0 = mx_max,
       y1 = mx_min,
       angle = 90,
       code = 3,
       length = 0.1)
dev.off()

########## Figure 4f ###########
dat<-read.csv("./Experimental/Metabolome_C18.csv",sep=",",dec=".",h=T)
extreme_rem <- readRDS("./Experimental/extreme_rec.rds") # list of previously analyzed metabolites to reject
extreme = "bacon"
Name = dat$Name
Formula = dat$Formula
Group_Area_NaCl_1 = dat$Group_Area_NaCl_1
Group_Area_NaCl_25 = dat$Group_Area_NaCl_25
Group_Area_NaCl_50 = dat$Group_Area_NaCl_50
Group_Area_NaCl_100 = dat$Group_Area_NaCl_100
Group_Area_NaCl_400 = dat$Group_Area_NaCl_400
Ratio_NaCl_1 = dat$Ratio_NaCl_1
Ratio_NaCl_25 = dat$Ratio_NaCl_25
Ratio_NaCl_50 = dat$Ratio_NaCl_50
Ratio_NaCl_100 = dat$Ratio_NaCl_100
Ratio_NaCl_400 = dat$Ratio_NaCl_400
salt_1 = 1
salt_25 = 1
salt_50 = 1
salt_100 = 1
salt_400 = 1
'%!in%' <- function(x,y)!('%in%'(x,y))
for(i in 1:length(Name)){
  temp <- c(Ratio_NaCl_1[i],Ratio_NaCl_25[i],Ratio_NaCl_50[i],Ratio_NaCl_100[i],Ratio_NaCl_400[i])
  if(temp[which.max(temp)]/temp[which.min(temp)]>=5 && Name[i]!="" && Name[i]%!in%extreme_rem){
    extreme[length(extreme)+1] = Name[i]
    salt_1[length(salt_1)+1] = Ratio_NaCl_1[i]
    salt_25[length(salt_25)+1] = Ratio_NaCl_25[i]
    salt_50[length(salt_50)+1] = Ratio_NaCl_50[i]
    salt_100[length(salt_100)+1] = Ratio_NaCl_100[i]
    salt_400[length(salt_400)+1] = Ratio_NaCl_400[i]
  }
}
extreme = extreme[-1]
salt_1 = salt_1[-1]
salt_25 = salt_25[-1]
salt_50 = salt_50[-1]
salt_100 = salt_100[-1]
salt_400 = salt_400[-1]
remove <-c(1,3,4,9,14,15,16,17,18,19,20,21,22,23,24,26,27,28,29,39,31,32,34,35,36,37,38,40,41,42,43,44,45,46,47,48,50,51,52,53,54,55,56,57,58,59,61,62,63,64,65,66,68,69,70,72,76,77,79,81,82,83,84,86,87,88,89,93,94,95,96,97,08,99,100,101,102,103,104,105,107,110,111,113,114,116,117,118,119,120,122,123,124,125,126,128,129,130,131,132,133,134,138,140,141,144,145)
extreme_C18 = extreme[-remove]
salt_1_C18 = salt_1[-remove]
salt_25_C18 = salt_25[-remove]
salt_50_C18 = salt_50[-remove]
salt_100_C18 = salt_100[-remove]
salt_400_C18 = salt_400[-remove]

dat<-read.csv("./Experimental/Metabolome_HILIQ.csv",sep=",",dec=".",h=T)
extreme = "bacon"
Name = dat$Name
Formula = dat$Formula
Group_Area_NaCl_1 = dat$Group_Area_NaCl_1
Group_Area_NaCl_25 = dat$Group_Area_NaCl_25
Group_Area_NaCl_50 = dat$Group_Area_NaCl_50
Group_Area_NaCl_100 = dat$Group_Area_NaCl_100
Group_Area_NaCl_400 = dat$Group_Area_NaCl_400
Ratio_NaCl_1 = dat$Ratio_NaCl_1
Ratio_NaCl_25 = dat$Ratio_NaCl_25
Ratio_NaCl_50 = dat$Ratio_NaCl_50
Ratio_NaCl_100 = dat$Ratio_NaCl_100
Ratio_NaCl_400 = dat$Ratio_NaCl_400
salt_1 = 1
salt_25 = 1
salt_50 = 1
salt_100 = 1
salt_400 = 1
for(i in 1:length(Name)){
  temp <- c(Ratio_NaCl_1[i],Ratio_NaCl_25[i],Ratio_NaCl_50[i],Ratio_NaCl_100[i],Ratio_NaCl_400[i])
  if(temp[which.max(temp)]/temp[which.min(temp)]>=5 && Name[i]!="" && Name[i]%!in%extreme_rem){
    extreme[length(extreme)+1] = Name[i]
    salt_1[length(salt_1)+1] = Ratio_NaCl_1[i]
    salt_25[length(salt_25)+1] = Ratio_NaCl_25[i]
    salt_50[length(salt_50)+1] = Ratio_NaCl_50[i]
    salt_100[length(salt_100)+1] = Ratio_NaCl_100[i]
    salt_400[length(salt_400)+1] = Ratio_NaCl_400[i]
  }
}
extreme = extreme[-1]
salt_1 = salt_1[-1]
salt_25 = salt_25[-1]
salt_50 = salt_50[-1]
salt_100 = salt_100[-1]
salt_400 = salt_400[-1]
keep_HILIQ <-c(4,7,16,19,19,23,35,36,38,45,47,52,53,54,62,80,89,90,95,106,107,111,114,116,119)
extreme_HILIQ = extreme[keep_HILIQ]
salt_1_HILIQ = salt_1[keep_HILIQ]
salt_25_HILIQ = salt_25[keep_HILIQ]
salt_50_HILIQ = salt_50[keep_HILIQ]
salt_100_HILIQ = salt_100[keep_HILIQ]
salt_400_HILIQ = salt_400[keep_HILIQ]
extreme<-c(extreme_C18,extreme_HILIQ)
salt_1<-c(salt_1_C18,salt_1_HILIQ)
salt_25<-c(salt_25_C18,salt_25_HILIQ)
salt_50<-c(salt_50_C18,salt_50_HILIQ)
salt_100<-c(salt_100_C18,salt_100_HILIQ)
salt_400<-c(salt_400_C18,salt_400_HILIQ)
remove <-c(5,6,9,11,14,15,17,19,21,22,23,24,26,28,32,33,35,36,38,40,45,46,52,59,63)
extreme = extreme[-remove]
salt_1 = salt_1[-remove]
salt_25 = salt_25[-remove]
salt_50 = salt_50[-remove]
salt_100 = salt_100[-remove]
salt_400 = salt_400[-remove]
remove <-c(5,12,13,36,17,15,10,8,19,23,18)
extreme = extreme[-remove]
salt_1 = salt_1[-remove]
salt_25 = salt_25[-remove]
salt_50 = salt_50[-remove]
salt_100 = salt_100[-remove]
salt_400 = salt_400[-remove]
datex <- as.matrix(data.frame(extreme,salt_1,salt_25,salt_50,salt_100,salt_400))

# Re-order for ease of read
highest = 1
for(i in 1:length(extreme)){
  temp <- c(salt_1[i],salt_25[i],salt_50[i],salt_100[i],salt_400[i])
  highest[i] = which.max(temp)
}
# Final all top with BiGG name
datex <- as.matrix(data.frame(salt_1[order(highest)],salt_25[order(highest)],
                              salt_50[order(highest)],salt_100[order(highest)],
                              salt_400[order(highest)]))
colnames(datex) <- c("1","25","50","100","400")
Heat_mets <- c("Alanine","Phosphoric acid","Thiosulfate","Ribulose 5-phosphate","Adenosine",
               "6-phosphogluconate","Malate","Phenylalanine","Guanosine","Fructose 1,6-bisphosphate",
               "Glutamate","Uracil","Proline","Histidine","Mannitol 1-phosphate",
               "ADP","Glutatione","N2-Succinyl-L-glutamate","Pyruvate","Valine","Xanthine",
               "Phophoenolpyruvate","Tyrosine","MEP","Maltose/Trehalose","Reduced glutathione",
               "Histidinol")
rownames(datex) = Heat_mets
# Generate heat map
hm <- heatmap(datex,Rowv = NA, Colv = NA)
labs <- colnames(datex)[hm$colInd]
pdf(file = "./Figures/Figure4f.pdf",
    width = 8, height = 6)
heatmap(datex,Rowv = NA, Colv = NA,labCol = "",xlab="NaCl concentration (mM)",
        cex.lab=4,add.expr = text(x = seq_along(labs), y = -.5, labels = colnames(datex), xpd = TRUE,cex=1.75))
dev.off()

