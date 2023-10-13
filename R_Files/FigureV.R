# This code is used to generate all the figures from the new metabolome data
setwd("C:/Users/alexa/OneDrive - University of Toronto/Documents/PhD_UofT/Paper/Osm/Exp")
########## Figure 6a ###########
dat<-read.csv("Metabolome_HILIQ.csv",sep=",",dec=".",h=T)
dat2<-read.csv("Metabolome_C18.csv",sep=",",dec=".",h=T)
## filter what is too low ##
too_low = 1
for(i in 1:length(dat$Name)){
  temp <- c(dat$Ratio_NaCl_1[i],dat$Ratio_NaCl_25[i],dat$Ratio_NaCl_50[i],dat$Ratio_NaCl_100[i],dat$Ratio_NaCl_400[i])
  if(max(temp)<5){
    too_low[length(too_low)+1] = i
  }
}
too_low_dat = too_low[-1]
too_low = 1
for(i in 1:length(dat2$Name)){
  temp <- c(dat2$Ratio_NaCl_1[i],dat2$Ratio_NaCl_25[i],dat2$Ratio_NaCl_50[i],dat2$Ratio_NaCl_100[i],dat2$Ratio_NaCl_400[i])
  if(max(temp)<5){
    too_low[length(too_low)+1] = i
  }
}
too_low_dat2 = too_low[-1]
Name = c(dat$Name[-too_low_dat],dat2$Name[-too_low_dat2])
Formula = c(dat$Formula[-too_low_dat],dat2$Formula[-too_low_dat2])
Group_Area_NaCl_1 = c(dat$Group_Area_NaCl_1[-too_low_dat],dat2$Group_Area_NaCl_1[-too_low_dat2])
Group_Area_NaCl_25 = c(dat$Group_Area_NaCl_25[-too_low_dat],dat2$Group_Area_NaCl_25[-too_low_dat2])
Group_Area_NaCl_50 = c(dat$Group_Area_NaCl_50[-too_low_dat],dat2$Group_Area_NaCl_50[-too_low_dat2])
Group_Area_NaCl_100 = c(dat$Group_Area_NaCl_100[-too_low_dat],dat2$Group_Area_NaCl_100[-too_low_dat2])
Group_Area_NaCl_400 = c(dat$Group_Area_NaCl_400[-too_low_dat],dat2$Group_Area_NaCl_400[-too_low_dat2])

highest = 1
second = 1
third = 1
last = 1
for(i in 1:length(Name)){
  temp <- c(Group_Area_NaCl_1[i],Group_Area_NaCl_25[i],Group_Area_NaCl_50[i],Group_Area_NaCl_100[i],Group_Area_NaCl_400[i])
  highest[i] = which.max(temp)
  last[i] = which.min(temp)
  temp[which.max(temp)] = 0
  second[i] = which.max(temp)
  temp[which.max(temp)] = 0
  third[i] =  which.max(temp)
}
# factor most abundant by a weight
NaCl_1_sum_HILIQ = (length(which(highest[1:length(dat$Name[-too_low_dat])]==1))*5 + length(which(second[1:length(dat$Name[-too_low_dat])]==1))*3 + length(which(third[1:length(dat$Name[-too_low_dat])]==1))*1)/(length(dat$Name[-too_low_dat])*5)
NaCl_25_sum_HILIQ = (length(which(highest[1:length(dat$Name[-too_low_dat])]==2))*5 + length(which(second[1:length(dat$Name[-too_low_dat])]==2))*3 + length(which(third[1:length(dat$Name[-too_low_dat])]==2))*1)/(length(dat$Name[-too_low_dat])*5)
NaCl_50_sum_HILIQ = (length(which(highest[1:length(dat$Name[-too_low_dat])]==3))*5 + length(which(second[1:length(dat$Name[-too_low_dat])]==3))*3 + length(which(third[1:length(dat$Name[-too_low_dat])]==3))*1)/(length(dat$Name[-too_low_dat])*5)
NaCl_100_sum_HILIQ = (length(which(highest[1:length(dat$Name[-too_low_dat])]==4))*5 + length(which(second[1:length(dat$Name[-too_low_dat])]==4))*3 + length(which(third[1:length(dat$Name[-too_low_dat])]==4))*1)/(length(dat$Name[-too_low_dat])*5)
NaCl_400_sum_HILIQ = (length(which(highest[1:length(dat$Name[-too_low_dat])]==5))*5 + length(which(second[1:length(dat$Name[-too_low_dat])]==5))*3 + length(which(third[1:length(dat$Name[-too_low_dat])]==5))*1)/(length(dat$Name[-too_low_dat])*5)

NaCl_1_sum_C18 = (length(which(highest[length(dat$Name[-too_low_dat])+1:length(highest)]==1))*5 + length(which(second[length(dat$Name[-too_low_dat])+1:length(highest)]==1))*3 + length(which(third[length(dat$Name[-too_low_dat])+1:length(highest)]==1))*1)/(length(dat2$Name[-too_low_dat2])*5)
NaCl_25_sum_C18 = (length(which(highest[length(dat$Name[-too_low_dat])+1:length(highest)]==2))*5 + length(which(second[length(dat$Name[-too_low_dat])+1:length(highest)]==2))*3 + length(which(third[length(dat$Name[-too_low_dat])+1:length(highest)]==2))*1)/(length(dat2$Name[-too_low_dat2])*5)
NaCl_50_sum_C18 = (length(which(highest[length(dat$Name[-too_low_dat])+1:length(highest)]==3))*5 + length(which(second[length(dat$Name[-too_low_dat])+1:length(highest)]==3))*3 + length(which(third[length(dat$Name[-too_low_dat])+1:length(highest)]==3))*1)/(length(dat2$Name[-too_low_dat2])*5)
NaCl_100_sum_C18 = (length(which(highest[length(dat$Name[-too_low_dat])+1:length(highest)]==4))*5 + length(which(second[length(dat$Name[-too_low_dat])+1:length(highest)]==4))*3 + length(which(third[length(dat$Name[-too_low_dat])+1:length(highest)]==4))*1)/(length(dat2$Name[-too_low_dat2])*5)
NaCl_400_sum_C18 = (length(which(highest[length(dat$Name[-too_low_dat])+1:length(highest)]==5))*5 + length(which(second[length(dat$Name[-too_low_dat])+1:length(highest)]==5))*3 + length(which(third[length(dat$Name[-too_low_dat])+1:length(highest)]==5))*1)/(length(dat2$Name[-too_low_dat2])*5)


score_max <- c(max(NaCl_1_sum_C18,NaCl_1_sum_HILIQ),max(NaCl_25_sum_C18,NaCl_25_sum_HILIQ),
               max(NaCl_50_sum_C18,NaCl_50_sum_HILIQ),max(NaCl_100_sum_C18,NaCl_100_sum_HILIQ),
               max(NaCl_400_sum_C18,NaCl_400_sum_HILIQ))
score_min <- c(min(NaCl_1_sum_C18,NaCl_1_sum_HILIQ),min(NaCl_25_sum_C18,NaCl_25_sum_HILIQ),
               min(NaCl_50_sum_C18,NaCl_50_sum_HILIQ),min(NaCl_100_sum_C18,NaCl_100_sum_HILIQ),
               min(NaCl_400_sum_C18,NaCl_400_sum_HILIQ))
score <-c(mean(c(NaCl_1_sum_C18,NaCl_1_sum_HILIQ)),mean(c(NaCl_25_sum_C18,NaCl_25_sum_HILIQ)),
             mean(c(NaCl_50_sum_C18,NaCl_50_sum_HILIQ)),mean(c(NaCl_100_sum_C18,NaCl_100_sum_HILIQ)),
             mean(c(NaCl_400_sum_C18,NaCl_400_sum_HILIQ)))
mx_max <- as.matrix(score_max)
mx_min <- as.matrix(score_min)
mx <- as.matrix(score)
pdf(file = "C:/Users/alexa/OneDrive - University of Toronto/Documents/PhD_UofT/Paper/Osm/Figure/Figure6a_col.pdf",
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

########## Figure 6b ###########
datHILIQ<-read.csv("Metabolome_HILIQ.csv",sep=",",dec=".",h=T)
datC18<-read.csv("Metabolome_C18.csv",sep=",",dec=".",h=T)
std <-c("Acetyl-CoA","Glutamic Acid","F6bp","Glycerol 3-phosphate","6-phosphogluconate",
        "Aspartic Acid","Valine","Glutamine","Alanine","Fumaric Acid","Citric Acid",
        "Malic Acid","G6P")
conc <-c(1,200,25,5,1,10,10,10,10,1,5,5,10)

# C18 #
data_conc_c18 <- data.frame(Name = "none",Nacl_1_1 =0,Nacl_1_2 = 0, Nacl_25_1 =0,Nacl_25_2 =0,
                            Nacl_50_1 =0, Nacl_50_2 =0, Nacl_100_1 =0, Nacl_100_2 =0,
                            Nacl_400_1 =0, Nacl_400_2 =0)
for(i in 1:length(std)){
  if(sum(datC18$Name==std[i])>0){
    for(j in 1:sum(datC18$Name==std[i])){
      y = c(1/conc[i],1/(conc[i]*0.1),1/(conc[i]*0.01))
      x = c(1/datC18$Area_Std1[datC18$Name==std[i]][j],
            1/datC18$Area_Std2[datC18$Name==std[i]][j],
            1/datC18$Area_Std3[datC18$Name==std[i]][j])
      cf = lm(y ~ x)
      if(summary(lm(y~x))$r.squared >= 0.90){
        S1_1 = 1/(datC18$Area_1_1[datC18$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S1_2 = 1/(datC18$Area_1_2[datC18$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S25_1 = 1/(datC18$Area_25_1[datC18$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S25_2 = 1/(datC18$Area_25_2[datC18$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S50_1 = 1/(datC18$Area_50_1[datC18$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S50_2 = 1/(datC18$Area_50_2[datC18$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S100_1 = 1/(datC18$Area_100_1[datC18$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S100_2 = 1/(datC18$Area_100_2[datC18$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S400_1 = 1/(datC18$Area_873_1[datC18$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S400_2 = 1/(datC18$Area_873_2[datC18$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        Std_val = data.frame(std[i],1/S1_1,1/S1_2,1/S25_1,1/S25_2,1/S50_1,
                             1/S50_2,1/S100_1,1/S100_2,1/S400_1,1/S400_2)
        data_conc_c18[nrow(data_conc_c18)+1,] = Std_val
      }
    }
  }
}

# HILIQ #
data_conc_HILIQ <- data.frame(Name = "none",Nacl_1_1 =0,Nacl_1_2 = 0, Nacl_25_1 =0,Nacl_25_2 =0,
                              Nacl_50_1 =0, Nacl_50_2 =0, Nacl_100_1 =0, Nacl_100_2 =0,
                              Nacl_400_1 =0, Nacl_400_2 =0)
for(i in 1:length(std)){
  if(sum(datHILIQ$Name==std[i])>0){
    for(j in 1:sum(datHILIQ$Name==std[i])){
      y = c(1/conc[i],1/(conc[i]*0.1),1/(conc[i]*0.01))
      x = c(1/datHILIQ$Area_Std1[datHILIQ$Name==std[i]][j],
            1/datHILIQ$Area_Std2[datHILIQ$Name==std[i]][j],
            1/datHILIQ$Area_Std3[datHILIQ$Name==std[i]][j])
      cf = lm(y ~ x)
      if(summary(lm(y~x))$r.squared >= 0.90){
        S1_1 = 1/(datHILIQ$Area_1_1[datHILIQ$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S1_2 = 1/(datHILIQ$Area_1_2[datHILIQ$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S25_1 = 1/(datHILIQ$Area_25_1[datHILIQ$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S25_2 = 1/(datHILIQ$Area_25_2[datHILIQ$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S50_1 = 1/(datHILIQ$Area_50_1[datHILIQ$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S50_2 = 1/(datHILIQ$Area_50_2[datHILIQ$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S100_1 = 1/(datHILIQ$Area_100_1[datHILIQ$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S100_2 = 1/(datHILIQ$Area_100_2[datHILIQ$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S400_1 = 1/(datHILIQ$Area_873_1[datHILIQ$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        S400_2 = 1/(datHILIQ$Area_873_2[datHILIQ$Name==std[i]][j])*cf$coefficients[2]+cf$coefficients[1]
        Std_val = data.frame(std[i],1/S1_1,1/S1_2,1/S25_1,1/S25_2,1/S50_1,
                             1/S50_2,1/S100_1,1/S100_2,1/S400_1,1/S400_2)
        data_conc_HILIQ[nrow(data_conc_HILIQ)+1,] = Std_val
      }
    }
  }
}

data_conc_HILIQ[3,1] = "Acetyl-CoA 2"
data_conc_HILIQ[4,1] = "F6bp 2"
data_conc_c18[data_conc_c18<0] = 0
data_conc_HILIQ[data_conc_HILIQ<0] = 0

C18_keep <- c("Alanine")
HILIQ_keep <-c("Acetyl-CoA")
Ala <-c(mean(c(data_conc_c18$Nacl_1_1[data_conc_c18$Name==C18_keep[1]],data_conc_c18$Nacl_1_2[data_conc_c18$Name==C18_keep[1]])),
           mean(c(data_conc_c18$Nacl_25_1[data_conc_c18$Name==C18_keep[1]],data_conc_c18$Nacl_25_2[data_conc_c18$Name==C18_keep[1]])),
           mean(c(data_conc_c18$Nacl_50_1[data_conc_c18$Name==C18_keep[1]],data_conc_c18$Nacl_50_2[data_conc_c18$Name==C18_keep[1]])),
           mean(c(data_conc_c18$Nacl_100_1[data_conc_c18$Name==C18_keep[1]],data_conc_c18$Nacl_100_2[data_conc_c18$Name==C18_keep[1]])),
           mean(c(data_conc_c18$Nacl_400_1[data_conc_c18$Name==C18_keep[1]],data_conc_c18$Nacl_400_2[data_conc_c18$Name==C18_keep[1]])))

Accoa <-c(mean(c(data_conc_HILIQ$Nacl_1_1[data_conc_HILIQ$Name==HILIQ_keep[1]],data_conc_HILIQ$Nacl_1_2[data_conc_HILIQ$Name==HILIQ_keep[1]])),
             mean(c(data_conc_HILIQ$Nacl_25_1[data_conc_HILIQ$Name==HILIQ_keep[1]],data_conc_HILIQ$Nacl_25_2[data_conc_HILIQ$Name==HILIQ_keep[1]])),
             mean(c(data_conc_HILIQ$Nacl_50_1[data_conc_HILIQ$Name==HILIQ_keep[1]],data_conc_HILIQ$Nacl_50_2[data_conc_HILIQ$Name==HILIQ_keep[1]])),
             mean(c(data_conc_HILIQ$Nacl_100_1[data_conc_HILIQ$Name==HILIQ_keep[1]],data_conc_HILIQ$Nacl_100_2[data_conc_HILIQ$Name==HILIQ_keep[1]])),
             mean(c(data_conc_HILIQ$Nacl_400_1[data_conc_HILIQ$Name==HILIQ_keep[1]],data_conc_HILIQ$Nacl_400_2[data_conc_HILIQ$Name==HILIQ_keep[1]])))

Ala_max <-c(max(c(data_conc_c18$Nacl_1_1[data_conc_c18$Name==C18_keep[1]],data_conc_c18$Nacl_1_2[data_conc_c18$Name==C18_keep[1]])),
           max(c(data_conc_c18$Nacl_25_1[data_conc_c18$Name==C18_keep[1]],data_conc_c18$Nacl_25_2[data_conc_c18$Name==C18_keep[1]])),
           max(c(data_conc_c18$Nacl_50_1[data_conc_c18$Name==C18_keep[1]],data_conc_c18$Nacl_50_2[data_conc_c18$Name==C18_keep[1]])),
           max(c(data_conc_c18$Nacl_100_1[data_conc_c18$Name==C18_keep[1]],data_conc_c18$Nacl_100_2[data_conc_c18$Name==C18_keep[1]])),
           max(c(data_conc_c18$Nacl_400_1[data_conc_c18$Name==C18_keep[1]],data_conc_c18$Nacl_400_2[data_conc_c18$Name==C18_keep[1]])))
Ala_min <-c(min(c(data_conc_c18$Nacl_1_1[data_conc_c18$Name==C18_keep[1]],data_conc_c18$Nacl_1_2[data_conc_c18$Name==C18_keep[1]])),
            min(c(data_conc_c18$Nacl_25_1[data_conc_c18$Name==C18_keep[1]],data_conc_c18$Nacl_25_2[data_conc_c18$Name==C18_keep[1]])),
            min(c(data_conc_c18$Nacl_50_1[data_conc_c18$Name==C18_keep[1]],data_conc_c18$Nacl_50_2[data_conc_c18$Name==C18_keep[1]])),
            min(c(data_conc_c18$Nacl_100_1[data_conc_c18$Name==C18_keep[1]],data_conc_c18$Nacl_100_2[data_conc_c18$Name==C18_keep[1]])),
            min(c(data_conc_c18$Nacl_400_1[data_conc_c18$Name==C18_keep[1]],data_conc_c18$Nacl_400_2[data_conc_c18$Name==C18_keep[1]])))

Accoa_max <-c(max(c(data_conc_HILIQ$Nacl_1_1[data_conc_HILIQ$Name==HILIQ_keep[1]],data_conc_HILIQ$Nacl_1_2[data_conc_HILIQ$Name==HILIQ_keep[1]])),
             max(c(data_conc_HILIQ$Nacl_25_1[data_conc_HILIQ$Name==HILIQ_keep[1]],data_conc_HILIQ$Nacl_25_2[data_conc_HILIQ$Name==HILIQ_keep[1]])),
             max(c(data_conc_HILIQ$Nacl_50_1[data_conc_HILIQ$Name==HILIQ_keep[1]],data_conc_HILIQ$Nacl_50_2[data_conc_HILIQ$Name==HILIQ_keep[1]])),
             max(c(data_conc_HILIQ$Nacl_100_1[data_conc_HILIQ$Name==HILIQ_keep[1]],data_conc_HILIQ$Nacl_100_2[data_conc_HILIQ$Name==HILIQ_keep[1]])),
             max(c(data_conc_HILIQ$Nacl_400_1[data_conc_HILIQ$Name==HILIQ_keep[1]],data_conc_HILIQ$Nacl_400_2[data_conc_HILIQ$Name==HILIQ_keep[1]])))
Accoa_min <-c(min(c(data_conc_HILIQ$Nacl_1_1[data_conc_HILIQ$Name==HILIQ_keep[1]],data_conc_HILIQ$Nacl_1_2[data_conc_HILIQ$Name==HILIQ_keep[1]])),
              min(c(data_conc_HILIQ$Nacl_25_1[data_conc_HILIQ$Name==HILIQ_keep[1]],data_conc_HILIQ$Nacl_25_2[data_conc_HILIQ$Name==HILIQ_keep[1]])),
              min(c(data_conc_HILIQ$Nacl_50_1[data_conc_HILIQ$Name==HILIQ_keep[1]],data_conc_HILIQ$Nacl_50_2[data_conc_HILIQ$Name==HILIQ_keep[1]])),
              min(c(data_conc_HILIQ$Nacl_100_1[data_conc_HILIQ$Name==HILIQ_keep[1]],data_conc_HILIQ$Nacl_100_2[data_conc_HILIQ$Name==HILIQ_keep[1]])),
              min(c(data_conc_HILIQ$Nacl_400_1[data_conc_HILIQ$Name==HILIQ_keep[1]],data_conc_HILIQ$Nacl_400_2[data_conc_HILIQ$Name==HILIQ_keep[1]])))

met_name <- c("Alanine","Acetyl-CoA")
met_mean = data.frame(Ala,Accoa)
met_min = data.frame(Ala_min,Accoa_min)
met_max = data.frame(Ala_max,Accoa_max)
mx_min <- as.matrix(met_min)
mx_max <- as.matrix(met_max)
mx <- as.matrix(met_mean)
colnames(mx) <- met_name

sim<-read.csv("C:/Users/alexa/OneDrive - University of Toronto/Documents/PhD_UofT/Paper/Osm/Analysis/IV/CVA_with_k/FigureIIb.csv",sep=",",dec=".",h=T)

pdf(file = "C:/Users/alexa/OneDrive - University of Toronto/Documents/PhD_UofT/Paper/Osm/Figure/Figure6b.pdf",
    width = 10, height = 6)
par(mar = c(5.1, 5.1, 4.1, 2.1))  
base_r_barplot <- barplot(mx,ylim = c(0, 1.5),col = c("red","forestgreen","#1b98e0","purple","grey"),
                          names.arg = met_name,beside=TRUE,ylab="Metabolite concentration (mM)",
                          cex.lab=1.75,cex.names = 1.75,cex.axis=1.75)
legend(x = "top",inset = c(0.25, -0.2), legend = c("1", "25","50","100","400"),
       cex = 1.5, title = "NaCl Concentration (mM)",
       fill = c("red","forestgreen","#1b98e0","purple","grey"),xpd = TRUE, horiz = TRUE,bty = 'n')
arrows(x0 = base_r_barplot,                           # Add error bars
       y0 = mx_max,
       y1 = mx_min,
       angle = 90,
       code = 3,
       length = 0.1)
dev.off()
########## Figure 6c ###########
dat<-read.csv("Metabolome_C18.csv",sep=",",dec=".",h=T)
extreme_rem <- readRDS("C:/Users/alexa/OneDrive - University of Toronto/Documents/PhD_UofT/Paper/Osm/Exp/extreme_rec.rds")
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

dat<-read.csv("Metabolome_HILIQ.csv",sep=",",dec=".",h=T)
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
remove <-c(4,11,12,35,16,14,9,7,18,22,17)+1
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
               "Phophoenolpyruvate","Tyrosine","MEP","Maltose","Reduced glutathione",
               "Histidinol")
rownames(datex) = Heat_mets
# Generate heat map
hm <- heatmap(datex,Rowv = NA, Colv = NA)
labs <- colnames(datex)[hm$colInd]
pdf(file = "C:/Users/alexa/OneDrive - University of Toronto/Documents/PhD_UofT/Paper/Osm/Figure/Figure6c.pdf",
    width = 8, height = 6)
heatmap(datex,Rowv = NA, Colv = NA,labCol = "",xlab="NaCl concentration (mM)",
        cex.lab=4,add.expr = text(x = seq_along(labs), y = -.5, labels = colnames(datex), xpd = TRUE,cex=1.75))
dev.off()

###### Figure 6d ######
setwd("C:/Users/alexa/OneDrive - University of Toronto/Documents/PhD_UofT/Paper/Osm/Analysis/IV/CVA_with_k")
dat<-read.csv("accoa_bioreac.csv",sep=",",dec=".",h=T)
pdf(file = "C:/Users/alexa/OneDrive - University of Toronto/Documents/PhD_UofT/Paper/Osm/Figure/Figure6d.pdf",
    width = 10, height = 6)
par(mar = c(5.1, 5.1, 4.1, 2.1))  
barplot<-barplot(log(dat$mean[dat$mets=="x_accoa_c"]*1000*1000),col = c("red","forestgreen","#1b98e0","purple","grey"),
        names.arg = c("0.636","0.754","0.872","0.902","0.444"),beside=TRUE,ylab="log of basal acetyl-coa concentration (µM)",
        ylim=c(0,8),cex.lab=1.75,cex.names = 1.75,cex.axis=1.75,xlab = "Growth rate (1/h)")

# hack: we draw arrows but with very special "arrowheads"
arrows(x0 = barplot,                           
       y0 = log(dat$max[dat$mets=="x_accoa_c"]*1000*1000),
       y1 = log(dat$min[dat$mets=="x_accoa_c"]*1000*1000),
       angle = 90,
       code = 3,
       length = 0.1)
dev.off()
