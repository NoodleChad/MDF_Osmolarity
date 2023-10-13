#### FigureS1 - With charge ####
setwd("C:/Users/alexa/OneDrive - University of Toronto/Documents/PhD_UofT/Paper/Osm/Analysis/II")
dat<-read.csv("analysis_II_sup_charge.csv",sep=",",dec=".",h=T)
dat1<-read.csv("analysis_II_last.csv",sep=",",dec=".",h=T)

# Draw first plot using axis y1
pdf(file = "C:/Users/alexa/OneDrive - University of Toronto/Documents/PhD_UofT/Paper/Osm/Figure/Figure5a.pdf",
    width = 8.27, height = 5.83)
par(mar = c(7, 4, 1, 4) + 0.3)  
plot(dat$K, dat$max_growth, pch = 16, col = 1, 
     xlab = "K (mOsm)",
     ylab = "Growth rate (1/h)",cex.axis=1.5, cex.lab = 1.5) 
points(dat1$K, dat1$max_growth,pch=18, col = "red")
legend(x=0.016, y=0.8,legend = c("With Charge Balance", "Without Charge Balance"), pch = c(16,18),col = c("black","red"), cex=1.5)
dev.off()

#### FigureIV b) - With charge ####
setwd("C:/Users/alexa/OneDrive - University of Toronto/Documents/PhD_UofT/Paper/Osm/Analysis/IV/CVA_alex")
dat<-read.csv("CVA_charge_compare.csv",sep=",",dec=".",h=T)
mets_4_fig <-c("ac","accoa","gln__L","glu__L","orot","thdp")
name <-c("Acetate","Acetyl-CoA","Glutamine","Glutamate",
         "Orotate","THDP")
mets = 1
mets_name = 1
min = 1
max = 1
mean = 1
Charge = "test"
for(i in 1:length(mets_4_fig)){
  mets[length(mets)+1] = mets_4_fig[i]
  mets[length(mets)+1] = mets_4_fig[i]
  mets_name[length(mets_name)+1] = name[i]
  mets_name[length(mets_name)+1] = name[i]
  min[length(min)+1] = dat$min[dat$mets==mets_4_fig[i]][1]
  min[length(min)+1] = dat$min[dat$mets==mets_4_fig[i]][2]
  max[length(max)+1] = dat$max[dat$mets==mets_4_fig[i]][1]
  max[length(max)+1] = dat$max[dat$mets==mets_4_fig[i]][2]
  mean[length(mean)+1] = dat$mean[dat$mets==mets_4_fig[i]][1]
  mean[length(mean)+1] = dat$mean[dat$mets==mets_4_fig[i]][2]
  Charge[length(Charge)+1] = dat$charge[dat$mets==mets_4_fig[i]][1]
  Charge[length(Charge)+1] = dat$charge[dat$mets==mets_4_fig[i]][2]
}
mets = mets[-1]
mets_name = mets_name[-1]
min = min[-1]
max = max[-1]
mean = mean[-1]
Charge = Charge[-1]
figure = data.frame(mets,mets_name,min,max,mean,Charge)
library(ggplot2)
pdf(file = "C:/Users/alexa/OneDrive - University of Toronto/Documents/PhD_UofT/Paper/Osm/Figure/Figure5b.pdf",
    width = 12, height = 7)
par(mar = c(7, 4, 1, 4) + 0.3)
pd = position_dodge(.35)    ### How much to jitter the points on the plot
ggplot(figure,                ### The data frame to use.
       aes(x     = mets_name,
           y     = mean*1000,
           color = Charge)) +
  geom_point(shape = 15,
             size  = 4,
             position = pd) +
  geom_errorbar(aes(ymin  = min*1000,
                    ymax  = max*1000),
                width = 0.2,
                size  = 0.7,
                position = pd) +
  ylab("Concentration (mM)") +
  scale_color_manual(values = c("red", "black"))+
  xlab("Metabolites") +
  labs(fill='Growth rate')+
  theme_classic(base_size = 23)

ylab("Mean steps")
dev.off()