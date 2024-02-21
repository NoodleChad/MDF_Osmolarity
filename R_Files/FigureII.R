#### Figure II b) ####
setwd("./Results")
dat<-read.csv("Figure2b.csv",sep=",",dec=".",h=T)

# Draw first plot using axis y1
pdf(file = "./Figures/Figure2b.pdf", 
    width = 18, height = 7)
par(mar = c(7, 4, 1, 4) + 0.3)  
plot(dat$k_val*1000, dat$max_growth, pch = 16, col = 1, xlab = "K (mOsm)", xlim=c(0,30),
     ylab = "Growth rate (1/h)",cex.axis=1.5, cex.lab = 1.5, xaxt = 'n',) 
axis(1, at = seq(0, 30, by = 1), las=1,cex.axis = 1.5)
# set parameter new=True for a new axis
par(new = TRUE)         

# Draw second plot using axis y2
plot(dat$k_val*1000, dat$num_open, pch = 18, col = 4, axes = FALSE, xlab = "", 
     ylab = "",ylim = c(985,1015),xlim=c(0,30))


axis(side = 4, at = pretty(range(dat$num_open[1:length(dat$k_val)])),cex.axis=1.35)      
mtext("Total number of feasible reactions", side = 4, line = 3,cex=1.5)
legend(x=24, y=1010,legend = c("Growth rate", "Active reactions"), pch = c(16,18),col = c(1,4),cex = 1.5)
dev.off()

#### Figure II c) ####
dat100<-read.csv("Figure_2c_100.csv",sep=",",dec=".",h=T)
dat200<-read.csv("Figure_2c_200.csv",sep=",",dec=".",h=T)
dat500<-read.csv("Figure_2c_500.csv",sep=",",dec=".",h=T)
dat1000<-read.csv("Figure_2c_1000.csv",sep=",",dec=".",h=T)
pdf(file = "./Figures/Figure_2c.pdf", 
    width = 12, height = 8)
par(mar = c(5.1, 5.1, 4.1, 2.1)) 
plot(dat100$k_val*1000, dat100$max_growth, pch = 18, col = 1,ylim=c(0,1),
     cex.axis=2, cex.lab = 2,xlab = expression(paste(Phi, " (mmol/L)")),cex=1.5, 
     ylab = expression(paste("Growth (h"^"-1",")")))
points(dat200$k_val*1000, dat200$max_growth, pch = 15, col = 2,cex=1.5)
points(dat500$k_val*1000, dat500$max_growth, pch = 17, col = 3,cex=1.5)
points(dat1000$k_val*1000, dat1000$max_growth, pch = 19, col = 4,cex=1.5)
abline(v=200, col="blue",lwd = 2)
text(x=225,y=0.95,expression(paste(Phi["threshold"])),col="blue",cex=2)
legend(x=205, y=0.35,legend = c(expression(paste(psi, " = 100")), 
                              expression(paste(psi, " = 200")),
                              expression(paste(psi, " = 500")),
                              expression(paste(psi, " = 1000"))),
       pch = c(18,15,17,19),col = c(1,2,3,4),cex=2)
dev.off()

#### Figure II d) ####
library(ggplot2)
library(ggnewscale)
dat<-read.csv("Figure_2d.csv",sep=",",dec=".",h=T)
pdf(file = "./Figures/Figure2d.pdf",
    width = 20, height = 7)
par(mar = c(5.1, 5.1, 4.1, 5.1))  
pd = position_dodge(.55)
ggplot() +                              
  geom_point(data = dat, size = 2.5, aes(k_val*1000, mean, color = mets, shape=mets), position = pd,show.legend = FALSE)+
  geom_errorbar(data = dat, size = 0.5, aes(x = k_val*1000, ymin  = min, ymax  = max, color = mets), 
                width = 1, position = pd) +
  xlab(expression(paste(Phi, " (mmol/L)"))) +
  scale_x_continuous(breaks=seq(145,171,2),limits=c(145, 171))+
  scale_y_continuous("Concentration (mmol/L)")+
  scale_color_manual(values = c("red", "blue",
                                "darkgreen","black"))+
  guides(fill="none", color=FALSE)+

  theme_classic(base_size = 26)
dev.off()


### Figure II e) ###
dat<-read.csv("Figure_2e.csv",sep=",",dec=".",h=T)
mets.name<-c("3pg","accoa","asp__L","atp","dhap","g3p","glu__L",'imp',"mal__L","nad")

library(plotrix)
par(mar=c(5,4,4,2))
# Each values for all values of Phi
same_in_all = "bacon"
val_902 = 1
val_9008 = 1
val_8966 = 1
val_8727 = 1
val_7432 = 1
val_902_max = 1
val_9008_max = 1
val_8966_max = 1
val_8727_max = 1
val_7432_max = 1
for(i in 1:length(mets.name)){
    val_902[i] = (dat$min[dat$mets==mets.name[i]][1]-min(dat$min))/(max(dat$min[dat$mets==mets.name[i]])-min(dat$min))
    val_9008[i] = (dat$min[dat$mets==mets.name[i]][2]-min(dat$min))/(max(dat$min[dat$mets==mets.name[i]])-min(dat$min))
    val_8966[i] = (dat$min[dat$mets==mets.name[i]][3]-min(dat$min))/(max(dat$min[dat$mets==mets.name[i]])-min(dat$min))
    val_8727[i] = (dat$min[dat$mets==mets.name[i]][4]-min(dat$min))/(max(dat$min[dat$mets==mets.name[i]])-min(dat$min))
    val_7432[i] = (dat$min[dat$mets==mets.name[i]][5]-min(dat$min))/(max(dat$min[dat$mets==mets.name[i]])-min(dat$min))
    val_902_max[i] = (dat$max[dat$mets==mets.name[i]][1]-min(dat$max))/(max(dat$max[dat$mets==mets.name[i]])-min(dat$max))
    val_9008_max[i] = (dat$max[dat$mets==mets.name[i]][2]-min(dat$max))/(max(dat$max[dat$mets==mets.name[i]])-min(dat$max))
    val_8966_max[i] = (dat$max[dat$mets==mets.name[i]][3]-min(dat$max))/(max(dat$max[dat$mets==mets.name[i]])-min(dat$max))
    val_8727_max[i] = (dat$max[dat$mets==mets.name[i]][4]-min(dat$max))/(max(dat$max[dat$mets==mets.name[i]])-min(dat$max))
    val_7432_max[i] = (dat$max[dat$mets==mets.name[i]][5]-min(dat$max))/(max(dat$max[dat$mets==mets.name[i]])-min(dat$max))
}
ions<-c(1,1,1,1,1,1,1,1,1,1)
ion.names<-c("3pg","Acetyl-CoA","Aspartate","ATP","DHAP","g3p","Glutamate","IMP","Malate","NAD")
pdf(file = "./Figures/Figure2e.pdf",
    width = 8.27, height = 5.83)
par(mar = c(7, 4, 1, 4) + 0.3)
radial.plot(ions,labels=ion.names,main="",
            grid.unit="meq/l",radial.lim=c(0,1),lwd=3,line.col=1,
            show.grid.labels=0)
# add points inside the polygon - radial.lim is supplied by plotrix_env
radial.plot(val_8966_max,rp.type="r",lwd=3, line.col="forestgreen",add=TRUE,cex=1.5)
radial.plot(val_8966,rp.type="r",lwd=3, line.col=1,add=TRUE,cex=1.5)
radial.plot(val_9008,rp.type="r",lwd=3, line.col="blue",add=TRUE,cex=1.5)
radial.plot(val_9008_max,rp.type="r",lwd=3, line.col=1,add=TRUE,cex=1.5)
radial.plot(val_902,rp.type="r",lwd=3, line.col="red",add=TRUE,cex=1.5)
radial.plot(val_902_max,rp.type="r",lwd=3, line.col=1,add=TRUE,cex=1.5)

radial.plot(val_8966_max,rp.type="s",point.symbols=15,point.col="forestgreen",add=TRUE,cex=1.5)
radial.plot(val_8966,rp.type="s",point.symbols=15,point.col="forestgreen",add=TRUE,cex=1.5)
radial.plot(val_902,rp.type="s",point.symbols=16,point.col="red",add=TRUE,cex=1.5)
radial.plot(val_902_max,rp.type="s",point.symbols=16,point.col="red",add=TRUE,cex=1.5)
radial.plot(val_9008,rp.type="s",point.symbols=17,point.col="blue",add=TRUE,cex=1.5)
radial.plot(val_9008_max,rp.type="s",point.symbols=17,point.col="blue",add=TRUE,cex=1.5)
radial.plot(val_8727,rp.type="s",point.symbols=3,point.col="black",add=TRUE,cex=1.5)
radial.plot(val_8727_max,rp.type="s",point.symbols=3,point.col="black",add=TRUE,cex=1.5)
radial.plot(val_7432,rp.type="s",point.symbols=18,point.col="yellow",add=TRUE,cex=1.5)
radial.plot(val_7432_max,rp.type="s",point.symbols=18,point.col="yellow",add=TRUE,cex=1.5)
legend(1.1,0.3,title="Growth rate (1/h)",legend=c("0.8770","0.8751","0.8478","0.8045","0.7174"),
       pch = c(16,17,15,3,18),col=c("red","blue","forestgreen","black","yellow"),cex=1,bty="n")
dev.off()