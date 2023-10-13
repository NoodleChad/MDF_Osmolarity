#### Figure II a) ####
setwd("./Results")
dat<-read.csv("Figure2a.csv",sep=",",dec=".",h=T)

# Draw first plot using axis y1
pdf(file = "./Figures/Figure2a.pdf", 
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

#### Figure II b) ####
library(ggplot2)
library(ggnewscale)
setwd("C:/Users/alexa/OneDrive - University of Toronto/Documents/PhD_UofT/Paper/Osm/Analysis/IV/CVA_with_k")
dat<-read.csv("FigureIIb.csv",sep=",",dec=".",h=T)

growth<-read.csv("Figure2b.csv",sep=",",dec=".",h=T)
pdf(file = "./Figures/Figure2b.pdf",
    width = 20, height = 7)
par(mar = c(7, 4, 1, 4) + 0.3)
pd = position_dodge(.75)

#data_breaks <- data.frame(start = c(0.0023,0.0035,0.0095, 0.0145, 0.0205),  # Create data with breaks
                          #end = c(0.0031,0.0045, 0.0115, 0.0195, 0.0255),
                          #colors = factor(1:5))
ggplot() +                                        # Add background colors to plot
  #geom_rect(data = data_breaks,
   #         aes(xmin = start,
    #            xmax = end,
     #           ymin = - Inf,
      #          ymax = Inf,
       #         fill = colors),
        #    alpha = 0.5) +
  geom_point(data = dat, size = 2.5, aes(k_val*1000, mean, color = Metabolites, shape=Metabolites), position = pd)+
  geom_errorbar(data = dat, size = 0.5, aes(x = k_val*1000, ymin  = min, ymax  = max, color = Metabolites), 
                width = 1, position = pd) +
  #scale_fill_manual(values = c("grey80","grey65", "grey50", "grey35","grey15"),labels=c("0.7432","0.8727","0.8966","0.9008","0.9024")) +
  xlab("K (mOsm)") +
  scale_x_continuous(breaks=seq(0,30,1),limits=c(0, 31))+
  #labs(fill='Growth rate')+
  #scale_y_continuous("Concentration (mol/L)",sec.axis = sec_axis(~ . * 100, name = "Growth Rate (1/h)"))+
  scale_y_continuous("Concentration (mol/L)")+
  scale_color_manual(values = c("Acetyl-CoA"="red", "Acetate"="blue",
                                "ATP"="darkgreen","NAD"="black","Growth Rate"="black"))+
  #new_scale_colour() +
  #geom_line(aes(color = "Growth Rate",growth$K[1:2745]*1000, growth$max_growth[1:2745]/100),lwd=1.5) +
  #scale_colour_discrete(labels = "Growth Rate",name = "",type = "black")+
  theme_classic(base_size = 23)
dev.off()

### Figure II c) ###
setwd("./Results/CVA/")
dat<-read.csv("CVA.csv",sep=",",dec=".",h=T)
mets.name<-c("3pg","accoa","asp__L","atp","dhap","g3p","glu__L",'imp',"mal__L","nad","succoa")

library(plotrix)
par(mar=c(5,4,4,2))
# Each values for all values of k
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
ions<-c(1,1,1,1,1,1,1,1,1,1,1)
ion.names<-c("3pg","Acetyl-CoA","Aspartate","ATP","DHAP","g3p","Glutamate","IMP","Malate","NAD","Succinyl-CoA")
pdf(file = "./Figures/Figure2c.pdf", width = 8.27, height = 5.83)
par(mar = c(7, 4, 1, 4) + 0.3)
radial.plot(ions,labels=ion.names,main="",
            grid.unit="meq/l",radial.lim=c(0,1),lwd=3,line.col=1,
            show.grid.labels=0)
# add points inside the polygon - radial.lim is supplied by plotrix_env
radial.plot(val_8966_max,rp.type="r",lwd=3, line.col="forestgreen",add=TRUE,cex=1.5)
radial.plot(val_8966,rp.type="r",lwd=3, line.col=1,add=TRUE,cex=1.5)

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
legend(1.1,0.3,title="Growth rate (1/h)",legend=c("0.9028","0.9008","0.8966","0.8727","0.7432"),
       pch = c(16,17,15,3,18),col=c("red","blue","forestgreen","black","yellow"),cex=1,bty="n")
dev.off()
