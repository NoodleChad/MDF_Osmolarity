## Figure S2 ##
dat<-read.csv("bioreactor_2.csv",sep=",",dec=".",h=T)
pdf(file = "C:/Users/alexa/OneDrive - University of Toronto/Documents/PhD_UofT/Paper/Osm/Figure/FigureS2.pdf",
    width = 10, height = 6)
plot(dat$time,dat$Bio5_1,xlab="Time (h)",ylab="OD600nm", pch =17,,type="b",
     main="",cex = 1.25,cex.axis=1.5, cex.lab = 1.5,col="black",ylim=c(0,0.55))
points(dat$time,dat$Bio5_2,pch=18,col="black",cex=1.5,type="b")
points(dat$time,dat$Bio1_1,pch=17,col="red",cex=1.5,type="b")
points(dat$time,dat$Bio1_2,pch=18,col="red",cex=1.5,type="b")
points(dat$time,dat$Bio2_1,pch=17,col="forestgreen",cex=1.5,type="b")
points(dat$time,dat$Bio2_2,pch=18,col="forestgreen",cex=1.5,type="b")
points(dat$time,dat$Bio3_1,pch=17,col="#1b98e0",cex=1.5,type="b")
points(dat$time,dat$Bio3_2,pch=18,col="#1b98e0",cex=1.5,type="b")
points(dat$time,dat$Bio4_1,pch=17,col="blueviolet",cex=1.5,type="b")
points(dat$time,dat$Bio4_2,pch=18,col="blueviolet",cex=1.5,type="b")
legend("topleft", inset=c(0.025,0.025),
       legend = c("1","25", "50","100","400"),
       title = "NaCl concentration (mM)",
       fill=c("red","forestgreen","#1b98e0","blueviolet","black"),cex=1.5,bty = "n")
dev.off()