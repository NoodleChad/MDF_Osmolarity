#### FigureS1 - With charge ####
setwd("./Results/CVA/")
dat<-read.csv("CVA_charge.csv",sep=",",dec=".",h=T)
dat1<-read.csv("CVA.csv",sep=",",dec=".",h=T)

# Draw first plot using axis y1
pdf(file = "./Figures/FigureS1.pdf",
    width = 8.27, height = 5.83)
par(mar = c(7, 4, 1, 4) + 0.3)  
plot(dat$K, dat$max_growth, pch = 16, col = 1, 
     xlab = "K (mOsm)",
     ylab = "Growth rate (1/h)",cex.axis=1.5, cex.lab = 1.5) 
points(dat1$K, dat1$max_growth,pch=18, col = "red")
legend(x=0.016, y=0.8,legend = c("With Charge Balance", "Without Charge Balance"), pch = c(16,18),col = c("black","red"), cex=1.5)
dev.off()