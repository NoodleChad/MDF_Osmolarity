#### FigureS1
dat<-read.csv("./Results/FigureS1.csv",sep=",",dec=".",h=T)
dat$ln_met = dat$ln_met*1000
dat$exp_met = dat$exp_met*1000
pdf(file = "./Figures/FigureS1.pdf",
    width = 12, height = 8)
par(mfrow = c(2, 2),mar = c(5.1, 5.1, 4.1, 2.1))
plot(dat$ln_met[dat$k_val==0.0026],dat$exp_met[dat$k_val==0.0026],xlab = expression(paste(e^{'Xj'},' (mM)')),ylab=expression(paste(Cj^{'*'},' (mM)')),
     pch = 16,cex = 1.25,cex.axis=2, cex.lab = 1.75, main="K = 0.0026 M",cex.main = 2)
plot(dat$ln_met[dat$k_val==0.02],dat$exp_met[dat$k_val==0.02],xlab = expression(paste(e^{'Xj'},' (mM)')),ylab=expression(paste(Cj^{'*'},' (mM)')),
     pch = 16,col = 'red',cex = 1.25,cex.axis=2, cex.lab = 1.75, main="K = 0.02 M",cex.main = 2)
plot(dat$ln_met[dat$k_val==0.2],dat$exp_met[dat$k_val==0.2],xlab = expression(paste(e^{'Xj'},' (mM)')),ylab=expression(paste(Cj^{'*'},' (mM)')),
     pch = 16,col = 'blue',cex = 1.25,cex.axis=2, cex.lab = 1.75, main="K = 0.2 M",cex.main = 2)
plot(dat$ln_met[dat$k_val==2],dat$exp_met[dat$k_val==2],xlab = expression(paste(e^{'Xj'},' (mM)')),ylab=expression(paste(Cj^{'*'},' (mM)')),
     pch = 16,col = 'forestgreen',cex = 1.25,cex.axis=2, cex.lab = 1.75, main="K = 2 M",cex.main = 2)
dev.off()