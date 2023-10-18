#### Figure IV  ####
setwd("./Experimental/")
dat<-read.csv("Max_growth_output.csv",sep=",",dec=".",h=T)
dat1<-read.csv("./Results/analysis_II_fitted_kcl.csv",sep=",",dec=".",h=T)
dat2<-read.csv("./Results/analysis_II_fitted_nacl.csv",sep=",",dec=".",h=T)
dat1$k_val = dat1$k_val*1000
dat2$k_val = dat2$k_val*1000
pdf(file = "./Figures/Figure4.pdf",
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
    dat_nacl<-read.csv("./Experimental/NaCl.csv",sep=",",dec=".",h=T)
    # linear
      linear_salt_effect.lm <- lm(dat_nacl$mu_norm ~ dat_nacl$NaCl)
      summary(linear_salt_effect.lm)
    # quadratic
      quad_salt_effect.lm <- lm(dat_nacl$mu[31:length(dat_nacl$mu)] ~ poly(dat_nacl$NaCl[31:length(dat_nacl$mu)], 2, raw=TRUE))
      summary(quad_salt_effect.lm)
    # log all
      log_salt_effect.lm <- lm(dat_nacl$mu[1:24] ~ log(dat_nacl$NaCl)[1:24])
      summary(log_salt_effect.lm)
    # linear bottom
      lin_salt_effect.lm <- lm(dat_nacl$mu[31:length(dat_nacl$mu)] ~ dat_nacl$NaCl[31:length(dat_nacl$mu)])
      summary(lin_salt_effect.lm)
      cor.test(dat_nacl$NaCl[31:length(dat_nacl$mu)], dat_nacl$mu[31:length(dat_nacl$mu)],method = 'spearman')
      cor.test(dat_nacl$NaCl[31:length(dat_nacl$mu)], dat_nacl$mu[31:length(dat_nacl$mu)],method = 'pearson')
    # AOV
      bottom_nacl <- aov(mu[1:24] ~ NaCl[1:24], data = dat_nacl)
      summary(bottom_nacl)
      top_nacl <- aov(mu[25:length(dat_nacl$mu)] ~ NaCl[25:length(dat_nacl$mu)], data = dat_nacl)
      summary(top_nacl)
    # Spearman
      cor.test(dat_nacl$NaCl[1:30], dat_nacl$mu_norm[1:30], method = 'spearman')
      cor.test(dat_nacl$NaCl[1:30], dat_nacl$mu_norm[1:30], method = 'pearson')
    # Distance correlation
      dcor(dat_nacl$NaCl, dat_nacl$mu_norm,index=1)
      
  # KCl
    dat_kcl<-read.csv("./Experimental/KCl.csv",sep=",",dec=".",h=T)
    # Quadratic
      quad_salt_effect.lm <- lm(dat_kcl$mu_norm ~ poly(dat_kcl$KCl, 2, raw=TRUE))
      summary(quad_salt_effect.lm)
    # Log
      log_salt_effect.lm <- lm(dat_kcl$mu[1:24] ~ log(dat_kcl$KCl[1:24]))
      summary(log_salt_effect.lm)
     # linear bottom
      lin_salt_effect.lm <- lm(dat_kcl$mu[31:length(dat_kcl$mu)] ~ dat_kcl$KCl[31:length(dat_kcl$mu)])
      summary(lin_salt_effect.lm)
      cor.test(dat_kcl$KCl[31:length(dat_kcl$mu)], dat_kcl$mu[31:length(dat_kcl$mu)],method = 'spearman')
      cor.test(dat_kcl$KCl[31:length(dat_kcl$mu)], dat_kcl$mu[31:length(dat_kcl$mu)],method = 'pearson')
    # Spearman
      cor.test(dat_kcl$KCl[1:30], dat_kcl$mu_norm[1:30], method = 'spearman')
    # Distance correlation
      dcor(dat_kcl$KCl, dat_kcl$mu_norm,index=1)


