setwd("C://Users//Ken//Desktop//School//MSC1_Winter//cpsc546//Project")

library(mlbench)
data(Ozone)

dat1 = airquality[complete.cases(airquality),]
fit1 = lm(Ozone~., dat1)
            
dat2 = Ozone[,-9]            
dat2 = dat2[complete.cases(dat2),]
fit2 = lm(V4~., data=dat2)            

par(mfrow=c(2,1))
plot(dat1$Ozone, fit1$fitted.values, xlim=c(0,170), ylim=c(-30,120))
abline(a=0, b=1)
plot(dat2$V4, fit2$fitted.values, xlim=c(0,50), ylim=c(-8,40))
abline(a=0, b=1)


