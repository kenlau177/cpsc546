setwd("C://Users//Ken//Desktop//School//MSC1_Winter//cpsc546//Project//code")

rm(list=ls())

source("predict2.rpart.R")

library(mlbench)
library(rpart)
library(optimx)
library(plyr)
library(gbm)
library(snow)
library(snowfall)
library(parallel)
library(reshape2)
data(Ozone)

#sfInit(parallel=T, cpus=5, type="SOCK", socketHosts=list("localhost", 
#  "localhost", "localhost", "localhost", "localhost"))
#sfInit(parallel=T, cpus=2)
#sfLibrary(snowfall)

process_Ozone = function(x){
  out = x[,-9]
  names(out) = c("Month", "Month_Day", "Week_Day", "Max_O3",
                "Pressure_Height", "Wind_Speed", "Humidity", 
                "Temperature", "Inversion_Base_Height", 
                "Pressure_Gradient", "Inversion_Base_Temperature", 
                "Visibility")
  temp = out[,-match("Max_O3",colnames(out))]
  out = data.frame("Max_O3"=out$Max_O3, temp)  
  out = out[complete.cases(out),]
  return(out)
}
split_train_test = function(x){
  id_q = with(x, cut(Max_O3, breaks=c(0,quantile(Max_O3, probs=.75), 
          max(Max_O3)), labels=c("low","up"), include.lowest=T))
  id_low = which(id_q=="low")
  id_up = which(id_q=="up")
  
  num_low_trn = round(length(id_low)*.75)
  id_low_trn = sample(id_low, size=num_low_trn)
  id_low_test = id_low[!id_low%in%id_low_trn]
  
  num_up_trn = round(length(id_up)*.75)
  id_up_trn = sample(id_up, size=num_up_trn)
  id_up_test = id_up[!id_up%in%id_up_trn]

  train = x[c(id_low_trn,id_up_trn),]
  test = x[c(id_low_test,id_up_test),]
  
  return(list("train"=train, "test"=test))
}
loss.sqep = function(x, ep){
  out = pmax(0, x^2-ep)
  return(out)
}
loss.ls = function(x){
  out = x^2
  return(out)
}
rss.sqep = function(rho, y, ep){
  out = sum(loss.sqep((y-rho), ep=ep))
  return(out)
}
rss.ls = function(rho, y){
  out = sum(loss.ls((y-rho)))
  return(out)
}
neg_grad.sqep = function(y, fj, ep){
 #indc = as.numeric(((y-fj)^2-ep)>0)
 #out = (2*(y-fj))*indc
  D = diag(as.numeric(((y-fj)^2-ep)>0))
  out = 2*c(D%*%(y-fj))
  return(out)
}
neg_grad.ls = function(y, fj){
  out = 2*(y-fj)
  return(out)
}
grad.rss.sqep = function(y, fj, ep){
  D = diag(as.numeric(((y-fj)^2-ep)>0))
  out = 2*c(D%*%(y-fj))
  return(out)
}
hess.rss.sqep = function(y, fj, ep){
  D = diag(as.numeric(((y-fj)^2-ep)>0))
  out = 2*D
  return(out)
}
subsample = function(y, x, zj){
  sub_size = round(nrow(x)*.5)
  isub = sample.int(nrow(x), size=sub_size)
  out.x = x[isub,]
  out.y = y[isub]
  out.zj = zj[isub]
  return(list("x"=out.x, "y"=out.y, "zj"=out.zj))
}
brt_update = function(y, x, fj, fj.tst, ep=NULL, rpart.control, y.tst, x.tst, shrinkage){
  if(is.null(ep)){
    zj = neg_grad.ls(y=y, fj=fj)
  } else {
    zj = neg_grad.sqep(y=y, fj=fj, ep=ep)
  }
  out_subsample = subsample(y=y, x=x, zj=zj)
  ysub = out_subsample$y
  xsub = out_subsample$x
  zsub = out_subsample$zj
  
  zdf = data.frame("z"=zsub, xsub)
  out_rpart = rpart(z~., zdf, method="anova", y=T, 
                control=rpart.control)
  
  if(is.null(ep)){
    rhoj = predict(out_rpart, newdata=x, type="vector")
    rhoj.tst = predict(out_rpart, newdata=x.tst, type="vector")
  } else {
    rhoj = predict2.rpart(out_rpart, newdata=x)
    rhoj.tst = predict2.rpart(out_rpart, newdata=x.tst)
  }
  
  fj = fj + shrinkage*rhoj
  fj.tst = fj.tst + shrinkage*rhoj.tst
  return(list("fj"=fj, "fj.tst"=fj.tst))
}
cverr.calc = function(fj.tst, y.tst, ep=NULL){
  if(is.null(ep)){
    out = rss.ls(fj.tst, y.tst)
  } else {
    out = rss.sqep(fj.tst, y.tst, ep=ep)
  }
  return(out)
}
brt_fit = function(ifold=NULL, y.trn, x.trn, ep=NULL, rpart.control, y.tst=NULL, 
                    x.tst=NULL, shrinkage, n.trees){
  if(!is.null(ifold)){
    cverr = numeric(n.trees)
    ifold_trn = which(!1:nrow(x.trn)%in%ifold)
    y.tst = y.trn[ifold]
    x.tst = x.trn[ifold,]
    y.trn = y.trn[ifold_trn]
    x.trn = x.trn[ifold_trn,]
  } else {
    cverr = NULL
  }
  
  rhoj = mean(y.trn)
  fj = rhoj
  rhoj.tst = mean(y.tst)
  fj.tst = rhoj.tst

  for(j in 1:n.trees){
    out_brt_update = brt_update(y.trn, x.trn, fj, fj.tst, ep, 
                      rpart.control, y.tst, x.tst, shrinkage)            
    fj = out_brt_update$fj
    fj.tst = out_brt_update$fj.tst
    if(!is.null(ifold)){
      cverr[j] = cverr.calc(fj.tst, y.tst, ep=ep)
    }
    print(j)
  }
  if(!is.null(ifold)){
    return(cverr)
  } else {
    return(list("fj"=fj, "fj.tst"=fj.tst))
  }
}

dat = process_Ozone(x=Ozone)
out_split_train_test = split_train_test(x=dat)
dat_trn = out_split_train_test$train
dat_tst = out_split_train_test$test
dat_trn_x = dat_trn[,-match("Max_O3",names(dat_trn))]
dat_trn_y = dat_trn[,match("Max_O3",names(dat_trn))]
dat_tst_x = dat_tst[,-match("Max_O3",names(dat_tst))]
dat_tst_y = dat_tst[,match("Max_O3",names(dat_tst))]
ep = 12
n.trees = 1400
interaction.depth = 4
shrinkage = .002
cv.folds = 5
#bag.fraction = .5
minobsinnode = 10
rpart.control = list(minbucket=minobsinnode, cp=0,
                  maxcompete=0, maxsurrogate=0, usesurrogate=0, 
                  xval=0, maxdepth=interaction.depth)

idx = sample(1:nrow(dat_trn_x))
maxfold = ceiling(nrow(dat_trn_x)/cv.folds)
ifold.list = split(idx, ceiling(seq_along(idx)/maxfold))

cverr.list = lapply(ifold.list, FUN=brt_fit, y.trn=dat_trn_y, x.trn=dat_trn_x, 
              rpart.control=rpart.control, shrinkage=shrinkage, 
              n.trees=n.trees, ep=ep)
cverr.mat = do.call("cbind", cverr.list)
cverr.df = data.frame(cverr.mat)
cverr.df$cverr = apply(cverr.mat, FUN=mean, MARGIN=1)
n2.trees = which.min(cverr.df$cverr)

out_brt_fit = brt_fit(y.trn=dat_trn_y, x.trn=dat_trn_x, ep=ep, 
                rpart.control=rpart.control, y.tst=dat_tst_y, 
                x.tst=dat_tst_x, shrinkage=shrinkage, 
                n.trees=n2.trees)

fit_lm = lm(Max_O3~., data=dat_trn)
fj_tst_lm = predict(fit_lm, newdata=dat_tst_x)

outf = "C://Users//Ken//Desktop//School//MSC1_Winter//cpsc546//Project//Visualizations"
jpeg(paste0(outf, "//plot3.jpg"))
plot(dat_tst_y, fj_tst_lm, xlim=c(0,55), ylim=c(-15,55))
points(dat_tst_y, out_brt_fit$fj.tst, col="red")
lines(-5:55, -5:55, lty=3)
title(paste0("depth=",interaction.depth, ", trees=",n2.trees, 
  ", shrinkage=", shrinkage, ", ep=", ep))
dev.off()

plot(dat_trn_y, fit_lm$fitted.values, xlim=c(0,55), ylim=c(-15,55))
points(dat_trn_y, out_brt_fit$fj, col="red")
lines(-5:55, -5:55, lty=3)

a = brt_fit(ifold.list[[1]], y.trn=dat_trn_y, x.trn=dat_trn_x, 
      ep=ep, rpart.control, y.tst=NULL, x.tst=NULL, 
      shrinkage, n.trees)

#scatterplot training
gg1.df = data.frame("lm"=fit_lm$fitted.values, 
          "brt"=out_brt_fit$fj, "trn"=dat_trn_y)
gg1.df = melt(gg1.df, measure.vars=c("lm","brt"), variable.name="model", 
          value.name="fitted.values")
gg1 = ggplot(gg1.df, aes(x=trn, y=fitted.values)) + 
  geom_point(aes(colour=model), size=2) + 
  geom_abline(intercept = 0, slope = 1) +
  xlab("Observed") + ylab("Fitted Value") + 
  ggtitle("Observed vs Fitted Values") + 
  theme(axis.title=element_text(size=16), 
    axis.text=element_text(size=16), title=element_text(size=16))  
ggsave("scatter_train.pdf", gg1, width=6, height=6)  
  
#scatterplot test
gg2.df = data.frame("lm"=predict(fit_lm, newdata=dat_tst_x), 
          "brt"=out_brt_fit$fj.tst, "tst"=dat_tst_y)
gg2.df = melt(gg2.df, measure.vars=c("lm","brt"), variable.name="model",
          value.name="prediction")
gg2 = ggplot(gg2.df, aes(x=tst, y=prediction)) + 
  geom_point(aes(colour=model), size=2) + 
  geom_abline(intercept=0, slope=1) + 
  xlab("Observed") + ylab("Predicted") + 
  ggtitle("Observed vs Predicted") + 
  theme(axis.title=element_text(size=16), 
    axis.text=element_text(size=16), title=element_text(size=16))  
  ggsave("scatter_test.pdf", gg2, width=6, height=6)  







