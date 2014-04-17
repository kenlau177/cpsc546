setwd("C://Users//Ken//Desktop//School//MSC1_Winter//cpsc546//Project//code")

rm(list=ls()) #clear all R objects from environment

source("process_data_fn.R") 
source("loss_fn.R")
source("predict2.rpart.R")
source("verif_fn.R")
source("brt_fit.R")
source("subsample.R")

library(mlbench) #ozone data
library(rpart)  #regression tree
#library(optimx) #optimization
library(plyr) #vectorization functions by Hadley
#library(gbm)  #gbm function as implemented by Greg Ridgeway
library(reshape2) #data reshaping functions by Hadley
library(foreach)  #foreach function for parallization
library(doSNOW) #parallelization functions
library(ggplot2)
library(gridExtra)
data(Ozone) #initialize ozone data

dat = process_Ozone(x=Ozone)
#seeds to reproduce results
seeds = c(1, 264, 551, 834, 896)

s = 8
set.seed(seeds[s])

out_split_train_test = split_train_test(x=dat)
dat_trn = out_split_train_test$train
dat_tst = out_split_train_test$test
#training data for explantory variables
dat_trn_x = dat_trn[,-match("Max_O3",names(dat_trn))]
#training data for response variable
dat_trn_y = dat_trn[,match("Max_O3",names(dat_trn))]
#testing data for explantory variables
dat_tst_x = dat_tst[,-match("Max_O3",names(dat_tst))]
#testing data for response variable
dat_tst_y = dat_tst[,match("Max_O3",names(dat_tst))]

#number of folds for cross-validation
cv.folds = 5
#minimum obs in each terminal node of any regression tree fitted
minobsinnode = 10

### Stage 1 ###
#At this point, I want to fit the models once through to compute 
#the optimal number of trees I want to fit for each model.
#Cross-validation on mspe is used to pick out number of optimal trees

#sample sequence 1 to number of obs in traning data
idx = sample(1:nrow(dat_trn_x))
#round up number of points in each fold
maxfold = ceiling(nrow(dat_trn_x)/cv.folds)
#split into folds for cross-validation
ifold.list = split(idx, ceiling(seq_along(idx)/maxfold))

#Try different parameters for brt
#set up parameters for brt
#number of trees (this is conservative)
n.trees.vec = 1200
#epsilon for SEEI Loss
ep.vec = c(0, 4, 8, 12)
#shrinkage (or step-size)
shrinkage.vec = c(.002, .004)
#interaction depth (depth of each tree)
depth.vec = c(4, 5)

#generate combinations of each parameter given above
params.df = expand.grid(n.trees.vec, ep.vec, shrinkage.vec, depth.vec)
names(params.df) = c("n.trees", "ep", "shrinkage", "depth")
#create result dataframe with 2 additional columns 
#n2.trees is the optimal number of trees found through CV
#cverr is the cross-validation error obtained based on the 
#optimal number of trees fitted to the model
params_res.df = matrix(NA, nrow=nrow(params.df), ncol=(ncol(params.df)+2))
params_res.df = data.frame(params_res.df)
names(params_res.df) = c("n.trees", "ep", "shrinkage", "depth", 
                         "n2.trees", "cverr")

#set-up 5 cores
cl<-makeCluster(5)
registerDoSNOW(cl)

#loop through each combination of parameters input into brt
for(p in 1:nrow(params.df)){
  
  n.trees = params.df[p,"n.trees"]
  ep = params.df[p,"ep"]
  shrinkage = params.df[p,"shrinkage"]
  interaction.depth = params.df[p,"depth"]
  #control parameters for regression tree
  #cp is criterion to continue splitting
  #maxcompete is to keep track of next best split variable
  #maxsurrogate is surrogate variables to use for missing data
  #xval is number of CV to choose split depth
  #maxdepth is depth of tree
  #all these variables are set to 0 except for maxdepth is 
  #because thre tree should be grown fully for each iteration
  #of brt
  rpart.control = list(minbucket=minobsinnode, cp=0,
                  maxcompete=0, maxsurrogate=0, usesurrogate=0, 
                  xval=0, maxdepth=interaction.depth)
  #split the job of cross-validation in parallel for each fold
  #fit the brt for each fold using input parameters set above
  cverr.list = foreach(i=1:cv.folds, 
                .export=c("rpart", "ddply", ".")) %dopar% {
    brt_fit(ifold.list[[i]], y.trn=dat_trn_y, x.trn=dat_trn_x,
     ep=ep, rpart.control=rpart.control, shrinkage=shrinkage, 
      n.trees=n.trees)
  }
  #combine CV err output for each fold
  cverr.mat = do.call("cbind", cverr.list)
  cverr.df = data.frame(cverr.mat)
  #calculate average of CV error based on each iteration 
  #of brt on the cv folds
  cverr.df$cverr = apply(cverr.mat, FUN=mean, MARGIN=1)
  #find mininum CV error and the number of trees corresponding
  #to the minimum
  cverr = min(cverr.df$cverr)
  n2.trees = which.min(cverr.df$cverr)
  params_res.df[p,] = c(n.trees, ep, shrinkage, interaction.depth, 
                        n2.trees, cverr)
  
  print(p)

}

stopCluster(cl)


### Stage 2 ###
#After we obtained the optimal number of trees for each parameter
#combination, we can fit the full model to the training data
#The parallelization is on the parameters this time

cl<-makeCluster(4)
registerDoSNOW(cl)

#Choose optimal model for standard brt and SEEI boosting

i1 = which(with(params_res.df, ep==0))
i1 = i1[which.min(params_res.df$cverr[i1])]
i2 = which(with(params_res.df, ep!=0))
i2 = i2[which.min(params_res.df$cverr[i2])]
pr2.df = params_res.df[c(i1,i2),]

#make copy of parameter results from Stage 1
#pr2.df = params_res.df
out_brt_par = foreach(i=1:nrow(pr2.df), 
                .export=c("rpart", "ddply", ".")) %dopar% {
  rpart.control = list(minbucket=minobsinnode, cp=0,
                   maxcompete=0, maxsurrogate=0, usesurrogate=0, 
                   xval=0, maxdepth=pr2.df$depth[i])
  brt_fit(y.trn=dat_trn_y, x.trn=dat_trn_x, ep=pr2.df$ep[i], 
    rpart.control=rpart.control, y.tst=dat_tst_y, 
    x.tst=dat_tst_x, shrinkage=pr2.df$shrinkage[i], 
    n.trees=pr2.df$n2.trees[i], 
    interaction.depth=pr2.df$depth[i])
}

stopCluster(cl)

## Stage 3 ##
#Verification and Plotting
#############

#multiple linear regression fit on training data
fit_lm = lm(Max_O3~., data=dat_trn)
#prediction on testing data with mlr
pred_lm = predict(fit_lm, newdata=dat_tst_x)
#output directory for verification plots based on seed
outf.vis = paste0("C://Users//Ken//Desktop//School//MSC1_Winter//cpsc546//Project//Visualizations//ObsPred//Seed", 
            seeds[s])
pred_brt_ls = out_brt_par[[1]]$fj.tst
pred_brt_sqep = out_brt_par[[2]]$fj.tst
gg.df = data.frame("lm"=pred_lm, "brt_ls"=pred_brt_ls, 
          "brt_sqep"=pred_brt_sqep, "tst"=dat_tst_y)
gg.df = melt(gg.df, measure.vars=c("lm","brt_ls","brt_sqep"), 
          variable.name="model", value.name="prediction")

err.lm.up75 = rss.ls.quant(rho=pred_lm, y=dat_tst_y, quant=c(.75,1))
err.lm.low60 = rss.ls.quant(rho=pred_lm, y=dat_tst_y, quant=c(0,.6))
err.lm.all = rss.ls.quant(rho=pred_lm, y=dat_tst_y, quant=c(.1,.75))
err.lm.up75 = round(err.lm.up75, 2)
err.lm.low60 = round(err.lm.low60, 2)
err.lm.all = round(err.lm.all, 2)

err.brt1.up75 = rss.ls.quant(rho=pred_brt_ls, y=dat_tst_y, quant=c(.75,1))
err.brt1.low60 = rss.ls.quant(rho=pred_brt_ls, y=dat_tst_y, quant=c(0,.6))
err.brt1.all = rss.ls.quant(rho=pred_brt_ls, y=dat_tst_y, quant=c(.1,.75))
err.brt1.up75 = round(err.brt1.up75, 2)
err.brt1.low60 = round(err.brt1.low60, 2)
err.brt1.all = round(err.brt1.all, 2)

err.brt2.up75 = rss.ls.quant(rho=pred_brt_sqep, y=dat_tst_y, quant=c(.75,1))
err.brt2.low60 = rss.ls.quant(rho=pred_brt_sqep, y=dat_tst_y, quant=c(0,.6))
err.brt2.all = rss.ls.quant(rho=pred_brt_sqep, y=dat_tst_y, quant=c(.1,.75))
err.brt2.up75 = round(err.brt2.up75, 2)
err.brt2.low60 = round(err.brt2.low60, 2)
err.brt2.all = round(err.brt2.all, 2)

gg1 = ggplot(subset(gg.df, model %in% c("lm","brt_sqep")), 
  aes(x=tst, y=prediction)) + 
  geom_point(aes(colour=model), size=2) + 
  geom_abline(intercept=0, slope=1, linetype="dashed") + 
  xlab("Observed") + ylab("Predicted") + 
  xlim(c(0,40)) + ylim(c(-10,40)) +
  ggtitle(paste0("Observation vs Prediction of MLR \n and BRT with SEEI ep=",
    pr2.df$ep[2])) +  
  theme(axis.title=element_text(size=16), 
    axis.text=element_text(size=16), title=element_text(size=16), 
    legend.text=element_text(size=16))
ggsave(paste0(outf.vis, "//Obs_Pred_MLR_BRTSEEI_", err.lm.up75, "_",
  err.lm.low60, "_", err.lm.all, "_", err.brt2.up75, "_", 
  err.brt2.low60, "_", err.brt2.all, ".pdf"), gg1,
  width=6, height=6)

gg2 = ggplot(subset(gg.df, model %in% c("brt_ls","brt_sqep")), 
  aes(x=tst, y=prediction)) + 
  geom_point(aes(colour=model), size=2) + 
  geom_abline(intercept=0, slope=1) + 
  xlab("Observed") + ylab("Predicted") + 
  xlim(c(0,40)) + ylim(c(-10,40)) +
  ggtitle(paste0("Observation vs Prediction of BRT with LS \n and BRT with SEEI ep=", 
    pr2.df$ep[2])) +  
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=16), title=element_text(size=16), 
        legend.text=element_text(size=16))
ggsave(paste0(outf.vis, "//Obs_Pred_BRTLS_BRTSEEI_", err.brt1.up75, "_",
 err.brt1.low60, "_", err.brt1.all, "_", err.brt2.up75, "_", 
 err.brt2.low60, "_", err.brt2.all, ".pdf"), gg2,
 width=6, height=6)

pdf(paste0(outf.vis,"//scatter.pdf"), width=12, height=6, compress=F)
grid.arrange(gg1, gg2, ncol=2)
dev.off()


  





# ## Stage 3 ##
# #multiple linear regression fit on training data
# fit_lm = lm(Max_O3~., data=dat_trn)
# #prediction on testing data with mlr
# pred_lm = predict(fit_lm, newdata=dat_tst_x)
# 
# #output directory for verification plots based on seed
# outf.vis = paste0("C://Users//Ken//Desktop//School//MSC1_Winter//cpsc546//Project//Visualizations//ObsPred//Seed", 
#             seeds[s])
# for(r in 1:nrow(pr2.df)){
#   
#   out_brt_fit = out_brt_par[[r]]
#   pred_brt = out_brt_fit$fj.tst
#   
#   outf.pre = paste0(outf.vis, "//trees=",pr2.df$n2.trees[r], 
#               "_shrinkage=", pr2.df$shrinkage[r], "_ep=", pr2.df$ep[r], 
#               "_depth=", pr2.df$depth[r])
#   outf.gg = paste0(outf.pre, "_gg.jpg")
# 
#   err.30.lm = rss.ls.quant(rho=pred_lm, y=dat_tst_y, quant=.3)
#   err.30.brt = rss.ls.quant(rho=pred_brt, y=dat_tst_y, quant=.3)
#   #rss.ls.quant(rho=pred_lm, y=dat_tst_y, quant=0)
#   #rss.ls.quant(rho=pred_brt, y=dat_tst_y, quant=0)
#   err.30.lm = round(err.30.lm, 1)
#   err.30.brt = round(err.30.brt, 1)
#   
#   gg.df = data.frame("lm"=pred_lm, 
#                       "brt"=pred_brt, "tst"=dat_tst_y)
#   gg.df = melt(gg.df, measure.vars=c("lm","brt"), variable.name="model",
#                 value.name="prediction")
#   gg = ggplot(gg.df, aes(x=tst, y=prediction)) + 
#     geom_point(aes(colour=model), size=2) + 
#     geom_abline(intercept=0, slope=1) + 
#     xlab("Observed") + ylab("Predicted") + 
#     xlim(c(0,45)) + ylim(c(-15,45)) +
#     ggtitle(paste0("trees=",pr2.df$n2.trees[r], 
#       ", shrinkage=", pr2.df$shrinkage[r], ", ep=", pr2.df$ep[r], 
#       ", depth=", pr2.df$depth[r], ",\n",  
#       "err.lm=", err.30.lm, ", err.brt=", err.30.brt)) +  
#     theme(axis.title=element_text(size=16), 
#           axis.text=element_text(size=16), title=element_text(size=16))
#   
#   ggsave(outf.gg, gg, width=6.5, height=6)  
#   
# }
  