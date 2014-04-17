
#update estimates of current fitted values by fitting an additional tree
#Fits the negative gradient of loss function
#performs subsampling step in brt algorithm
#adjusts terminal node estimates of fitted regression tree
#
#There are two scenarios again:
#1. if ep is null then fit quartic loss function
#2. if ep is not null then fit SEEI loss
#
#input: similar to brt_fit
#output: updated estimates for training and testing
brt_update = function(y.trn, x.trn, fj, fj.tst, ep=NULL, rpart.control, y.tst, x.tst, shrinkage){
  if(is.null(ep) | is.na(ep)){
    #compute negative gradients on quartic loss function
    zj = neg_grad.quar(y=y.trn, fj=fj)
    #adjust shrinkage
    shrinkage = shrinkage^2
  } else {
    #compute negative gradients on SEEI loss function
    zj = neg_grad.sqep(y=y.trn, fj=fj, ep=ep)
  }
  #perform subsampling step
  out_subsample = subsample(y=y.trn, x=x.trn, zj=zj)
  ysub = out_subsample$y
  xsub = out_subsample$x
  zsub = out_subsample$zj
  zdf = data.frame("z"=zsub, xsub)
  
  #fit a regression tree on subsampled negative gradient
  out_rpart = rpart(z~., zdf, method="anova", y=T, 
                    control=rpart.control)
  
  #adjust terminal node estimates in regression tree and return
  #regression tree fit
  out_rpart_update = predict2.rpart(out_rpart, newdata=x.trn, ep=ep)
  #training estimates with adjusted terminal node estimates
  rhoj = predict(out_rpart_update, newdata=x.trn, type="vector")
  #rhoj = predict2.rpart(out_rpart, newdata=x.trn, ep=ep)
  rhoj.tst = predict(out_rpart_update, newdata=x.tst, type="vector")
    
#   if(is.null(ep) | is.na(ep)){
#     out_rpart_update = predict2.rpart(out_rpart, newdata=x.trn, ep=ep)
#     rhoj = predict(out_rpart_update, newdata=x.trn)
#     rhoj.tst = predict(out_rpart_update, newdata=x.tst)
#     #rhoj = predict2.rpart(out_rpart, newdata=x.trn, ep=ep)
#     #rhoj.tst = predict2.rpart(out_rpart, newdata=x.tst, ep=ep)
#   } else {
#     out_rpart_update = predict2.rpart(out_rpart, newdata=x.trn, ep=ep)
#     rhoj = predict(out_rpart_update, newdata=x.trn)
#     rhoj.tst = predict(out_rpart_update, newdata=x.tst)
#     #rhoj = predict2.rpart(out_rpart, newdata=x.trn, ep=ep)
#     #rhoj.tst = predict2.rpart(out_rpart, newdata=x.tst, ep=ep)
#   }
  
  #update training estimates
  fj = fj + shrinkage*rhoj
  #update testing estimates
  fj.tst = fj.tst + shrinkage*rhoj.tst
  return(list("fj"=fj, "fj.tst"=fj.tst))
}
#
cverr.calc = function(fj.tst, y.tst, ep=NULL){
  #if(cvls){
  out = rss.ls(fj.tst, y.tst)
  #} else {
  #  out = rss.sqep(fj.tst, y.tst, ep=ep)
  #}
  return(out)
}
#fit boosted regression tree. There are two scenarios.
#1. ifold is not null, then we are interested in using CV
#to calculate the minimal CV error to optimize number of trees
#2. ifold is null, then we y.tst and x.tst is also null, and 
#we fit model to all the data.

#input:
#ifold: fold left out to calculate CV error
#y.trn: response variable on training data
#x.trn: explantory variables on training data
#ep: epsilon parameter for SEEI
#rpart.control: regression tree control parameters for rpart
#y.tst: response variable on test data
#x.tst: explantory variables on test data
#shrinkage, n.trees, int.depth: 
#
#output:
#based on the two scenarios
#1. CV error calculate at each iteration of the tree
#2. fj: training fitted values
#fj.tst: testing fitted values
brt_fit = function(ifold=NULL, y.trn, x.trn, ep=NULL, rpart.control, 
                   y.tst=NULL, x.tst=NULL, shrinkage, n.trees, 
                   interaction.depth){
  #if(!is.null(interaction.depth)){
  #  rpart.control$interaction.depth = interaction.depth
  #}
  
  #Scenario 1
  #split the data according to CV, leave ifold to calcul CV error
  if(!is.null(ifold)){
    cverr = numeric(n.trees)
    ifold_trn = which(!1:nrow(x.trn)%in%ifold)
    y.tst = y.trn[ifold]
    x.tst = x.trn[ifold,]
    y.trn = y.trn[ifold_trn]
    x.trn = x.trn[ifold_trn,]
  } else {
    #Scenario 2
    #No need to keep CV error as we're fitting to all the data
    cverr = NULL
  }
  
  #calculate initial estimate which is the mean of all 
  #observed values
  rhoj = mean(y.trn)
  fj = rhoj
  #calculate initial estimate of test which is the mean of 
  #observed values from training where we have data for
  rhoj.tst = mean(y.trn)
  fj.tst = rhoj.tst
  
  for(j in 1:n.trees){
    out_brt_update = brt_update(y.trn, x.trn, fj, fj.tst, ep=ep, 
                      rpart.control, y.tst, x.tst, shrinkage)            
    fj = out_brt_update$fj
    fj.tst = out_brt_update$fj.tst
    if(any(is.nan(fj.tst))){
      break
    }
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



