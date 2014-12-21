
#initial data cleansing: remove variable, data completion

#input: x, raw data from Los Angeles Pollution ozone data
#output: data after cleaning
process_Ozone = function(x){
  #remove temperature variable at Saunders
  out = x[,-9]
  #rename variables
  names(out) = c("Month", "Month_Day", "Week_Day", "Max_O3",
                 "Pressure_Height", "Wind_Speed", "Humidity", 
                 "Temperature", "Inversion_Base_Height", 
                 "Pressure_Gradient", "Inversion_Base_Temperature", 
                 "Visibility")
  #temp = out[,-match("Max_O3",colnames(out))]
  #out = data.frame("Max_O3"=out$Max_O3, temp)  
  
  #remove any incomplete rows or data
  out = out[complete.cases(out),]
  return(out)
}

#given dataframe ozone data, split the data into lower 75% and 
#upper 25% quantiles. Do the following 2 sampling steps:
#1. In the lower quantile, sample 75% for training and 25% 
#for testing
#2. In the upper quantile, sample 75% for training and 25% 
#for testing
#combine the data respectively for training and testing
#
#input: dataframe x corresponding to ozone data
#output: list containing dataframes for training and testing
split_train_test = function(x){
  #create vector of factors corresponding to lower or upper quantile
  #of x
  id_q = with(x, cut(Max_O3, breaks=c(0,quantile(Max_O3, probs=.75), 
          max(Max_O3)), labels=c("low","up"), include.lowest=T))
  #index of x corresponding to lower and upper quantile
  id_low = which(id_q=="low")
  id_up = which(id_q=="up")
  
  num_low_trn = round(length(id_low)*.75)
  #sample from indices of low quantile for training data
  id_low_trn = sample(id_low, size=num_low_trn)
  #leave the rest for testing data
  id_low_test = id_low[!id_low%in%id_low_trn]
  
  num_up_trn = round(length(id_up)*.75)
  #sample from indices of upper quantile for training data
  id_up_trn = sample(id_up, size=num_up_trn)
  #leave the rest for testing data
  id_up_test = id_up[!id_up%in%id_up_trn]
  
  #concantenate indices for upper and lower for training and assign
  #to x
  train = x[c(id_low_trn,id_up_trn),]
  #concantenate indices for upper and lower for testing and assign
  #to x
  test = x[c(id_low_test,id_up_test),]
  
  return(list("train"=train, "test"=test))
}




# add_outlier = function(x, k, mean=8, sd=2){
#   id_q = with(x, cut(Max_O3, breaks=c(0,quantile(Max_O3, probs=.9), 
#                                       max(Max_O3)), labels=c("low","up"), include.lowest=T))
#   id_up = which(id_q=="up")
#   x_up = x[id_up,]
#   newx = x_up[sample(nrow(x_up), size=k),]
#   newx = newx + matrix(rnorm(length(newx),sd=2), nrow=nrow(newx), 
#                        ncol=ncol(newx))
#   newx$Max_O3 = newx$Max_O3 + rnorm(k, mean=mean, sd=sd)
#   out = rbind(x, newx)
#   return(out)
# }

