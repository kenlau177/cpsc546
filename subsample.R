
#perform the stochastic part of boosting step of subsampling
#50% of data prior to fitting regression tree to the negative
#gradients

#input: y response variable 1 dimension
#x explnatory variables p dimensions
#zj negative gradient of current iteration 1 dimension
#
#output: 50% subsampled data of y, x, and zj
subsample = function(y, x, zj){
  sub_size = round(nrow(x)*.5)
  isub = sample.int(nrow(x), size=sub_size)
  out.x = x[isub,]
  out.y = y[isub]
  out.zj = zj[isub]
  return(list("y"=out.y, "x"=out.x, "zj"=out.zj))
}

