
rss.ls.quant = function(rho, y, quant){
  q = quantile(y, probs=quant)
  idx = which(y>q[1] & y<q[2])
  
  y = y[idx]
  rho = rho[idx]
  out = rss.ls(rho=rho, y=y)
  return(out)
}





