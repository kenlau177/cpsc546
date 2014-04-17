
#SEEI loss function
loss.sqep = function(x, ep){
  out = pmax(0, x^2-ep)
  return(out)
}
#quartic loss function
loss.quar = function(x){
  out = x^4
  return(out)
}
#least squares loss function
loss.ls = function(x){
  out = x^2
  return(out)
}
#RS (residula sum) of SEEI loss
rss.sqep = function(rho, y, ep){
  out = sum(loss.sqep((y-rho), ep=ep))
#  D = diag(as.numeric(((y-rho)^2-ep)>0))
#  sum(D%*%loss.sqep((y-rho), ep=ep))
#  I = rep(1, length(y))
#  out = c(t(y-I*rho)%*%D%*%(y-I*rho)) - sum(D*ep)
  return(out)
}
#RS of quartic loss
rss.quar = function(rho, y){
  out = sum(loss.quar(y-rho))
  return(out)
}
#RSS of LS
rss.ls = function(rho, y){
  out = sum(loss.ls(y-rho))
  return(out)
}
#negative gradient of SEEI loss
neg_grad.sqep = function(y, fj, ep){
  D = diag(as.numeric(((y-fj)^2-ep)>0))
  out = 2*c(D%*%(y-fj))
  return(out)
}
#negative gradient of quartic loss
neg_grad.quar = function(y, fj){
  out = 4*(y-fj)^3
  return(out)
}
#negative gradient of LS
neg_grad.ls = function(y, fj){
  out = 2*(y-fj)
  return(out)
}
#gradient function for SEEI loss
grad.rss.sqep = function(fj, y, ep){
  D = diag(as.numeric(((y-fj)^2-ep)>0))
  I = rep(1, length(y))
  out = c(-2*(t(y) - fj*t(I))%*%(D%*%I))
  #out = -t(y)%*%D%*%I - t(I)%*%D%*%y + 2*fj*t(I)%*%D%*%I
  return(out)
}
#gradient function for quartic loss
grad.rss.quar = function(fj, y){
  sy2 = sum(y^2)
  sy = sum(y)
  n = length(y)
  out = -4*(sy2*sy - 3*fj*(sy^2) + 3*(fj^2)*n*sy - (fj^3)*(n^2))
  return(out)
}
#hessian function for SEEI loss
hess.rss.sqep = function(fj, y, ep){
  D = diag(as.numeric(((y-fj)^2-ep)>0))
  I = rep(1, length(y))
  out = c(2*(t(I)%*%D%*%I))
  return(out)
}
#hessian function for quartic loss
hess.rss.quar = function(fj, y){
  sy = sum(y)
  n = length(y)
  out = 12*(sy^2 - 2*fj*n*sy + (fj^2)*(n^2))
  return(out)
}




