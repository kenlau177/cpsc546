#newton-raphson optimization for SEEI
#direction computation optimization for quartic loss
#
#input: current terminal node estimate
#training data
#
#output: new terminal node estimate
newton_rap = function(rho0, y, ep=NULL, tol=1e-2){
  if(is.null(ep) | is.na(ep)){
    #for(nr in 1:500){
    #  rho.tmp = rho
    #  H = hess.rss.quar(fj=rho, y=y)
    #  g = grad.rss.quar(fj=rho, y=y)
      #print(H)
    #  if(g==0 | H==0){
    #    break
    #  }
    #  rho = rho - .01*g/H
    #  if(abs(rho-rho.tmp) < tol){
    #    break
    #  }
    
    #direct compuation to solve optimization for
    #quartic loss function
    sy2 = sum(y^2)
    sy = sum(y)
    n = length(y)
    rho.tmp = (sy2*sy*n^4 - (sy^3)*(n^3))
    rho = (sign(rho.tmp)*(abs(rho.tmp))^(1/3))/(n^2) + sy/n
    if(is.nan(rho)){
      return(rho0)
    }
    if(hess.rss.quar(fj=rho, y=y) < 0){
      return(rho0)
    } else {
      return(rho)
    }
  
  } else {
    rho = rho0
    #max iterations set to 500
    for(nr in 1:500){
      rho.tmp = rho
      #calculate second derivative of objective function
      H = hess.rss.sqep(fj=rho, y=y, ep=ep)
      #calculate first derivative of objective function
      g = grad.rss.sqep(fj=rho, y=y, ep=ep)
      if(g==0 | H==0){
        break
      }
      rho = rho - .1*g/H
      if(abs(rho-rho.tmp) < tol){
        break
      }
    }
  }
  return(rho)  
}

#input function for ddply vectorization
#takes current terminal node estimate, rho0, and calculates
#new terminal node estimate based on:
#1. SEEI loss then performs newton-raphson
#2. quartic loss then performs direct calculation
#
#input: tol for newton-raphson, generally we do not need something
#very close because it will most likely be corrected on the 
#next iteration of boosting
#
#output: new terminal node estimate
inply_sqep = function(x, ep=NULL, tol=1e-2){
  #current terminal node estimate and response variable training data
  rho0 = unique(x$rho); y = x$y;
  if(length(rho0)>1){
    #for whatever reason, I get two terminal node estimates
    #then quit program
    stop("rho0 is wrong")
  }
  #out = try(optimx(par=rho0, fn=rss.sqep, gr=NULL, 
  #        hess=NULL, method="BFGS", 
  #        itnmax=2, control=list("trace"=0), y=y, ep=ep)$p1)
  
  #function to calculate new terminal node estimates
  #I called this newton_rap, but in fact newton raphson is 
  #done for SEEI loss, and direct computation is done for 
  #quartic loss.
  rho = newton_rap(rho0=rho0, y=y, ep=ep, tol=tol)
  return(rho)
}

#modifies fitted regression tree object by adjusting terminal node
#estimates. Handles 2 loss functions:
#1. quartic loss
#2. SEEI loss
#
#input:
#object: fitted regression tree object
#newdata: training data
predict2.rpart = function(object, newdata, ep=NULL){
  #indices of terminal node estimates in object's frame
  where = object$where
  frame = object$frame
  frame$where = 1:nrow(frame)
  y = object$y
  
  #dataframe containing the index of terminal node, the 
  #current terminal node estimates, and training data, y
  rho_df = data.frame("where"=where, "rho"=frame$yval[where], 
            "y"=y)
  #for each terminal node, apply function to adjust terminal node
  #estimate based on loss function and using newton-rapson for 
  #SEEI and direct computation for quartic loss
  rho_df_update = ddply(rho_df, .(where), .fun=inply_sqep, ep=ep)
  names(rho_df_update) = c("where", "yval")
  #replace current terminal node estimates with adjusted
  frame$yval = replace(frame$yval, frame$where %in% rho_df_update$where, 
                rho_df_update$yval)
  frame = frame[,-match("where", names(frame))]
  
  object$frame = frame
  
  return(object)
  
  #out = predict(object, newdata=newdata, type=type)
  #return(out)
}


