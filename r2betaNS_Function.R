

r2beta.NS <- function(model){
  
  # model: an lmer object
  # Note: As of now, this function only works for 2 level models
  
  # This function calculates a simplified version of R2_Beta 
  # using the average variance from \hat{\Sigma}.
  # The resulting R^2 statistic is approximately equivalent
  # to the marginal R^2 proposed by Nakagawa and Schielzeth (2014) 
  
  # Get fixed effects
  beta = fixef(model)
  p <- length(beta)
  
  # Get random effects
  Sig.b = raneff.covmat(model)
  
  # Get model matrices
  X = model.matrix(model)
  n <- nrow(X)
  Z = getME(model, 'Z')
  
  # C matrix defines the Wald Test for Fixed Effects
  C = cbind(rep(0, p-1),diag(p-1))
  
  # Counts of level 1 units within level 2 
  N = as.numeric(table(eval(model@call$data)[, names(model@flist)[1]]))

  # Stacked (diagonal) random effects covariance matrix
  Sig.B = kronecker(diag(length(N)), Sig.b)
  
  # Covariance Matrix from the model
  SigHat = Z%*%Sig.B%*%t(Z) + sigma(model)^2*diag(nrow(Z))
  
  # Mean of the trace of the matrix above
  sighat = mean(diag(SigHat))
  
  # Wald F stat using the mean covariance estimate, 
  WaldF = t(C %*% beta) %*% ginv(C %*% solve(t(X)%*%X * (sighat)^(-1)) %*% t(C)) %*% C%*%beta / (ncol(C)-1)
  
  ddf = n-1
  ndf = ncol(C)-1
  ss = ndf/ddf * WaldF
  R2 = ss/(1+ss) 
  
  return(R2)
   
}
