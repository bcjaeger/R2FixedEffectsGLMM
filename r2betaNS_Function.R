

r2beta.NS <- function(model = m, slope = orthodont$age){
  
  # This function calculates a simplified version of R2_Beta 
  # using the Wald F statistic from a linear mixed model.
  # The resulting R^2 statistic is approximately equivalent
  # to the marginal R^2 proposed by Nakagawa and Schielzeth (2014) 
  
  # Note: this function is sloppily coded and only works properly
  # for the dental simulation study. 
  
  beta = fixef(model)
  X = model.matrix(model)
  n <- nrow(X)
  if (ncol(VarCorr(model)[1]$Subject) > 1){Z = cbind(1, slope)}
  else {Z = matrix(1, nrow = n)}
  p <- length(beta)
  C = cbind(rep(0, p-1),diag(p-1))
  
  # Actual Covariance Matrix from the model
  SigHat = Z%*%VarCorr(model)[1]$Subject%*%t(Z) + sigma(model)^2*diag(nrow(Z))
  
  # NS mean of the trace of the matrix above
  sighat = mean(diag(SigHat))
  
  # Wald F stat using the NS matrix
  WaldF = t(C %*% beta) %*% ginv(C %*% solve(t(X)%*%X * (sighat)^(-1)) %*% t(C)) %*% C%*%beta / (ncol(C)-1)
  
  # Fstat using the real matrix
  WaldF.r = t(C %*% beta) %*% ginv(C %*% solve(t(X)%*%solve(SigHat)%*%X ) %*% t(C)) %*% C%*%beta / (ncol(C)-1)
  
  ddf = n-p
  ndf = p-1
  ss = ndf/ddf * WaldF
  
  R2 = ss/(1+ss) 
  
  return(R2)
}

