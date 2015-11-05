

r2beta.glmm <- function(model, niter = 20, verbose = F, partial = F, ...){
  
  # Computes R^2 statistic for the GLMM.
  #
  # Args:
  #    formula: A formula that can be used as an lmer() formula
  #    family : The GLM family used for modelling non-normal outcomes
  #    data   : The data used to fit the GLMM
  #    niter  : maximum number of iterations for PQL algorithm
  #    verbose: if TRUE, R will print commentary on computations
  #    partial: if TRUE, includes computation of partial R2
  #             for all fixed effects in the model
  #
  # Returns:
  #   a dataframe containing the model F statistic, numerator  
  #   and denominator degrees of freedom, and R^2 statistic. If
  #   partial = TRUE, then the dataframe also contains partial
  #   R^2 statistics for all fixed effects in the model
  # 
  # Programmer:  
  #   Byron C. Jaeger
  #   bcjaeger@live.unc.edu
  
  require(lme4)
  require(pbkrtest) 
  require(afex)
  require(dplyr)
  
  formula = model@call$formula
  family  = model@call$family
  data    = eval(model@call$data)
  
  if(is.null(family)) {
    R2 = r2beta(model=model, partial = partial)
  }
  else{    
    ## get family
    if(is.character(family)) family <- get(family)
    if(is.function(family)) family <- family()
    if(is.null(family$family)) {print(family); stop("'family' not recognized")}
    
    ## get fixed effects 
    fixed <- nobars(formula)
    
    ## get formula for the null model
    
    ran.eff = paste('(',unlist(findbars(formula)), ')', sep = '', collapse = '+')
    null.form <- as.formula(paste('. ~ 1 + ', ran.eff))                           
    
    ## get initial weights and pseudo outcomes
    fit0 <- glm(formula=fixed, family=family, data=data)
    eta <- fit0$linear.predictors
    data$zz <- eta + fit0$residuals # initial pseudo outcome values
    data$wz <- fit0$weights # GLM weights from initial fit
    formula[[2L]] <- quote(zz) # Change response to pseudo-variables in model formula
    
    ## PQL algorithm : iterate between updating weights/outcomes and fitting a lmm 
    for(i in 1:niter) {
      if(verbose) message(gettextf("iteration %d", i), domain = NA)
      fit <- lmer(formula = formula, data = data, weights = wz)
      etaold <- eta
      eta <- fitted(fit) 
      if(sum((eta-etaold)^2) < 1e-10*sum(eta^2)) break;
      mu <- family$linkinv(eta) # mu = h(XB + Zb), where h is inverse link function
      mu.eta.val <- family$mu.eta(eta) # delta = derivative of mu with respect to eta
      data$zz <- eta + (fit0$y - mu)/mu.eta.val 
      data$wz <- mu.eta.val^2 / family$variance(mu)
    }
    
      mc <- pbkrtest::KRmodcomp(fit, update(fit, null.form))$stat 
      
      ### Compute the R2beta statistic for the full model
      ### Store results in a dataframe R2
      ### ndf and ddf are the numer
      ### the non-centrality parameter (ncp) is used to test hypotheses regarding covariance selection.

      
      
      R2 = data.frame(Effect = 'Model', F = mc$Fstat, v1 = mc$ndf, v2 = mc$ddf, ncp = mc$Fstat * mc$ndf,
                      Rsq       = with(mc, (ndf * Fstat  / ddf   ) / (1 + ndf * Fstat  / ddf  )))
      
      ### use mixed function to conduct approximate KR F tests for each
      ### fixed effect in the full model
      
      if (partial == T){

        ### Compute partial R2beta statistics for all fixed effects
        
        partials <- afex::mixed(formula = formula, data = data, progress = verbose)$anova.table
        
        r2part = mutate(data.frame(partials[c('Effect', 'F', 'ndf', 'ddf')]), 
                        Rsq = (ndf * F / ddf) / (1 + ndf*F/ddf), 
                        ncp = F * ndf) %>% rename(v1 = ndf, v2 = ddf)
        
        R2 = rbind(R2, r2part)    
        }
    }
  return(R2)
}



