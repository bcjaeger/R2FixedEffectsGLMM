

r2beta <- function(model = NULL, partial = F){
  
  # Computes R^2 statistic from edwards et al., 2008.
  #
  # Args:
  #    model: A linear mixed model fit with the lmer() function
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
  
  ### Get random effects from model call
  
  random = paste('(',unlist(lme4::findbars(model@call)), ')', sep = '', collapse = '+')
  
  ### Null model formula: 
  ### exactly the same covariance structure with all fixed effects removed except the intercept 
  
  null.form <- as.formula(paste('. ~ 1 +', random))
  
  ### Compute Kenward Roger approximate F test using the null model defined above

  mc <- pbkrtest::KRmodcomp(model, update(model, null.form))$stat 

  ### Compute the R2beta statistic for the full model
  ### Store results in a dataframe r2
  
  R2 = data.frame(Effect = 'Model', F = mc$Fstat, v1 = mc$ndf, v2 = mc$ddf, ncp = mc$Fstat * mc$ndf,
                  Rsq       = with(mc, (ndf * Fstat  / ddf   ) / (1 + ndf * Fstat  / ddf  )))
  
  ### use mixed function to conduct approximate KR F tests for each
  ### fixed effect in the full model
  
  if (partial == T){
    
    partials <- afex::mixed(model@call$formula, data = model@frame, progress = F)$anova.table
    
    ### Compute partial R2beta statistics for all fixed effects
    
    r2part = mutate(data.frame(partials[c('Effect', 'F', 'ndf', 'ddf')]), 
                    Rsq = (ndf * F / ddf) / (1 + ndf*F/ddf), 
                    ncp = F * ndf) %>% rename(v1 = ndf, v2 = ddf)
    
    R2 = rbind(R2, r2part) 
    
    }
  
  return(R2)

}

