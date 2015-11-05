

library(lme4)
library(MuMIn)
library(ggplot2)
library(nlme) 
library(MASS)
library(plyr)
library(dplyr)
library(tables)

# My functions

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
r2beta.NS <- function(model = m, slope = Orthodont$age){
  
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

# Some functions for tables
myfun=function(x){c(Mean = mean(x), Sd = sd(x))} 
tabfun <- function(x,y, dig=2){paste(round(x,dig), ' (', round(y, dig), ')', sep = '')} 


# Set the number of sims to run
nsims = 10000

# Names for the output matrix
nms <- c('R2m', 'R2c', 'R2betaNS', 'R2betaKR')

# Get some dental data
data(Orthodont)

# Make three matrices to hold results for the three outcome types
# R2.Norm1 holds the R2 values from the correct model
# R2.Norm2 holds R2 values from an underspecified model
# R2.Norm3 holds R2 values from model with noise parameters

R2.Norm1 = R2.Norm2 = R2.Norm3 = matrix(0, nrow = nsims, ncol = length(nms), dimnames = list(1:nsims,  nms))
Tabs = R2.Dat = list()

m <- lmer(distance~age*Sex+(1+age|Subject), data = Orthodont)

subs = unique(Orthodont$Subject)
n = length(subs)

for(cov in c('Covariance 1', 'Covariance 2')){

  model.cov = ifelse(cov == 'Covariance 1', '(1|Subject)', '(1+age|Subject)')
  null.form = paste('sim_1 ~ 1 +', model.cov)  

# Dental data simulation
for (sim in 1:nsims){

  # So I can watch something
  print(sim)
  
  # Generate the data by simulating data from m
  dat = cbind(Orthodont,simulate(m, 1, seed = sim))
  dat = merge(dat, data.frame(Subject = subs, grp = rbinom(length(subs), 1, 0.5)))
  dat$noise = rnorm(nrow(dat))
  
  # Model 0 is the null model
  m0 <- lmer(as.formula(null.form), data = dat)

  # Model 1 is the correct model '
  m1 <- update(m0, .~. + age*Sex)

  # Model 2 is an underspecified model
  m2 <- update(m0, .~. + age)
  
  # Model 3 is a misspecified model
  m3 <- update(m0, .~. + noise*grp)
             
  R2.Norm1[sim, ] = unlist(c(r.squaredGLMM(m1), r2beta.NS(m1), r2beta(m1, partial = F)[c('Rsq')]))
  R2.Norm2[sim, ] = unlist(c(r.squaredGLMM(m2), r2beta.NS(m2), r2beta(m2, partial = F)[c('Rsq')]))
  R2.Norm3[sim, ] = unlist(c(r.squaredGLMM(m3), r2beta.NS(m3), r2beta(m3, partial = F)[c('Rsq')]))
  
}

# Get summary statistics
r2.norm1 = data.frame(t(apply(R2.Norm1, 2, myfun))) %>% mutate(Model = 'Full', Statistic = colnames(R2.Norm1))
r2.norm2 = data.frame(t(apply(R2.Norm2, 2, myfun))) %>% mutate(Model = 'Reduced', Statistic = colnames(R2.Norm2))
r2.norm3 = data.frame(t(apply(R2.Norm3, 2, myfun))) %>% mutate(Model = 'Noise', Statistic = colnames(R2.Norm3))

r2.df <- rbind(r2.norm1, r2.norm2, r2.norm3) %>% 
  dplyr::mutate(Model = factor(Model, levels = c('Full', 'Reduced', 'Noise')),
                Statistic = factor(Statistic),
                tabvar = tabfun(Mean, Sd)) %>%
  plyr::arrange(Model, Statistic)

R2.Dat[[cov]] = r2.df
Tabs[[cov]] = tabular(Model ~ Statistic * Heading() * (tabvar * Heading()*paste), data = r2.df)

}

rbind(Tabs[[1]], Tabs[[2]])

# Some scratchwork used to show that R^2_beta reduces to R^2_NS in the normal case

m0 <- lmer(distance~1+(1+age|Subject), data = Orthodont)
m <- lmer(distance~age+(1+age|Subject), data = Orthodont)
mm=summary(m); mm0 = summary(m0)

# Check to make sure my r2beta.NS function agrees with Naka/Schielz
r2beta.NS(m)
r.squaredGLMM(m)

# pred is the predicted values of y from the model m
# here we only use the population values - no random effect
pred = model.matrix(m)%*%fixef(m)
pred0 = model.matrix(m0)%*%fixef(m0)

# y is the response variable
y = Orthodont$distance

# sigma^2 from the full model
s2 = var(y - pred)

# sigma^2 from the null model
s20 = var(y - pred0)

# This is sigma^2_f
s20 - s2
# Variance of fixed effects (Snijders and Bosker, 2012)
s2f = var(pred) 









