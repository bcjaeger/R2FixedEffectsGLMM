
library(lme4)
library(MuMIn)
library(ggplot2)
library(nlme) 
library(MASS)
library(plyr)
library(dplyr)
library(tables)

# Note: r2beta and r2beta.NS functions are needed to run the simulation.

# Some functions for tables
myfun=function(x){c(Mean = mean(x), Sd = sd(x))} 
tabfun <- function(x,y, dig=2){paste(round(x,dig), ' (', round(y, dig), ')', sep = '')} 

# Set the number of sims to run
nsims = 10000

# Names for the output matrix
nms <- c('R2m', 'R2c', 'R2betaNS', 'R2betaKR')

# Get the dental data
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
    
    # Something to watch
    # print(sim)
    
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
  
  r2.df$Cov = factor(cov)
  R2.Dat[[cov]] = r2.df
  Tabs[[cov]] = tabular(Cov * Model ~ Statistic * Heading() * (tabvar * Heading()*paste), data = r2.df)
  
}

rbind(Tabs[[1]], Tabs[[2]])
