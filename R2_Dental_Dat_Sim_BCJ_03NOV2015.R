
rm(list=ls())

library(lme4)
library(MuMIn)
library(r2glmm)
library(nlme) 
library(MASS)
library(plyr)
library(dplyr)
library(tables)

#  Summary function for tables
msd <- function(x,dig=2){
  res = format(round(c(mean(x), sd(x)) , dig), nsmall=dig)
  return(paste0(res[1], ' (', res[2], ')'))
}

# Function to save results in a dataframe
turn.to.df <- function(lst){
  df = lapply(lst, data.frame)
  for(i in names(df)){ df[[i]]=mutate(df[[i]], Model = i, cov = cov) }
  
  df <- do.call(rbind,df)%>% 
    reshape2::melt(value.name = 'R2', id = c('Model', 'cov'))%>%
    dplyr::rename(Statistic = variable)
  
  return(df)
}

# Set the number of sims to run
nsims = 10

# Names for the output matrix

# Get the dental data
data(Orthodont)

m <- lmer(distance~age*Sex+(1+age|Subject), data = Orthodont)

subs = unique(Orthodont$Subject)
n = length(subs)
r2dat=data.frame()

for(cov in c('CS', 'GC')){
  
  r2=list()
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
  
    modlist=list('Full'=m1, 'Reduced'=m2, 'Noise'=m3)

    new.r2 = lapply(modlist, function(x){
      res = unlist(c(
        r.squaredGLMM(x), 
        'r2bKR'=r2beta(x, data = dat, ddf='kr')['Rsq'],
        'r2bNS'=r2beta(x, data = dat, ddf='res', NS.adj = T)['Rsq']
      ))
      names(res)=gsub('.Rsq', '', names(res))
      return(res)
    })
    
    for(i in names(new.r2)) r2[[i]] = rbind(r2[[i]], new.r2[[i]])
    
  }
  
  r2.df = turn.to.df(r2)
  
  r2dat = rbind(r2dat, r2.df)
  
}

ch.vars <- lapply(r2dat, class) == "character"
r2dat[, ch.vars] <- lapply(r2dat[, ch.vars], as.factor)
r2dat$Model = factor(r2dat$Model, 
                     levels = c('Full', 'Reduced', 'Noise'),
                     ordered = T)

tab = tabular(cov*Model~Statistic*Heading()*R2*Heading()*(msd), data = r2dat)
