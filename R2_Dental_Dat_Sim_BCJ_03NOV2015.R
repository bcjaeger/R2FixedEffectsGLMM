

library(lme4)
library(MuMIn)
library(ggplot2)
library(nlme) 
library(MASS)
library(plyr)
library(dplyr)

# Load my babies
source('My_Functions.R')

# Set the number of sims to run
nsims = 10000

# Names for the output matrix
nms <- c('R2m', 'R2c', 'R2betaNS', 'R2betaKR')

# Get some dental data
data(Orthodont)

# Make three matrices to hold results for the three outcome types
# R2.Norm1 holds the R2 values from the correct model
# R2.Norm2 holds R2 values from an underspecified model
# R2.Norm3 holds R2 values from a misspecified model

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
             
  R2.Norm1[sim, ] = unlist(c(r.squaredGLMM(m1), r2beta.NS(m1), r2beta.glmm(m1, partial = F)[c('Rsq')]))
  R2.Norm2[sim, ] = unlist(c(r.squaredGLMM(m2), r2beta.NS(m2), r2beta.glmm(m2, partial = F)[c('Rsq')]))
  R2.Norm3[sim, ] = unlist(c(r.squaredGLMM(m3), r2beta.NS(m3), r2beta.glmm(m3, partial = F)[c('Rsq')]))
  
}

# Look at summary stats for each matrix
myfun=function(x){c(Mean = mean(x), Sd = sd(x))}
tabfun <- function(x,y, dig=2){paste(round(x,dig), ' (', round(y, dig), ')', sep = '')}

r2.norm1 = data.frame(t(apply(R2.Norm1, 2, myfun))) %>% mutate(Model = 'Full', Statistic = colnames(R2.Norm1))
r2.norm2 = data.frame(t(apply(R2.Norm2, 2, myfun))) %>% mutate(Model = 'Reduced', Statistic = colnames(R2.Norm2))
r2.norm3 = data.frame(t(apply(R2.Norm3, 2, myfun))) %>% mutate(Model = 'Noise', Statistic = colnames(R2.Norm3))

r2.df <- rbind(r2.norm1, r2.norm2, r2.norm3) %>% 
  dplyr::mutate(Model = factor(Model, levels = c('Full', 'Reduced', 'Noise')),
                Statistic = factor(Statistic),
                tabvar = tabfun(Mean, Sd)) %>%
  plyr::arrange(Model, Statistic)

library(tables)

R2.Dat[[cov]] = r2.df
Tabs[[cov]] = tabular(Model ~ Statistic * Heading() * (tabvar * Heading()*paste), data = r2.df)

}

gdat = rbind(R2.Dat$`Covariance 1` %>% mutate(Cov = 'Covariance 1'),
             R2.Dat$`Covariance 2` %>% mutate(Cov = 'Covariance 2')) %>% 
  dplyr::select(-Mean, -Sd) %>% reshape2::melt(value.name = 'tabvar')

latex(rbind(Tabs[[1]], Tabs[[2]]))

# Graphs

gdat <- reshape2::melt(gdat, value.name = 'R2')
gdat <- dplyr::rename(gdat, Statistic = variable)

ggplot(gdat, aes(x=R2, fill=Statistic)) +
  geom_histogram(binwidth=0.01, colour="white", position = position_stack())

r2sim.plot(r2.df)



# y is the response variable
y = Orthodont$distance

m0 <- lmer(distance~1+(1+age|Subject), data = Orthodont)
m00 <- lm(distance~1, data = Orthodont)
mm=summary(m); mm0 = summary(m0)

m <- lmer(distance~age+(1+age|Subject), data = Orthodont)

r2beta.NS(m)

r.squaredGLMM(m)

# pred is the predicted values of y from the model m
# here we only use the population values - no random effect
pred = model.matrix(m)%*%fixef(m)
pred0 = model.matrix(m0)%*%fixef(m0)

# sigma^2 from the full model
s2 = var(y - pred)
 (t(y) %*% y - t(y)%*%M%*% y) / (nrow(Orthodont)-1)

# sigma^2 from the null model
s20 = var(y - pred0)
( t(y) %*% y - t(y)%*%M0%*% y) / (nrow(Orthodont)-1)

# This is sigma^2_f
s20 - s2

s2f = var(pred) # Variance of fixed effects (Snijders and Bosker, 2012)
s2e = mm$sigma^2 # Variance of residual error
s2b0 = mm$varcor$Subject[1,1] # variance of random intercept
s2b1 = mm$varcor$Subject[2,2] # variance of random slopes for age
s2b01 = mm$varcor$Subject[2,1] # covariance of intercept and clope

# This is what Naka/Schielz do  
# Calculate sigma^2_l, the variance explained by random effects
s2l = s2b0 + mean(Orthodont$age^2*s2b1 + 2*Orthodont$age*s2b01)
# Marginal R2
s2f / (s2f + s2e + s2l)
# Conditional R2
(s2f + s2l) / (s2f + s2e + s2l)

# The MLE for s2e from a lm is close to s2l + s2e
summary(lm(distance ~ age*Sex, data = Orthodont))$sigma^2








