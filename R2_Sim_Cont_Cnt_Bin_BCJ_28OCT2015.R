

library(lme4)
library(MuMIn)
library(stringr)
library(tables)
library(plyr)
library(dplyr)

# Load my babies
source('My_Functions.R')

# Set the number of sims to run
nsims = 10

# Set location of data output
file.loc <- 'C:/Users/Byron/Desktop/School/Dissertation/R Code/Data'

# Names for the output matrix
nms <- c('True R2', 'True R2m', 'True R2c', 'R2m', 'R2c')

# Make three matrices to hold results for the three outcome types
R2.Norm = R2.Pois = R2.Bin = matrix(0, nrow = nsims, ncol = length(nms), dimnames = list(1:nsims,  nms))
#R2.Norm1 = R2.Norm2 = R2.Norm3 = matrix(0, nrow = nsims, ncol = length(nms), dimnames = list(1:nsims,  nms))

Tabs = R2.Dat = list()

for (sim in 1:nsims){
  # Generate a dataset using one of my functions
  dat <- Sim.LongDat(seed = sim, subs = 100, beta = c(0,-1,1/2,1/2),
                     minobs = 10, maxobs = 20, 
                     sigma = 2.5, tau0 = 5, tau1 = 0.2, tau01 = 0.3)
  # Store it in my archive
  write.csv(dat, file = paste(file.loc, '/LongDat', sim, sep = ''), row.names = FALSE, quote = FALSE)
}

y <- c(Normal = 'yij', Poisson = 'yij_pois', Binomial = 'yij_bin')

for (fix in c('Full', 'Reduced')){
  for(cov in c('Covariance 1', 'Covariance 2')){
    
    model.cov = ifelse(cov == 'Covariance 1', '(1|id)', '(1+obstime|id)')
    model.fix = ifelse(fix == 'Full', '(X1i * obstime)', 'X1i + obstime')
    mod.form = paste(y, '~', model.fix, '+', model.cov) ; names(mod.form) = names(y)
    
    # Not too complicated of a loop
    
    for (sim in 1:nsims){
      
      # So I can watch
      print(sim)
      
      # read the dataset from stored file
      dat <- read.csv(file = paste(file.loc, '/LongDat', sim, sep = ''))
      # make three models (Normal, Poisson, and Binomial)
      lmm <- lmer(as.formula(mod.form['Normal']), data = dat)
      glmm.cnt <- glmer(as.formula(mod.form['Poisson']), family = 'poisson', data = dat)
      glmm.bin <- glmer(as.formula(mod.form['Binomial']), family = 'binomial', data = dat)
      
      # Make R^2_{Likelihood} using r2.ll function
      
      null.form = paste(y, '~ 1 + (', findbars(as.formula(mod.form['Normal'])), ')'); names(null.form) = names(y)
      
      subs = length(unique(dat$id))

      # Likelihood Ratio R^2      
      # r2.L <- c(nrm = r2.ll(lmm, update(lmm, as.formula(null.form['Normal'])), n = subs), 
      # cnt = r2.ll(glmm.cnt, update(glmm.cnt, as.formula(null.form['Poisson'])), n = subs),
      # bin = r2.ll(glmm.bin, update(glmm.bin, as.formula(null.form['Binomial'])), n = subs))
      
      # true marginal R2 
      r2.m = with(dat, var(eta) / var(yij))
      
      # Estimate of true R^2
      r2.t = sum(dat$exp) / sum(dat$exp + dat$unexp)
      
      # true conditional R2 
      r2.c = with(dat, var(etac) / var(yij))
      
      # Put all the R^2 stats into their matrices 
      R2.Norm[sim, ] = unlist(c(r2.t, r2.m, r2.c, r.squaredGLMM(lmm)))
      R2.Pois[sim, ] = unlist(c(r2.t, r2.m, r2.c, r.squaredGLMM(glmm.cnt)))
      R2.Bin[sim,  ] = unlist(c(r2.t, r2.m, r2.c, r.squaredGLMM(glmm.bin)))
    }  
    
    r2.norm = data.frame(R2.Norm) %>% mutate(Outcome = names(y)[1], Model = fix, Cov = cov)
    r2.pois = data.frame(R2.Pois) %>% mutate(Outcome = names(y)[2], Model = fix, Cov = cov)
    r2.bin = data.frame(R2.Bin) %>% mutate(Outcome = names(y)[3], Model = fix, Cov = cov)
    
    r2.df <- rbind(r2.norm, r2.pois, r2.bin) %>% 
      dplyr::mutate(Outcome = factor(Outcome, levels = names(y)),
                    Cov = factor(Cov), Model = factor(Model))%>%
      reshape2::melt(value.name = 'R2')
    r2.df <- dplyr::rename(r2.df, Statistic = variable)
    R2.Dat[[paste(fix, cov, collapse = '_')]] <- r2.df
    
  }  
}

R2.Dat = rbind(R2.Dat[[1]], R2.Dat[[2]], R2.Dat[[3]], R2.Dat[[4]])

# Read data from SAS: 
sasdat = read.csv(paste(file.loc, '/R2_SIM_13OCT2015.csv', sep = '')) %>%
  mutate(Outcome = factor(Outcome, ordered=T), levels = 'Normal', 'Poisson', 'Binomia',
         Statistic = 'R2b') %>%
  dplyr::rename(R2 = Rsq_B, Model = Pred) %>%
  dplyr::select(R2, NumDF, DenDF, Outcome, Model, Cov, Statistic)
  sasdat$Outcome = revalue(sasdat$Outcome, c("Binomia"="Binomial"))
  sasdat$Model = revalue(sasdat$Model, c('Full Model' = 'Full', 'Reduced Model' = 'Reduced'))


tabular(Cov * Model ~ Outcome * R2 * (mean + sd), data = sasdat)

simdat = rbind(subset(R2.Dat, Statistic %in% c('R2m', 'R2c')), sasdat) %>% 
  mutate(Statistic = factor(Statistic, order = T))

levels(simdat$Statistic) <- c('R2m', 'R2c', 'R2_B')

#tabular(Cov*Model*Outcome~Statistic*Heading()*R2*Heading()*function(x,dig=2){paste(round(mean(x),dig),' (',round(sd(x),dig),')',sep='')},data=simdat)

gdat <- reshape2::melt(data.frame(simdat), value.name = 'R2')
gdat <- dplyr::rename(gdat, Statistic = variable)

ggplot(gdat, aes(x=R2, fill=Statistic)) +
  geom_histogram(binwidth=0.01, colour="white", position = position_stack()) + facet_wrap(Cov ~ Model + Outcome, ncol = 3)


# Descriptive Graphs

dat$id = factor(dat$id)
set.seed(260689)

sdat1 <- subset(dat, id %in% sample(unique(subset(dat, grp == 'Treatment')$id), 5))
sdat2 <- subset(dat, id %in% sample(unique(subset(dat, grp == 'Control')$id), 5))

sdat <- rbind(sdat1, sdat2)

# Black/White plot

require(ggplot2)

pnorm <- ggplot(data = dat, aes(x = obstime, y = yij)) + 
  geom_point(data = sdat, size = 0.1) + 
  geom_line(data = sdat, aes(group = id, linetype = grp)) + 
  geom_smooth(aes(group = grp, linetype = grp), method = 'loess', size = 1.2, se = F, colour = 'black')+ 
  labs(x = "Time", y = "", title = 'Normal Outcomes') +
  theme_classic() + theme(legend.position = "none")

pbin <- ggplot(data = dat, aes(x = obstime, y = yij_bin)) + 
  geom_point(data = sdat, size = 0.1) + 
  geom_line(data = sdat,aes(group = id, linetype = grp)) + 
  geom_smooth(aes(group = grp, linetype = grp), method = 'loess', size = 1.2, se = F, colour = 'black' )+ 
  labs(x = "Time", y = "", title = 'Binary Outcomes', linetype = 'Group' ) +
  theme_classic() + theme(legend.position = "top") 

ppois <- ggplot(data = dat, aes(x = obstime, y = yij_pois)) + 
  geom_point(data = sdat, size = 0.1) + 
  geom_line(data = sdat,aes(group = id, linetype = grp)) + 
  geom_smooth(aes(group = grp, linetype = grp), method = 'loess', size = 1, se = F, colour = 'black')+ 
  labs(x = "Time", y = "", title = 'Count Outcomes') +
  theme_classic()  + theme(legend.position = "none") 
  
multiplot(pnorm, pbin, ppois, cols = 3)

