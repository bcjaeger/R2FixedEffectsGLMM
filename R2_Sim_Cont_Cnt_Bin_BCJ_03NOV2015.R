

library(lme4)
library(MuMIn)
library(stringr)
library(tables)
library(plyr)
library(dplyr)

Sim.LongDat <- function(seed = 260689, subs = 100, beta = c(0,-1,1/2,1/2),
                        minobs = 10, maxobs = 20, 
                        sigma = 2.5, tau0 = 2.5, tau1 = 0.2, tau01 = 0.3){
  
  # Simulated longitudinal data from a simple specification
  # with two fixed effects and two random effects (intercept and slope)
  #
  # Args:
  #    seed     : A number for random seed
  #    subs     : The number of subjects
  #    minobs   : Minimum number of observations per subject
  #    maxobs   : Maximum number of observations per subject
  #    beta     : fixed effects vector, ordered as follows
  #    
  #     (intercept, X1i, , time, X1i*time)
  #    
  #    sigma    : true random error SD
  #    tau0     : true random intercept SD
  #    tau1     : true random slope SD 
  #    tau01    : true random intercept-slope correlation
  #
  # Returns:
  #   a dataframe containing simulated longitudinal data
  # 
  # Programmer:  
  #   Byron C. Jaeger
  #   bcjaeger@live.unc.edu
  
  set.seed(seed)
  
  require(MASS)
  
  # Grouping variable
  X1i <- rbinom(subs,size=1,prob=.5)
  
  # simulate (correlated) random effects for intercepts and slopes
  S   <- matrix(c(1, tau01, tau01, 1), nrow=2)
  tau <- c(tau0, tau1)
  G   <- diag(tau) %*% S %*% diag(tau)
  
  bi  <- mvrnorm(subs,mu=rep(0, nrow(G)),Sigma=G)
  
  # vector of observations per subject
  j   <- round(runif(subs,minobs,maxobs)) 
  
  # vector of observation times
  tij <- as.vector(unlist(sapply(j, function(x) (1*(1:x-1)))))  
  
  dat <- data.frame(id=factor(rep(1:subs,times=j)),obstime=tij)
  dat$X1i <- rep(X1i,times=j)
  dat$grp <- factor(dat$X1i, labels = c('Control', 'Treatment'))
  dat$eij <- as.vector(unlist(sapply(j, function(x) rnorm(x,mean=0,sd=sigma))))
  dat$b1i <- rep(bi[,1],times=j)
  dat$b2i <- rep(bi[,2],times=j)
  
  # Make a continuous normally distributed outcome
  dat$eta <- with(dat,beta[1]+beta[2]*X1i+beta[3]*obstime+beta[4]*X1i*obstime)
  dat$etac <- with(dat, eta + b1i + b2i*obstime)
  dat$yij <- with(dat, eta + eij + b1i + b2i*obstime)
  
  # Make a binary outcome
  dat$yij_bin <- factor(dat$yij > median(dat$yij), labels = c('0', '1'))
  
  # Make a count outcome
  dat$yij_pois <- pmax(0, round( (dat$yij - min(dat$yij))/ 15  ,0))
  return(dat)
}

# Set the number of sims to run
nsims = 1000

# Set a location for simulated data to be stored
file.loc <- 'C:/Users/Byron/Desktop/School/Dissertation/R Code/Data'

# Names for the results matrix
nms <- c('True R2m', 'True R2c', 'R2m', 'R2c')

# Make three matrices to hold results for the three outcome types
R2.Norm = R2.Pois = R2.Bin = matrix(0, nrow = nsims, ncol = length(nms), dimnames = list(1:nsims,  nms))

Tabs = R2.Dat = list()

for (sim in 1:nsims){
  # Generate a dataset using one of my functions
  dat <- Sim.LongDat(seed = sim)
  # Store it 
  write.csv(dat, file = paste(file.loc, '/LongDat', sim, sep = ''), row.names = FALSE, quote = FALSE)
}

# y is the vector of outcomes for the three models
y <- c(Normal = 'yij', Poisson = 'yij_pois', Binomial = 'yij_bin')

for (fix in c('Full', 'Reduced')){
  for(cov in c('Covariance 1', 'Covariance 2')){
    
    # Set up the equations for the normal, log linear, and logit models
    model.cov = ifelse(cov == 'Covariance 1', '(1|id)', '(1+obstime|id)')
    model.fix = ifelse(fix == 'Full', '(X1i * obstime)', 'X1i + obstime')
    mod.form = paste(y, '~', model.fix, '+', model.cov) ; names(mod.form) = names(y)
    
    for (sim in 1:nsims){
      
      # Something to watch while you wait
      print(sim)
      
      # read the dataset from stored file
      dat <- read.csv(file = paste(file.loc, '/LongDat', sim, sep = ''))
      
      # make three models (Normal, Poisson, and Binomial)
      lmm <- lmer(as.formula(mod.form['Normal']), data = dat)
      glmm.cnt <- glmer(as.formula(mod.form['Poisson']), family = 'poisson', data = dat)
      glmm.bin <- glmer(as.formula(mod.form['Binomial']), family = 'binomial', data = dat)
      
      
      # Put all the R^2 stats into their matrices 
      R2.Norm[sim, ] = unlist(c(r.squaredGLMM(lmm)))
      R2.Pois[sim, ] = unlist(c(r.squaredGLMM(glmm.cnt)))
      R2.Bin[sim,  ] = unlist(c(r.squaredGLMM(glmm.bin)))
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

# Put it all together
R2.Dat = rbind(R2.Dat[[1]], R2.Dat[[2]], R2.Dat[[3]], R2.Dat[[4]])

# The KR DDF method in SAS is more reliable, so we use the estimated R^2_\beta from SAS 
# Note: The SAS code must be run and results stored before this code is run.

sasdat = read.csv(paste(file.loc, '/R2_SIM_13OCT2015.csv', sep = '')) %>%
  mutate(Outcome = factor(Outcome, ordered=T), levels = 'Normal', 'Poisson', 'Binomia',
         Statistic = 'R2b') %>%
  dplyr::rename(R2 = Rsq_B, Model = Pred) %>%
  dplyr::select(R2, NumDF, DenDF, Outcome, Model, Cov, Statistic)
  sasdat$Outcome = revalue(sasdat$Outcome, c("Binomia"="Binomial"))
  sasdat$Model = revalue(sasdat$Model, c('Full Model' = 'Full', 'Reduced Model' = 'Reduced'))

# a quick look at the results from SAS
tabular(Cov * Model ~ Outcome * R2 * (mean + sd), data = sasdat)

# Put all the R^2 statistics together
simdat = rbind(subset(R2.Dat, Statistic %in% c('R2m', 'R2c')), sasdat) %>% 
  mutate(Statistic = factor(Statistic, order = T))
levels(simdat$Statistic) <- c('R2m', 'R2c', 'R2_B')

# Melt the data for graphing
gdat <- reshape2::melt(data.frame(simdat), value.name = 'R2')
gdat <- dplyr::rename(gdat, Statistic = variable)

# A helpful way to look at the simulation results
ggplot(gdat, aes(x=R2, fill=Statistic)) +
  geom_histogram(binwidth=0.01, colour="white", position = position_stack()) + facet_wrap(Cov ~ Model + Outcome, ncol = 3)

# Descriptive Graphs of simulated data

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
  
pnorm
pbin
ppois

