# Chapter 3 preliminary analysis

# Occupancy and/or abundance models
# accounting for:
# 1. phenology
# 2. subtransect habitat use
# 3. clustering (for abundance)
# 4. detection probability (if possible)

library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(mgcv)


# Simulate data
# Possibly use Reto's code here, or Matechou adaptation I made for phenology


# Organize monitoring data

# By week or by specific day?
# Week: could be specified in wide format
# Day: long format needed, but long might be necessary for subtransects


# Need untrimmed data with subsection counts
# Remember issues with duplicate species on same survey

rawcounts <- read_csv("data/rawcounts.csv")
rawcovs <- read_csv("data/rawcovariates.csv")
rawcovs <- rawcovs %>%
  mutate(year = year(SiteDate))

# get all surveys at sites/years matching sites where species has ever been seen

species <- "West Virginia White"
counts <- rawcounts %>% filter(CommonName == species)
allsurvs <- rawcovs %>% 
  filter(SiteID %in% unique(counts$SiteID))
allcounts <- rawcounts %>%
  filter(SiteID %in% unique(counts$SiteID))

# losing empty years at this step
test <- allcounts %>%
  group_by(SiteID, year) %>%
  complete(SeqID, subtransect) %>%
  select(SiteID, SeqID, subtransect, year) %>%
  distinct()
surveys <- merge(test, allsurvs, all.x = TRUE)

test <- merge(surveys, counts, by = c("SeqID", "SiteID", "SiteDate", "Week", "year", "subtransect"), all.x = TRUE)

test$Total <- plyr::mapvalues(test[, "Total"], from = NA, to = 0)
test$count <- plyr::mapvalues(test[, "count"], from = NA, to = 0)
finaldat <- test %>%
  group_by(SiteID, year) %>%
  mutate(uniqsubtran = length(unique(subtransect)),
         present = ifelse(count > 0, 1, 0),
         CommonName = species)

# put my data into Stebel format
# 1st try lumping subtransects together
# this would make assumption that spatial replicates are independent
# equal probability of habitat use
datsimple <- finaldat %>%
  group_by(SeqID, year) %>%
  mutate(nrep = length(unique(subtransect)),
         ndet = length(which(count > 0))) %>%
  ungroup() %>%
  mutate(day = yday(SiteDate)) %>%
  select(nrep, ndet, SiteID, day, year) %>%
  rename(site = SiteID) %>%
  distinct()

dat <- datsimple

# Specify model in BUGS language
sink("splinesSiteOcc_simple.txt")
cat("
    model { 
    ### Define seasonal and annual patterns in detectability
    for (m in 1:nyear) {  
    for (i in 1:n) {
    logit(p[m,i]) <- lp[m,i]
    lp[m,i] <- mfe[m,i]+mre[m,i]
    mfe[m,i] <- a[m]*X[i,1]+b[m]*X[i,2]+c[m]*X[i,3]
    mre[m,i]<-sum(n.mre[m,i,1:nknots])
    for (k in 1:nknots) {
    n.mre[m,i,k]<-b.k[m,k]*Z[i,k]
    }
    }
    
    ### Random regression coefficients corresponding to the truncated polynomial functions
    for (k in 1:nknots) {
    b.k[m,k] ~ dnorm(0,taub)
    }
    
    ### Fixed regression coefficients corresponding to the 'plus' functions
    
    a[m] ~ dnorm(0,0.01)
    b[m] ~ dnorm(0,0.01)
    c[m] ~ dnorm(0,0.01)
    }
    
    ### precision for random regression coefficients corresponding to the truncated polynomial functions
    taub~dgamma(1.0E-6,1.0E-6)      
    
    # Specify priors
    for (k in 1:nyear) {
    psi[k] ~ dunif(0, 1)
    }
    
    # Ecological submodel: Define state conditional on parameters
    for (i in 1:nsite){
    for (k in 1:nyear){
    z[i,k] ~ dbern(psi[k])
    }
    }
    
    # Observation model
    for (i in 1:nobs){
    muy[site[i],survey[i],year[i]] <- z[site[i],year[i]]*p[year[i],survey[i]]
    y[i] ~ dbin(muy[site[i],survey[i],year[i]], nrep[i])
    }
    
    }
    ",fill = TRUE)
sink()

### Read observation data from Acrocephalus arundinaceus
# dat<-read.table(file="dat S4.txt",header=T)

### The following procedure is based on the models presented in Crainiceanu et al. 2005 and in Gimenez et al. 2006 
# Degree of splines
degree <- 2

# covariate
covariate<-as.numeric(scale(range(dat$day)[1]:range(dat$day)[2]))

# covariate length
n <- length(covariate)

# location of knots
nk<-round((max(dat$day)-min(dat$day)+1)/4)
nknots<-ifelse(nk<35,nk,35)
knots<-quantile(unique(covariate),seq(0,1,length=(nknots+2))[-c(1,(nknots+2))])

# fixed effects matrix
X<-NULL
for (l in 0:degree) {
  X<-cbind(X,covariate^l)  
}

# random coefficients matrix
Z_K<-(abs(outer(covariate,knots,"-")))^3
OMEGA_all<-(abs(outer(knots,knots,"-")))^3
svd.OMEGA_all<-svd(OMEGA_all)
sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*% (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d)))
Z<-t(solve(sqrt.OMEGA_all,t(Z_K)))

# Input data
site <- as.numeric(factor(dat$site))
survey <- dat$day-min(dat$day)+1
nobs <- length(unique(paste(dat$site,dat$day,dat$year)))
nrep <- dat$nrep
nsite <- length(unique(dat$site))
nyear <- length(unique(dat$year))
year <- as.numeric(factor(dat$year))
zst <- array(1, dim=c(nsite,nyear))  
y <- dat$ndet

# Simulation parameters
ni=500; nc=2; nb=250; nt=5

# List input data
# jags.data <- list("site","survey","nobs","nrep","nsite","nyear","year","nknots","n","X","Z","nc", "nb", "ni", "nt","zst","y")
jags.data <- list(site = site,
                  survey = survey,
                  nobs = nobs,
                  nrep = nrep,
                  nsite = nsite,
                  nyear = nyear,
                  year = year,
                  nknots = nknots,
                  n = n,
                  X = X,
                  Z = Z,
                  zst = zst,
                  y = y)


# Inits function
f.inits <- function(){list(a=rep(0,nyear), b=rep(0,nyear), c=rep(0,nyear), z=zst)}

f.inits.par <- function(chain){ 
  list(a=rep(0,nyear), b=rep(0,nyear), c=rep(0,nyear), z=zst,
            .RNG.seed=chain, .RNG.name = c("base::Wichmann-Hill","base::Marsaglia-Multicarry",
                                           "base::Super-Duper","base::Mersenne-Twister")[chain])
}
# specify the parameters to be monitored
parameters <- c("a","b","c","b.k","lp","psi","taub")

### Run MCMC Analysis using jags
library(runjags)
jags.out<- autorun.jags(model = "splinesSiteOcc_simple.txt", 
                        data = jags.data, 
                        inits = f.inits.par, 
                        monitor = parameters,
                        n.chains = 3,
                        # burnin = 20000,
                        # sample = 50000,
                        # adapt = 1000,
                        thin = 10,
                        modules = "glm",
                        method = "parallel",
                        summarise = TRUE)
jags.out<-jags(jags.data,f.inits,parameters,"splinesSiteOcc_simple.txt",nc,ni,nb,nt)
out<-jags.out$BUGSoutput



#############################################################
# Try simple model, but with jagam instead of Strebel spline
#############################################################
datsimple <- finaldat %>%
  group_by(SeqID, year) %>%
  mutate(nrep = length(unique(subtransect)),
         ndet = length(which(count > 0))) %>%
  ungroup() %>%
  mutate(day = yday(SiteDate)) %>%
  select(nrep, ndet, SiteID, day, year) %>%
  distinct()

dat <- merge(datsimple, distinct(rawcovs[, c("SiteID", "lat", "lon")]), all.x = TRUE, all.y = FALSE)
dat <- dat %>% filter(year == 2010) %>% data.frame()



dat$failure <- dat$nrep - dat$ndet
dat$propor <- dat$ndet / dat$nrep
# test <- list()
# test$response <- cbind(dat$ndet, dat$nrep - dat$ndet)
# test$day <- dat$day
# test$site <- dat$site

jd <- jagam(propor ~ te(day, lat, k = c(15,5)), family = binomial, 
                 weights = dat$nrep,
                 data = dat, file = "jagam_unmod.txt")

# modify text in model file to account for other things

jd$jags.data$site <- as.numeric(factor(dat$SiteID))
jd$jags.data$nsite <- length(unique(dat$SiteID))
jd$jags.data$survey <- dat$day-min(dat$day)+1
# jd$jags.data$zst <- rep(1, nsite)  

# jd$jags.ini$zst <- jd$jags.data$zst
# jd$jags.ini$y <- 
# 
require(rjags); load.module("glm")

jm <- jags.model("jagamtest.txt", data = jd$jags.data,
                    inits=NULL, n.chains = 2)
jm2 <- update(jm, n.iter = 100)
jm3 <- jags.samples(jm, variable.names = c("b", "rho", "mu", "psi", "z"),
                                        n.iter = 1000, thin = 10)
# sam <- jags.samples(jm, c("b", "rho", "mu", "psi", "z"),
#                     n.iter = 1000, thin = 10)
# jam <- sim2jam(sam, jd$pregam)
# 
# 
# newdata <- data.frame(day = c(min(testdat$day), max(testdat$day)), lat = 39)
# plot(predict.gam(jam, newdata, type = "response"))
# 
# #plot not working
# Xp <- predict(jam, newdata, type = "lpmatrix")
# ii <- 1:25 * 20 + 500
# for (i in 1:25) {
#   fv <- Xp[1, ] %*% sam$b[, 1, 1]
#   if (i==1) plot(days, fv, type = "l", ylim = c(4, 7)) else
#     lines(days, fv)
#  }
library(rjags)
library(jagsUI); load.module("glm"); load.module('dic')


out <- jagsUI::jags(data = jd$jags.data,
            inits = NULL,
            parameters.to.save = c("b", "rho", "mu", "psi", "z"),
            model.file = "jagamtest.txt",
            n.chains = 3,
            n.adapt = 10,
            n.iter = 100,
            n.burnin = 50,
            n.thin = 2)

# using this jags approach works, but output strange, lots of zeros in 'out'
library(R2jags); load.module("glm")

# Simulation parameters
ni=100; nc=2; nb=50; nt=1
jags.out<-R2jags::jags(jd$jags.data, 
               inits = NULL, parameters.to.save = c("b", "rho", "mu", "psi", "z"),
               model.file = "jagamtest.txt", nc,ni,nb,nt)
out<-jags.out$BUGSoutput

saveRDS(out, "jagsmtest_output.rds") # moved to desktop, too big for github

# plot with strebel's code
# Input data
site <- as.numeric(factor(dat$site))
survey <- dat$day-min(dat$day)+1
nobs <- length(unique(paste(dat$site,dat$day,dat$year)))
nrep <- dat$nrep
nsite <- length(unique(dat$site))
nyear <- length(unique(dat$year))
year <- as.numeric(factor(dat$year))
zst <- array(1, dim=c(nsite,nyear))  
y <- dat$ndet



library(jagstools)
jags.out <- jm3

# show all results except for the many N nodes
round(jagsresults(x=jags.out, params='mu', invert=TRUE), 2)

mus <- jagsresults(x=jags.out, param='mu')
# N.ests <- array(data = round(jags.Ns[, "mean"]), dim = c(nsites, nyears, nspecies))
mu_array <- rearray(x = jags.out, param = "mu", fields = "mean")[[1]]

# this shows that it might be working!
plot(colMeans(mu_array, na.rm = TRUE))

rownames(mu_array) <- unique(jd$jags.data$site)
colnames(mu_array) <- min(jd$jags.data$survey):max(jd$jags.data$survey)
mu_df <- as.data.frame(mu_array)
mu_df$site <- unique(jd$jags.data$site) 
mu_df <- tidyr::gather(mu_df, key = "survey", value = "pheno", 1:214)

preddat <- dat
preddat$site <- as.numeric(factor(dat$SiteID))
preddat$survey <- dat$day-min(dat$day)+1

plotdat <- merge(preddat, mu_df, by = c("site", "survey"))


library(ggplot2)

p <- ggplot(data = plotdat, aes(x = day, y = pheno, group = SiteID, color = lat)) +
  geom_line()
p



