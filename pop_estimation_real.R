#script combining:
#counts filtering
#gam fitting and summary stats
#mclust broods
source("funcs_pop.R")
# packages to load
library(mclust)
library(lubridate)
library(MASS)
library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(mixsmsn)
library(caret)
library(zoo)

# for parallel simulations with control over seed for reproducibility
# need different packages for windows computers
library(doRNG)
library(foreach) # for parallelized loops
if(.Platform$OS.type == "unix"){
  library(doMC)    # parallel backend for foreach, only for linux/mac
}else if(.Platform$OS.type == "windows"){
  library(doSNOW)
}


data <- readr::read_csv("data/data.trim.csv") %>% 
  mutate(SiteID = formatC(SiteID.x, width = 3, format = "d", flag = "0"),
         SiteDate = lubridate::ymd(SiteDate))
data$CommonName[which(data$CommonName == "Spring/Summer Azure")] <- "Azures"

# filter unidentified species
species <- data %>% 
  filter(CommonName %in% unique(CommonName)[1:122]) %>% 
  group_by(CommonName) %>% 
  summarise(n = sum(Total)) %>% 
  arrange(n)

surveys <- distinct(data[, c("SeqID", "SiteID", "SiteDate", "Week")])

covdata <- data %>%
  group_by(SeqID) %>%
  summarise(listlength = length(which(unique(CommonName) %in% species$CommonName)),
            temperature = mean(c(StartTemp, EndTemp), na.rm = TRUE),
            duration = duration[1]) %>%
  distinct()

sites <- read.csv("data/OHsites_reconciled_update2016.csv") %>% 
  mutate(SiteID = formatC(Name, width = 3, format = "d", flag = "0"))
# gdd <- readRDS("data/dailyDD.rds")
gdd <- readRDS("../ohiogdd/dailyDD.rds")


gdd <- left_join(gdd, sites) %>% 
  dplyr::select(SiteID, SiteDate, degday530, lat, lon, maxT, minT) %>% 
  mutate(Year = year(SiteDate),
         DOY = yday(SiteDate)) %>% 
  group_by(SiteID, Year) %>% 
  arrange(DOY) %>% 
  mutate(AccumDD = cumsum(degday530))

siteGDD <- gdd %>%
  group_by(SiteID, lat, lon) %>% 
  filter(DOY == 365) %>%
  summarise(meanGDD = mean(AccumDD))
sitemod <- densityMclust(siteGDD[,c(2:3)], G = 1:4, modelNames = "EEV")
siteGDD$region <- as.character(sitemod$classification)
# # visualize regional clusters
# a <- ggplot(data = siteGDD, aes(x = lon, y = lat, group = region, color = region)) + geom_point()
# a
siteGDD$region <- plyr::mapvalues(siteGDD$region, from = c("1", "2", "3", "4"), 
                                  to = c("NE", "NW", "CN", "SW"))
gdd <- gdd %>% 
  left_join(siteGDD[, c("SiteID", "region")])

# add spring/fall frost to gdd as covariate
gdd <- gdd %>% 
  group_by(SiteID, Year) %>% 
  arrange(DOY) %>% 
  mutate(frost5day = rollmean(x = minT, 5, align = "right", fill = NA)) %>% 
  mutate(fallfrost = ifelse(DOY > 200, frost5day, NA),
         springfrost = ifelse(DOY < 200, frost5day, NA),
         photo = geosphere::daylength(lat, DOY))

# test <- gdd %>% filter(Year > 2010) 
# plt <- ggplot(test, aes(x = SiteDate, y = springfrost)) +
#   geom_point() +
#   facet_wrap(~region, ncol = 2)
# plt

models <- c("doy", "gdd")
cutoff <- "loose" #c("adapt","strict", "loose")
params <- expand.grid(species$CommonName[100:122], models, cutoff,
                      stringsAsFactors = FALSE)
names(params) <- c("species", "model", "cutoff")
# params <- params %>% arrange(species)
# params <- params[c(5),]

ncores <- 3

if(.Platform$OS.type == "unix"){
  registerDoMC(cores = ncores)
}else if(.Platform$OS.type == "windows"){
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
}

mcoptions <- list(preschedule = FALSE)

# foreach loop
outfiles <- foreach(sim = 1:nrow(params),
                    .combine='c',
                    .packages= c("mgcv", "dplyr", "tidyr", "purrr",
                                 "lubridate", "mclust", "mixsmsn"),
                    .export = c("data", "surveys", "siteGDD", "covdata", "params", "gdd"),
                    .inorder = FALSE,
                    .options.multicore = mcoptions) %dopar% {
                      
                      species <- params$species[sim]
                      model <- params$model[sim]
                      cutoff <- params$cutoff[sim]
                      
                      # print(species)
                      # print(Sys.time())
                      
                      reduced <- NA
                      pars <- data.frame(species, model, cutoff, reduced)
                      
                      counts <- data %>% 
                        filter(CommonName == species) %>% 
                        mutate(DOY = yday(SiteDate),
                               Year = year(SiteDate))
                      
                      #get unique surveys, including those where species not counted
                      survs <- surveys %>% 
                        filter(year(SiteDate) %in% unique(counts$Year)) %>% 
                        filter(SiteID %in% unique(counts$SiteID)) %>% 
                        mutate(Year = year(SiteDate))
                      
                      #Add zeros to surveys when species not counted during a survey
                      test <- left_join(survs, counts, by = c("SiteID", "SiteDate", "Week", "SeqID", "Year"))
                      counts <- test[, c("SiteID", "SiteDate", "Week", "SeqID", "Total", "Year")]
                      counts$Total <- plyr::mapvalues(counts$Total, from = NA, to = 0)
                      counts <- left_join(counts, siteGDD, by = "SiteID")
                      counts <- left_join(counts, covdata, by = "SeqID")
                      counts$temperature[which(counts$temperature < 50)] <- NA
                      counts$duration[which(counts$duration == 0)] <- NA
                      
                      #scaling covariates
                      #scaled by region9/week
                      counts <- counts %>% 
                        group_by(month(SiteDate)) %>% 
                        mutate(zlistlength = as.numeric(scale(log(listlength+1))),
                               ztemperature = as.numeric(scale(temperature))) %>% 
                        group_by(SiteID) %>% 
                        mutate(zduration = as.numeric(scale(duration)))
                      
                      
                      # trying to add GDD instead of ordinal date
                      counts <- left_join(counts, gdd, by = c("SiteID", "SiteDate", "Year", "lat", "lon", "region")) %>% 
                        group_by(SiteID, Year) %>% 
                        mutate(SurvPerYear = length(unique(SeqID)),
                               YearTotal = sum(Total))
                      
                      
                      if(cutoff == "strict"){
                        datGAM <- counts[YearTotal >= 3]
                        datGAM <- datGAM[SurvPerYear >= 15]
                      }
                      if(cutoff == "loose"){
                        datGAM <- counts %>% filter(YearTotal >= 1, SurvPerYear >= 10)
                      }
                      # #to decrease computation of common species, adapt filter removes low data siteyears
                      # if(cutoff == "adapt"){
                      #   datGAM <- counts[YearTotal >= 1]
                      #   datGAM <- datGAM[SurvPerYear >= 11] 
                      #   #this SurvPerYear cutoff chosen to include 95% of survey effort
                      #   #16 surveys would be including 90% of surveys 
                      #   counts1 <- datGAM %>% 
                      #     group_by(SiteID, Year) %>% 
                      #     mutate(SiteYearSurvs = length(unique(SeqID)))
                      #   counts1 <- counts1 %>% 
                      #     group_by(SiteID) %>% 
                      #     mutate(SpeciesSiteYears = length(unique(Year)))
                      #   counts2 <- counts1 %>% 
                      #     group_by(SiteID, Year, SiteYearSurvs, SpeciesSiteYears) %>% 
                      #     summarise(PresenceSiteYear = length(unique(SeqID[Total>=1])),
                      #               SiteYearTotal = sum(Total))
                      #   
                      #   filtered <- counts2 %>% 
                      #     filter(SpeciesSiteYears >= quantile(counts2$SpeciesSiteYears, 0.1),
                      #            PresenceSiteYear >= quantile(counts2$PresenceSiteYear, 0.1),
                      #            SiteYearTotal >= quantile(counts2$SiteYearTotal, 0.1)) %>% 
                      #     mutate(SiteYear = paste(SiteID, Year, sep = "_"))
                      #   
                      #   datGAM <- datGAM %>% 
                      #     filter(SiteYear %in% unique(filtered$SiteYear))
                      # }
                      
                      if(nrow(datGAM) == 0){
                        mod <- NA
                      }else{
                        # datGAM[, SitesObserved := length(unique(SiteID)), by = list(Year)]
                        # datGAM <- datGAM[SitesObserved >= 1]
                        
                        # without leading zeros
                        dat <- datGAM %>% 
                          group_by(SiteID, Year) %>%
                          mutate(YearTotal = sum(Total),
                                 SurvSeen = length(which(Total > 0))) %>%
                          filter(YearTotal > 0,
                                 SurvSeen > 0)
                        
                        dat$Year <- as.factor(as.character(dat$Year))
                        dat$region <- as.factor(as.character(dat$region))
                        dat$SiteID <- as.factor(as.character(dat$SiteID))
                        dat$SiteYear <- as.factor(paste(dat$SiteID, dat$Year, sep = "_"))
                        dat$zlistlength[which(is.na(dat$zlistlength))] <- 0
                        dat$ztemperature[which(is.na(dat$ztemperature))] <- 0
                        dat$zduration[which(is.na(dat$zduration))] <- 0
                        dat$RegYear <- as.factor(paste(dat$region, dat$Year, sep = "_"))
                        dat <- as.data.frame(dat)
                        #silly filters for univoltine species with outliers
                        
                        if(species == "Baltimore"){
                          dat <- dat[-which(dat$DOY > 220 & dat$Total >= 1), ]
                        }
                        if(species == "Leonard's Skipper"){
                          dat <- dat[-which(dat$DOY < 220 & dat$Total >= 1), ]
                        }
                        
                        starttime <- Sys.time()
                        
                        temp <- dat
                        if(sum(temp$Total) < 20|length(unique(temp$SiteID)) < 2|length(unique(temp$Year)) < 2|
                           length(unique(temp$SiteYear))<5|length(unique(temp$RegYear))<2) {
                          mod <- NA
                        }else{
                          modtime <- system.time({ 
                            safe_gam <- purrr::safely(gam)
                            
                            if(model == "gdd"){
                              mod <- safe_gam(Total ~ 
                                                s(zlistlength) +
                                                s(ztemperature) +
                                                te(lat, lon, AccumDD, bs = c("tp", "cr"), k = c(5, 30), d = c(2, 1)) +
                                                s(SiteID, bs = "re") +
                                                s(RegYear, AccumDD, bs = "fs", k = 5, m = 1),
                                              family = nb(theta = NULL, link = "log"),
                                              # family = poisson(link = "log"),
                                              data = temp,
                                              method = "REML", 
                                              optimizer = c("outer", "newton"), 
                                              # gamma = 1.4, 
                                              control = list(maxit = 500))
                              
                            }
                            
                            if(model == "doy"){
                              mod <- safe_gam(Total ~ 
                                                s(zlistlength) +
                                                s(ztemperature) +
                                                te(lat, lon, DOY, bs = c("tp", "cr"), k = c(5, 30), d = c(2, 1)) +
                                                s(SiteID, bs = "re") +
                                                s(RegYear, DOY, bs = "fs", k = 5, m = 1),
                                              family = nb(theta = NULL, link = "log"),
                                              # family = poisson(link = "log"),
                                              data = temp,
                                              method = "REML", 
                                              optimizer = c("outer", "newton"), 
                                              # gamma = 1.4, 
                                              control = list(maxit = 500))
                            }
                            
                          })
                        }
                      }
                      
                      if(is.na(mod)[1]){
                        return(sim)
                      }else{
                        
                        if(is.null(mod$error)){
                          pars$modtime <- modtime[1]
                          pars$AIC <- AIC(mod$result)
                          summod <- summary(mod$result)
                          pars$N <- summod$n
                          pars$dev.expl <- summod$dev.expl
                          pars$negbin <- mod$result$family$getTheta(TRUE)
                          outlist <- list()
                          outlist[["params"]] <- pars
                          outlist[["gammod"]] <- mod$result
                          outlist[["datGAM"]] <- temp
                          saveRDS(outlist, paste(species, model, cutoff, "rds", sep = "."))
                          
                        }else{
                          outlist <- list()
                          outlist[["params"]] <- pars
                          outlist[["gammod"]] <- mod$error
                          outlist[["datGAM"]] <- temp
                          saveRDS(outlist, paste("gamerr", species, model, cutoff, "rds", sep = "."))
                        }
                        return(sim)
                      }
                    }
if(.Platform$OS.type == "windows"){
  stopCluster(cl)
}



# Now that GAMs fit, run through mixture models as in simulation

traits <- read.csv("data/speciesphenology.csv", header = TRUE) %>% 
  filter(UseMismatch == "y") %>% 
  select(CommonName, BroodsGAMmin, BroodsGAMmax, UseMV, SyncedBroods, UseMismatch, Model)

fs <- list.files("C:/Users/Tyson/Desktop/OHGAMS", full.names = TRUE)
fs <- fs[grep(pattern = "loose", x = fs, fixed = TRUE)]

# mixmod parameters to run for each species
mixmod <- c("hom", "het", "skew")
mod_region <- c("ALL", "reg4")
params <- list(mixmod = mixmod, mod_region = mod_region)
params <- expand.grid(params)

#for each fs

f <- fs[1]
spp <- stringr::str_split(string = fs, pattern = coll("/"), 6)
spp <- stringr::str_split(string = spp, pattern = coll("."), 4)[[1]][1]
if(spp %in% traits$CommonName){
  
  tmp <- readRDS(f)
  pars <- tmp$params
  mod <- tmp$gammod
  counts <- tmp$datGAM
  
  params$gam_scale <- toupper(x = pars$model)
  params$ngen <- traits$BroodsGAMmax[which(traits$CommonName == spp)]
  
  adjcounts <- gdd %>% 
    mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>% 
    filter(SiteYear %in% unique(counts$SiteYear)) %>%
    group_by(SiteID, Year) %>% 
    filter(DOY %in% (seq(90, 305, 4) + sample.int(n=3, size=54, replace=TRUE))) %>% 
    ungroup() %>% 
    mutate(zlistlength = 0,
           ztemperature = 0,
           zduration = 0,
           RegYear = paste(region, Year, sep = "_" )) 
  
  # prediction with simulated counts, stochastic but integers (if n is odd)
  Xp <- predict.gam(object = mod, newdata = adjcounts, type="lpmatrix") ## map coefs to fitted curves
  beta <- coef(mod)
  Vb   <- vcov(mod) ## posterior mean and cov of coefs
  n <- 5 # choose number of simulations
  mrand <- MASS::mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
  ilink <- family(mod)$linkinv
  linklppreds <- Xp %*% t(mrand)
  nbpreds <- apply(X = linklppreds,
                   MARGIN = 1,
                   FUN = function(x){
                     # temp <- sort(x)
                     # bounds <- quantile(1:n, probs = c(0.025, 0.975))
                     # x <- temp[bounds[1]:bounds[2]]
                     x <- ilink(x)
                     x <- rnbinom(n = length(x),
                                  mu = x,
                                  size = mod$family$getTheta(TRUE))
                     x <- quantile(x, .5)
                     return(x)
                   })
  adjcounts$adjY <- nbpreds
  
  if (param$mod_region == "ALL"){
    adjcounts$region <- "ALL"
  }
  
  results <- adjcounts %>%
    group_by(region) %>% 
    do(mixmods = CompareMixMods(dat = ., param = param))
  
}else{
  next
}

# This df could be used to see if BIC/AIC choose right number of gen
genright <- results %>% 
  do(rightgen = RightNumGen(.$mixmods, param = param, reg = .$region)) %>% 
  unnest()

# regroup generation assignments by SiteYear
siteresults <- results %>%
  do(GenClass = AssignGeneration(mixmod = .$mixmods, 
                                 dat = adjcounts, 
                                 param = param,
                                 reg = .$region)) %>% 
  unnest()


outlist <- list(param, counts, adjcounts, truth, genright, siteresults)   
saveRDS(object = outlist, file = paste("popest", param$index, sep = "_"))

f <- fs[sim]
tmp <- readRDS(f)
param <- tmp[[1]]
counts <- tmp[[2]]
adjcounts <- tmp[[3]]
truth <- tmp[[4]]
genright <- tmp[[5]]
siteresults <- tmp[[6]]

# 1. when is number of generations correct compared to params? 
# What others errors, like degenerate modes with redundant mu?
# Compare across models for best fit to true ngen, compare within models to see if correct ngen selected
ngen <- genright %>% 
  group_by(region) %>% 
  mutate(aicpick = maxgen[which(within_model_aic == 0)][1],
         bicpick = maxgen[which(within_model_bic == 0)][1]) %>% 
  filter(right_ngen == "yes") %>% 
  slice(1L) %>% 
  dplyr::select(maxgen, mixmod_flag, badmixmod:bicpick)

trueweight <- truth %>% 
  mutate(Year = as.character(Year)) %>% 
  left_join(adjcounts[, c("SiteID", "Year", "region")]) %>% 
  filter(Gen == max(Gen)) %>% 
  group_by(region) %>% 
  summarise(true_weight = mean(gen_weight))

ngen_score <- bind_cols(param[rep(seq_len(nrow(param)), each=nrow(ngen)),], ngen) %>% 
  left_join(trueweight)

tmp[[7]] <- ngen_score

# 2. what is error in generation weights estimated compared to simulated "truth"?
# 3. what is error in phenology for each site/year compared to "truth" (quantiles/max/mean)?


# Do curve summaries in terms of GDD for comparison
siteresults$Gen <- siteresults$gen
siteresults$gen <- NULL
phen_est <- siteresults %>% 
  dplyr::select(-c(curve_mean, curve_max, curve_q0.1, curve_q0.5, curve_q0.9, n)) %>% 
  filter(complete.cases(.)) %>% 
  filter(maxgen == param$ngen) %>% 
  ungroup() %>% 
  group_by(SiteID, Year, Gen, AccumDD) %>% 
  summarise(Total = sum(count)) %>% 
  group_by(SiteID, Year, Gen) %>% 
  do(Summ_curve(t = .$AccumDD, y = .$Total)) %>% 
  group_by(SiteID, Year) %>%
  mutate(gen_weight = estN / sum(estN))

truth <- truth %>% 
  ungroup() %>% 
  mutate(PopIndex = gen_weight * M,
         Data = "truth") %>% 
  dplyr::select(-n, -M)
# 3a. regression of random effects for site/year compared to estimated variation?
# maybe not necessary/interesting...


# 4. what is error in site/year population size compared to "truth" using trapezoid rule?
# do population index in terms of DOY like UKBMS
N_est <- siteresults %>% 
  dplyr::select(-c(curve_mean, curve_max, curve_q0.1, curve_q0.5, curve_q0.9, n)) %>% 
  filter(complete.cases(.)) %>% 
  filter(maxgen == param$ngen) %>% 
  ungroup() %>% 
  group_by(SiteID, Year, Gen, DOY) %>% 
  summarise(Total = sum(count)) %>% 
  group_by(SiteID, Year, Gen) %>% 
  summarise(PopIndex = TrapezoidIndex(DOY, Total)) %>% 
  right_join(phen_est) %>% 
  group_by(Gen) %>% 
  mutate(meanmu = mean(curve_mean)) %>% 
  group_by(SiteID, Gen) %>% 
  mutate(Site_RE = mean(curve_mean) - meanmu) %>% 
  group_by(Year, Gen) %>% 
  mutate(Year_RE = mean(curve_mean) - meanmu,
         Data = "estimate") %>% 
  dplyr::select(-estN, -meanmu) %>% 
  ungroup() %>% 
  mutate(Year = as.numeric(as.character(Year)))

alldat <- bind_rows(truth, N_est)

allscore <- alldat %>% 
  dplyr::select(-Site_RE, -Year_RE) %>% 
  tidyr::gather(key = "metric", value = "value", gen_weight:PopIndex) %>% 
  group_by(SiteID, Year, Gen, metric) %>% 
  summarise(error = value[1] - value[2])

# NAs if a generation not counted at a particulr site due to low numbers
sumscore <- allscore %>% 
  ungroup() %>% 
  dplyr::filter(complete.cases(.)) %>% 
  group_by(metric) %>% 
  summarise(rmse = sqrt(mean(error^2)), 
            mae = mean(abs(error))) %>% 
  ungroup() %>% 
  mutate(Gen = "ALL")

genscore <- allscore %>% 
  ungroup() %>% 
  filter(complete.cases(.)) %>% 
  group_by(Gen, metric) %>% 
  summarise(rmse = sqrt(mean(error^2)), 
            mae = mean(abs(error))) %>% 
  ungroup() %>% 
  mutate(Gen = as.character(Gen)) %>% 
  bind_rows(sumscore)

score <- bind_cols(param[rep(seq_len(nrow(param)), each=nrow(genscore)),], genscore)
tmp[[8]] <- score

# complete combos of data for last generations missing
relpop <- alldat %>% 
  tidyr::complete(SiteID, Year, Gen, Data, fill = list(PopIndex = 0)) %>% 
  group_by(Gen) %>% 
  arrange(SiteID, Year) %>% 
  summarise(corrpop = cor(PopIndex[Data == "truth"], PopIndex[Data == "estimate"]))
popscore <- bind_cols(param[rep(seq_len(nrow(param)), each=nrow(relpop)),], relpop)
tmp[[9]] <- popscore

saveRDS(tmp, file = f)

