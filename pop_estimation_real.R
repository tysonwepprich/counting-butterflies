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
gdd <- readRDS("data/dailyDD.rds")

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
sitemod <- densityMclust(scale(siteGDD[,c(2:3)]), G = 1:4, modelNames = "EEV")
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
params <- expand.grid(species$CommonName, models, cutoff,
                      stringsAsFactors = FALSE)
names(params) <- c("species", "model", "cutoff")
# params <- params %>% arrange(species)
# params <- params[c(5),]

ncores <- 20

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
                    .export = c("data", "surveys", "siteGDD", "covdata", "params", "gdd", "Abund_Curve", 
                                "Adjust_Counts", "CompareMixMods",
                                "Simulate_Counts", "Summ_curve",
                                "Summ_mixmod", "mainDir", "RightNumGen", "AssignGeneration"),
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
                          dat <- dat[-which(dat$cumdegday > 1000 & dat$Total >= 1), ]
                        }
                        if(species == "Leonard's Skipper"){
                          dat <- dat[-which(dat$cumdegday < 1000 & dat$Total >= 1), ]
                        }
                        
                        starttime <- Sys.time()
                        
                        temp <- dat
                        if(sum(temp$Total) < 20|length(unique(temp$SiteID)) < 2|length(unique(temp$Year)) < 2|
                           length(unique(temp$SiteYear))<5|length(unique(temp$RegYear))<2) {
                          mod <- NA
                        }else{
                          modtime <- system.time({ 
                            if(model == "gdd"){
                              mod <- try(gam(Total ~ 
                                               s(zlistlength) +
                                               s(ztemperature) +
                                               s(zduration) +
                                               te(lat, lon, AccumDD, bs = c("tp", "cr"), k = c(5, 30), d = c(2, 1)) +
                                               s(SiteYear, bs = "re"),
                                             family = nb(theta = NULL, link = "log"),
                                             # family = poisson(link = "log"),
                                             data = temp,
                                             method = "REML", 
                                             optimizer = c("outer", "newton"), 
                                             # gamma = 1.4, 
                                             control = list(maxit = 500)))
                            }
                            
                            if(model == "doy"){
                              mod <- try(gam(Total ~ 
                                               s(zlistlength) +
                                               s(ztemperature) +
                                               s(zduration) +
                                               te(lat, lon, DOY, bs = c("tp", "cr"), k = c(5, 30), d = c(2, 1)) +
                                               s(SiteYear, bs = "re"),
                                             family = nb(theta = NULL, link = "log"),
                                             # family = poisson(link = "log"),
                                             data = temp,
                                             method = "REML", 
                                             optimizer = c("outer", "newton"), 
                                             # gamma = 1.4, 
                                             control = list(maxit = 500)))
                            }
                            
                          })
                        }
                      }
                      
                      
                      pars$modtime <- modtime
                      if(is.na(mod) == FALSE){
                        
                        pars$AIC <- AIC(mod)
                        summod <- summary(mod)
                        pars$N <- summod$n
                        pars$dev.expl <- summod$dev.expl
                        pars$negbin <- mod$family$getTheta(TRUE)
                        outlist <- list()
                        outlist[["params"]] <- pars
                        outlist[["gammod"]] <- mod
                        outlist[["datGAM"]] <- temp
                        saveRDS(outlist, paste(species, model, cutoff, "rds", sep = "."))
                        
                      }
                      return(sim)
                    }

