# ON hold until I get new data through 2016
# then rerun GAMs and do population index all at once

#script combining:
#counts filtering
#gam predictions for a collated index
#differs from poptrends.R by accounting for zero counts at sites
# use sites where sum(total) >= 5, YearSeen at site > 1
# otherwise Site's population indices very low and bring down collated index

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
library(stringr)

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

# NB: these regions are wrong, but GAM's fit with them
siteGDD$region <- plyr::mapvalues(siteGDD$region, from = c("1", "2", "3", "4"), 
                                  to = c("NE", "NW", "CN", "SW"))
gdd <- gdd %>% 
  left_join(siteGDD[, c("SiteID", "region")])


fs <- list.files("OHGAMS/all", full.names = TRUE)
# for species not using mixture models
# still want population index using GAMs
outN <- list()
for (sim in 1:length(fs)){
  f <- fs[sim]
  print(f)
  spp <- stringr::str_split(string = f, pattern = coll("/"), 3) %>% map(3)
  spp <- stringr::str_split(string = spp, pattern = coll("."), 4) %>% map(1) %>% unlist()
  tmp <- readRDS(f)
  pars <- tmp$params
  mod <- tmp$gammod
  # counts <- tmp$datGAM # this doesn't include zeros
  
  if("error" %in% class(mod)){
    next
  }
  
  counts <- data %>% 
    filter(CommonName == spp) %>% 
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
  
  dat <- counts %>% filter(YearTotal >= 0, SurvPerYear >= 10)
  
  dat$Year <- as.factor(as.character(dat$Year))
  dat$region <- as.factor(as.character(dat$region))
  dat$SiteID <- as.factor(as.character(dat$SiteID))
  dat$SiteYear <- as.factor(paste(dat$SiteID, dat$Year, sep = "_"))
  dat$zlistlength[which(is.na(dat$zlistlength))] <- 0
  dat$ztemperature[which(is.na(dat$ztemperature))] <- 0
  dat$zduration[which(is.na(dat$zduration))] <- 0
  dat$RegYear <- as.factor(paste(dat$region, dat$Year, sep = "_"))
  dat <- as.data.frame(dat)
  
  # silly filters for univoltine species with outliers
  # remove if I'm not doing this systematically for others???
  if(spp == "Baltimore"){
    dat <- dat[-which(dat$DOY > 220 & dat$Total >= 1), ]
  }
  if(spp == "Leonard's Skipper"){
    dat <- dat[-which(dat$DOY < 220 & dat$Total >= 1), ]
  }
  
  # simplified gam with fixed effects of Site and Year like UKBMS models
  # one option, fit gam with SiteID and Year
  mod <- gam(Total ~ 
                SiteID + Year +
                     s(zlistlength) +
                     s(ztemperature) +
                     te(lat, lon, AccumDD, bs = c("tp", "cr"), k = c(5, 30), d = c(2, 1)) +
                     # s(SiteID, bs = "re") +
                     s(DOY) +
                     s(RegYear, AccumDD, bs = "fs", k = 5, m = 1),
                   family = nb(theta = NULL, link = "log"),
                   # family = poisson(link = "log"),
                   data = dat,
                   method = "REML", 
                   optimizer = c("outer", "newton"), 
                   # gamma = 1.4, 
                   control = list(maxit = 500))
  
  
  
  
  counts <- tmp$datGAM # this doesn't include zeros
  # GAM predictions to use as offset
  preds <- gdd %>%
    mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>%
    filter(SiteYear %in% unique(dat$SiteYear)) %>%
    group_by(SiteID, Year) %>%
    filter(DOY %in% (seq(90, 305, 4) + sample.int(n=3, size=54, replace=TRUE))) %>%
    ungroup() %>%
    mutate(zlistlength = 0,
           ztemperature = 0,
           zduration = 0,
           RegYear = paste(region, Year, sep = "_" ))
  
  preds$adjY <- predict.gam(object = mod2, newdata = preds, type="response")
  # tmp[["preds"]] <- preds
  
  N2 <- preds %>%
    group_by(SiteID, Year) %>%
    summarise(PopIndex = TrapezoidIndex(DOY, adjY)) %>%
    mutate(Year = as.numeric(as.character(Year)),
           species = spp) %>% 
    
    ggplot(data = N, aes(x = Year, y = PopIndex)) +
    geom_point() +
    facet_wrap(~SiteID, scales = "free")
    
  
  outN[[length(outN)+1]] <- N