# bind all species results together for gen and N
# for N, measure precision as in Cayton et al. 
# errors for mixmod and mod_region to determine best options
# tendency to select ngen for mixmod and mod_region (without judgement on whether right?)
# could plot GAM predictions
# could plot generation weights/means to check for separation

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
library(stringr)
library(viridis)
theme_set(theme_bw(base_size = 16)) 

# gdd <- readRDS("data/dailyDD.rds")
gdd <- readRDS("../ohiogdd/dailyDD.rds")

sites <- read.csv("data/OHsites_reconciled_update2016.csv") %>% 
  mutate(SiteID = formatC(Name, width = 3, format = "d", flag = "0"))

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

###NB: these regions are right
siteGDD$region <- plyr::mapvalues(siteGDD$region, from = c("1", "2", "3", "4"), 
                                  to = c("CN", "NE", "NW", "SW"))
gdd <- gdd %>% 
  left_join(siteGDD[, c("SiteID", "region")])

gddyear <- gdd %>% 
  mutate(season = ifelse(month(SiteDate) %in% c(12, 1, 2), "Winter",
                         ifelse(month(SiteDate) %in% c(3, 4, 5), "Spring",
                                ifelse(month(SiteDate) %in% c(6, 7, 8), "Summer", "Fall")))) %>% 
  group_by(season, Year) %>% 
  summarise(meanT = mean((maxT + minT)/ 2))

gddyear$season <- factor(gddyear$season, levels = c("Winter", "Spring", "Summer", "Fall"))

tempts <- ggplot(gddyear, aes(Year, meanT)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~season, ncol = 2, scales = "free_y") +
  labs(y = "Mean daily temperature (Celsius)")

tempts

gddgdd <- gdd %>% 
  filter(DOY == 365) %>% 
  group_by(region, Year) %>% 
  summarise(meanGDD = mean(AccumDD))
gddpts <- ggplot(gddgdd, aes(Year, meanGDD)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~region, ncol = 2)
gddpts

# messy, but looks like spring/early summer have increasing trend in GDD
# late summer and fall might show less GDD lately
gddmonth <- gdd %>% 
  mutate(Month = month(SiteDate)) %>% 
  group_by(region, SiteID, Year, Month) %>% 
  summarise(totgdd = sum(degday530)) %>% 
  group_by(region, Year, Month) %>% 
  summarise(meangdd = mean(totgdd))
gddpts <- ggplot(gddmonth, aes(Year, meangdd)) +
  geom_point() +
  geom_smooth() +
  facet_grid(Month~region, scales = "free_y")
gddpts


