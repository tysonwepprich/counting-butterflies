# take results of pop_estimation_real.R
# bind all species results together for gen and N
# for N, measure precision as in Cayton et al. 
# errors for mixmod and mod_region to determine best options
# tendency to select ngen for mixmod and mod_region (without judgement on whether right?)
# could plot GAM predictions
# could plot generation weights/means to check for separation
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
library(stringr)
library(viridis)
theme_set(theme_bw(base_size = 12)) 


fs <- list.files("OHGAMS/separates", full.names = TRUE)
pp <- stringr::str_split(string = fs, pattern = coll("/"), 3) %>% map(3)
spp <- stringr::str_split(string = pp, pattern = coll("."), 4) %>% 
  map(1) %>% 
  unlist() %>% 
  unique()

gdd <- readRDS("data/dailyDD.rds")
# gdd <- readRDS("../ohiogdd/dailyDD.rds")

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
siteGDD$region <- plyr::mapvalues(siteGDD$region, from = c("1", "2", "3", "4"), 
                                  to = c("NE", "NW", "CN", "SW"))
gdd <- gdd %>% 
  left_join(siteGDD[, c("SiteID", "region")])

gdd$region <- factor(gdd$region, levels = c("NW", "NE", "SW", "CN"))



# test out PV metric

PropVariation <- function(numvect){
  reldiff <- sapply(numvect, function(x) sapply(numvect, function(y) 1 - min(x, y) / max(x, y)))
  pv <- mean(as.numeric(reldiff[upper.tri(reldiff)]))
  return(pv)
}

test2 <- tmp$N %>%
  filter(Data == gam_scale) %>% 
  group_by(Gen, mixmod, region, gam_scale) %>% 
  summarise(pv = PropVariation(curve_mean),
            cv = sd(curve_mean) / mean(curve_mean))

system.time({
  genscore <- vector("list", length(fs))
  siteyrscore <- vector("list", length(fs))
  relpopscore <- vector("list", length(fs))
  for (i in 1:length(fs)){
    f <- fs[i]
    tmp <- readRDS(f)
    mod <- tmp$gammod
    
    counts <- tmp$datGAM %>% 
      filter(YearTotal > 1, SurvSeen > 1) %>% 
      droplevels()
    
    ##### 
    # plot GAM predictions for species/model
    preds <- gdd %>% 
      mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>% 
      filter(SiteYear %in% unique(counts$SiteYear)) %>%
      group_by(SiteID, Year) %>% 
      mutate(SiteYearGDD = max(AccumDD)) %>% 
      filter(DOY %in% seq(90, 305, 1)) %>% 
      ungroup() %>% 
      mutate(zlistlength = 0,
             ztemperature = 0,
             zduration = 0,
             RegYear = paste(region, Year, sep = "_" )) 
    preds$adjY <- predict.gam(object = mod, newdata = preds, type="response")
    
    preds <- preds %>% 
      group_by(SiteYear) %>%
      mutate(Gamma = adjY / as.vector(sum(adjY)),
             SiteYearTotal = sum(adjY)) %>%
      filter(SiteYearTotal >= 20)
    
    if(tmp$params$model == "doy"){
      gamplt <- ggplot(preds, aes(x = DOY, y = Gamma, group = SiteYear, color = SiteYearGDD)) +
        geom_path(alpha = .5) + 
        scale_color_viridis() + 
        facet_wrap(~region, scales = "free_y") +
        ggtitle(tmp$params$species)
    }else{
      gamplt <- ggplot(preds, aes(x = AccumDD, y = Gamma, group = SiteYear, color = SiteYearGDD)) +
        geom_path(alpha = .5) + 
        scale_color_viridis() + 
        facet_wrap(~region, scales = "free_y") +
        ggtitle(tmp$params$species)
    }
    ggsave(filename = paste(tmp$params$species, tmp$params$model, "GAM", "png", sep = "."), 
           plot = gamplt, device = "png", path = "plots", width = 8, height = 6, units = "in")
    
    
    ##### 
    # Brood weights
    # 
    # do( complete(., SiteYear, nsim, brood, 
    #              fill = list(num = 0, weight = 0, mu = NA, sigma = NA, sim = 0))) %>% 
    # 
    weights <- tmp$N %>% 
      filter(Data == gam_scale,
             mod_region == "ALL",
             mixmod == "hom") %>% 
      dplyr::select(SiteID, Year, Gen, PopIndex, region, ngen, species)
      group_by(SiteID, Year) %>% 
      mutate(EstPerGen = ifelse(max(Gen) > length(which(PopIndex > 0)),
                                "no", "yes"),
             DevFromNgen = max(Gen) - ngen,
             SiteYearTotal = sum(PopIndex)) %>% 
      filter(EstPerGen == "yes",
             SiteYearTotal >= 10)
      
    
    # plots brood weights by region for each SiteYear
    plt <- ggplot(N, aes(x = curve_q0.5, y = gen_weight, color = as.factor(Gen))) +
      geom_point() +
      facet_grid(region~mixmod)
    plt
    
    
    
    genscore[[i]] <- tmp[[7]]
    siteyrscore[[i]] <- tmp[[8]]
    relpopscore[[i]] <- tmp[[9]]
  }
  gendf <- bind_rows(genscore)
  phendf <- bind_rows(siteyrscore)
  popdf <- bind_rows(relpopscore)
})
saveRDS(gendf, "gendf.rds")
saveRDS(phendf, "phendf.rds")
saveRDS(popdf, "popdf.rds")







