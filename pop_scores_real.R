# take results of pop_estimation_real.R
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
theme_set(theme_bw(base_size = 14)) 


fs <- list.files("OHGAMS/all", full.names = TRUE, recursive = TRUE)

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

###NB: these regions are wrong, but were also modeled with wrong names
# correct below after GAM predictions from models
siteGDD$region <- plyr::mapvalues(siteGDD$region, from = c("1", "2", "3", "4"), 
                                  to = c("NE", "NW", "CN", "SW"))
gdd <- gdd %>% 
  left_join(siteGDD[, c("SiteID", "region")])


traits <- read.csv("data/speciesphenology.csv", header = TRUE) %>% 
  filter(UseMV == "y") %>% 
  select(CommonName, BroodsGAMmin, BroodsGAMmax, UseMV, SyncedBroods, UseMismatch, Model)


outgen <- list()
outweight <- list()
outpops <- list()
outprec <- list()
for (i in 1:length(fs)){
  f <- fs[i]
  pp <- stringr::str_split(string = f, pattern = coll("/"), 3) %>% map(3)
  spp <- stringr::str_split(string = pp, pattern = coll("."), 4) %>% 
    map(1) %>% unlist()
  
  if(spp %in% traits$CommonName){
    tmp <- readRDS(f)
    mod <- tmp$gammod
    if("error" %in% class(tmp$gammod)){
      next
    }
    
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
      # filter(SiteYearTotal >= 20)
      group_by(region) %>% 
      filter(SiteYear %in% sample(unique(SiteYear), 20, replace = TRUE))
    
    preds$region <- plyr::mapvalues(preds$region, from = c("NE", "NW", "CN", "SW"), 
                                    to = c("CN", "NE", "NW", "SW"))
    preds$region <- factor(preds$region, levels = c("NW", "NE", "SW", "CN"))
    
    # outliers
    outs <- tmp$datGAM %>% 
      filter(Total > 0) %>% 
      dplyr::select(AccumDD, DOY, region) %>% 
      filter(complete.cases(.))
    
    outs$region <- plyr::mapvalues(outs$region, from = c("NE", "NW", "CN", "SW"), 
                                   to = c("CN", "NE", "NW", "SW"))
    outs$region <- factor(outs$region, levels = c("NW", "NE", "SW", "CN"))
    
    
    if(tmp$params$model == "doy"){
      gamplt <- ggplot(preds, aes(x = DOY, y = Gamma, group = SiteYear, color = SiteYearGDD)) +
        geom_path(alpha = .5) + 
        scale_color_viridis() + 
        facet_wrap(~region, ncol = 2) +
        geom_rug(data = outs, aes(x = DOY), sides="b", inherit.aes = FALSE, alpha = 0.5) +
        ggtitle(paste0(tmp$params$species, " seasonal phenology"),
                subtitle = "modeled on calendar scale") +
        labs(color = "Total degree-days\n for site and year") +
        labs(x = "Day of year") +
        labs(y = "Scaled phenology (model predictions)")
    }else{
      gamplt <- ggplot(preds, aes(x = AccumDD, y = Gamma, group = SiteYear, color = SiteYearGDD)) +
        geom_path(alpha = .5) + 
        scale_color_viridis() + 
        facet_wrap(~region, ncol = 2) +
        geom_rug(data = outs, aes(x = AccumDD), sides="b", inherit.aes = FALSE,  alpha = 0.5) +
        ggtitle(tmp$params$species,
                subtitle = "modeled on degree-day scale") +
        labs(color = "Total degree-days\n for site and year") +
        labs(x = "Degree-days accumulated (5/30C thresholds)") +
        labs(y = "Scaled phenology (model predictions)")
    }
    ggsave(filename = paste(tmp$params$species, tmp$params$model, "GAM", "png", sep = "."), 
           plot = gamplt, device = "png", path = "plots", width = 8, height = 6, units = "in")
    
    # }
    
    if("N" %in% names(tmp)){
      # regions <- tmp$datGAM %>% 
      #   ungroup() %>% 
      #   dplyr::select(SiteID, region) %>% 
      #   filter(complete.cases(.)) %>% 
      #   distinct()
      # if("region" %in% names(tmp$N)){
      #   N <- tmp$N %>% 
      #     ungroup() %>% 
      #     dplyr::select(-region) %>% 
      #     left_join(regions)
      # }else{
      #   N <- tmp$N %>% 
      #     ungroup() %>% 
      #     left_join(regions)
      # }
      N <- tmp$N %>% ungroup()
      
      N$region <- plyr::mapvalues(N$region, from = c("NE", "NW", "CN", "SW", "ALL"), 
                                  to = c("CN", "NE", "NW", "SW", "ALL"))
      N$region <- factor(N$region, levels = c("NW", "NE", "SW", "CN", "ALL"))
      
      ##### 
      
      # mixture model performed
      if("gen" %in% names(tmp)){
        outgen[[length(outgen)+1]] <- tmp$gen
        # Plot of brood separation
        weights <- N %>% 
          spread(key = metric, value = value) %>% 
          group_by(region, SiteID, Year) %>% 
          mutate(AnnPopIndex = sum(PopIndex),
                 gen_weight = PopIndex / AnnPopIndex) %>% 
          ungroup() %>% 
          mutate(zpopindex = AnnPopIndex / max(AnnPopIndex))
        
        if(tmp$params$model == "doy"){
          genplt <- ggplot(weights, aes(x = curve_q0.5, y = gen_weight, color = as.factor(Gen))) +
            geom_point(aes(alpha = zpopindex)) + 
            # scale_color_viridis(discrete = TRUE) + 
            facet_wrap(~region, ncol = 2) +
            geom_rug(data = outs, aes(x = DOY), sides="b", inherit.aes = FALSE,  alpha = 0.5) +
            ggtitle(paste0(tmp$params$species, " generation size and timing"),
                    subtitle = "modeled on calendar scale") +
            labs(color = "Generation") +
            labs(x = "Day of year") +
            labs(y = "Scaled generation size")
        }else{
          genplt <- ggplot(weights, aes(x = curve_q0.5, y = gen_weight, color = as.factor(Gen))) +
            geom_point(aes(alpha = zpopindex)) + 
            # scale_color_viridis(discrete = TRUE) + 
            facet_wrap(~region, ncol = 2) +
            geom_rug(data = outs, aes(x = AccumDD), sides="b", inherit.aes = FALSE,  alpha = 0.5) +
            ggtitle(paste0(tmp$params$species, " generation size and timing"),
                    subtitle = "modeled on degree-day scale") +
            labs(color = "Generation") +
            labs(x = "Degree-days accumulated (5/30C thresholds)") +
            labs(y = "Scaled generation size")
        }
        ggsave(filename = paste(tmp$params$species, tmp$params$model, "gens", "png", sep = "."), 
               plot = genplt, device = "png", path = "plots", width = 8, height = 6, units = "in")
        
        # Brood weights with zero counts added if missing estimates
        pops <- N %>% 
          ungroup() %>% 
          dplyr::select(region, SiteID, Year, Gen, metric, value) %>%
          filter(metric == "PopIndex") %>% 
          complete(nesting(region, SiteID, Year), Gen, fill = list(metric = "PopIndex", value = 0)) %>% 
          mutate(maxgen = max(Gen)) %>% 
          group_by(SiteID, Year) %>% 
          mutate(AnnPopIndex = sum(value),
                 gen_weight = value / AnnPopIndex,
                 gam_scale = toupper(tmp$params$model),
                 species = tmp$params$species)
        
        
        # CV/PV metric
        
        test1 <- N %>%
          filter(region != "ALL") %>% 
          group_by(Gen, metric, region) %>% 
          summarise(pv = PropVariation(value),
                    cv = sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE),
                    avg = mean(value, na.rm = TRUE))
        
        test2 <- tmp$N %>%
          filter(region == "ALL") %>% 
          group_by(Gen, metric, region) %>% 
          summarise(pv = PropVariation(value),
                    cv = sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE),
                    avg = mean(value, na.rm = TRUE))
        
        precision <- bind_rows(test1, test2) %>% 
          mutate(gam_scale = toupper(tmp$params$model),
                 species = tmp$params$species)
        
        outweight[[length(outweight)+1]] <- weights
        outpops[[length(outpops)+1]] <- pops
        outprec[[length(outprec)+1]] <- precision
        
        
      }else{
        # no mixture model, extract Popindex only
        pops <- N %>% 
          mutate(gam_scale = toupper(tmp$params$model),
                 AnnPopIndex = PopIndex)
        outpops[[length(outpops)+1]] <- pops
        
      } # close if gen exists in tmp
    } # close if N exists in tmp
  }else{
    next
  }
} # close fs loop

allgen <- bind_rows(outgen)
allweight <- bind_rows(outweight)
allpops <- bind_rows(outpops)
allprec <- bind_rows(outprec)
saveRDS(allgen, "MVallgen.rds")
saveRDS(allweight, "MVallweight.rds")
saveRDS(allpops, "MVallpops.rds")
saveRDS(allprec, "MVallprec.rds")


