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

fs <- list.files("OHGAMS/separates", full.names = TRUE)
pp <- stringr::str_split(string = fs, pattern = coll("/"), 3) %>% map(3)
spp <- stringr::str_split(string = pp, pattern = coll("."), 4) %>% 
  map(1) %>% 
  unlist() %>% 
  unique()




system.time({
  genscore <- vector("list", length(fs))
  siteyrscore <- vector("list", length(fs))
  relpopscore <- vector("list", length(fs))
  for (i in 1:length(fs)){
    f <- fs[i]
    tmp <- readRDS(f)
    
    preds <- tmp$preds #%>% 
      # group_by(SiteYear) %>% 
      # mutate(Gamma = adjY / sum(adjY),
      #        SiteYearTotal = sum(adjY)) %>% 
      # group_by(RegYear, DOY) %>% 
      # mutate(meanGamma = mean(Gamma))
    
    
    
    gamplt <- ggplot(preds, aes(x = DOY, y = adjY, group = SiteYear, color = region)) +
      geom_point(alpha = .5) + 
      facet_wrap(~Year, scales = "free_y")
    gamplt
    
    N <- tmp$N %>% 
      filter(Data == "DOY") %>% 
      group_by(SiteID, Year, region, gam_scale, mixmod) %>% 
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







