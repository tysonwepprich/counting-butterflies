# plot GAM to check voltinism fit
# for gdd and doy


fs <- list.files("C:/Users/Tyson/Desktop/OHGAMS", full.names = TRUE)
fs <- fs[grep(pattern = "loose", x = fs, fixed = TRUE)]
f <- fs[171]

# fits <- fs %>% 
#   map_df(function(x){tmp <- readRDS(x); return(tmp$params)})
# 
# fits2 <- fits %>% 
#   filter(years == "all") %>% 
#   group_by(species) %>% 
#   arrange(model) %>% 
#   summarise(diffAIC = AIC[1] - AIC[2],
#             diffdev = dev.expl[1] - dev.expl[2])

tmp <- readRDS(f)
params <- tmp$params
mod <- tmp$gammod
counts <- tmp$datGAM

# ran front stuff from pop_estimation_real.R for gdd and sites
# Add sites for predicting, centroids of 4 regions
library(geosphere)

modsites <- counts %>% 
  select(SiteID, lon, lat) %>% 
  distinct()
sitedist <- distm (cbind(modsites$lon, modsites$lat), t(sitemod$parameters$mean[c(2, 1), ]), fun = distHaversine)
pltsite <- modsites[apply(X = sitedist, MARGIN = 2, FUN = function(x) which(x == min(x, na.rm = TRUE))[1]), ]

pred <- gdd %>% 
  mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>% 
  filter(SiteID %in% pltsite$SiteID) %>%
  filter(SiteYear %in% unique(counts$SiteYear)) %>%
  group_by(SiteID, Year) %>% 
  filter(DOY %in% (seq(90, 305, 4) + sample.int(n=3, size=54, replace=TRUE))) %>% 
  ungroup() %>% 
  mutate(zlistlength = 0,
         ztemperature = 0,
         zduration = 0,
         # SiteYear = counts$SiteYear[1],
         RegYear = paste(region, Year, sep = "_" )) #dummy placeholder for prediction

pred$pred <- predict(mod, newdata = pred, type = "response", exclude = "s(SiteYear)")

plt <- ggplot(pred, aes(x = DOY, y = pred, group = SiteYear, color = Year)) +
  geom_path() + 
  facet_wrap(~region, scales = "free_y")
plt

plt <- ggplot(counts, aes(x = AccumDD, y = Total, group = SiteYear, color = Year)) +
  geom_point() + 
  facet_wrap(~region, scales = "free_y")
plt

  
