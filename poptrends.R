# population trends with GAM prediction based indices

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
  counts <- tmp$datGAM
  
  if("error" %in% class(mod)){
    next
  }
  # issue with near infinite gam predictions,
  # try removing siteyears with very few counts
  # happened at beginning of season, might need zero anchoring
  counts <- counts %>%
    filter(YearTotal > 1, SurvSeen > 1) %>%
    droplevels()
  
  preds <- gdd %>%
    mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>%
    filter(SiteYear %in% unique(counts$SiteYear)) %>%
    group_by(SiteID, Year) %>%
    filter(DOY %in% (seq(90, 305, 4) + sample.int(n=3, size=54, replace=TRUE))) %>%
    ungroup() %>%
    mutate(zlistlength = 0,
           ztemperature = 0,
           zduration = 0,
           RegYear = paste(region, Year, sep = "_" ))
  
  preds$adjY <- predict.gam(object = mod, newdata = preds, type="response")
  # tmp[["preds"]] <- preds
  
  N <- preds %>%
    group_by(SiteID, Year) %>%
    summarise(PopIndex = TrapezoidIndex(DOY, adjY)) %>%
    mutate(Year = as.numeric(as.character(Year)),
           species = spp)
  
  outN[[length(outN)+1]] <- N
  
}


pops <- bind_rows(outN)
saveRDS(pops, "allpops.rds")

pops <- readRDS("allpops.rds")
pops$YearFact <- as.factor(as.character(pops$Year))
pops$SiteID <- as.factor(as.character(pops$SiteID))
pops$zyear <- scale(pops$Year)[,1]

# a <- ggplot(pops, aes(x = Year, y = PopIndex)) +
  # geom_point()+
  # geom_smooth() +
  # facet_wrap(~species, scales = "free_y")
# a

CollInd <- function(temp){
  temp <- temp %>% droplevels()
  mod <- glm(round(PopIndex) ~ YearFact + SiteID - 1, 
             data = temp, family = poisson(link = "log"))
  out <- data.frame(Year = levels(temp$YearFact), 
                    Index = coef(mod)[1:length(levels(temp$YearFact))],
                    zIndex = scale(coef(mod)[1:length(levels(temp$YearFact))])[,1])
}


popmod <- pops %>%  
  filter(Year != 1995) %>%
  group_by(species, SiteID) %>% 
  mutate(yrpersite = length(unique(Year)),
         yrsincestart = Year - min(Year)) %>% 
  group_by(species, Year) %>% 
  mutate(siteperyr = length(unique(SiteID))) %>% 
  filter(yrpersite > 1,
         siteperyr > 1) %>%
  droplevels() %>% 
  group_by(species) %>% 
  mutate(uniqyr = length(unique(Year))) %>%
  filter(uniqyr >= 5) %>%
  do(., CollInd(.)) %>% 
  mutate(Year = as.numeric(as.character(Year)))

poptrend <- popmod %>% 
  group_by(species) %>% 
  do(fit = lm(Index ~ Year, data = .)) %>% 
  broom::tidy(fit) %>% 
  filter(term == "Year") %>% 
  arrange(estimate)

perctrend <- function(data){
  data <- arrange(data, Year)
  pred <-  predict(lm(Index ~ Year, data = data))
  last <- length(pred)
  out <- (exp(pred[last]) - exp(pred[1]))/exp(pred[1])
}

poptrendperc <- popmod %>% 
  group_by(species) %>% 
  do(trend = perctrend(.)) %>% 
  unnest()


a <- ggplot(popmod, aes(x = Year, y = Index)) +
  geom_point()+
  geom_smooth(method = "lm") +
  facet_wrap(~species, scales = "free_y")
a


traits <- read.csv("data/speciesphenology.csv", header = TRUE)

poptrend <- poptrend %>% 
  left_join(traits, by = c("species" = "CommonName")) %>% 
  dplyr::select(species:p.value, CombinedLatin:HostCategory, ResStatus, WinterStage) %>% 
  left_join(poptrendperc)
write.csv(poptrend, file = "signifpoptrends.csv", row.names = FALSE)

summary(lm(estimate~BroodsGAMmax + HostCategory + ResStatus + WinterStage, data = poptrend))

# all individuals?

pops <- readRDS("allpops.rds")
pops$YearFact <- as.factor(as.character(pops$Year))
pops$SiteID <- as.factor(as.character(pops$SiteID))
pops$zyear <- scale(pops$Year)[,1]


popmod <- pops %>%  
  filter(Year != 1995) %>%
  filter(species != "Cabbage White") %>%  # with or without Cabbage White
  # group_by(species) %>% 
  # mutate(uniqyr = length(unique(Year))) %>%
  # filter(uniqyr >= 10) %>%
  group_by(SiteID) %>% 
  mutate(yrpersite = length(unique(Year))) %>% 
  group_by(Year) %>% 
  mutate(siteperyr = length(unique(SiteID))) %>% 
  # filter(yrpersite > 1) %>% 
         # siteperyr > 20) %>%
  group_by(SiteID, YearFact, Year) %>%
  summarise(PopIndex = sum(PopIndex)) %>% 
  ungroup() %>% 
  do(., CollInd(.)) %>% 
  mutate(Year = as.numeric(as.character(Year)))
        #Pieris = "yes")
# popmod <- rbind(popmod, popmod2)

a <- ggplot(popmod, aes(x = Year, y = exp(Index))) +
  geom_point()+
  geom_smooth(method = "lm") +
  # scale_color_viridis(discrete = TRUE, begin = .2, end = .8) +
  ggtitle("Trend of total # of butterflies counted") +
  labs(y = "Collated Index (butterfly-days)")
a


poptrendperc <- popmod %>% 
  do(trend = perctrend(.)) %>% 
  unnest()

# is there an effect of years since Site initiation?
popmod <- pops %>%  
  # filter(Year != 1995,
  #        species != "Cabbage White") %>%  # with or without Cabbage White
  # group_by(species) %>% 
  # mutate(uniqyr = length(unique(Year))) %>%
  # filter(uniqyr >= 10) %>%
  group_by(SiteID) %>% 
  mutate(yrpersite = length(unique(Year)),
         yrsincestart = Year - min(Year)) %>% 
  group_by(Year) %>% 
  mutate(siteperyr = length(unique(SiteID))) %>% 
  # filter(yrpersite > 1,
  #        siteperyr > 20) %>%
  group_by(SiteID, YearFact, Year, yrsincestart) %>%
  summarise(PopIndex = sum(PopIndex))

# inconclusive, need to use random effects to get more df
mod <- glm(round(PopIndex) ~ SiteID + Year + yrsincestart + I(yrsincestart^2), 
           data = popmod, family = poisson(link = "log"))
summary(mod)

library(lme4)

pops <- readRDS("allpops.rds")
pops$YearFact <- as.factor(as.character(pops$Year))
pops$SiteID <- as.factor(as.character(pops$SiteID))
pops$zyear <- scale(pops$Year)[,1]

popmod <- pops %>%  
  filter(Year != 1995) %>%
  group_by(species, SiteID) %>% 
  mutate(yrpersite = length(unique(Year)),
         yrsincestart = Year - min(Year)) %>% 
  group_by(species, Year) %>% 
  mutate(siteperyr = length(unique(SiteID))) %>% 
  filter(yrpersite > 1,
         siteperyr > 1) %>%
  droplevels() %>% 
  group_by(species) %>% 
  mutate(uniqyr = length(unique(Year))) %>%
  filter(uniqyr >= 5) %>% 
  # ungroup() %>% 
  # mutate(zyrsince = scale(yrsincestart)) %>% 
  filter(species == "Cabbage White")

gddsite <- gdd %>% 
  filter(Year >= 1995) %>%
  group_by(SiteID, Year, lat, lon) %>% 
  summarise(accumgdd = max(AccumDD)) %>% 
  group_by(SiteID) %>% 
  mutate(meangdd = mean(accumgdd),
         vargdd = as.numeric(scale(accumgdd - meangdd))) %>% 
  ungroup() %>% 
  mutate(zlat = as.numeric(scale(lat)),
         zmeangdd = as.numeric(scale(meangdd)))
  

popdat <- left_join(popmod, gddsite) %>% 
  mutate(rowid = as.factor(row_number()))

moddat <- popdat %>% 
  filter(species == "Least Skipper")
mod <- glmer(round(PopIndex) ~ zyear 
             + zmeangdd + zyear:zmeangdd
             + vargdd + zmeangdd:vargdd
             #+ zyrsince 
               + (1 + zyear + zmeangdd + zyear:zmeangdd
                  + vargdd + zmeangdd:vargdd|species)
             + (1 + zyear|SiteID) + (1|rowid) + (1|YearFact),
             data = moddat, family = poisson(link = "log"))
summary(mod)
saveRDS(mod, "allspeciesclimtrends.rds")

test <- data.frame(row.names(coef(mod)$SiteID), coef(mod)$SiteID$zyear)
names(test) <- c("site", "trend")
test$site <- as.numeric(as.character(test$site))



mod <- glmer(round(PopIndex) ~ zyear 
             #+ zyrsince 
             # + (1 + zyear + zyrsince|species) 
             + (1 + zyear|SiteID) + (1|rowid),
             data = popdat, family = poisson(link = "log"))
summary(mod)

# pesticides
pests <- read.csv("../NCEAS-RENCI_2014/Pesticides/pest_bfly_buff_overtime_MOA.csv") %>% 
  mutate(Year = YEAR,
         YearFact = as.factor(as.character(Year)),
         SiteID = formatC(site, width = 3, format = "d", flag = "0")) %>% 
  filter(buffer == 2000) %>% 
  select(Year, SiteID, ag_km2_buffer, MOAuse, meanbuff_kgkm2ag) %>% 
  group_by(Year, SiteID, ag_km2_buffer) %>% 
  spread(MOAuse, meanbuff_kgkm2ag, fill = 0) %>% 
  select(Year:ag_km2_buffer, her_G, ins_4)



moddat2 <- left_join(moddat, pests) %>% 
  filter(complete.cases(.)) %>% 
  mutate(rowid = as.factor(row_number()))



mod <- glmer(round(PopIndex) ~ zyear 
             + zmeangdd + zyear:zmeangdd
             + vargdd + zmeangdd:vargdd
             + log(her_G)
             # + her_G:zyear
             # + ins_4
             # + ins_4:zyear
             + log(ag_km2_buffer)
             + log(ag_km2_buffer):zyear
             + log(ag_km2_buffer):log(her_G)
             # + ag_km2_buffer:ins_4
             #+ zyrsince 
             # + (1 + zyear + zyrsince|species) 
             + (1 + zyear|SiteID) + (1|rowid) + (1|YearFact),
             data = moddat2, family = poisson(link = "log"))
summary(mod)

# moddat2$SiteID <- as.factor(moddat2$SiteID)
# moddat2$rowid <- as.factor(moddat2$rowid)
# 
# modg <- gam(round(PopIndex) ~
#               s(zyear) 
#             + s(zmeangdd) + ti(zyear, zmeangdd)
#             + s(vargdd) + ti(zmeangdd, vargdd)
#             # + her_G
#             # + her_G:zyear
#             # + s(ins_4)
#             # + ins_4:zyear
#             # + s(ag_km2_buffer)
#             # + ag_km2_buffer:zyear
#             # + ag_km2_buffer:her_G
#             # + ti(ag_km2_buffer, ins_4)
#             + s(SiteID, bs = "re")
#             + s(YearFact, bs = "re"),
#             data = moddat2, family = poisson(link = "log"))


# pesticides
pests <- read.csv("../NCEAS-RENCI_2014/Pesticides/pest_bfly_buff_overtime_MOA.csv") %>% 
  mutate(Year = YEAR,
         YearFact = as.factor(as.character(Year)),
         SiteID = formatC(site, width = 3, format = "d", flag = "0"),
         PopIndex = kg_buff) %>% 
  filter(buffer == 2000,
         PopIndex > 0) %>% 
  group_by(MOAuse, SiteID) %>% 
  mutate(yrpersite = length(unique(Year))) %>% 
  group_by(MOAuse, Year) %>% 
  mutate(siteperyr = length(unique(SiteID))) %>% 
  filter(yrpersite > 5,
         siteperyr > 5) %>%
  group_by(MOAuse) %>% 
  do(CollInd(.)) %>% 
  mutate(numzero = length(which(Index < 0)) / length(Index)) %>% 
  filter(numzero < .6) %>% 
  droplevels()

plt <- ggplot(pests, aes(x = Year, y = Index, group = MOAuse)) +
  geom_point() + 
  geom_smooth() +
  facet_wrap(~MOAuse, scales = "free")
plt




pests <- read.csv("../NCEAS-RENCI_2014/Pesticides/pest_bfly_buff_overtime_MOA.csv") %>% 
  mutate(Year = YEAR,
         YearFact = as.factor(as.character(Year)),
         SiteID = formatC(site, width = 3, format = "d", flag = "0")) %>% 
  filter(buffer == 2000) %>% 
  select(Year, SiteID, ag_km2_buffer, MOAuse, meanbuff_kgkm2ag) %>% 
  group_by(Year, SiteID, ag_km2_buffer) %>% 
  spread(MOAuse, meanbuff_kgkm2ag, fill = 0) %>% 
  select(Year:ag_km2_buffer, her_G, ins_4)



popdat <- left_join(popdat, pests) %>% 
  filter(complete.cases(.)) %>% 
  mutate(rowid = as.factor(row_number()))
  


mod <- glmer(round(PopIndex) ~ zyear 
             + zmeangdd + zyear:zmeangdd
             + vargdd + zmeangdd:vargdd
             # + her_G 
             # + her_G:zyear
             + ins_4
             + ins_4:zyear
             # + ag_km2_buffer
             # + ag_km2_buffer:zyear
             # + ag_km2_buffer:her_G
             # + ag_km2_buffer:ins_4:zyear
             #+ zyrsince 
             # + (1 + zyear + zyrsince|species) 
             + (1 + zyear|SiteID) + (1|rowid) + (1|YearFact),
             data = popdat, family = poisson(link = "log"))
summary(mod)

pestdat <- popdat %>%
  # filter(complete.cases(.)) %>% 
  ungroup() %>% 
  select(fun_A:ins_U)
df <- pestdat[, sapply(pestdat, function(x) { sd(x, na.rm = TRUE) != 0} )]
df <- df[complete.cases(df), ]


pestpca <- prcomp(df, scale. = TRUE, center = TRUE)


