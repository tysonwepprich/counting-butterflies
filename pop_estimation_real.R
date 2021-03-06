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
library(stringr)

# for parallel simulations with control over seed for reproducibility
# need different packages for windows computers
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
sitemod <- densityMclust(scale(siteGDD[,c(2:3)]), G = 1:4, modelNames = "EEV")
siteGDD$region <- as.character(sitemod$classification)
# visualize regional clusters
a <- ggplot(data = siteGDD, aes(x = lon, y = lat, group = region, color = region)) + geom_point()
a
siteGDD$region <- plyr::mapvalues(siteGDD$region, from = c("1", "2", "3", "4"), 
                                  to = c("NE", "NW", "CN", "SW"))
gdd <- gdd %>% 
  left_join(siteGDD[, c("SiteID", "region")])

# # add spring/fall frost to gdd as covariate
# gdd <- gdd %>% 
#   group_by(SiteID, Year) %>% 
#   arrange(DOY) %>% 
#   mutate(frost5day = rollmean(x = minT, 5, align = "right", fill = NA)) %>% 
#   mutate(fallfrost = ifelse(DOY > 200, frost5day, NA),
#          springfrost = ifelse(DOY < 200, frost5day, NA),
#          photo = geosphere::daylength(lat, DOY))

# # what is range of gdd experienced within a week?
# # can be over 100 within same site and year!
# test <- gdd %>% 
#   mutate(week = week(SiteDate)) %>% 
#   group_by(week) %>% 
#   mutate(maxgddall = max(AccumDD),
#          mingddall = min(AccumDD),
#          rangeall = maxgddall - mingddall) %>% 
#   group_by(Year, week) %>% 
#   mutate(maxgddyear = max(AccumDD),
#          mingddyear = min(AccumDD),
#          rangeyear = maxgddyear - mingddyear) %>% 
#   group_by(Year, week, SiteID) %>% 
#   mutate(maxgddsite = max(AccumDD),
#          mingddsite = min(AccumDD),
#          rangesite = maxgddsite - mingddsite) %>% 
#   dplyr::select(Year, week, SiteID, rangeall, rangeyear, rangesite) %>% 
#   distinct()



# test <- gdd %>% filter(Year > 2010) 
# plt <- ggplot(test, aes(x = SiteDate, y = springfrost)) +
#   geom_point() +
#   facet_wrap(~region, ncol = 2)
# plt

# CHOOSE 10 species for methods paper
traits <- read.csv("data/speciesphenology.csv", header = TRUE) %>% 
  # filter(UseMV == "y") %>% 
  select(CommonName, BroodsGAMmin, BroodsGAMmax, UseMV, SyncedBroods, UseMismatch, Model) %>% 
  filter(CommonName %in% c("American Copper", "Zabulon Skipper", "Leonard's Skipper", "Cabbage White",
                           "Pearl Crescent", "Black Swallowtail", "Least Skipper", "Meadow Fritillary",
                           "Viceroy", "Red-spotted Purple", "Great Spangled Fritillary", 
                           "Hackberry Emperor", "Horace's Duskywing", "Peck's Skipper",
                           "Northern Pearly-Eye"))






models <- c("doy", "gdd")
cutoff <- "loose" #c("adapt","strict", "loose")
years <- c("train", "test") #"all"  #c("all", "train", "test")
params <- expand.grid(traits$CommonName, models, cutoff, years,
                      stringsAsFactors = FALSE)
names(params) <- c("species", "model", "cutoff", "years")
# params <- params[c(5),]


ncores <- 30
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
                      years <- params$years[sim]
                      # print(species)
                      # print(Sys.time())
                      
                      pars <- data.frame(species, model, cutoff, years)
                      
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

                      # years to use in GAM
                      if(years == "train"){
                        datGAM <- datGAM %>%
                          filter(Year < 2010) %>%
                          droplevels()
                      }
                      if(years == "test"){
                        datGAM <- datGAM %>%
                          filter(Year >= 2010) %>%
                          droplevels()
                      }

                      
                      if(nrow(datGAM) == 0){
                        mod <- NA
                      }else{
                        # dat <- datGAM
                        
                        # leading zeros
                        # outlist <- list()
                        # for (y in sort(unique(dat$Year))){
                        #   temp <- dat %>% filter(Year == y) %>%
                        #     # mutate(DOY = yday(SiteDate)) %>% 
                        #     group_by(SiteID, Year) %>%
                        #     mutate(YearTotal = sum(Total),
                        #            SurvSeen = length(which(Total > 0))) %>%
                        #     filter(YearTotal > 0,
                        #            SurvSeen > 0) %>%
                        #     dplyr::select(SiteID, Total, zlistlength, lat, lon, region, DOY, 
                        #                   Year, ztemperature, zduration, AccumDD) %>% 
                        #     data.frame()
                        #   
                        #   # pad zeros at beginning and end for GAM fitting
                        #   zeros <- c(60, 70, 320, 330)
                        #   tempdf <- expand.grid(unique(temp$SiteID), zeros)
                        #   names(tempdf) <- c("SiteID", "DOY")
                        #   tempdf$SiteID <- as.character(tempdf$SiteID)
                        #   tempdf$Year <- y
                        #   tempdf$zlistlength <- 0
                        #   tempdf$ztemperature <- 0
                        #   tempdf$zduration <- 0
                        #   tempdf$Total <- 0
                        #   tempdf <- left_join(tempdf, gdd[,c("SiteID", "Year", "DOY", "AccumDD")], 
                        #                       by = c("SiteID", "Year", "DOY"))
                        #   names(tempdf)[2] <- "DOY"
                        #   tempdf$Year <- NULL
                        #   # tempdf$cumdegday <- c(rep(0, nrow(tempdf)/2), rep(max(temp$cumdegday), nrow(tempdf)/2))
                        #   
                        #   outdf <- temp %>% 
                        #     dplyr::select(SiteID, region, lat, lon, Year) %>% 
                        #     distinct() %>% 
                        #     right_join(tempdf) %>% 
                        #     bind_rows(temp)
                        # 
                        #   outlist[[y]] <- outdf
                        # }
                        # 
                        # dat <- bind_rows(outlist)
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
                                                # s(zlistlength) +
                                                # s(ztemperature) +
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
                                                # s(zlistlength) +
                                                # s(ztemperature) +
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
                            
                            # if(model == "both"){
                            #   mod <- safe_gam(Total ~ 
                            #                     s(zlistlength) +
                            #                     s(ztemperature) +
                            #                     te(lat, lon, AccumDD, bs = c("tp", "cr"), k = c(5, 30), d = c(2, 1)) +
                            #                     s(SiteYear, bs = "re") +
                            #                     s(Year, DOY, bs = "fs", k = 5, m = 1) +
                            #                     s(RegYear, AccumDD, bs = "fs", k = 5, m = 1),
                            #                   family = nb(theta = NULL, link = "log"),
                            #                   # family = poisson(link = "log"),
                            #                   data = temp,
                            #                   method = "REML", 
                            #                   optimizer = c("outer", "newton"), 
                            #                   # gamma = 1.4, 
                            #                   control = list(maxit = 500))
                            # }
                            
                          })
                        }
                      }
                      
                        
                        if(is.null(mod$error)){
                          pars$modtime <- as.numeric(modtime)[1]
                          pars$AIC <- AIC(mod$result)
                          summod <- summary(mod$result)
                          pars$N <- summod$n
                          pars$dev.expl <- summod$dev.expl
                          pars$negbin <- mod$result$family$getTheta(TRUE)
                          outlist <- list()
                          outlist[["params"]] <- pars
                          outlist[["gammod"]] <- mod$result
                          outlist[["datGAM"]] <- temp
                          saveRDS(outlist, paste(species, model, years, "noDP", "rds", sep = "."))
                          
                        }else{
                          outlist <- list()
                          outlist[["params"]] <- pars
                          outlist[["gammod"]] <- mod$error
                          outlist[["datGAM"]] <- temp
                          saveRDS(outlist, paste("gamerr", species, model, cutoff, "rds", sep = "."))
                        }
                        return(sim)
                      }
                    
if(.Platform$OS.type == "windows"){
  stopCluster(cl)
}



# Now that GAMs fit, run through mixture models as in simulation
# Check errors in baltimore and leonard's

fs <- list.files("OHGAMS/manuscript", full.names = TRUE)
# fs <- fs[grep(pattern = "loose", x = fs, fixed = TRUE)]
# fs <- fs[c(15, 16, 111, 112)]
# mixmod parameters to run for each species
mixmod <- c("hom", "het", "skew")
mod_region <- "ALL" # c("ALL", "reg4")
modyear <- "ALL"
gam_smooth <- c("none", "interpolate", "preds_8day", "preds_4day")

modparams <- list(mixmod = mixmod, gam_smooth = gam_smooth, mod_region = mod_region, modyear = modyear)
modparams <- expand.grid(modparams)



ncores <- 30
if(.Platform$OS.type == "unix"){
  registerDoMC(cores = ncores)
}else if(.Platform$OS.type == "windows"){
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
}

mcoptions <- list(preschedule = FALSE)
# foreach loop
outfiles <- foreach(sim = 1:length(fs),
                    .combine='c',
                    .packages= c("mgcv", "dplyr", "tidyr", "purrr",
                                 "lubridate", "mclust", "mixsmsn"),
                    .export = c("data", "surveys", "siteGDD", "covdata", "params", 
                                "gdd", "fs", "modparams", "traits"),
                    .inorder = FALSE,
                    .options.multicore = mcoptions) %dopar%{
                      
                      f <- fs[sim]
                      print(f)
                      spstr <- stringr::str_split(string = f, pattern = coll("/"), 3) %>% map(3)
                      spp <- stringr::str_split(string = spstr, pattern = coll("."), 5) %>% map(1) %>% unlist()
                      detprob <- stringr::str_split(string = spstr, pattern = coll("."), 5) %>% map(4) %>% unlist()
                      
                      # if(spp %in% traits$CommonName){
                        
                        tmp <- readRDS(f)
                        pars <- tmp$params
                        mod <- tmp$gammod
                        counts <- tmp$datGAM
                        
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
                          filter(DOY %in% (seq(90, 294, 4) + sample.int(n=3, size=52, replace=TRUE))) %>% 
                          ungroup() %>% 
                          mutate(zlistlength = 0,
                                 ztemperature = 0,
                                 zduration = 0,
                                 RegYear = paste(region, Year, sep = "_" ),
                                 Week = 1 + ((yday(SiteDate) - yday(ymd(paste0(year(SiteDate), "-04-01")))) %/% 7)) %>% 
                          filter(Week > 0)
                        
                        # Simulation from GAM parameter posterior probs
                        # filter results with high predictions off the chart
                        # 100 used for estimating phenology/abundance
                        Xp <- predict.gam(object = mod, newdata = preds, type="lpmatrix") ## map coefs to fitted curves
                        beta <- coef(mod)
                        Vb   <- vcov(mod) ## posterior mean and cov of coefs
                        n <- 250  #527 # choose number of simulations
                        mrand <- MASS::mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
                        ilink <- family(mod)$linkinv
                        linklppreds <- Xp %*% t(mrand)
                        nbpreds <- apply(X = linklppreds,
                                         MARGIN = 1,
                                         FUN = function(x){
                                           x <- ilink(x)
                                           x <- rnbinom(n = length(x),
                                                        mu = x,
                                                        size = mod$family$getTheta(TRUE))
                                           # x <- quantile(x, .5)
                                           x[x > 10*max(counts$Total)] <- NA
                                           return(x)
                                         })
                        
                        # check out if super high values, remove those predictions
                        test <- apply(X = nbpreds,
                                      MARGIN = 1,
                                      FUN = function(x) max(x))
                        nbpreds <- nbpreds[which(is.finite(test)), ]
                        nbpreds <- nbpreds[sample(x = 1:nrow(nbpreds), size = 100, replace = TRUE), ]
                        
                        
                        outgen <- list()
                        outN <- list()
                        # for(predsim in 1:nrow(nbpreds)){
                          for(predsim in 1:5){
                            
                          # # old way if only want median prediction
                          # adjY <- predict(mod, newdata = preds, type = "link", se.fit = TRUE)
                          # ilink <- family(mod)$linkinv
                          # preds$adjY <-  map2_dbl(.x = adjY$fit, .y = adjY$se.fit, .f = function(x, y) {
                          #   pred <- ilink(rnorm(n = 1, mean = x, sd = y))
                          #   nbpred <- qnbinom(.5, size = mod$family$getTheta(TRUE), mu = pred)
                          #   return(nbpred)
                          # })
                          
                          preds$adjY <- nbpreds[predsim, ]
                          
                          # predictions still giving ridiculous numbers at start/end of season
                          countstart <- min(counts$DOY[which(counts$Total > 0)]) - 14
                          countend <- max(counts$DOY[which(counts$Total > 0)]) + 14
                          
                          preds$adjY[which(preds$DOY < countstart | preds$DOY > countend)] <- 0
                        
                      
                          
                        tmpgen <- list()
                        tmpN <- list()
                        for (p in 1:nrow(modparams)){
                          # for (p in 1:3){
                            
                          # NEED LOOP FOR MODPARAMS
                          param <- modparams[p, ]
                          param$gam_scale <- toupper(x = pars$model)
                          param$ngen <- traits$BroodsGAMmax[which(traits$CommonName == spp)]
                          param$seed <- sim
                          param$index <- sim
                          param$species <- spp
                          param$detprob <- detprob
                          
                          adjcounts <- Adjust_Counts_AfterGAM(param, counts, preds)
                          
                          if (param$mod_region == "ALL"){
                            adjcounts$region <- "ALL"
                          }
                          
                          if (param$modyear == "ALL"){
                            adjcounts$modyear <-  "ALL"
                          }
                          
                          if (param$gam_scale == "GDD"){
                            adjcounts$Timescale <- adjcounts$AccumDD
                          } else if (param$gam_scale == "DOY"){
                            adjcounts$Timescale <- adjcounts$DOY
                          }
                        
                          
                          results <- adjcounts %>%
                            group_by(region, modyear) %>% 
                            do(mixmods = CompareMixMods(dat = ., param = param))
                          
                          genright <- results %>% 
                            do(rightgen = RightNumGen(.$mixmods, param = param, reg = .$region, yr = .$modyear)) %>% 
                            unnest()
                          
                          siteresults <- results %>%
                            do(GenClass = AssignGeneration(mixmod = .$mixmods, 
                                                           dat = adjcounts, 
                                                           param = param,
                                                           reg = .$region,
                                                           yr = .$modyear)) %>% 
                            unnest()
                          
                          ngen <- genright %>% 
                            group_by(region, modyear) %>% 
                            mutate(aicpick = maxgen[which(within_model_aic == 0)][1],
                                   bicpick = maxgen[which(within_model_bic == 0)][1]) %>% 
                            filter(bicpick == maxgen) %>% 
                            slice(1L) %>% 
                            dplyr::select(maxgen, mixmod_flag, aic, bic, badmixmod:bicpick)
                          ngen <- ngen %>% 
                            bind_cols(param[rep(seq_len(nrow(param)), each=nrow(ngen)),]) %>% 
                            mutate(predsim = predsim)
                          
                          # get phenology for each SiteYear and generation for "best" maxgen
                          siteresults$Gen <- siteresults$gen
                          siteresults$gen <- NULL
                          
                          phen_est <- siteresults %>% 
                            select(region, SiteID, Year, Gen, Timescale, count, maxgen) %>% 
                            filter(complete.cases(.)) %>% 
                            left_join(ngen[, c("region", "bicpick")]) %>% 
                            filter(maxgen == bicpick) %>% 
                            ungroup() %>% 
                            group_by(region, SiteID, Year, Gen, Timescale) %>% 
                            summarise(Total = sum(count)) %>% 
                            group_by(region, SiteID, Year, Gen) %>% 
                            do(Summ_curve(t = .$Timescale, y = .$Total)) %>% 
                            group_by(SiteID, Year) %>%
                            mutate(gen_weight = estN / sum(estN))
                          
                          # do population index in terms of DOY like UKBMS
                          N_est <- siteresults %>% 
                            select(region, SiteID, Year, Gen, DOY, count, maxgen) %>% 
                            filter(complete.cases(.)) %>% 
                            left_join(ngen[, c("region", "bicpick")]) %>% 
                            filter(maxgen == bicpick) %>% 
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
                            mutate(Year_RE = mean(curve_mean) - meanmu) %>% 
                            # dplyr::select(-estN, -meanmu) %>%
                            ungroup() %>% 
                            mutate(Year = as.numeric(as.character(Year)))
                          
                          N_est <- N_est %>% 
                            bind_cols(param[rep(seq_len(nrow(param)), each= nrow(N_est)),]) %>% 
                            mutate(predsim = predsim)
                          
                          tmpgen[[p]] <- ngen
                          tmpN[[p]] <- N_est

                        } # close p
                        tmpgendf <- bind_rows(tmpgen)
                        tmpNdf <- bind_rows(tmpN)
                        
                        outgen[[predsim]] <- tmpgendf
                        outN[[predsim]] <- tmpNdf
                        } # close predsim
                        
                        Q05 <- function(x){
                          y <- quantile(x, probs = .05, na.rm = TRUE)
                        }
                        
                        Q50 <- function(x){
                          y <- quantile(x, probs = .5, na.rm = TRUE)
                        }
                        
                        Q95 <- function(x){
                          y <- quantile(x, probs = .95, na.rm = TRUE)
                        }
                        
                        gen <- bind_rows(outgen)
                        N <- bind_rows(outN) %>% 
                          tidyr::gather(key = "metric", value = "value", PopIndex, curve_mean:Year_RE) %>% 
                          group_by(SiteID, Year, Gen, region, gam_scale, species, mixmod, gam_smooth, metric) %>% 
                          summarise_at("value", funs(q05 = Q05, q5 = Q50, q95 = Q95))
                        
                        
                        tmp[["params"]] <- left_join(tmp$params, param)
                        tmp[["gen"]] <- gen
                        tmp[["N"]] <- N
                        saveRDS(tmp, file = f)
                        
                        
                        return(sim)
                      # }else{
                        # return(sim)
                      # }
                    }
