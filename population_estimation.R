# Population Estimation with parallel processing
# Similar to RMD, but hopefully without Rstudio bugs


#########
#setup
source("funcs_pop.R")
# packages to load
library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(mclust)
library(mixsmsn)
library(caret)
library(viridis)
theme_set(theme_bw(base_size = 14)) 


# for parallel simulations with control over seed for reproducibility
# need different packages for windows computers
library(doRNG)
library(foreach) # for parallelized loops
if(.Platform$OS.type == "unix"){
  library(doMC)    # parallel backend for foreach, only for linux/mac
}else if(.Platform$OS.type == "windows"){
  library(doSNOW)
}


#####
# load gdd and site data
sites <- read.csv("data/OHsites_reconciled_update2016.csv")
sites$SiteID <- formatC(as.numeric(sites$Name), width = 3, format = "d", flag = "0")
sites$Name <- NULL
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
# sites$region9 <- as.character(sitemod$classification)
# sitemod <- densityMclust(scale(siteGDD[,c(2:3)]), G = 4)
siteGDD$region <- as.character(sitemod$classification)
# # visualize regional clusters
a <- ggplot(data = siteGDD, aes(x = lon, y = lat, group = region, color = region)) + geom_point()
a

siteGDD$region <- plyr::mapvalues(siteGDD$region, from = c("1", "2", "3", "4"), 
                                 to = c("NE", "NW", "CN", "SW"))
gdd <- gdd %>% 
  left_join(siteGDD[, c("SiteID", "region")])

# remove some sites that are clustered close together
uniqsites <- gdd %>% ungroup() %>%  dplyr::select(SiteID, lat, lon) %>% distinct()

points_matrix <- as.matrix(dist(uniqsites[, c("lat", "lon")], diag = TRUE, upper = TRUE))
points_matrix[lower.tri(points_matrix, diag=TRUE)] <- NA
points_matrix <- points_matrix <= 0.07

v <- colSums(points_matrix, na.rm=TRUE) == 0
uniqsites <- uniqsites[v, ]

gdd <- gdd %>% 
  filter(SiteID %in% uniqsites$SiteID)

# Point out that DOY and GDD are highly correlated in spring/summer, less so in fall 
# gdd %>% group_by(month(SiteDate), region) %>% 
# summarise(corr = cor(DOY, AccumDD)) %>% data.frame()



#####
# define parameters
nsite <- 40 # up to 140 OH sites
nyear <- 5 # up to 35 years in Daymet data
nsurv <- 30
surv_missing <- c(0, .2, .4) # remove surveys from all counts
# site_missing <- .2 # turnover across years
group_struct <- c("all", "year")
gam_scale <- c("DOY", "GDD")
gam_smooth <- c("none", "interpolate", "preds_8day", "preds_4day")

ngen <- c(1:4)
gen_size <- "equal"  #c("equal", "inc", "dec")
# volt_flex <- "Y"
# gen_ddreq <- 600 # depends on ngen
peak_sd <- 25 # 75 might be too much... c(25, 75) # low and high overlap?
death_rate <- c(.005, .0075) # on degree day rate
# peak_sd * death_rate needs to be < 1 for positive beta parameters

# site total population dispersion
# negbin_mu <- 100
# negbin_disp <- 1
# pois_lam <- c(100, 500, 1000)
# detection probability parameters (logit scale)
detprob_b0 <- c(-7,-9)
detprob_b1 <- .6
detprob_b2 <- -.01
detprob_model <- c("known", "covariate", "none")

# site/year variation in mu
site_mu_sd <- 50
year_mu_sd <- 50

# mixmod choices
mixmod <- c("hom", "het", "skew")
mod_region <- c("ALL", "reg4")

# # nonlinear det prob
x <- seq(-5, 45, .1)
y <- plogis(-9 + x * .6 + -.01 * x^2)
plot(x, y)

params <- list(nsite = nsite, nyear = nyear, nsurv = nsurv, surv_missing = surv_missing, 
               group_struct = group_struct, gam_scale = gam_scale,
               gam_smooth = gam_smooth, ngen = ngen, gen_size = gen_size, peak_sd = peak_sd, 
               death_rate = death_rate, 
               detprob_b0 = detprob_b0, detprob_b1 = detprob_b1, detprob_b2 = detprob_b2, 
               detprob_model = detprob_model, site_mu_sd = site_mu_sd, year_mu_sd = year_mu_sd,
               mixmod = mixmod, mod_region = mod_region)
params <- expand.grid(params)

# some parameter combinations don't make sense, try to exclude to save time
params <- params[-which(params$ngen == 1 & params$mixmod == "het"), ]
params <- params[-which(params$surv_missing == 0 & params$gam_smooth == "interpolate"), ]

# for random numbers
params$seed <- 1:nrow(params)
params$index <- formatC(params$seed, width=5, flag="0")


# params <- params[c(1, 2001, 2501, 3001, 3501, 4001), ]
# 
# fs <- list.files(pattern = "popest_", recursive = TRUE, full.names = TRUE)
# details <- file.info(fs)
# done <- fs[which((ymd_hms(details$mtime) - ymd_hms(Sys.time())) > -100)]
# done <- c(list.files('results_0'), list.files('results_1'))
# 
# done <- as.numeric(stringr::str_split_fixed(string = done, pattern = "_", n = 2)[,2])
# 
# params <- params[-which(params$seed %in% done), ]


######
# analysis workflow


# test <- params[5, ]
# # test$site_mu_sd <- 50
# # test$year_mu_sd <- 50
# # test$gam_smooth <- "preds_8day"
# # test$detprob_model <- "covariate"
# 
# counts <- Simulate_Counts(data = test, gdd = gdd)
# 
# # checking counts
# pops <- counts %>%
#   group_by(SiteID, Year, M) %>%
#   summarise(bflydays = sum(Y))
plt <- ggplot(counts, aes(x = AccumDD, y = Y, group = as.factor(Year), color = as.factor(Year))) +
  geom_path() +
  facet_wrap(~SiteID, ncol = 4, scales = "free_y")
plt
# 
# # GAM interpolation/prediction
# adjcounts <- Adjust_Counts(data = test, counts)
plt <- ggplot(adjcounts, aes(x = DOY, y = adjY, group = Year, color = Year)) +
  geom_point() +
  facet_wrap(~SiteID, ncol = 4, scales = "free_y")
plt
# 
# system.time({
#   results <- adjcounts %>%
#     # mutate(mvmin = 1, mvmax = 3) %>%
#     group_by(SiteID, Year) %>%
#     mutate(weeks_obs = length(which(adjY > 0)),
#            total_obs = sum(round(adjY))) %>%
#     filter(weeks_obs >= 2, total_obs >= 5) %>%
#     do(mixmods = CompareMixMods(dat = ., mvmax = 3))
# })


# # test run with a few rows
# params <- params %>% 
#   filter(nsite == 4, nyear == 4, ngen == 2, gen_size == "equal",
#          peak_sd == 10, death_rate == 0.6, site_mu_sd == 25, 
#          year_mu_sd == 50, surv_missing != 0.2)
#####
# Simulation
#####
# ncores <- 6
ncores <- 40

if(.Platform$OS.type == "unix"){
  registerDoMC(cores = ncores)
}else if(.Platform$OS.type == "windows"){
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
}

mainDir <- getwd()
mcoptions <- list(preschedule = FALSE)

# foreach loop
outfiles <- foreach(sim = 1:nrow(params),
                    .combine='c',
                    .packages= c("mgcv", "dplyr", "tidyr", "purrr",
                                 "lubridate", "mclust", "mixsmsn"),
                    .export = c("params", "gdd", "Abund_Curve", 
                                "Adjust_Counts", "CompareMixMods",
                                "Simulate_Counts", "Summ_curve",
                                "Summ_mixmod", "mainDir", "RightNumGen", "AssignGeneration"),
                    .inorder = FALSE,
                    .options.multicore = mcoptions) %dopar% {
                      
                      param <- params[sim, ]
                      
                      # make directories to store results
                      subDir <- paste("results", 
                                      formatC(floor(param$seed / 1000), width = 2, flag = "0"), sep = "_")

                      setwd(mainDir)
                      if (file.exists(subDir)){
                        setwd(file.path(mainDir, subDir))
                      } else {
                        dir.create(file.path(mainDir, subDir))
                        setwd(file.path(mainDir, subDir))
                      }
                      
                      
                      
                      simpop <- Simulate_Counts(data = param, gdd = gdd)
                      counts <- simpop[[1]]
                      truth <- simpop[[2]]
                      # GAM interpolation/prediction
                      adjcounts <- Adjust_Counts(data = param, counts)
                      
                      # Do regional mixture models, pooling years/sites
                      # Also try statewide
                      if (param$mod_region == "ALL"){
                        adjcounts$region <- "ALL"
                      }
                      if (param$group_struct == "all"){
                        adjcounts$modyear <- "all"
                      }else{
                        adjcounts$modyear <- adjcounts$Year
                      }
                      
                        results <- adjcounts %>%
                          group_by(region, modyear) %>% 
                          # mutate(weeks_obs = length(which(round(adjY) > 0)),
                          # total_obs = sum(round(adjY))) %>% 
                          # filter(weeks_obs >= 2, total_obs >= 5) %>% 
                          do(mixmods = CompareMixMods(dat = ., param = param))
                        
                        # This df could be used to see if BIC/AIC choose right number of gen
                        genright <- results %>% 
                          do(rightgen = RightNumGen(.$mixmods, param = param, reg = .$region, yr = .$modyear)) %>% 
                          unnest()
                        
                        # regroup generation assignments by SiteYear
                        siteresults <- results %>%
                          do(GenClass = AssignGeneration(mixmod = .$mixmods, 
                                                         dat = adjcounts, 
                                                         param = param,
                                                         reg = .$region,
                                                         yr = .$modyear)) %>% 
                          unnest()
                 
                      
                      # # what is siteresults?
                      # # plot shows generations estimated for different maxgen
                      # plt <- ggplot(siteresults, aes(x = Timescale, group = gen, color = gen)) +
                      #   geom_density() +
                      #   facet_wrap(~maxgen)
                      # plt
                      
                      # # Extract summary stats for each generation's distribution
                      # summ_mods <- results %>% 
                      #   ungroup() %>% 
                      #   group_by(region) %>% 
                      #   do(Summ_mixmod(.))
                      
                      outlist <- list(param, counts, adjcounts, truth, genright, siteresults)   
                      saveRDS(object = outlist, file = paste("popest", param$index, sep = "_"))
                      
                      rm(counts, adjcounts, summ_mods, results, outlist, simpop, truth, genright, siteresults)
                      gc()
                      
                      # hope this is the only thing returned
                      return(param$seed)
                    } # close dopar


if(.Platform$OS.type == "windows"){
  stopCluster(cl)
}
#####
# Summarize results
#####
#new notes:
# if mixmods don't work at all, maxgen skipped in bestmods/genright which could skew best model based on AIC if some missing
# the right maxgen may not work, which threw off siteresults
# now, if this happens, in siteresults all generations/counts assigned as NA, which is cue that mixture model didn't work


#####

# Expect ngen to be wrong when last generation weight falls lower.
# Many mixture models fail to fit, how to include or account for these?

# with simulation results in hand, next:
# 1. when is number of generations correct compared to params?
# For each SiteYear: model, maxgen, AIC/BIC, fit fails, other errors, last gen weight, # counts, smoothing
# Correct maxgen, mixture models "work" as responses


# 2. what is error in generation weights estimated compared to simulated "truth"?
# 3. what is error in phenology for each site/year compared to "truth"?
# 4. what is error in site/year population size compared to "truth" using trapezoid rule?

# NB: not all mixture model problems filtered out, like sigma2 close to zero to create a peak for one count

fs <- list.files(pattern = "popest_", recursive = TRUE, full.names = TRUE)
fs_index <- stringr::str_split_fixed(string = fs, pattern = "_", n = 3)[, 3]

md <- file.info(fs) %>% 
  mutate(name = row.names(file.info(fs)))
missed <- md %>% filter(mtime < as.Date("2018-02-13"))

# errors, just ignored because only 3
params[-which(params$index %in% fs_index), ]

# outlist <- list()
# outscore <- list()
# outpop <- list()
# for (f in fs){
ncores <- 40

if(.Platform$OS.type == "unix"){
  registerDoMC(cores = ncores)
}else if(.Platform$OS.type == "windows"){
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
}

mcoptions <- list(preschedule = FALSE)
outfiles <- foreach(sim = 1:length(fs),
                    .combine='c',
                    .packages= c("mgcv", "dplyr", "tidyr", "purrr",
                                 "lubridate", "mclust", "mixsmsn"),
                    .export = c("Abund_Curve", "fs",
                                "Adjust_Counts", "CompareMixMods",
                                "Simulate_Counts", "Summ_curve",
                                "Summ_mixmod", "RightNumGen", "AssignGeneration", "TrapezoidIndex"),
                    .inorder = FALSE,
                    .options.multicore = mcoptions) %dopar% {
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
                      trueweight <- truth %>%
                        mutate(Year = as.character(Year)) %>%
                        left_join(adjcounts[, c("SiteID", "Year", "region", "modyear")]) %>%
                        filter(Gen == max(Gen)) %>%
                        group_by(region, modyear) %>%
                        summarise(true_weight = mean(gen_weight))
                      
                      
                      ngen <- genright %>%
                        left_join(trueweight) %>% 
                        group_by(region, modyear) %>%
                        mutate(aicpick = maxgen[which(within_model_aic == 0)][1],
                               bicpick = maxgen[which(within_model_bic == 0)][1],
                               right_ngen = ifelse(true_weight == 0,
                                                   ifelse(maxgen == (param$ngen - 1),
                                                          "yes", "no"),
                                                   ifelse(maxgen == param$ngen,
                                                          "yes", "no"))) %>% 
                        filter(right_ngen == "yes") %>%
                        slice(1L) %>%
                        dplyr::select(maxgen, mixmod_flag, badmixmod:bicpick)

                      ngen_score <- bind_cols(param[rep(seq_len(nrow(param)), each=nrow(ngen)),], ngen)

                      tmp[[7]] <- ngen_score

                      # 2. what is error in generation weights estimated compared to simulated "truth"?
                      # 3. what is error in phenology for each site/year compared to "truth" (quantiles/max/mean)?
                      right_gen <- ngen %>% 
                        mutate(rightgen = maxgen) %>% 
                        dplyr::select(region, modyear, rightgen)
                      
                      # Do curve summaries in terms of GDD for comparison
                      siteresults$Gen <- siteresults$gen
                      siteresults$gen <- NULL
                      phen_est <- siteresults %>% 
                        left_join(right_gen) %>% 
                        # dplyr::select(-c(curve_mean, curve_max, curve_q0.1, curve_q0.5, curve_q0.9, n)) %>% 
                        filter(complete.cases(.)) %>% 
                        filter(maxgen == rightgen) %>% 
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
                        dplyr::select(-estN, -M)
                      # 3a. regression of random effects for site/year compared to estimated variation?
                      # maybe not necessary/interesting...
                      
                      
                      # 4. what is error in site/year population size compared to "truth" using trapezoid rule?
                      # do population index in terms of DOY like UKBMS
                      N_est <- siteresults %>% 
                        # dplyr::select(-c(curve_mean, curve_max, curve_q0.1, curve_q0.5, curve_q0.9, n)) %>% 
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
                        dplyr::select(-Site_RE, -Year_RE, -Site_M, -Year_M) %>% 
                        tidyr::gather(key = "metric", value = "value", gen_weight:PopIndex) %>% 
                        group_by(SiteID, Year, Gen, metric) %>% 
                        summarise(error = value[1] - value[2],
                                  value = value[2])
                      
                      # NAs if a generation not counted at a particulr site due to low numbers
                      sumscore <- allscore %>% 
                        ungroup() %>% 
                        dplyr::filter(complete.cases(.)) %>% 
                        group_by(metric) %>% 
                        summarise(rmse = sqrt(mean(error^2)), 
                                  mae = mean(abs(error)),
                                  pv = PropVariation(value),
                                  cv = sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE),
                                  avg = mean(value, na.rm = TRUE)) %>% 
                        ungroup() %>% 
                        mutate(Gen = "ALL")
                      
                      genscore <- allscore %>% 
                        ungroup() %>% 
                        filter(complete.cases(.)) %>% 
                        group_by(Gen, metric) %>% 
                        summarise(rmse = sqrt(mean(error^2)), 
                                  mae = mean(abs(error)),
                                  pv = PropVariation(value),
                                  cv = sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE),
                                  avg = mean(value, na.rm = TRUE)) %>%                         
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
                        summarise(corrpop = cor(PopIndex[Data == "truth"], PopIndex[Data == "estimate"], 
                                                method = "kendall"))
                      popscore <- bind_cols(param[rep(seq_len(nrow(param)), each=nrow(relpop)),], relpop)
                      tmp[[9]] <- popscore
                      
                      saveRDS(tmp, file = f)
                      
                    }
# Didnt work with purrr
# system.time({
# genscore <- fs %>% 
#   map_df(~ readRDS(.)[[7]])
# siteyrscore <- fs %>% 
#   map_df(~ readRDS(.)[[8]])
# relpopscore <- fs %>% 
#   map_df(~ readRDS(.)[[9]])
# })
 



# 
# ### Another issue
# # with GAM, predictions at ends of year for GDD near upper limit go really high
# # fixing this by removing lpmatrix prediction uncertainty
# system.time({
#  pops <- vector("list", length(fs))
#   for (i in 1:length(fs)){
#     f <- fs[1970]
#     tmp <- readRDS(f)
#    counts <- tmp[[2]]
#    adjcounts <- tmp[[3]]
#     maxy <- max(counts$Y)
#     maxadjy <- max(adjcounts$adjY)
#     ratio <- maxadjy / maxy
#     pops[[i]] <- ratio
#   }
# })
# ratdf <- as.numeric(flatten(pops))
# 
# counts$Year <- as.character(counts$Year)
# adjcounts$Year <- as.character(adjcounts$Year)
# ggplot(data = counts, aes(x = AccumDD, y = Y, group = SiteID)) +
#          geom_point(color = "red", alpha = .5) +
#          geom_point(data = adjcounts, aes(x = AccumDD, y = adjY, group = SiteID), alpha = .5) +
#          facet_wrap(~Year)
# 



# scoring results
gendf <- readRDS("gendf.rds")
phendf <- readRDS("phendf.rds")
popdf <- readRDS("popdf.rds")


# 1. Model pathologies
# Mixmod errors, GAM errors, SD of mixmod really small
test <- gendf %>% 
  group_by(index) %>% 
  summarise(badmixmod = max(badmixmod),
            redundant = min(redundant, na.rm = TRUE),
            nearzerosigma = min(nearzerosigma)) %>% 
  left_join(params) %>% 
  mutate(redundant = ifelse(redundant == Inf, NA, redundant))

# this shows that mixture models failed to fit when group_struct == year & mod_region == reg4 & mixmod == skew
# 21 cases of this error
test %>% filter(badmixmod == 1) %>% data.frame()

# from confusion matrix, we know overfitting of voltinism
# just state that it's due to redundant modes / near zero sigma for modes
test %>% filter(gam_scale == "DOY", redundant <= 10) %>% data.frame()

# redundant modes correlated with sigma being very low
ggplot(test %>% filter(gam_scale == "DOY"), aes(x = sqrt(nearzerosigma), y = redundant, group = gam_scale)) +
  geom_point(alpha = .1) +
  facet_grid(gam_scale~ngen, scale = "free")

# could compare previous graph with what's expected from the true counts

# do reg4 or year groupings help at all?
df <- gendf %>% 
  filter(group_struct == "year") %>%
  filter(mod_region == "reg4") %>%
  mutate(mismatch = bicpick - maxgen)

# versus
df <- gendf %>% 
  filter(group_struct == "year") %>%
  filter(mod_region == "reg4") %>%
  group_by(index) %>% 
  mutate(bicpick = max(bicpick),
         maxgen = max(maxgen),
         mismatch = bicpick - maxgen)

df <- gendf %>% 
  filter(group_struct != "year") %>%
  filter(mod_region != "reg4")
# 
# mod <- lm(mismatch ~ (true_weight + ngen + region)^2, data = df)
# summary(mod)

# confusion matrix for generations classified by mixmod
# AIC chooses too many gen, BIC better but not great especially for ngen == 1
true <- factor(df$maxgen, levels = c("1", "2", "3", "4", "5"))
pred <- factor(df$bicpick, levels = c("1", "2", "3", "4", "5"))
cfmat <- caret::confusionMatrix(pred, true)
cfmat



# 2. 


df <- gendf %>% 
  filter(group_struct == "all", mod_region == "ALL")

# NOTE: group_struct with year makes 5 rows for each simulation replicate for gendf
# how to account for this and compare?
# Also same issue for reg4 grouping

ggplot(df, aes(x = sqrt(nearzerosigma), group = gam_scale)) +
         geom_density() +
         facet_grid(ngen~gam_scale, scale = "free")



# first gauge gen classification as correctly selecting the max # gen in each simulation
# really bad accuracy over all replicates!
df <- gendf %>% 
  group_by(index) %>% 
  summarise(ngen = ngen[1],
            aicpick = max(aicpick),
            bicpick = max(bicpick)) %>% 
  left_join(params)

df <- df %>% 
  filter(group_struct == "all")
  filter(group_struct == "all", mod_region == "ALL", mixmod != "het", gam_scale == "GDD", gam_smooth == "preds_4day")


# confusion matrix for generations classified by mixmod
# AIC chooses too many gen, BIC better but not great especially for ngen == 1
true <- factor(df$ngen, levels = c("1", "2", "3", "4", "5"))
pred <- factor(df$bicpick, levels = c("1", "2", "3", "4", "5"))
cfmat <- caret::confusionMatrix(pred, true)
cfmat

# scores of parameter sets, one at a time
# could use confusion matrix accuracy CI for significant differences
outlist <- list()
groupset <- expand.grid(group_struct, mod_region)
names(groupset) <- c("group_struct", "mod_region")
for (j in 1:nrow(groupset)){
  df <- gendf %>% 
    filter(group_struct == groupset[j, 1],
           mod_region == groupset[j, 2])
  for(i in names(params)[c(4, 6:8, 11, 12, 15, 18)]){
    p <- pull(.data = df, var = i)
    # mixmod accuracy
    tmpdf <- df %>% 
      mutate(true = factor(ngen, levels = c("1", "2", "3", "4", "5")),
             pred = factor(bicpick, levels = c("1", "2", "3", "4", "5")),
             varsplit = i,
             varfact = as.factor(as.character(p))) %>% 
      group_by(varsplit, varfact) %>% 
      do(data.frame(t(caret::confusionMatrix(.$pred, .$true)$overall))) %>% 
      mutate(        group_struct = groupset[j, 1],
                     mod_region = groupset[j, 2])
    
    outlist[[length(outlist)+1]] <- tmpdf
    
    # df1 <- df %>% 
    #   filter(ngen != 1) 
    # 
    # p <- pull(df1, "mixmod")
    # tmpdf1 <- df1 %>% 
    #   mutate(true = factor(ngen, levels = c("1", "2", "3", "4", "5")),
    #          pred = factor(bicpick, levels = c("1", "2", "3", "4", "5")),
    #          varsplit = "mixmod",
    #          varfact = as.factor(as.character(p))) %>% 
    #   group_by(varsplit, varfact) %>% 
    #   do(data.frame(t(caret::confusionMatrix(.$pred, .$true)$overall))) %>% 
    #   mutate(        group_struct = groupset[j, 1],
    #                  mod_region = groupset[j, 2])
    # outlist[[length(outlist)+1]] <- tmpdf1
  }
}

generr <- bind_rows(outlist)
generr$varsplit <- factor(generr$varsplit, 
                          levels = c("ngen", "death_rate", "detprob_b0", "surv_missing", "detprob_model",
                                     "gam_scale", "gam_smooth", "mixmod"))
generr$varsplit <- plyr::mapvalues(generr$varsplit, 
                                   from = c("ngen", "death_rate", "detprob_b0", "surv_missing", "detprob_model",
                                            "gam_scale", "gam_smooth", "mixmod"),
                                   to = c("# Generations", "Death rate", "Detection Probability", "Missing surveys", "Detection Modeling",
                                          "GAM timescale", "Use of GAM predictions", "Mixture model distribution"))
generr$varfact <- factor(generr$varfact,
                         levels = c("none", "0", "0.2", "0.4", "DOY", "GDD", "interpolate", "preds_8day", "preds_4day",
                                    "1", "2", "3", "4", "0.005", "0.0075", "-9", "-7", "covariate", "known", "hom", "het", "skew"))
generr$varfact <- plyr::mapvalues(generr$varfact, 
                                   from = c("none", "0", "0.2", "0.4", "DOY", "GDD", "interpolate", "preds_8day", "preds_4day",
                                            "0.005", "0.0075", "-9", "-7", "covariate", "known", "hom", "het", "skew"),
                                   to = c("None", "0", "20%", "40%", "Day of year", "Degree-day", "Missing", 
                                          "Replace\n8-day", "Replace\n4-day", "Low", "High", "Low", "High", 
                                          "Covariate", "Known", "Equal\nvariance", "Unequal\nvariance", "Skew\nnormal"))

genplt <- ggplot(generr, aes(x = varfact, y = Accuracy, group = interaction(group_struct, mod_region),
                             color = interaction(group_struct, mod_region))) +
  geom_point(position = position_dodge(width=0.4)) +
  geom_linerange(aes(ymin=AccuracyLower, ymax=AccuracyUpper), position=position_dodge(width=0.4)) +
  facet_wrap(~varsplit, scales = "free_x") +
  expand_limits(x = 0, y = 0) +
  scale_color_viridis(name="Grouping of data\nfor mixture models",
                      breaks=c("all.ALL", "year.ALL", "all.reg4", "year.reg4"),
                      labels=c("All together", "By year", "By region", "By year x region"),
                      discrete = TRUE, begin = 0, end = .8) +
  xlab("Variables in simulation study") +
  ylab("Voltinism prediction accuracy\nfrom confusion matrix") +
  theme(legend.position = c(0.85, 0.15))
  
genplt 

ggsave(filename = paste("ConfMatrixAccuracy", "png", sep = "."), 
       plot = genplt, device = "png", path = "plots", width = 10, height = 8, units = "in")





# scores of parameter sets, one at a time
outlist1 <- list()
df <- phendf
# can choose to split by ngen or not, keep it simple for results
# for (j in 1:4){
  # df <- phendf %>% 
    # filter(ngen == j)
  for(i in names(params)[c(4:8, 11, 12, 15, 19)]){
    
    p <- pull(df, i)
    tmpdf <- df %>% 
      mutate( varsplit = i,
              varfact = as.factor(as.character(p)))  %>% 
      group_by(varsplit, varfact, metric, Gen) %>% # add ngen here if wanted
      summarise(RMSE = mean(pv),
                sem = sd(pv)/sqrt(length(pv)))
    outlist1[[length(outlist1)+1]] <- tmpdf
  }
  
df <- phendf %>% 
  filter(ngen != 1) 
  
  p <- pull(df, "mixmod")
  tmpdf1 <- df %>% 
    mutate( varsplit = "mixmod",
            varfact = as.factor(as.character(p)))  %>% 
    group_by(varsplit, varfact, metric, Gen) %>% 
    summarise(RMSE = mean(pv),
              sem = sd(pv)/sqrt(length(pv)))
  outlist1[[length(outlist1)+1]] <- tmpdf1
  
# }


# scores of parameter sets, one at a time
outlist1 <- list()
df <- popdf
# for (j in 1:4){
  # df <- popdf %>% 
    # filter(ngen == j)
  for(i in names(params)[c(4:8, 11, 12, 15, 18, 19)]){
    
    p <- pull(df, i)
    tmpdf <- df %>% 
      mutate( varsplit = i,
              varfact = as.factor(as.character(p)))  %>% 
      group_by(varsplit, varfact, Gen) %>% 
      summarise(RMSE = mean(corrpop),
                sem = sd(corrpop)/sqrt(length(corrpop))) %>% 
      data.frame()
    outlist1[[length(outlist1)+1]] <- tmpdf
  }
# }



phenerr <- bind_rows(outlist1) %>% 
  rowwise() %>% 
  mutate(RMSElower = RMSE - 2 * sem,
         RMSEupper = RMSE + 2 * sem) %>% 
  filter(metric == "curve_mean", Gen != "ALL")
  # filter(Gen == "ALL", 
         # metric %in% c("curve_q0.1", "curve_q0.5", "curve_q0.9", "curve_max", "curve_mean"))
         filter(metric %in% c("PopIndex"))


phenerr$varsplit <- factor(phenerr$varsplit, 
                          levels = c("ngen", "death_rate", "detprob_b0", "surv_missing", "detprob_model",
                                     "gam_scale", "gam_smooth", "mixmod", "group_struct", "mod_region"))
phenerr$varsplit <- plyr::mapvalues(phenerr$varsplit, 
                                   from = c("ngen", "death_rate", "detprob_b0", "surv_missing", "detprob_model",
                                            "gam_scale", "gam_smooth", "mixmod", "group_struct", "mod_region"),
                                   to = c("Number of generations", "Death rate", "Detection Probability", "Missing surveys", "Detection Modeling",
                                          "GAM timescale", "Use of GAM predictions", "Mixture model distribution",
                                          "Mixture model time", "Mixture model space"))
phenerr$varfact <- factor(phenerr$varfact,
                         levels = c("none", "0", "0.2", "0.4", "DOY", "GDD", "interpolate", "preds_8day", "preds_4day",
                                    "1", "2", "3", "4", "0.005", "0.0075", "-9", "-7", "covariate", "known", "hom", "het", "skew",
                                    "all", "ALL", "year", "reg4"))
phenerr$varfact <- plyr::mapvalues(phenerr$varfact, 
                                  from = c("none", "0", "0.2", "0.4", "DOY", "GDD", "interpolate", "preds_8day", "preds_4day",
                                           "0.005", "0.0075", "-9", "-7", "covariate", "known", "hom", "het", "skew", "all", "ALL", "year", "reg4"),
                                  to = c("None", "0", "20%", "40%", "Day of year", "Degree-day", "Missing", 
                                         "Replace\n8-day", "Replace\n4-day", "Low", "High", "Low", "High", 
                                         "Covariate", "Known", "Equal\nvariance", "Unequal\nvariance", "Skew\nnormal",
                                         "All years", "Statewide", "By year", "By region"))

genplt <- ggplot(phenerr, aes(x = varfact, y = RMSE, group = as.factor(Gen), color = as.factor(Gen))) +
# genplt <- ggplot(phenerr, aes(x = varfact, y = RMSE, group = metric, color = metric)) +
  geom_point( position = position_dodge(width=0.4)) +
  # geom_point()+
# geom_linerange(aes(ymin=RMSElower, ymax=RMSEupper)) +
  geom_linerange(aes(ymin=RMSElower, ymax=RMSEupper), position=position_dodge(width=0.4)) +
  facet_wrap(~varsplit, scales = "free_x") +
  # expand_limits(x = 0, y = 0) +
  scale_color_viridis(
    name="By generation",
                      # breaks=c("curve_max", "curve_mean", "curve_q0.1", "curve_q0.5", "curve_q0.9"),
                      # labels=c("Peak", "Weighted mean", "10th percentile", "50th percentile", "90th percentile"),
                      discrete = TRUE, begin = 0, end = .8) +
  xlab("Variables in simulation study") +
  ylab("Precision of estimated phenology\nProportional Variability (PV) of mean phenology") +
  theme(legend.position = c(0.75, 0.15))

genplt 

ggsave(filename = paste("PVPhenology", "png", sep = "."), 
       plot = genplt, device = "png", path = "plots", width = 12, height = 8, units = "in")




# for (splt in unique(phenerr$varsplit)){
#   plt <- phenerr %>% 
#   filter(varsplit == splt, Gen == "ALL") %>% 
#   ggplot(aes(x = varfact, y = RMSE)) +
#   geom_point() +
#   scale_y_continuous(limits = c(0, NA)) +
#   geom_linerange(aes(ymin=RMSElower, ymax=RMSEupper)) +
#   facet_wrap(~metric, scales = "free_y")
#   print(plt)
# }


# what about best combination of parameters?
# is accuracy good?


df <- gendf %>% 
  filter(group_struct != "year", mod_region != "reg4", surv_missing != 0.4, gam_scale == "GDD",
         gam_smooth != "none", detprob_model != "none", mixmod == "hom")

# confusion matrix for generations classified by mixmod
# AIC chooses too many gen, BIC better but not great especially for ngen == 1
true <- factor(df$maxgen, levels = c("1", "2", "3", "4", "5"))
pred <- factor(df$bicpick, levels = c("1", "2", "3", "4", "5"))
cfmat <- caret::confusionMatrix(pred, true)
cfmat


df <- phendf %>% 
  filter(group_struct != "year", mod_region != "reg4", surv_missing == 0.4,
         gam_smooth != "none", metric == "curve_mean")
df %>% group_by(gam_smooth, detprob_model, gam_scale, mixmod) %>%
  summarise(rmse = mean(rmse),
            pv = mean(pv)) %>% 
  arrange(rmse) %>% 
  data.frame()


df <- popdf %>% 
  filter(group_struct != "year", mod_region != "reg4", surv_missing != 0.4, gam_scale == "GDD",
         gam_smooth != "none", detprob_model != "none", mixmod == "hom")

ggplot(df, aes(x = corrpop, group = Gen, color = Gen)) +
  geom_density() +
  facet_wrap(~ngen)

# rmse comparison for different parameters
phendf %>% 
  # filter(gam_scale == "GDD", gam_smooth == "preds_4day", mixmod == "hom") %>% 
  group_by(metric, ngen, Gen) %>% 
  summarise(RMSE = mean(rmse)) %>% 
  ggplot(aes(x = ngen, y = RMSE, group = Gen, color = Gen)) +
  geom_point() +
  facet_wrap(~metric, scales = "free_y")

# correlation of relative pop size
df <- popdf %>% 
  filter(gam_scale == "GDD", gam_smooth == "preds_4day", mixmod == "hom")
plt <- ggplot(df, aes(x = corrpop, group = Gen, color = Gen)) +
  geom_density() +
  facet_wrap(~ngen)
plt

# correlation of populations bad for earlier Generations (not simulated to be different)
# what about if using only last generation that varies between sites with latitude?
df <- popdf %>% 
  mutate(death_rate = as.factor(as.character(death_rate)),
         detprob_b0 = as.factor(as.character(detprob_b0)))

mod <- lm(corrpop ~ surv_missing + gam_scale + gam_smooth + ngen + death_rate + detprob_b0 + detprob_model + mixmod + mod_region,
          data = df)
summary(mod)



df <- phendf %>% 
  filter(metric == "curve_q0.5") %>% 
  filter(group_struct == "all", mod_region == "ALL") %>% 
  filter(Gen != "ALL") %>% 
  mutate(death_rate = as.factor(as.character(death_rate)),
         detprob_b0 = as.factor(as.character(detprob_b0)))

mod <- lm(rmse ~ surv_missing + gam_scale + gam_smooth + ngen + death_rate + detprob_b0 + detprob_model + mixmod,
          data = df)
summary(mod)

test <- step(object = lm(pv ~ 1, data = df),
             scope = pv ~ (surv_missing + gam_scale + gam_smooth + detprob_model + mixmod)^2,
             # scope = rmse ~ (surv_missing + gam_scale + gam_smooth + ngen + death_rate + detprob_b0 + detprob_model + mixmod + mod_region + group_struct)^2,
             direction = "forward", k = log(nrow(df)))

round_df(broom::tidy(test), 3)


# switch RMSE and PV to gauge combined effects of GDD/DOY and Gen
df <- phendf %>% 
  filter(metric == "curve_q0.5") %>% 
  group_by(gam_scale, ngen, Gen, metric) %>% 
  summarise(RMSE = mean(pv),
            sem = sd(pv)/sqrt(length(pv))) %>% 
  rowwise() %>% 
  mutate(RMSElower = RMSE - 2 * sem,
         RMSEupper = RMSE + 2 * sem) %>% 
  data.frame()

ggplot(df, aes(x = ngen, y = RMSE, color = Gen)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  geom_linerange(aes(ymin=RMSElower, ymax=RMSEupper)) +
  facet_wrap(~gam_scale)

outdf <- bind_rows(outlist)
saveRDS(outdf, file = "genright.rds")

outdf <- readRDS("genright.rds")

# more than 1000 errors (no results returned) of unknown cause, DOY gam_scale in common and 3 or 4 ngen
errparam <- params[-which(params$index %in% unique(outdf$index)), ]
  
outdf <- outdf %>% 
  mutate(corrAIC = ifelse(within_model_aic > 0, 0, 1),
         corrBIC = ifelse(within_model_bic > 0, 0, 1))

mod <- glm(corrBIC ~ surv_missing + gam_scale + gam_smooth + ngen + death_rate + pois_lam + detprob_model + mixmod + region,
           family = binomial(link = "logit"), data = outdf)
mod <- glm(corrAIC ~ surv_missing + gam_scale + gam_smooth + ngen + death_rate + pois_lam + detprob_model + mixmod + region,
           family = binomial(link = "logit"), data = outdf)


modsumm <- outdf %>% 
  group_by(ngen, death_rate, gam_scale, mixmod) %>% 
  summarise(meanrignt = mean(corrBIC))

# hom, GDD, short-lived increase accuracy
modsumm <- outdf %>% 
  filter(death_rate == 0.0075, mixmod == "hom") %>%
  group_by(ngen, gam_scale, pois_lam) %>% 
  summarise(meanrignt = mean(corrBIC))
modsumm %>% data.frame()

# library(randomForest)
#   
# ?randomForest
# datrf <- as.data.frame(unclass(outdf)) %>% 
#   dplyr::select(corrAIC, surv_missing, gam_scale, gam_smooth, ngen, death_rate, pois_lam, detprob_model, mixmod, region)
# datrf$corrAIC <- as.factor(datrf$corrAIC)
# 
# modrf <- randomForest(x = datrf[, -1], y = datrf[, 1], ntree = 500, importance = TRUE)
# summary(modrf)

# graveyard
  
  mods_work <- summ_mods %>% 
    ungroup() %>% 
    select(SiteID, Year, maxgen, model) %>% 
    distinct() %>% 
    dplyr::filter(complete.cases(.))
  
  mods_all <- expand.grid(1:(true_ngen + 1), c("Skew.normal", "Normal", "Skew.t", "t", "Mclust_E", "Mclust_V"))
  names(mods_all) <- c("maxgen", "model")
  mods_all <- mods_all[-which(mods_all$maxgen == 1 & mods_all$model == "Mclust_V"), ] 
  
  mods_skip <- summ_mods %>% 
    ungroup() %>% 
    select(SiteID, Year) %>% 
    distinct() %>% 
    group_by(SiteID, Year) %>% 
    do(., newcols = mods_all) %>% 
    unnest() %>% 
    anti_join(mods_work)

  best_mods <- summ_mods %>% 
    ungroup() %>% 
    filter(is.na(bic) == FALSE) %>% 
    mutate(bic = ifelse(bic < 0 & model %in% c("Mclust_E", "Mclust_V"), -bic, bic)) %>% 
    group_by(SiteID, Year, model, maxgen) %>%
    mutate(badmixmod = sum(mixmod_flag, na.rm = TRUE),
           redundant = ifelse(maxgen > 1, min(diff(mu, lag = 1), na.rm = TRUE), NA),
           nearzerosigma = min(sigma2, na.rm = TRUE)) %>%
    filter(badmixmod == 0, nearzerosigma > 100) %>% 
    group_by(SiteID, Year, model) %>%
    mutate(within_model_bic = bic - min(bic, na.rm = TRUE),
           within_model_aic = aic - min(aic, na.rm = TRUE)) %>%
    group_by(SiteID, Year, maxgen) %>%
    mutate(within_maxgen_bic = bic - min(bic, na.rm = TRUE),
           within_maxgen_aic = aic - min(aic, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(right_ngen = ifelse(maxgen == true_ngen, "yes", "no")) %>% 
    distinct() # removes duplicated mclust_E for ngen == 1

  bestest <- best_mods %>% 
    group_by(model) %>% 
    filter(right_ngen == "yes") %>% 
    dplyr::select(SiteID, Year, model, within_maxgen_aic, within_maxgen_bic, within_model_aic,
                  within_model_bic) %>% 
    distinct() %>% 
    summarise(meanaic = mean(within_maxgen_aic),
              meanbic = mean(within_maxgen_bic),
              meanaicngen = mean(within_model_aic),
              meanbicngen = mean(within_model_bic),
              n = length(unique(paste0(SiteID, Year))),
              rightaic = length(which(within_model_aic == 0)),
              rightbic = length(which(within_model_bic == 0)),
              bestaic = length(which(within_maxgen_aic == 0)),
              bestbic = length(which(within_maxgen_bic == 0)))

  # plot to assess model selection of ngen vs siteYear N and weight of last gen
  plot_best <- truth %>% 
    filter(Gen == true_ngen) %>% 
    mutate(SiteID = as.factor(SiteID),
           Year = as.factor(as.character(Year))) %>% 
    dplyr::select(SiteID, Year, gen_weight, M) %>% 
    right_join(best_mods) %>% 
    group_by(SiteID, Year, model) %>% 
    filter(within_model_bic == 0) %>% 
    dplyr::select(SiteID, Year, model, gen, gen_weight, M, maxgen, w) %>% 
    filter(gen == maxgen) %>% 
    distinct()
  
  plt <- ggplot(plot_best, aes(x = M, y = gen_weight, color = as.factor(maxgen))) +
    geom_point() +
    theme_bw() +
    facet_wrap(~model, ncol = 3)
  plt
  
  
  summ_mods$seed <- stringr::str_split_fixed(string = f, 
                                             pattern = "_", n = 2)[2]
  
  
  
  
  outlist[[length(outlist)+1]] <- summ_mods
}
res <- bind_rows(outlist)

right_ngen <- res %>% 
  mutate(seed = as.numeric(seed)) %>% 
  left_join(params, by = "seed") %>% 
  group_by(seed, SiteID, Year, model, surv_missing, gam_smooth, detprob_model) %>% 
  summarise(correct = length(which(maxgen == 2))/2) %>% 
  group_by(model) %>% 
  summarise(propright = sum(correct, na.rm = TRUE) / length(correct))
right_ngen



# Extract summary stats for each generation's distribution
summ_mods <- results %>% 
  ungroup() %>% 
  group_by(SiteID, Year) %>% 
  do(Summ_mixmod(.))
# model errors lead to blank rows of NA
# anti_join to learn what these are?
mods_work <- summ_mods %>% 
  select(SiteID, Year, maxgen, model) %>% 
  distinct()
mods_all <- expand.grid(1:3, unique(summ_mods$model))
names(mods_all) <- c("maxgen", "model")
mods_all <- mods_all[complete.cases(mods_all), ]

mods_skip <- summ_mods %>% 
  ungroup() %>% 
  select(SiteID, Year) %>% 
  distinct() %>% 
  group_by(SiteID, Year) %>% 
  do(., newcols = mods_all) %>% 
  unnest() %>% 
  anti_join(mods_work)

mods_skip <- test[-which(test$maxgen == 1 & test$model == "Mclust_V"), ]



# chooses the first generation of each group, not optimal
# might also want to summarise bic delta for each worse model in each grouping
best_mods <- summ_mods %>% 
  ungroup() %>% 
  filter(is.na(bic) == FALSE) %>% 
  mutate(gen = ifelse(is.na(gen), 1, gen),
         maxgen = ifelse(is.na(maxgen), 1, maxgen),
         bic = ifelse(bic < 0, -bic, bic)) %>% 
  group_by(SiteID, Year, model, maxgen) %>%
  mutate(badmixmod = sum(mixmod_flag)) %>% 
  filter(badmixmod == 0) %>% 
  group_by(SiteID, Year) %>% 
  # arrange(bic)
  filter(bic == min(bic))


# right # gen?
res1 <- best_mods %>% 
  group_by(SiteID, Year, model) %>% 
  summarise(maxgen = maxgen[1])


size <- adjcounts %>% 
  group_by(SiteID, Year) %>% 
  summarise(popsize = sum(adjY)) %>% 
  merge(res1)

mutest <- adjcounts %>% 
  select(SiteID, Year, Site_RE, Year_RE, M) %>% 
  distinct() %>% 
  left_join(best_mods) %>% 
  mutate(Total_RE = Site_RE + Year_RE)

plt <- ggplot(mutest, aes(x = Total_RE, y = curve_max)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~gen, scales = "free_y")
plt
