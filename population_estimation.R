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
gdd <- readRDS("data/dailyDD.rds")
# gdd <- readRDS("../ohiogdd/dailyDD.rds")

gdd <- left_join(gdd, sites) %>% 
  dplyr::select(SiteID, SiteDate, degday530, lat, lon, maxT) %>% 
  mutate(Year = year(SiteDate),
         DOY = yday(SiteDate)) %>% 
  group_by(SiteID, Year) %>% 
  arrange(DOY) %>% 
  mutate(AccumDD = cumsum(degday530))

# remove some sites that are clustered close together
uniqsites <- gdd %>% ungroup() %>%  dplyr::select(SiteID, lat, lon) %>% distinct()

points_matrix <- as.matrix(dist(uniqsites[, c("lat", "lon")], diag = TRUE, upper = TRUE))
points_matrix[lower.tri(points_matrix, diag=TRUE)] <- NA
points_matrix <- points_matrix <= 0.07

v <- colSums(points_matrix, na.rm=TRUE) == 0
uniqsites <- uniqsites[v, ]

gdd <- gdd %>% 
  filter(SiteID %in% uniqsites$SiteID)

#####
# define parameters
nsite <- c(4, 16) # up to 140 OH sites
nyear <- c(4, 16) # up to 35 years in Daymet data
nsurv <- 30
surv_missing <- c(0, .2, .4) # remove surveys from all counts
# site_missing <- .2 # turnover across years
# group_struct <- c("all", "site", "year", "siteyear")
gam_smooth <- c("none", "interpolate", "preds_8day", "preds_4day")

ngen <- c(1:4)
gen_size <- c("equal", "inc", "dec")
# volt_flex <- "Y"
# gen_ddreq <- 600 # depends on ngen
peak_sd <- c(25, 100) # low and high overlap?
death_rate <- c(.6, .8) # on weekly rate

# site total population dispersion
negbin_mu <- 100
negbin_disp <- 1
# detection probability parameters (logit scale)
detprob_b0 <- -7
detprob_b1 <- .6
detprob_b2 <- -.01
detprob_model <- c("known", "covariate", "none")

# site/year variation in mu
site_mu_sd <- c(10, 50)
year_mu_sd <- c(10, 50)

# # nonlinear det prob
# x <- seq(-5, 45, .1)
# y <- plogis(-10 + x * .6 + -.01 * x^2)
# plot(x, y)

params <- list(nsite = nsite, nyear = nyear, nsurv = nsurv, surv_missing = surv_missing, 
               gam_smooth = gam_smooth, ngen = ngen, gen_size = gen_size, peak_sd = peak_sd, 
               death_rate = death_rate, negbin_mu = negbin_mu, negbin_disp = negbin_disp, 
               detprob_b0 = detprob_b0, detprob_b1 = detprob_b1, detprob_b2 = detprob_b2, 
               detprob_model = detprob_model, site_mu_sd = site_mu_sd, year_mu_sd = year_mu_sd)
params <- expand.grid(params)

# some parameter combinations don't make sense, try to exclude to save time
params <- params[-which(params$ngen == 1 & params$gen_size %in% c("inc", "dec")), ]
params <- params[-which(params$surv_missing == 0 & params$gam_smooth == "interpolate"), ]

# for random numbers
params$seed <- 1:nrow(params)
# params <- params[1:2000, ]


done <- c(list.files('results_0'), list.files('results_1'))

done <- as.numeric(stringr::str_split_fixed(string = done, pattern = "_", n = 2)[,2])

params <- params[-which(params$seed %in% done), ]
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
# plt <- ggplot(counts, aes(x = AccumDD, y = Y, group = as.factor(Year), color = as.factor(Year))) +
#   geom_path() +
#   facet_wrap(~SiteID, ncol = 4, scales = "free_y")
# plt
# 
# # GAM interpolation/prediction
# adjcounts <- Adjust_Counts(data = test, counts)
# plt <- ggplot(adjcounts, aes(x = AccumDD, y = adjY, group = Year, color = Year)) +
#   geom_point() +
#   facet_wrap(~SiteID, ncol = 4, scales = "free_y")
# plt
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


ncores <- 30

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
                                "Summ_mixmod", "mainDir"),
                    .inorder = FALSE,
                    .options.multicore = mcoptions) %dopar% {
                      
                      test <- params[sim, ]
                      
                      # make directories to store results
                      subDir <- paste("results", floor(test$seed / 1000), sep = "_")
                      
                      setwd(mainDir)
                      if (file.exists(subDir)){
                        setwd(file.path(mainDir, subDir))
                      } else {
                        dir.create(file.path(mainDir, subDir))
                        setwd(file.path(mainDir, subDir))
                      }
                      
                      
                      
                      counts <- Simulate_Counts(data = test, gdd = gdd)
                      
                      # GAM interpolation/prediction
                      adjcounts <- Adjust_Counts(data = test, counts)
                      
                      
                      results <- adjcounts %>%
                        group_by(SiteID, Year) %>% 
                        mutate(weeks_obs = length(which(round(adjY) > 0)),
                               total_obs = sum(round(adjY))) %>% 
                        filter(weeks_obs >= 2, total_obs >= 5) %>% 
                        do(mixmods = CompareMixMods(dat = ., mvmax = (test$ngen + 1), seed = test$seed))
                      
                      # Extract summary stats for each generation's distribution
                      summ_mods <- results %>% 
                        ungroup() %>% 
                        group_by(SiteID, Year) %>% 
                        do(Summ_mixmod(.))
                      
                      outlist <- list(test, counts, adjcounts, summ_mods)   
                      saveRDS(object = outlist, file = paste("popest", test$seed, sep = "_"))
                      
                      rm(counts, adjcounts, summ_mods, results, outlist)
                      gc()
                      
                      # hope this is the only thing returned
                      return(test$seed)
                    } # close dopar


if(.Platform$OS.type == "windows"){
  stopCluster(cl)
}


# with simulation results in hand, next:
# 1. when is number of generations correct compared to params?
# 2. what is error in generation weights estimated compared to simulated "truth"?
# 3. what is error in phenology for each site/year compared to "truth"?
# 4. what is error in site/year population size compared to "truth" using trapezoid rule?



fs <- list.files(pattern = "popest_", recursive = TRUE, full.names = TRUE)

outlist <- list()
for (f in fs){
  tmp <- readRDS(f)
  param <- tmp[[1]]
  counts <- tmp[[2]]
  adjcounts <- tmp[[3]]
  summ_mods <- tmp[[4]]
  # 1. when is number of generations correct compared to params? 
  # What others errors, like degenerate modes with redundant mu?
  # Compare across models for best fit to true ngen, compare within models to see if correct ngen selected
  true_ngen <- param$ngen
    
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
    mutate(bic = ifelse(bic < 0, -bic, bic)) %>% 
    group_by(SiteID, Year, model, maxgen) %>%
    mutate(badmixmod = sum(mixmod_flag, na.rm = TRUE),
           degenerate = ifelse(length(which((diff(mu, lag = 1) < 50) == TRUE)) == 0, 
                               rep(0, length(mu)), 
                               rep(1, length(mu)))) %>%
    mutate(degenerate = sum(degenerate, na.rm = TRUE)) %>% 
    filter(badmixmod == 0, degenerate == 0) %>% 
    group_by(SiteID, Year, model) %>% 
    # arrange(bic)
    mutate(within_delta_bic = bic - min(bic),
           within_delta_aic = aic - min(aic, na.rm = TRUE)) %>% 
    group_by(SiteID, Year, maxgen) %>% 
    mutate(among_delta_bic = bic - min(bic),
           among_delta_aic = aic - min(aic, na.rm = TRUE))

  
  truth <- Simulate_Truth(data = param, counts = counts, gdd = gdd)
  
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
