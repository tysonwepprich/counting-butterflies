# Functions for population estimation

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}


# calculate Proportional variation (similar to CV)
PropVariation <- function(numvect){
  reldiff <- sapply(numvect, function(x) sapply(numvect, function(y) 1 - min(x, y, na.rm = TRUE) / max(x, y, na.rm = TRUE)))
  pv <- mean(as.numeric(reldiff[upper.tri(reldiff)]), na.rm = TRUE)
  return(pv)
}

# UKBMS count index approximation
TrapezoidIndex <- function(timescale, counts){
  dat <- data.frame(t = timescale, y = counts)
  dat <- arrange(dat, timescale)
  if(nrow(dat) == 1){
    return(dat$y)
  }else{
    temp <- rep(NA, nrow(dat)-1)
    for (i in 2:length(dat$t)){
      temp[i-1] <- (dat$y[i] + dat$y[i-1]) * (dat$t[i] - dat$t[i-1]) / 2
    }
  }
  return(sum(temp))
}


# Uses generalized Zonneveld model from Calabrese to model emergence and total abundance of one generation
Abund_Curve <- function(t = t, alpha = .3, beta = 1, mu = 10, sig = 1){
  
  z <- exp((t - mu) / beta) / (1 + exp((t - mu) / beta))
  a <- 1 + alpha * beta
  b <- sig - alpha * beta
  
  # incomplete beta
  ibeta <- function(z, a, b){ pbeta(z, a, b) * beta(a, b) }
  
  e <- (sig * (exp((t - mu)/beta)))/(beta * (1 + exp((t - mu)/beta)))^(sig + 1)
  e <- e / sum(e)
  y <- alpha * sig * exp(-alpha * (t - mu)) * ibeta(z, a, b)
  y <- y / sum(y)
  return(data.frame(t, y, e))
}

# Summary function for a generation's count distribution
# To compare mixture model results when skewed
Summ_curve <- function(t, y, quants = c(.1, .25, .5, .75, .9)){
  n <- length(y)
  cdf <- cumsum(y) / sum(y)
  cmean <- weighted.mean(t, y)
  cmax <- t[which(y == max(y))][1]
  estN <- sum(y)
  
  cquant <- rep(NA, length(quants))
  for (i in seq_along(quants)){
    cquant[i] <- t[which(abs(cdf - quants[i]) == min(abs(cdf - quants[i])))][1]
  }
  
  df <- matrix(data = c(cmean, cmax, cquant, estN), nrow = 1)
  df <- data.frame(df)
  names(df) <- c("curve_mean", "curve_max", paste0("curve_q", quants), "estN")
  return(df)
}


# data argument is row of parameters
# output list of counts, truth
Simulate_Counts <- function(param, gdd){
  set.seed(param$seed) # use for all rows?
  # randomly select sites and years to use their coordinates and historical gdd accumulation
  # sites <- sample(x = unique(gdd$SiteID), size = param$nsite, replace = FALSE)
  
  years <- sample(x = unique(gdd$Year), size = param$nyear, replace = FALSE)
  
  # sample sites so each of 4 regions equally represented
  sites <- gdd %>% 
    ungroup() %>% 
    dplyr::select(SiteID, region) %>% 
    distinct() %>% 
    group_by(region) %>% 
    sample_n(param$nsite / 4) %>% 
    pull(SiteID)
  
  dddat <- gdd %>% 
    filter(SiteID %in% sites, Year %in% years) # %>% 
    # filter(DOY %in% (seq(90, 293, 7) + sample.int(n=6, size=30, replace=TRUE)))
    # filter(DOY %in% (seq(90, 293, 3)))
  
  
  # make some years have voltinism variation
  # originally did this by available gdd, but this makes last generations small for most years
  # NEW: split years into 3 treatments: 0, 1, 2 multipliers of site variation
  val <- length(years) / 3
  gddannvar <- dddat %>% 
    group_by(SiteID, Year) %>%
    summarise(maxgdd = max(AccumDD)) %>% 
    group_by(Year) %>% 
    summarise(anngdd = mean(maxgdd)) %>% 
    mutate(zanngdd = sample(x = rep.int(0:2, times = val), size = length(years), replace = FALSE))
    # mutate(zanngdd = runif(length(years), .5, 2))
    # mutate(zanngdd = (anngdd - min(anngdd)) / (max(anngdd) - min(anngdd)))
  
  
  
  # voltinism varying by site
  ngen <- param$ngen
  if(ngen == 1){
    gen_ddreq <- param$emerg_dd
    gen_relsize = 1
    gen_relsprd = 1
  }else{
    gen_ddreq <- c(param$emerg_dd, param$emerg_dd + param$gen_dd * 1:(ngen - 1))
    gen_relsize <- switch(as.character(param$gen_size), 
                          "equal" = rep(1, ngen),
                          "inc" = seq(1, ngen, length.out = ngen),
                          "dec" = seq(ngen, 1, length.out = ngen))
    gen_relsprd <- seq(1, 1.5, length.out = ngen)
  }
  # simulate phenology for each generation and combine
  year_var <- rnorm(param$nyear, mean = 0, sd = param$year_mu_sd)
  site_var <- rnorm(param$nsite, mean = 0, sd = param$site_mu_sd)
  year_Mvar <- seq(from = .5, to = 2, length.out = length(years))
  site_Mvar <- seq(from = 100, to = 1000, length.out =  length(sites))
  
  dflist <- list()
  for(y in seq_along(years)){
    for(z in seq_along(sites)){
      yr <- years[y]
      site <- sites[z]
      # site re includes dependence on latitude, lat of 40 assumed to be baseline
      site_lat <- gdd$lat[gdd$SiteID == site][1]
      site_re <- 1000 - site_lat * 25 + site_var[z]
      year_re <- year_var[y]
      site_M <- site_Mvar[z]
      year_M <- year_Mvar[y]
      if(param$ngen == 1){
        year_volt <- 1
      }else{
        year_volt <- gddannvar %>% 
          filter(Year == yr) %>% 
          select(zanngdd) %>% 
          as.numeric()
      }
      
      for (g in 1:param$ngen){
        simdat <- dddat %>% filter(SiteID == site, Year == yr)
        simt <-  simdat$AccumDD
        simalpha <- param$death_rate
        simbeta <- gen_relsprd[g] * param$peak_sd
        simmu <- (gen_ddreq[g] + site_re + year_re)
        df <- Abund_Curve(t = simt, alpha = simalpha, 
                          beta = simbeta, mu = simmu, sig = 1)
        df <- cbind(data.frame(simdat), df)
        df$Gen <- g
        df$Site_RE <- site_re
        df$Year_RE <- year_re
        df$Site_M <- site_M
        df$Year_M <- year_M
        df$gen_weight <- gen_relsize[g] * ifelse(g == param$ngen, year_volt * ((42 - df$lat) / 3.35), 1)
        dflist[[length(dflist)+1]] <- df
      }
    }
  }
  dfall <- bind_rows(dflist) %>% 
    arrange(SiteID, SiteDate)
  
  
  
  # counting process
  counts <- dfall %>% 
    group_by(SiteID, SiteDate, lat, lon, Year, DOY, AccumDD, region, maxT, minT, Site_RE, Year_RE, Site_M, Year_M) %>% 
    summarise(RelProp = sum(y * gen_weight),
              DP = plogis(param$detprob_b0 + param$detprob_b1 * maxT[1] + param$detprob_b2 * maxT[1]^2)) %>% 
    group_by(SiteID, Year) %>% 
    mutate(RelProp = RelProp / sum(RelProp),
           # M = rnbinom(1, mu = param$negbin_mu, size = param$negbin_disp),
           M = rpois(1, lambda = Site_M * Year_M * ngen), # should also multiply by ngen to standardize M in each gen, but overlooked
           N = rpois(length(RelProp), lambda = RelProp * M),
           Y = rbinom(length(N), size = N, prob = DP)) %>% 
    data.frame()
  
  # ggplot(counts, aes(x = DOY, y = N, color = lat)) + 
  #   geom_point(alpha = .3) + 
  #   facet_wrap(~Year)
  
  # truth
  true_weight <- dfall %>% 
    group_by(SiteID, SiteDate) %>% 
    mutate(gen_weight = gen_weight / sum(gen_weight))
  true_phen_gdd <- dfall %>% 
    ungroup() %>% 
    group_by(SiteID, Year, Gen) %>% 
    do(Summ_curve(t = .$AccumDD, y = .$y)) %>% 
    mutate(resultTimescale = "GDD")
  true_phen <- dfall %>% 
    ungroup() %>% 
    group_by(SiteID, Year, Gen) %>% 
    do(Summ_curve(t = .$DOY, y = .$y)) %>% 
    mutate(resultTimescale = "DOY") %>% 
    bind_rows(true_phen_gdd)
  true_N <- counts %>% 
    ungroup() %>% 
    dplyr::select(SiteID, Year, Site_RE, Year_RE, Site_M, Year_M, M) %>% 
    distinct()
  
  truthout <- true_weight %>% 
    ungroup() %>% 
    dplyr::select(SiteID, Year, Gen, gen_weight) %>% 
    distinct() %>% 
    right_join(true_phen) %>% 
    right_join(true_N)
  
  out <- list(counts, truthout)
  return(out)
}


# accounts for det prob, GAM smoothing
Adjust_Counts <- function(param, counts){
  set.seed(param$seed)
  # for GAM, need factors
  counts$SiteYearID <- as.factor(paste(counts$SiteID, counts$Year, sep = "_"))
  counts$SiteID <- as.factor(counts$SiteID)
  counts$region <- as.factor(counts$region)
  counts$Year <- as.factor(as.character(counts$Year))
  counts$RegYear <- as.factor(paste(counts$region, counts$Year, sep = "_"))
  if(param$gam_scale == "DOY"){
    counts$Timescale <- counts$DOY
  }
  if(param$gam_scale == "GDD"){
    counts$Timescale <- counts$AccumDD
  }
  
  
  # counts were simulated each day, reduce to every 4 days
  counts_4day <- counts %>% 
    group_by(SiteYearID) %>% 
    filter(DOY %in% (seq(90, 294, 4) + sample.int(n=3, size=52, replace=TRUE))) %>% 
    ungroup()
  
  counts_8day <- counts_4day %>% 
    group_by(SiteYearID) %>% 
    slice(seq(1, length(DOY), 2)) %>% 
    ungroup()

  adjcounts <- counts_8day %>% 
    sample_frac(size = 1 - param$surv_missing) %>% 
    droplevels()
  
  if(param$detprob_model == "known"){
    adjcounts <- mutate(adjcounts, adjY = Y / DP)
  }else{
    adjcounts <- mutate(adjcounts, adjY = Y)
  }
  
  if(param$gam_smooth == "none"){
    adjcounts$gam_flag <- 0
    summod <- NA
    return(list(adjcounts, summod))
  }else{
    # model GAM once for next two cases
    
    safe_gam <- purrr::safely(gam)
    
    if(param$detprob_model %in% c("none", "known")){
      gammod <- safe_gam(adjY ~ 
                           te(lat, lon, Timescale, bs = c("tp", "cr"), k = c(5, 30), d = c(2, 1)) +
                           s(SiteID, bs = "re") + 
                           s(RegYear, Timescale, bs = "fs", k = 5, m = 1),
                         family = nb(theta = NULL, link = "log"),
                         # family = poisson(link = "log"),
                         data = adjcounts,
                         method = "REML", 
                         optimizer = c("outer", "newton"), 
                         control = list(maxit = 500))
    }
    if(param$detprob_model == "covariate"){
      gammod <- safe_gam(adjY ~ 
                           s(maxT) +
                           te(lat, lon, Timescale, bs = c("tp", "cr"), k = c(5, 30), d = c(2, 1)) +
                           s(SiteID, bs = "re") +
                           s(RegYear, Timescale, bs = "fs", k = 5, m = 1),
                         family = nb(theta = NULL, link = "log"),
                         # family = poisson(link = "log"),
                         data = adjcounts,
                         method = "REML", 
                         optimizer = c("outer", "newton"), 
                         control = list(maxit = 500))
    }
  }
  
  if(is.null(gammod$error)){
    if(param$gam_smooth == "preds_4day"){
      newdata <- counts_4day
    }else{
      newdata <- counts_8day
    }
    
    # HOW TO DEAL WITH TEMP COVARIATE HERE?
    # If excluded, predicted counts are really low. Could scale and use mean like before
    # if(param$detprob_model == "covariate"){
    #   newdata$adjY <- predict(gammod$result, newdata = newdata, type = "response", exclude = "s(maxT)")
    # }else{
    #   newdata$adjY <- predict(gammod$result, newdata = newdata, type = "response")
    # }
    
    mod <- gammod$result
    summod_all <- summary(gammod$result)
    summod <- summod_all[c("residual.df", "r.sq", "n", "dev.expl", "p.table", "s.table", "sp.criterion")]
    summod$nb_theta <- gammod$result$family$getTheta(TRUE)
    
    # prediction with response
    adjY <- predict(gammod$result, newdata = newdata, type = "response")
    newdata$adjY <- qnbinom(.5, size = mod$family$getTheta(TRUE), mu = adjY)
    # # prediction with simulated counts, stochastic but integers (if n is odd)
    # Xp <- predict.gam(object = mod, newdata = newdata, type="lpmatrix") ## map coefs to fitted curves
    # beta <- coef(mod)
    # Vb   <- vcov(mod) ## posterior mean and cov of coefs
    # n <- 5 # choose number of simulations
    # mrand <- MASS::mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    # ilink <- family(mod)$linkinv
    # linklppreds <- Xp %*% t(mrand)
    # nbpreds <- apply(X = linklppreds,
    #                  MARGIN = 1,
    #                  FUN = function(x){
    #                    # temp <- sort(x)
    #                    # bounds <- quantile(1:n, probs = c(0.025, 0.975))
    #                    # x <- temp[bounds[1]:bounds[2]]
    #                    x <- ilink(x)
    #                    x <- rnbinom(n = length(x),
    #                                 mu = x,
    #                                 size = mod$family$getTheta(TRUE))
    #                    x <- quantile(x, .5)
    #                    return(x)
    #                  })
    # newdata$adjY <- nbpreds
    
    if(param$gam_smooth == "interpolate"){
      interps <- anti_join(newdata, adjcounts, by = c("SiteID", "SiteDate"))
      adjcounts <- bind_rows(adjcounts, interps)
      adjcounts$gam_flag <- 0
      return(list(adjcounts, summod))
    }
    
    if(param$gam_smooth %in% c("preds_4day", "preds_8day")){
      adjcounts <- newdata
      adjcounts$gam_flag <- 0
      return(list(adjcounts, summod))
      }
  }else{
    # error in gam fit
    adjcounts$gam_flag <- 1
    summod <- NA
    return(list(adjcounts, summod))
  }
}


Adjust_Counts_AfterGAM <- function(param, counts, preds){
  counts <- counts %>% 
    filter(SiteYear %in% preds$SiteYear) %>% 
    mutate(Year = as.numeric(as.character(Year)),
           adjY = Total)
  
  # counts were simulated each day, reduce to every 4 days
  counts_4day <- preds %>% 
    ungroup()
  
  counts_8day <- counts_4day %>% 
    group_by(SiteYear, Week) %>% 
    sample_n(size = 1) %>% 
    ungroup()
  
  if(param$gam_smooth == "none"){
    adjcounts <- counts
    return(adjcounts)  
  }else{
    
    if(param$gam_smooth == "interpolate"){
      countweeks <- counts %>% 
        select(SiteYear, Week) %>% 
        distinct()
      predweeks <- counts_8day %>% 
        select(SiteYear, Week) %>% 
        distinct()
      
      interps <- anti_join(predweeks, countweeks, by = c("SiteYear", "Week")) %>% 
        filter(Week <= 30) %>% 
        mutate(SYW = paste(SiteYear, Week, sep = "_"))
      
      interps2 <- counts_8day %>% 
        mutate(SYW = paste(SiteYear, as.character(Week), sep = "_")) %>% 
        filter(SYW %in% interps$SYW)
      
      adjcounts <- bind_rows(counts, interps2)
      
      return(adjcounts)
    }
    
    if(param$gam_smooth == "preds_4day"){
      adjcounts <- counts_4day
      return(adjcounts)
    }
    
    if(param$gam_smooth == "preds_8day"){
      adjcounts <- counts_8day
      return(adjcounts)
    }
  }
}


# fit different mixture models, output a list of model results
# adapt this so it will work with different groupings, from statewide to siteyear
CompareMixMods <- function(dat, param){
  set.seed(param$seed) # use for all rows?
  # print(c(as.character(dat$Year)[1], as.character(dat$SiteID)[1]))
  dd <- dat$Timescale
  y <- round(dat$adjY)
  dd_dist <- rep(dd, y)

  mvmin <- 1
  mvmax <- ifelse(param$mod_maxgen == "ngen", param$ngen, param$ngen + 1) # tried ngen + 1 in simulation and saw overfitting
  gens <- c(mvmin:mvmax)
  maxtry <- 5 # repeating smsn.mix function if errors
  # out <- as.list(mvmin:mvmax)
  # names(out) <- paste0("gen", c(mvmin:mvmax))
  out <- list()
  for (i in seq_along(gens)){
    safe_mix <- safely(smsn.mix)
    if(param$mixmod == "skew"){
      # skew normal
      worked <- 0
      tries <- 1
      while(worked == 0 & tries <= maxtry){
        mod <- safe_mix(dd_dist, nu = 5, g = i, get.init = TRUE, family = "Skew.normal", 
                        calc.im = FALSE, obs.prob = TRUE)
        if(is.null(mod$error)){
          worked <- 1
        }else{
          tries <- tries + 1
        }
      }
    }
    
    # # Normal het from smsn
    # worked <- 0
    # tries <- 1
    # while(worked == 0 & tries <= maxtry){
    #   mod_norm <- safe_mix(dd_dist, nu = 5, g = i, get.init = TRUE, family = "Normal", 
    #                        calc.im = FALSE, obs.prob = TRUE, kmeans.param = list(n.start = 5))
    #   if(is.null(mod_norm$error)){
    #     worked <- 1
    #   }else{
    #     tries <- tries + 1
    #   }
    # }
    # 
    
    # # T skew
    # worked <- 0
    # tries <- 1
    # while(worked == 0 & tries <= maxtry){
    #   mod_skt <- safe_mix(dd_dist, nu = 5, g = i, get.init = TRUE, family = "Skew.t", 
    #                       calc.im = FALSE, obs.prob = TRUE, kmeans.param = list(n.start = 5))
    #   if(is.null(mod_skt$error)){
    #     worked <- 1
    #   }else{
    #     tries <- tries + 1
    #   }
    # }
    # 
    # # T het from smsn
    # worked <- 0
    # tries <- 1
    # while(worked == 0 & tries <= maxtry){
    #   mod_t <- safe_mix(dd_dist, nu = 5, g = i, get.init = TRUE, family = "t", 
    #                     calc.im = FALSE, obs.prob = TRUE, kmeans.param = list(n.start = 5))
    #   if(is.null(mod_t$error)){
    #     worked <- 1
    #   }else{
    #     tries <- tries + 1
    #   }
    # }
    
    # mclust
    safe_mclust <- safely(Mclust)
    # hom
    if(param$mixmod == "hom"){
      worked <- 0
      tries <- 1
      while(worked == 0 & tries <= maxtry){
        mod <- safe_mclust(dd_dist, G=i, modelNames = "E")
        if(is.null(mod$error) & !is.null(mod$result)){
          worked <- 1
        }else{
          tries <- tries + 1
        }
      }
      if(is.null(mod$result) & is.null(mod$error)){
        mod$error <- "none fitted"
      }
    }
    
    # het
    if(param$mixmod == "het"){
      worked <- 0
      tries <- 1
      while(worked == 0 & tries <= maxtry){
        mod <- safe_mclust(dd_dist, G=i, modelNames = "V")
        if(is.null(mod$error) & !is.null(mod$result)){
          worked <- 1
        }else{
          tries <- tries + 1
        }
      }
      if(is.null(mod$result) & is.null(mod$error)){
        mod$error <- "none fitted"
      }
    }
    out[[i]] <- list(mod)
    # out[[i]] <- list(mod_sknorm, mod_norm, mod_skt, mod_t, mod_mc_hom, mod_mc_het)
  }
  return(out)
}
  
RightNumGen <- function(mixmods, param, reg, yr){
  test <- Summ_mixmod(mixmods)
  
  # problem with model classes
  test$model <- param$mixmod
  
  best_mods <- test %>% 
    ungroup() %>% 
    filter(is.na(bic) == FALSE) %>% 
    mutate(bic = ifelse(bic < 0 & model %in% c("hom", "het"), -bic, bic)) %>% 
    group_by(model, maxgen) %>%
    mutate(badmixmod = sum(mixmod_flag, na.rm = TRUE),
           redundant = ifelse(maxgen > 1, min(diff(mu, lag = 1), na.rm = TRUE), NA),
           nearzerosigma = min(sigma2, na.rm = TRUE)) %>%
    # filter(badmixmod == 0, nearzerosigma > 100) %>% 
    group_by(model) %>%
    mutate(within_model_bic = bic - min(bic, na.rm = TRUE),
           within_model_aic = aic - min(aic, na.rm = TRUE)) %>%
    ungroup() %>% 
    mutate(right_ngen = ifelse(maxgen == param$ngen, "yes", "no"),
           region = reg,
           modyear = yr) %>% 
    distinct()
  
  return(best_mods)
}



# Extract results from mixture model objects
Summ_mixmod <- function(df){
  # print(df$SiteID)
  # print(df$Year)
  allmods <- flatten(df)
  alldf <- list()
  for (m in seq_along(allmods)){
    mod <- allmods[[m]]
    if(is.null(mod$error) && !is.null(mod$result)){
      modtype <- attr(mod$result, "class")
      if(modtype == "Mclust"){
        # mclust param
        modtype <- paste(modtype, mod$result$modelName, sep = "_")
        if(modtype == "Mclust_X") modtype <- "Mclust_E"
        #calc AIC
        res <- data.frame(mu = mod$result$parameters$mean, sigma2 = mod$result$parameters$variance$sigmasq, shape = NA,
                          w = mod$result$parameters$pro, nu = NA, aic = 2*mod$result$df - 2*mod$result$loglik,
                          bic = mod$result$bic,
                          edc = NA, icl = NA, model = modtype)
        res <- res %>% arrange(mu)
        rownames(res) <- 1:nrow(res)
        
        # simulate each generation
        # extract curve summary
        gen_res <- list()
        for (i in seq_along(as.numeric(rownames(res)))){
          x <- seq(0, 3000, 10)
          y <- dnorm(x = x, mean = res$mu[i], sd = sqrt(res$sigma2[i]))
          if(max(y) == 0){
            gen_res[[i]] <- data.frame(curve_mean = NA, curve_max = NA, 
                                       curve_q0.1 = NA, curve_q0.5 = NA, curve_q0.9 = NA, mixmod_flag = 1)
          }else{
            gen_res[[i]] <- Summ_curve(t = x, y = y)
            gen_res[[i]]$mixmod_flag <- 0
          }
        }
        gen_res_df <- bind_rows(gen_res) %>% 
          mutate(gen = as.numeric(rownames(res)),
                 maxgen = max(gen)) %>% 
          bind_cols(res)
        
      }else{
        # mixsmsn
        
        res <- data.frame(mu = mod$result$mu, sigma2 = mod$result$sigma2, shape = mod$result$shape,
                          w = mod$result$pii, nu = mod$result$nu, aic = mod$result$aic, bic = mod$result$bic,
                          edc = mod$result$edc, icl = mod$result$icl, model = modtype)
        res <- res %>% arrange(mu)
        rownames(res) <- 1:nrow(res)
        
        # simulate each generation and
        # extract curve summary
        gen_res <- list()
        for (i in seq_along(as.numeric(rownames(res)))){
          x <- seq(0, 3000, 10)
          if (modtype %in% c("Skew.normal", "Normal")){
            y <- mixsmsn:::dSN(y = x, mu = res$mu[i], sigma2 = res$sigma2[i], shape = res$shape[i])
          }else if(modtype %in% c("t", "Skew.t")){
            y <- mixsmsn:::dt.ls(x = x, loc = res$mu[i], sigma2 = res$sigma2[i], shape = res$shape[i], nu = res$nu[i])
          }
          if(max(y) == 0 | length(which(is.finite(y))) == 0){
            gen_res[[i]] <- data.frame(curve_mean = NA, curve_max = NA, 
                                       curve_q0.1 = NA, curve_q0.5 = NA, curve_q0.9 = NA, mixmod_flag = 1)
          }else{
          gen_res[[i]] <- Summ_curve(t = x, y = y)
          gen_res[[i]]$mixmod_flag <- 0
          }
        }
        gen_res_df <- bind_rows(gen_res) %>% 
          mutate(gen = as.numeric(rownames(res)),
                 maxgen = max(gen)) %>% 
          bind_cols(res)
      }
    }else{
      # error statement
      gen_res_df <- data.frame(curve_mean = NA, curve_max = NA, curve_q0.1 = NA, curve_q0.5 = NA,
                               curve_q0.9 = NA, gen = NA, maxgen = NA, mu= NA, sigma2 = NA, shape = NA,
                               w = NA, nu = NA, aic = NA, bic = NA, edc = NA, icl = NA, model = NA, mixmod_flag = 1)
    }
    alldf[[m]] <- gen_res_df
  }
  
  # make data.frame for output values
  outdf <- bind_rows(alldf)
  return(outdf)
}

# from http://coleoguy.blogspot.com/2016/04/stochasticprobabilistic-rounding.html
StochRound = function(x){
  ## extract the decimal portion
  q <-  abs(x - trunc(x))
  
  ## draw a value 0 or 1 with probability
  ## based on how close we already are
  adj <- purrr::map_int(q, ~base::sample(0:1, size = length(.x), prob = c(1 - .x, .x)))
  
  ## make it negative if x is
  adj <- ifelse(x < 0, adj * -1, adj)
  
  ## return our new value
  trunc(x) + adj
}



# intending with region grouping
AssignGeneration <- function(mixmod, dat, param, reg, yr){
  
  set.seed(param$seed) # use for all rows?
  dat <- dat %>%
    filter(region == reg,
           modyear == yr)
  y <- StochRound(dat$adjY)
  
  # df to be combined with generation classifications at the end
  dat$row <- 1:nrow(dat)
  outdf <- dat[, c("SiteID", "SiteDate", "Timescale")]
  outdf <- outdf[rep(dat$row, y), ]
  
  outlist <- list()
  for (ngen in 1:(param$ngen)){
    if(ngen == 1){
      outclass <- outdf %>% 
        mutate(gen = 1,
               count = 1,
               year = lubridate::year(SiteDate)) %>% 
        filter(count > 0) %>% 
        mutate(region = reg, 
               index = param$index)
      
    }else{
      
      allmods <- flatten(mixmod)
      # only use model with correct ngen
      mod <- allmods[[ngen]]
      if(is.null(mod$error) && !is.null(mod$result)){
        modtype <- attr(mod$result, "class")
        if(modtype == "Mclust"){
          # assign generations
          res <- as.data.frame(mod$result$z)
          quickclass <- function(x){rmultinom(n = 1, size = 1, prob = x)}
          resclass <- apply(X = res, MARGIN = 1, FUN = quickclass)
          resclass <- data.frame(t(resclass))
          names(resclass) <- paste("gen", rank(mod$result$parameters$mean), sep = "_")
          outclass <- cbind(outdf, resclass) %>% 
            tidyr::gather(key = gen, value = count, contains("gen")) %>% 
            mutate(gen = as.numeric(stringr::str_split_fixed(gen, pattern = "_", n = 2)[, 2]),
                   year = lubridate::year(SiteDate)) %>% 
            filter(count > 0) %>% 
            mutate(region = reg, 
                   modyear = yr,
                   index = param$index)
          
        }else if(modtype == "Skew.normal"){
          # mixsmsn
          res <- as.data.frame(mod$result$obs.prob)
          quickclass <- function(x){rmultinom(n = 1, size = 1, prob = x)}
          resclass <- apply(X = res, MARGIN = 1, FUN = quickclass)
          resclass <- data.frame(t(resclass))
          names(resclass) <- paste("gen", rank(mod$result$mu), sep = "_")
          outclass <- cbind(outdf, resclass) %>% 
            tidyr::gather(key = gen, value = count, contains("gen")) %>% 
            mutate(gen = as.numeric(stringr::str_split_fixed(gen, pattern = "_", n = 2)[, 2]),
                   year = lubridate::year(SiteDate)) %>% 
            filter(count > 0) %>% 
            mutate(region = reg, 
                   modyear = yr,
                   index = param$index)
        }else{
          # error statement
          outclass <- data.frame(SiteID = "999", gen = ngen,
                                 region = reg, modyear = yr, index = param$index)
        }
      }else{
        # error statement
        outclass <- data.frame(SiteID = "999", gen = ngen,
                               region = reg, modyear = yr, index = param$index)
      }
    }
    outlist[[ngen]] <- dplyr::full_join(dat, outclass) %>% mutate(maxgen = ngen)
  }
  allout <- bind_rows(outlist)
  return(allout)
}