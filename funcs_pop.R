# Functions for population estimation

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
Summ_curve <- function(t, y, quants = c(.1, .5, .9)){
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
Simulate_Counts <- function(data, gdd){
  set.seed(data$seed) # use for all rows?
  # randomly select sites and years to use their coordinates and historical gdd accumulation
  sites <- sample(x = unique(gdd$SiteID), size = data$nsite, replace = FALSE)
  years <- sample(x = unique(gdd$Year), size = data$nyear, replace = FALSE)
  dddat <- gdd %>% 
    filter(SiteID %in% sites, Year %in% years) # %>% 
    # filter(DOY %in% (seq(90, 293, 7) + sample.int(n=6, size=30, replace=TRUE)))
    # filter(DOY %in% (seq(90, 293, 3)))
  
  # voltinism varying by site
  if(data$ngen == 1){
    gen_ddreq <- 800
    gen_relsize = 1
    gen_relsprd = 1
  }else{
    ngen <- data$ngen
    gen_ddreq <- seq(600, 1600, length.out = ngen)
    gen_relsize <- switch(as.character(data$gen_size), 
                          "equal" = rep(1, ngen),
                          "inc" = seq(1, ngen, length.out = ngen),
                          "dec" = seq(ngen, 1, length.out = ngen))
    gen_relsprd <- seq(1, 1.5, length.out = ngen)
  }
  # simulate phenology for each generation and combine
  year_var <- rnorm(data$nyear, mean = 0, sd = data$year_mu_sd)
  site_var <- rnorm(data$nsite, mean = 0, sd = data$site_mu_sd)
  
  dflist <- list()
  for(y in seq_along(years)){
    for(z in seq_along(sites)){
      yr <- years[y]
      site <- sites[z]
      # site re includes dependence on latitude, lat of 40 assumed to be baseline
      site_lat <- gdd$lat[gdd$SiteID == site][1]
      site_re <- 1000 - site_lat * 25 + site_var[z]
      year_re <- year_var[y]
      
      for (g in 1:data$ngen){
        simdat <- dddat %>% filter(SiteID == site, Year == yr)
        simt <-  simdat$AccumDD
        simalpha <- data$death_rate
        simbeta <- gen_relsprd[g] * data$peak_sd
        simmu <- (gen_ddreq[g] + site_re + year_re)
        df <- Abund_Curve(t = simt, alpha = simalpha, 
                          beta = simbeta, mu = simmu, sig = 1)
        df <- cbind(data.frame(simdat), df)
        df$Gen <- g
        df$Site_RE <- site_re
        df$Year_RE <- year_re
        df$gen_weight <- gen_relsize[g] * ifelse(g == data$ngen, ((42 - df$lat) / 3.35), 1)
        dflist[[length(dflist)+1]] <- df
      }
    }
  }
  dfall <- bind_rows(dflist) %>% 
    arrange(SiteID, SiteDate)
  
  # counting process
  counts <- dfall %>% 
    group_by(SiteID, SiteDate, lat, lon, Year, DOY, AccumDD, region, maxT, minT, Site_RE, Year_RE) %>% 
    summarise(RelProp = sum(y * gen_weight),
              DP = plogis(data$detprob_b0 + data$detprob_b1 * maxT[1] + data$detprob_b2 * maxT[1]^2)) %>% 
    group_by(SiteID, Year) %>% 
    mutate(RelProp = RelProp / sum(RelProp),
           # M = rnbinom(1, mu = data$negbin_mu, size = data$negbin_disp),
           M = rpois(1, lambda = data$pois_lam),
           N = rpois(length(RelProp), lambda = RelProp * M),
           Y = rbinom(length(N), size = N, prob = DP)) %>% 
    data.frame()
  
  # truth
  true_weight <- dfall %>% 
    group_by(SiteID, SiteDate) %>% 
    mutate(gen_weight = gen_weight / sum(gen_weight))
  true_phen <- dfall %>% 
    ungroup() %>% 
    group_by(SiteID, Year, Gen) %>% 
    do(Summ_curve(t = .$AccumDD, y = .$y))
  true_N <- counts %>% 
    ungroup() %>% 
    dplyr::select(SiteID, Year, Site_RE, Year_RE, M) %>% 
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
Adjust_Counts <- function(data, counts){
  set.seed(data$seed)
  # for GAM, need factors
  counts$SiteYearID <- as.factor(paste(counts$SiteID, counts$Year, sep = "_"))
  counts$SiteID <- as.factor(counts$SiteID)
  counts$region <- as.factor(counts$region)
  counts$Year <- as.factor(as.character(counts$Year))
  if(data$gam_scale == "DOY"){
    counts$Timescale <- counts$DOY
  }
  if(data$gam_scale == "GDD"){
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
    sample_frac(size = 1 - data$surv_missing) %>% 
    droplevels()
  
  if(data$detprob_model == "known"){
    adjcounts <- mutate(adjcounts, adjY = Y / DP)
  }else{
    adjcounts <- mutate(adjcounts, adjY = Y)
  }
  
  if(data$gam_smooth == "none"){
    adjcounts$gam_flag <- 0
    return(adjcounts)  
  }else{
    # model GAM once for next two cases
    
    safe_gam <- purrr::safely(gam)
    
    if(data$detprob_model %in% c("none", "known")){
      gammod <- safe_gam(adjY ~ 
                            # s(AccumDD, bs = "cc", k = 30) +
                           te(lat, lon, Timescale, bs = c("tp", "cc"), k = c(5, 30), d = c(2, 1)) +
                            s(SiteYearID, bs = "re"),
                         family = nb(theta = NULL, link = "log"),
                         # family = poisson(link = "log"),
                         data = adjcounts,
                         method = "REML", 
                         optimizer = c("outer", "newton"), 
                         control = list(maxit = 500))
    }
    if(data$detprob_model == "covariate"){
      gammod <- safe_gam(adjY ~ 
                           s(maxT) +
                           te(lat, lon, Timescale, bs = c("tp", "cc"), k = c(5, 30), d = c(2, 1)) +
                           s(SiteYearID, bs = "re"),
                         family = nb(theta = NULL, link = "log"),
                         # family = poisson(link = "log"),
                         data = adjcounts,
                         method = "REML", 
                         optimizer = c("outer", "newton"), 
                         control = list(maxit = 500))
    }
  }
  
  if(is.null(gammod$error)){
    if(data$gam_smooth == "preds_4day"){
      newdata <- counts_4day
    }else{
      newdata <- counts_8day
    }
    
    # HOW TO DEAL WITH TEMP COVARIATE HERE?
    # If excluded, predicted counts are really low. Could scale and use mean like before
    # if(data$detprob_model == "covariate"){
    #   newdata$adjY <- predict(gammod$result, newdata = newdata, type = "response", exclude = "s(maxT)")
    # }else{
    #   newdata$adjY <- predict(gammod$result, newdata = newdata, type = "response")
    # }
    
    mod <- gammod$result
    summod <- summary(gammod$result)

    # # prediction with response, gives non-integers which makes mixmod harder
    # adjY <- predict(gammod$result, newdata = newdata, type = "response")

    # prediction with simulated counts, stochastic but integers (if n is odd)
    Xp <- predict.gam(object = mod, newdata = newdata, type="lpmatrix") ## map coefs to fitted curves
    beta <- coef(mod)
    Vb   <- vcov(mod) ## posterior mean and cov of coefs
    n <- 5 # choose number of simulations
    mrand <- MASS::mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    ilink <- family(mod)$linkinv
    linklppreds <- Xp %*% t(mrand)
    nbpreds <- apply(X = linklppreds,
                     MARGIN = 1,
                     FUN = function(x){
                       # temp <- sort(x)
                       # bounds <- quantile(1:n, probs = c(0.025, 0.975))
                       # x <- temp[bounds[1]:bounds[2]]
                       x <- ilink(x)
                       x <- rnbinom(n = length(x),
                                    mu = x,
                                    size = mod$family$getTheta(TRUE))
                       x <- quantile(x, .5)
                       return(x)
                     })
    newdata$adjY <- nbpreds
    
    if(data$gam_smooth == "interpolate"){
      interps <- anti_join(newdata, adjcounts, by = c("SiteID", "SiteDate"))
      adjcounts <- bind_rows(adjcounts, interps)
      adjcounts$gam_flag <- 0
      adjcounts$gam_devexpl <- summod$dev.expl
      adjcounts$nb_theta <- gammod$result$family$getTheta(TRUE)
      
      return(adjcounts)
    }
    
    if(data$gam_smooth %in% c("preds_4day", "preds_8day")){
      adjcounts <- newdata
      adjcounts$gam_flag <- 0
      adjcounts$gam_devexpl <- summod$dev.expl
      adjcounts$nb_theta <- gammod$result$family$getTheta(TRUE)
      return(adjcounts)
    }
  }else{
    # error in gam fit
    adjcounts$gam_flag <- 1
    adjcounts$nb_disp <- NA
    return(adjcounts)  
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
  mvmax <- param$ngen
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
  
RightNumGen <- function(mixmods, param, reg){
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
           region = reg) %>% 
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



# intending with region grouping
AssignGeneration <- function(mixmod, dat, param, reg){
  
  set.seed(param$seed) # use for all rows?
  dat <- dat %>%
    filter(region == reg)
  y <- round(dat$adjY)
  
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
                   index = param$index)
        }else{
          # error statement
          outclass <- data.frame(SiteID = "999", gen = ngen,
                                 region = reg, index = param$index)
        }
      }else{
        # error statement
        outclass <- data.frame(SiteID = "999", gen = ngen,
                               region = reg, index = param$index)
      }
    }
    outlist[[ngen]] <- dplyr::full_join(dat, outclass) %>% mutate(maxgen = ngen)
  }
  allout <- bind_rows(outlist)
  return(allout)
}