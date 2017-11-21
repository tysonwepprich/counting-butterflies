# Functions for population estimation

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
  cdf <- cumsum(y) / sum(y)
  cmean <- weighted.mean(t, y)
  cmax <- t[which(y == max(y))]
  
  cquant <- rep(NA, length(quants))
  for (i in seq_along(quants)){
  cquant[i] <- t[which(abs(cdf - quants[i]) == min(abs(cdf - quants[i])))]
  }
  
  df <- matrix(data = c(cmean, cmax, cquant), nrow = 1)
  df <- data.frame(df)
  names(df) <- c("curve_mean", "curve_max", paste0("curve_q", quants))
  return(df)
}

# data argument is row of parameters
# this assumes all sites/years have common physiology
Simulate_Truth <- function(data){
  set.seed(111) # use for all rows?
  
  # voltinism varying by site
  if(data$ngen == 1){
    gen_ddreq <- 800
    gen_relsize = 1
    gen_relsprd = 1
  }else{
    ngen <- data$ngen
    gen_ddreq <- seq(600, 1600, length.out = ngen)
    gen_relsize <- switch(data$gen_size, 
                          "equal" = rep(1, ngen),
                          "inc" = seq(1, 2, length.out = ngen),
                          "dec" = seq(2, 1, length.out = ngen))
    gen_relsprd <- seq(1, 2, length.out = ngen)
  }
  # simulate phenology for each generation and combine
  dflist <- list()
  for (g in 1:data$ngen){
    df <- Abund_Curve(t = seq(0, 3000, 1)/100, alpha = data$death_rate, 
                      beta = gen_relsprd[g] * data$peak_sd / 100, mu = gen_ddreq[g] / 100, sig = 1)
    df$Gen <- g
    dflist[[g]] <- df
  }
  dfall <- bind_rows(dflist)
  return(dfall)
}


# data argument is row of parameters
Simulate_Counts <- function(data, gdd){
  # set.seed(111) # use for all rows?
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
                          "inc" = seq(1, 2, length.out = ngen),
                          "dec" = seq(2, 1, length.out = ngen))
    gen_relsprd <- seq(1, 2, length.out = ngen)
  }
  # simulate phenology for each generation and combine
  dflist <- list()
  for(site in sites){
    site_re <- rnorm(1, mean = 0, sd = data$site_mu_sd)
    
    for(yr in years){
      year_re <- rnorm(1, mean = 0, sd = data$year_mu_sd)
      
      for (g in 1:data$ngen){
        simdat <- dddat %>% filter(SiteID == site, Year == yr)
        simt <-  simdat$AccumDD / 100
        simalpha <- data$death_rate
        simbeta <- gen_relsprd[g] * data$peak_sd / 100
        simmu <- (gen_ddreq[g] + site_re + year_re) / 100
        df <- Abund_Curve(t = simt, alpha = simalpha, 
                          beta = simbeta, mu = simmu, sig = 1)
        df <- cbind(data.frame(simdat), df)
        df$Gen <- g
        df$Site_RE <- site_re
        df$Year_RE <- year_re
        dflist[[length(dflist)+1]] <- df
      }
    }
  }
  dfall <- bind_rows(dflist) %>% 
    arrange(SiteID, SiteDate)
  
  # counting process
  counts <- dfall %>% 
    group_by(SiteID, SiteDate, lat, lon, Year, DOY, AccumDD, maxT, Site_RE, Year_RE) %>% 
    summarise(RelProp = sum(y * gen_relsize[Gen]),
              DP = plogis(data$detprob_b0 + data$detprob_b1 * maxT[1] + data$detprob_b2 * maxT[1]^2)) %>% 
    group_by(SiteID, Year) %>% 
    mutate(RelProp = RelProp / sum(RelProp),
           M = rnbinom(1, mu = data$negbin_mu, size = data$negbin_disp),
           N = rpois(length(RelProp), lambda = RelProp * M),
           Y = rbinom(length(N), size = N, prob = DP)) %>% 
    data.frame()
  return(counts)
}


# accounts for det prob, GAM smoothing
Adjust_Counts <- function(data, counts){
  # for GAM, need factors
  counts$SiteYearID <- as.factor(paste(counts$SiteID, counts$Year, sep = "_"))
  counts$SiteID <- as.factor(counts$SiteID)
  counts$Year <- as.factor(as.character(counts$Year))
  # counts were simulated each day, reduce to every 4 days
  counts_4day <- counts %>% 
    filter(DOY %in% (seq(90, 294, 4) + sample.int(n=3, size=52, replace=TRUE)))
  
  counts_8day <- counts_4day[seq(1, nrow(counts_4day), 2), ]

  adjcounts <- counts_8day %>% 
    sample_frac(size = 1 - data$surv_missing)
  
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
                           # s(zlistlength)+
                           # s(maxT)+
                           # s(zduration)+
                           # s(SiteYearID, AccumDD, bs = "fs", k = 5, m = 1) +
                           ti(SiteID, AccumDD, bs = c("re", "cc"), k = c(4, 10)) +
                           ti(Year, AccumDD, bs = c("re", "cc"), k = c(4, 10)) +
                           s(AccumDD, bs = "cc", k = 20) +
                           # te(lat, lon, AccumDD, bs = c("tp", "cc"), k = c(4, 30), d = c(2, 1)) +
                           s(Year, bs = "re") +
                           s(SiteID, bs = "re"),
                         # s(Ordinal, bs = "cc", k = 10),
                         family = nb(theta = NULL, link = "log"),
                         # family = poisson(link = "log"),
                         data = adjcounts,
                         method = "REML", 
                         optimizer = c("outer", "newton"), 
                         # gamma = 1.4, 
                         control = list(maxit = 500))
    }
    if(data$detprob_model == "covariate"){
      gammod <- safe_gam(adjY ~ 
                           # s(zlistlength)+
                           s(maxT)+
                           # s(zduration)+
                           # s(SiteYearID, AccumDD, bs = "fs", k = 5, m = 1) +
                           ti(SiteID, AccumDD, bs = c("re", "cc"), k = c(4, 10)) +
                           ti(Year, AccumDD, bs = c("re", "cc"), k = c(4, 10)) +
                           s(AccumDD, bs = "cc", k = 20) +
                           # te(lat, lon, AccumDD, bs = c("tp", "cc"), k = c(5, 20), d = c(2, 1)) +
                           s(Year, bs = "re") +
                           s(SiteID, bs = "re"),
                         # s(Ordinal, bs = "cc", k = 10),
                         family = nb(theta = NULL, link = "log"),
                         # family = poisson(link = "log"),
                         data = adjcounts,
                         method = "REML", 
                         optimizer = c("outer", "newton"), 
                         # gamma = 1.4, 
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
    newdata$adjY <- predict(gammod$result, newdata = newdata, type = "response")
    
    
    if(data$gam_smooth == "interpolate"){
      interps <- anti_join(newdata, adjcounts, by = c("SiteID", "SiteDate"))
      adjcounts <- bind_rows(adjcounts, interps)
      adjcounts$gam_flag <- 0
      return(adjcounts)
    }
    
    if(data$gam_smooth %in% c("preds_4day", "preds_8day")){
      adjcounts <- newdata
      adjcounts$gam_flag <- 0
      return(adjcounts)
    }
  }else{
    # error in gam fit
    adjcounts$gam_flag <- 1
    return(adjcounts)  
  }
}


# fit different mixture models, output a list of model results
# 
CompareMixMods <- function(dat, mvmax){
  dd <- dat$AccumDD
  y <- round(dat$adjY)
  dd_dist <- rep(dd, y)
  mvmin <- 1
  mvmax <- mvmax
  gens <- c(mvmin:mvmax)
  maxtry <- 10 # repeating smsn.mix function if errors
  # out <- as.list(mvmin:mvmax)
  # names(out) <- paste0("gen", c(mvmin:mvmax))
  out <- list()
  for (i in seq_along(gens)){
    safe_mix <- safely(smsn.mix)
    # skew normal
    worked <- 0
    tries <- 1
    while(worked == 0 & tries <= maxtry){
      mod_sknorm <- safe_mix(dd_dist, nu = 5, g = i, get.init = TRUE, family = "Skew.normal", 
                             calc.im = FALSE, obs.prob = FALSE, kmeans.param = list(n.start = 5))
      if(is.null(mod_sknorm$error)){
        worked <- 1
      }else{
        tries <- tries + 1
      }
    }
    
    
    # Normal het from smsn
    worked <- 0
    tries <- 1
    while(worked == 0 & tries <= maxtry){
      mod_norm <- safe_mix(dd_dist, nu = 5, g = i, get.init = TRUE, family = "Normal", 
                           calc.im = FALSE, obs.prob = FALSE, kmeans.param = list(n.start = 5))
      if(is.null(mod_norm$error)){
        worked <- 1
      }else{
        tries <- tries + 1
      }
    }
    
    
    # T skew
    worked <- 0
    tries <- 1
    while(worked == 0 & tries <= maxtry){
      mod_skt <- safe_mix(dd_dist, nu = 5, g = i, get.init = TRUE, family = "Skew.t", 
                          calc.im = FALSE, obs.prob = FALSE, kmeans.param = list(n.start = 5))
      if(is.null(mod_skt$error)){
        worked <- 1
      }else{
        tries <- tries + 1
      }
    }
   
    # T het from smsn
    worked <- 0
    tries <- 1
    while(worked == 0 & tries <= maxtry){
      mod_t <- safe_mix(dd_dist, nu = 5, g = i, get.init = TRUE, family = "t", 
                        calc.im = FALSE, obs.prob = FALSE, kmeans.param = list(n.start = 5))
      if(is.null(mod_t$error)){
        worked <- 1
      }else{
        tries <- tries + 1
      }
    }
    # mclust
    safe_mclust <- safely(Mclust)
    # hom
    worked <- 0
    tries <- 1
    while(worked == 0 & tries <= maxtry){
      mod_mc_hom <- safe_mclust(dd_dist, G=i, modelNames = "E")
      if(is.null(mod_mc_hom$error) & !is.null(mod_mc_hom$result)){
        worked <- 1
      }else{
        tries <- tries + 1
      }
    }
    if(is.null(mod_mc_hom$result) & is.null(mod_mc_hom$error)){
      mod_mc_hom$error <- "none fitted"
    }
    # het
    worked <- 0
    tries <- 1
    while(worked == 0 & tries <= maxtry){
      mod_mc_het <- safe_mclust(dd_dist, G=i, modelNames = "V")
      if(is.null(mod_mc_het$error) & !is.null(mod_mc_het$result)){
        worked <- 1
      }else{
        tries <- tries + 1
      }
    }
    if(is.null(mod_mc_het$result) & is.null(mod_mc_het$error)){
      mod_mc_het$error <- "none fitted"
    }
    
    out[[i]] <- list(mod_sknorm, mod_norm, mod_skt, mod_t, mod_mc_hom, mod_mc_het)
  }
  return(out)
}



# Extract results from mixture model objects
Summ_mixmod <- function(df){
  allmods <- flatten(flatten(df$mixmods))
  alldf <- list()
  for (m in seq_along(allmods)){
    mod <- allmods[[m]]
    if(is.null(mod$error) && !is.null(mod$result)){
      modtype <- attr(mod$result, "class")
      if(modtype == "Mclust"){
        # mclust param
        modtype <- paste(modtype, mod$result$modelName, sep = "_")
        if(modtype == "Mclust_X") modtype <- "Mclust_E"
        res <- data.frame(mu = mod$result$parameters$mean, sigma2 = mod$result$parameters$variance$sigmasq, shape = NA,
                          w = mod$result$parameters$pro, nu = NA, aic = NA, bic = mod$result$bic,
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
      }
    }else{
      # error statement
      gen_res_df <- data.frame(curve_mean = NA, curve_max = NA, curve_q0.1 = NA, curve_q0.5 = NA,
                               curve_q0.9 = NA, gen = NA, maxgen = NA, mu= NA, sigma2 = NA, shape = NA,
                               w = NA, nu = NA, aic = NA, bic = NA, edc = NA, icl = NA, model = NA)
    }
    alldf[[m]] <- gen_res_df
  }
  
  # make data.frame for output values
  outdf <- bind_rows(alldf)
  return(outdf)
}


