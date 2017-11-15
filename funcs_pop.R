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
  set.seed(111) # use for all rows?
  # randomly select sites and years to use their coordinates and historical gdd accumulation
  sites <- sample(x = unique(gdd$SiteID), size = data$nsite, replace = FALSE)
  years <- sample(x = unique(gdd$Year), size = data$nyear, replace = FALSE)
  dddat <- gdd %>% 
    filter(SiteID %in% sites, Year %in% years) %>% 
    filter(DOY %in% (seq(90, 293, 7) + sample.int(n=6, size=30, replace=TRUE))) %>% 
    sample_frac(size = 1 - data$surv_missing)
  
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
  for (g in 1:data$ngen){
    df <- Abund_Curve(t = dddat$AccumDD/100, alpha = data$death_rate, 
                      beta = gen_relsprd[g] * data$peak_sd / 100, mu = gen_ddreq[g] / 100, sig = 1)
    df <- cbind(data.frame(dddat), df)
    df$Gen <- g
    dflist[[g]] <- df
  }
  dfall <- bind_rows(dflist) %>% 
    arrange(SiteID, SiteDate)
  
  # counting process
  counts <- dfall %>% 
    group_by(SiteID, SiteDate, lat, lon, Year, DOY, AccumDD, maxT) %>% 
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





# fit different mixture models, output a list of model results
# 
CompareMixMods <- function(dat){
  dd <- dat$AccumDD
  y <- dat$Y
  dd_dist <- rep(dd, y)
  mvmin <- dat$mvmin[1]
  mvmax <- dat$mvmax[1]
  gens <- c(mvmin:mvmax)
  # out <- as.list(mvmin:mvmax)
  # names(out) <- paste0("gen", c(mvmin:mvmax))
  out <- list()
  for (i in seq_along(gens)){
    
    safe_mix <- safely(smsn.mix)
    # skew normal
    mod_sknorm <- safe_mix(dd_dist, nu = 5, g = i, get.init = TRUE, family = "Skew.normal", 
                                calc.im = TRUE, obs.prob = TRUE, kmeans.param = list(n.start = 5))
    
    # Normal het from smsn
    mod_norm <- safe_mix(dd_dist, nu = 5, g = i, get.init = TRUE, family = "Normal", 
                              calc.im = TRUE, obs.prob = TRUE, kmeans.param = list(n.start = 5))
    
    
    
    # T DIST
    # get initial vals from search for best num gen
    mod_skt <- safe_mix(dd_dist, nu = 5, g = i, get.init = TRUE, family = "Skew.t", 
                             calc.im = TRUE, obs.prob = TRUE, kmeans.param = list(n.start = 5))
    
    # T het from smsn
    mod_t <- safe_mix(dd_dist, nu = 5, g = i, get.init = TRUE, family = "t", 
                           calc.im = TRUE, obs.prob = TRUE, kmeans.param = list(n.start = 5))
    
    # mclust
    safe_mclust <- safely(Mclust)
    mod_mc_hom <- safe_mclust(dd_dist, G=i, modelNames = "E")
    mod_mc_het <- safe_mclust(dd_dist, G=i, modelNames = "V")
    out[[i]] <- list(mod_sknorm, mod_norm, mod_skt, mod_t, mod_mc_hom, mod_mc_het)
  }
  return(out)
}


# take mixture model probability for each observation's class,
# then assign generations to each count (possibly splitting them)
Assign_generation <- function(counts, results){
  yrs <- unique(results$Year)
  for (y in seq_along(yrs)){
    cts <- counts %>% filter(Year == y)
    res <- results %>% filter(Year == y)
    allmods <- flatten(flatten(res$mixmods))
    alldf <- list()
    for (m in seq_along(allmods)){
      mod <- allmods[[m]]
      if(is.null(mod$error)){
  }
  
  if(is.null(mod$error)){
    modtype <- attr(mod$result, "class")
    if(modtype == "Mclust"){
      # mclust param
      
      # simulate each generation
      
      # extract curve summary
    }else{
      # mixsmsn
      
      # simulate each generation
      
      # extract curve summary
    }
  }else{
    # error statement
  }
  
}

# Extract results from mixture model objects
Summ_mixmod <- function(df){
  allmods <- flatten(flatten(df$mixmods))
  alldf <- list()
  for (m in seq_along(allmods)){
    mod <- allmods[[m]]
    if(is.null(mod$error)){
      modtype <- attr(mod$result, "class")
      if(modtype == "Mclust"){
        # mclust param
        modtype <- paste(modtype, mod$result$modelName, sep = "_")
        res <- data.frame(mu = mod$result$parameters$mean, sigma2 = mod$result$parameters$variance$sigmasq, shape = NA,
                          w = mod$result$parameters$pro, nu = NA, aic = NA, bic = mod$result$bic,
                          edc = NA, icl = NA, model = modtype)
        
        # simulate each generation
        # extract curve summary
        gen_res <- list()
        for (i in seq_along(as.numeric(rownames(res)))){
          x <- seq(0, 3000, 10)
          y <- dnorm(x = x, mean = res$mu[i], sd = sqrt(res$sigma2[i]))
          gen_res[[i]] <- Summ_curve(t = x, y = y)
        }
        gen_res_df <- bind_rows(gen_res) %>% 
          mutate(gen = as.numeric(rownames(res))) %>% 
          bind_cols(res)
        
      }else{
        # mixsmsn
        
        res <- data.frame(mu = mod$result$mu, sigma2 = mod$result$sigma2, shape = mod$result$shape,
                          w = mod$result$pii, nu = mod$result$nu, aic = mod$result$aic, bic = mod$result$bic,
                          edc = mod$result$edc, icl = mod$result$icl, model = modtype)
        
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
          gen_res[[i]] <- Summ_curve(t = x, y = y)
        }
        gen_res_df <- bind_rows(gen_res) %>% 
          mutate(gen = as.numeric(rownames(res))) %>% 
          bind_cols(res)
      }
    }else{
      # error statement
      gen_res_df <- data.frame(curve_mean = NA, curve_max = NA, curve_q0.1 = NA, curve_q0.5 = NA,
                               curve_q0.9 = NA, gen = NA, mu= NA, sigma2 = NA, shape = NA,
                               w = NA, nu = NA, aic = NA, bic = NA, edc = NA, icl = NA, model = NA)
    }
    alldf[[m]] <- gen_res_df
  }
  
  # make data.frame for output values
  outdf <- bind_rows(alldf)
  return(outdf)
}


