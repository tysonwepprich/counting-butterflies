test1 <- system.time({moda <- gam(Total ~ 
                  s(zlistlength) +
                  s(ztemperature) +
                  te(lat, lon, AccumDD, bs = c("tp", "cr"), k = c(5, 30), d = c(2, 1)) +
                  s(SiteID, bs = "re") +
                  s(Year, DOY, bs = "fs", k = 5, m = 1) +
                  s(RegYear, AccumDD, bs = "fs", k = 5, m = 1),
                # family = nb(theta = NULL, link = "log"),
                family = poisson(link = "log"),
                data = temp,
                method = "REML", 
                optimizer = c("outer", "newton"), 
                # gamma = 1.4, 
                control = list(maxit = 500))
})


test2 <- system.time({mod2 <- gam(Total ~ 
                                    s(zlistlength) +
                                    s(ztemperature) +
                                    Year +
                                    te(lat, lon, AccumDD, bs = "cr", k = c(5, 30), d = c(2, 1), by = Year) +
                                    s(SiteID, bs = "re"),
                                  family = nb(theta = NULL, link = "log"),
                                  # family = poisson(link = "log"),
                                  data = temp,
                                  method = "REML", 
                                  optimizer = c("outer", "newton"), 
                                  # gamma = 1.4, 
                                  control = list(maxit = 500))
})

test3 <- system.time({mod3 <- gam(Total ~ 
                                    s(zlistlength) +
                                    s(ztemperature) +
                                    Year +
                                    s(AccumDD, bs = "cr", k = 30, by = Year) +
                                    s(SiteID, bs = "re"),
                                  family = nb(theta = NULL, link = "log"),
                                  # family = poisson(link = "log"),
                                  data = temp,
                                  method = "REML", 
                                  optimizer = c("outer", "newton"), 
                                  # gamma = 1.4, 
                                  control = list(maxit = 500))
})

test4 <- system.time({mod4 <- gam(Total ~ 
                                    s(zlistlength) +
                                    s(ztemperature) +
                                    s(RegYear, AccumDD, bs = "fs", k = 5, m = 1) +
                                    s(AccumDD, bs = "cr", k = 30) +
                                    s(SiteID, bs = "re"),
                                  family = nb(theta = NULL, link = "log"),
                                  # family = poisson(link = "log"),
                                  data = temp,
                                  method = "REML", 
                                  optimizer = c("outer", "newton"), 
                                  # gamma = 1.4, 
                                  control = list(maxit = 500))
})


test5 <- system.time({mod5 <- gam(Total ~ 
                                    s(zlistlength) +
                                    s(ztemperature) +
                                    s(RegYear, AccumDD, bs = "fs", k = 10, m = 1) +
                                    s(AccumDD, bs = "cr", k = 30) +
                                    s(SiteID, bs = "re"),
                                  family = nb(theta = NULL, link = "log"),
                                  # family = poisson(link = "log"),
                                  data = temp,
                                  method = "REML", 
                                  optimizer = c("outer", "newton"), 
                                  # gamma = 1.4, 
                                  control = list(maxit = 500))
})

test7 <- system.time({mod7 <- gam(Total ~ 
                                    s(zlistlength) +
                                    s(ztemperature) +
                                    s(RegYear, AccumDD, bs = "fs", k = 10, m = 1) +
                                    s(AccumDD, bs = "cr", k = 30) +
                                    s(SiteYear, bs = "re"),
                                  family = nb(theta = NULL, link = "log"),
                                  # family = poisson(link = "log"),
                                  data = temp,
                                  method = "REML", 
                                  optimizer = c("outer", "newton"), 
                                  # gamma = 1.4, 
                                  control = list(maxit = 500))
})

test8 <- system.time({mod8 <- gamm(Total ~ 
                                    s(zlistlength) +
                                    s(ztemperature) +
                                    s(RegYear, AccumDD, bs = "fs", k = 5, m = 1) +
                                    s(AccumDD, bs = "cr", k = 30),
                                  family = nb(theta = NULL, link = "log"),
                                  random=list(SiteYear=~1),
                                  # family = poisson(link = "log"),
                                  data = temp,
                                  method = "REML", 
                                  # optimizer = c("outer", "newton"), 
                                  # gamma = 1.4, 
                                  control = list(maxit = 500))
})


# if variables correlated, should fix smooth parameters if evaluating reduced model
# for deviance explained by one variable compared to full

## simulate some data
set.seed(0)
n<-400
x1 <- runif(n, 0, 1)
## to see problem with not fixing smoothing parameters
## remove the `##' from the next line, and the `sp'
## arguments from the `gam' calls generating b1 and b2.
x2 <- runif(n, 0, 1) ## *.1 + x1
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
f <- f1(x1) + f2(x2)
e <- rnorm(n, 0, 2)
y <- f + e
## fit full and reduced models...
b <- gam(y~s(x1)+s(x2))
b1 <- gam(y~s(x1),sp=b$sp[1])
b2 <- gam(y~s(x2),sp=b$sp[2])
b0 <- gam(y~1)
## calculate proportions deviance explained...
(deviance(b1)-deviance(b))/deviance(b0) ## prop explained by s(x2)
(deviance(b2)-deviance(b))/deviance(b0) ## prop explained by s(x1)
