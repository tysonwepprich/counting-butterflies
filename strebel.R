##################################################################
# Worked example to run the model presented in Strebel et al., 2014
#   (Study of phenology by flexible estimation and modeling of seasonal detectability peaks)
# Author:  Nicolas Strebel, nicolas_strebel@gmx.ch
# Date:	4.2.2014
# Title:	Run model
##################################################################
#-----------------------------------------------------------------
# Codes prepare data and run the analysis
#-----------------------------------------------------------------
### Set working directory
#setwd(...)

# Specify model in BUGS language
sink("splinesSiteOcc S4.txt")
cat("
    model { 
    ### Define seasonal and annual patterns in detectability
    for (m in 1:nyear) {  
    for (i in 1:n) {
    logit(p[m,i]) <- lp[m,i]
    lp[m,i] <- mfe[m,i]+mre[m,i]
    mfe[m,i] <- a[m]*X[i,1]+b[m]*X[i,2]+c[m]*X[i,3]
    mre[m,i]<-sum(n.mre[m,i,1:nknots])
    for (k in 1:nknots) {
    n.mre[m,i,k]<-b.k[m,k]*Z[i,k]
    }
    }
    
    ### Random regression coefficients corresponding to the truncated polynomial functions
    for (k in 1:nknots) {
    b.k[m,k] ~ dnorm(0,taub)
    }
    
    ### Fixed regression coefficients corresponding to the 'plus' functions
    
    a[m] ~ dnorm(0,0.01)
    b[m] ~ dnorm(0,0.01)
    c[m] ~ dnorm(0,0.01)
    }
    
    ### precision for random regression coefficients corresponding to the truncated polynomial functions
    taub~dgamma(1.0E-6,1.0E-6)      
    
    # Specify priors
    for (k in 1:nyear) {
    psi[k] ~ dunif(0, 1)
    }
    
    # Ecological submodel: Define state conditional on parameters
    for (i in 1:nsite){
    for (k in 1:nyear){
    z[i,k] ~ dbern(psi[k])
    }
    }
    
    # Observation model
    for (i in 1:nobs){
    muy[site[i],survey[i],year[i]] <- z[site[i],year[i]]*p[year[i],survey[i]]
    y[i] ~ dbin(muy[site[i],survey[i],year[i]], nrep[i])
    }
    
    }
    ",fill = TRUE)
sink()

### Read observation data from Acrocephalus arundinaceus
dat<-read.table(file="dat S4.txt",header=T)

### The following procedure is based on the models presented in Crainiceanu et al. 2005 and in Gimenez et al. 2006 
# Degree of splines
degree <- 2

# covariate
covariate<-as.numeric(scale(range(dat$day)[1]:range(dat$day)[2]))

# covariate length
n <- length(covariate)

# location of knots
nk<-round((max(dat$day)-min(dat$day)+1)/4)
nknots<-ifelse(nk<35,nk,35)
knots<-quantile(unique(covariate),seq(0,1,length=(nknots+2))[-c(1,(nknots+2))])

# fixed effects matrix
X<-NULL
for (l in 0:degree) {
  X<-cbind(X,covariate^l)  
}

# random coefficients matrix
Z_K<-(abs(outer(covariate,knots,"-")))^3
OMEGA_all<-(abs(outer(knots,knots,"-")))^3
svd.OMEGA_all<-svd(OMEGA_all)
sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*% (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d)))
Z<-t(solve(sqrt.OMEGA_all,t(Z_K)))

# Input data
site <- dat$site
survey <- dat$day-min(dat$day)+1
nobs <- length(unique(paste(dat$site,dat$day,dat$year)))
nrep <- dat$nrep
nsite <- length(unique(dat$site))
nyear <- length(unique(dat$year))
year <- as.numeric(factor(dat$year))
zst <- array(1, dim=c(nsite,nyear))  
y <- dat$ndet

# Simulation parameters
ni=500; nc=2; nb=25; nt=1

# List input data
jags.data <- list("site","survey","nobs","nrep","nsite","nyear","year","nknots","n","X","Z","nc", "nb", "ni", "nt","zst","y")

# Inits function
f.inits <- function(){list(a=rep(0,nyear), b=rep(0,nyear), c=rep(0,nyear), z=zst)}

# specify the parameters to be monitored
parameters <- c("a","b","c","b.k","lp","psi","taub")

### Run MCMC Analysis using jags
library(R2jags)
jags.out<-jags.parallel(jags.data,f.inits,parameters,"splinesSiteOcc S4.txt",nc,ni,nb,nt)
jags.out<-jags(jags.data,f.inits,parameters,"splinesSiteOcc S4.txt",nc,ni,nb,nt)

out<-jags.out$BUGSoutput


# Save model output
#save(out,file="out S4")

#I usually run two chains over 50'000 iterations, this takes several hours on my PC (3.4GHz, 4GB RAM)
#Usually convergence is reached within the first 10'000; I set burnin to 25'000
#To plot model output run the following codes:

#-----------------------------------------------------------------
# Codes to summarize the output
#-----------------------------------------------------------------
### If you want to skip running the bugs function, then load the output here:
#load(file="out S4")

### get estimated date of peak detectability based on posterior distribution
# get date of peak detectability in each simulation
findmax.fn<-function(x) {
  mean(which(x==max(x)))
}
lpmax<-array(data=NA,dim=c(out$n.sims,nyear))
dimnames(lpmax)<-list(c(1:out$n.sims),c(sort(unique(dat$year))))
for (xj in sort(unique(as.numeric(factor(dat$year))))) { 
  lpmax[,xj]<-apply(out$sims.array[,,paste("lp[",xj[1],",",1:(max(dat$day)-min(dat$day)+1),"]",sep="")],MARGIN=c(if(out$n.chains>1) 1:2 else 1),findmax.fn)
}
lpmax<-lpmax+min(dat$day)-1
lpmax[lpmax==max(dat$day)]<-NA
lpmax[lpmax==min(dat$day)]<-NA

# summarize estimates
ann.res<-array(NA, dim=c(max(dat$year)-min(dat$year)+1,3),dimnames=list(c(min(dat$year):max(dat$year)),c("mean","2.5%","97.5%")))
res<-apply(lpmax,c(2),mean,na.rm=T)
ann.res[names(res),"mean"]<-res
res<-apply(lpmax,c(2),quantile,probs=0.025,na.rm=T)
ann.res[names(res),"2.5%"]<-res
res<-apply(lpmax,c(2),quantile,probs=0.975,na.rm=T)
ann.res[names(res),"97.5%"]<-res

# get estimate of trend in date of peak detectability over years
do.lm<-function(x) {
  lmres<-lm(x~as.numeric(names(x)))$coefficients
  return(lmres)
}
r<-matrix(NA,dim(lpmax)[1],2)
for (o in 1:(dim(lpmax)[1])) {
  if(!is.na(sum(lpmax[o,]))) {
    lm(lpmax[o,]~as.numeric(colnames(lpmax)))$coefficients->r[o,]
  }    
}
slopevec<-as.vector(r[,2])
intercept<-mean(r[,1],na.rm=T)
slope<-mean(r[,2],na.rm=T)

### Write results (in console if argument file is not specified in function cat)
cat(paste("summary results","Acrocephalus arundinaceus"),"\n",
    paste("annual change of activity peak:", round(mean(slopevec,na.rm=T),digits=2),"days"),
    paste("confidence interval from", round(quantile(slopevec,0.025,na.rm=T),digits=2),
          "to",round(quantile(slopevec,0.975,na.rm=T),digits=2)),
    "\n","mean estimate of activity peak","as date",
    as.character(as.Date(x=c(ann.res[,colnames(ann.res)=="mean"]),origin=c(paste(row.names(ann.res),"-01-01",sep="")))),"\n",
    sep="\n","as day of the year",
    paste(rownames(ann.res),round(ann.res[,"mean"])))   

#-----------------------------------------------------------------
# Plot output
#-----------------------------------------------------------------
# save plotted results as pdf
pdf(file=paste("Graphical summary S4.pdf"),width=6,height=4)

### plot estimates of peak detectability over all years
par(mfrow=c(1,1))
par(mai=c(1,1,1,0.5))
x=rownames(ann.res)
y=ann.res[,"mean"]
plot(x,y,xlab="",ylab="",axes=F,main=paste("Peak Detection Probability","\n","Acrocephalus arundinaceus"),
     ylim=c(min(ann.res),max(ann.res)),pch=16,type="p", col="black")
lines(x,ann.res[,"2.5%"],col="grey",lwd=2)
lines(x,ann.res[,"97.5%"],col="grey",lwd=2)  
axis(side=1,at=x)
axis(side=2,at=c(121,135,152,166),
     labels=c("1May","15May","1Jun","15Jun"))
abline(a=intercept,b=slope,lty=2,col=colors()[200])

### Plot annual detectability pattern
# loop over all years
years<-sort(unique(as.numeric((dat$year))))
for (xj in 1:nyear) {
  j<-years[xj]
  
  # Get BUGS estimates
  res.chains<-out$sims.array[,,paste("lp[",xj[1],",",1:(max(dat$day)-min(dat$day)+1),"]",sep="")]
  res=plogis(apply(res.chains,MARGIN=c(length(dim(res.chains))),quantile,probs=c(.025,.5,.975)))
  
  ### Plot "naive" estimate of detectability
  # prepare bars to compare barplot of observation data (bars) with estimates (line); barheight represents weekly proportion of detection events divided by all surveys
  barwidth<-7
  z<-1
  for(m in seq(from=min(dat$day),to=max(dat$day),by=barwidth)) {
    dat$barpos[dat$day >= m & dat$day < m+barwidth]<-z
    z<-z+1
  }
  barheight<-rep(NA,times=max(dat$day)+7)
  names(barheight)<-1:max(dat$day)
  
  # height of the bars equals to a seven day successful obs to all obs ratio
  n<-(max(dat$day)-min(dat$day)+1)
  res.height<-tapply(dat$ndet[dat$year==j],dat$barpos[dat$year==j],mean)
  barheight[as.numeric(names(res.height))*7-3+min(dat$day)]<-res.height
  
  # plot bars
  barplot(as.numeric(barheight[min(dat$day):max(dat$day)]),
          width=1,space=0,ylim=c(0,max(res[3,])),xlab="", ylab="Detection Probability", 
          main=paste("Acrocephalus arundinaceus",j),border=NA,axes=F)
  
  ### Plot model estimates  
  # plot seasonal estimates of detectability p
  lines(res[3,],lty=3,col=1,lwd=2.5) # lower bound of the 95% CI
  lines(res[2,],lty=1,col=1,lwd=2) # median
  lines(res[1,],lty=3,col=1,lwd=2.5) # upper bound of the 95% CI
  axis(2)
  axis(side=1,at=1-min(dat$day)+c(105,121,135,152,166,182),
       labels=c("15Apr","1May","15May","1Jun","15Jun","1Jul"))     
}

dev.off()