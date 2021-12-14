##################################################################
################### Master Thesis in Biostatistics ###############
#####' @author ZINZINHEDO Mahoukpégo Luc
#####' @field  Species Distribution Modeling
#####' @TOPIC: Evaluation of predictive performance of 
#####' inhomogeneous point process models trained with
#####' presence-only data in species distribution modeling
#####'
#####' @Objective The purpose of this study is to assess the
#####' predictive performance of some Point Process models
#####'
#####' @details In this part we're going to do 
#####' *extensive sampling on a well defined subregion* 
#####' *and see check whether point process models perform better*

#######################' *The beginning* ################# 
rm(list= ls(all=T))

##' Load packages 
library(spatstat) # For point pattern analysis and point process models
library(raster)   # To create virtual Rasters
library(maptools) # To translate rasters into .im object needed by "spatstat" package
library(fields)   # Helps to make plots
library(dplyr)    # Facilitate database subset
library(RandomFieldsUtils) # Needed by spatsat package
library(RandomFields)      # Needed by spatsat package to fit cox process model
library(writexl)           # Export simulated Dataset into Excel
library(nleqslv)           # Required by 'spatstat' to run cox model with 'adapcl' estimation method

D <- c(10, 10)  # Square Domaine D of side 300
Win <- owin(xrange =c(0, D[1]), yrange =c(0,D[2]), unnitname =c("km")) 
spatstat.options(npixel=c(D[1],D[2]))


ext <- extent(Win$xrange, Win$yrange) # Extent of the rasters
par(mfrow=c(1,1))

# First raster : Solar radiation sim
ras1 <- raster()
extent(ras1) <- ext
res(ras1) <- 1       # This is the resolution. Each pixel will be a 10*10 square
names(ras1) <- 'Solar radiation'
crs(ras1) <- "+proj=longlat +datum=WGS84 +no_defs"   # Coordinated reference system
vx1 <- sqrt(seq(219.87, 98,  length.out = 10))
vy1 <- sqrt(seq(98, 219.87, length.out = 10))
mv1 <- outer(vx1, vy1, FUN = "*")
values(ras1) <- mv1
ras1
plot(ras1, asp=1)

# Second Raster: Precipitation sim
ras2 <- raster()
extent(ras2) <- ext
res(ras2) <- 1
names(ras2) <- 'Precipitation'
crs(ras2) <- "+proj=longlat +datum=WGS84 +no_defs"
vx2 <- sqrt(seq(33, 244, length.out = 10))
vy2 <- sqrt(seq(244, 33, length.out = 10))
mv2 <- outer(vx2, vy2, FUN = "*")
values(ras2) <- mv2
ras2
plot(ras2, asp=1)

#Third raster ras3: Temperature sim
ras3 <- raster()
extent(ras3) <- ext
res(ras3) <- 1
names(ras3) <- 'Temperature'
crs(ras3) <- "+proj=longlat +datum=WGS84 +no_defs"
vx3 <- seq(20.41389, 25.06494, length.out = 10)
vy3 <- seq(20.41389, 27.06494, length.out = 10)
mv3 <- outer(vx3, vy3, function(x,y) ((1.04)*x+sin(y)))
values(ras3) <- t(mv3)
ras3
plot(ras3, asp=1)

#Fourth raster ras4: Wind speed sim
ras4 <- raster()
extent(ras4) <- ext
res(ras4) <- 1
names(ras4) <- 'wind speed'
crs(ras4) <- "+proj=longlat +datum=WGS84 +no_defs"
values(ras4) <- matrix(c(seq(from =0.58875, to =0.6, length.out=156), seq(from=0.6, to=0.7, length.out = 190), seq(from=0.7, to=0.8, length.out = 190), seq(from=0.8, to=0.97, length.out = 368), seq(from=0.97, to=1, length.out = 294), seq(from=1, to=1.66, length.out = 186),seq(from=1.66, to=2, length.out = 56),seq(from=2, to=2.78924, length.out = 160)), nrow = 10, ncol = 10)
ras4
plot(ras4, asp=1)

#Fife raster ras5: Water vapour pressure sim
ras5 <- raster()
extent(ras5) <- ext
res(ras5) <- 1
names(ras5) <- 'Water vapour pressure'
crs(ras5) <- "+proj=longlat +datum=WGS84 +no_defs"
vx5 <- sqrt(seq(2.32385, 3.0751, length.out = 10))
vy5 <- sqrt(seq(2.32385, 3.0751, length.out = 10))
mv5 <- outer(vx5, vy5, FUN = "*")
values(ras5) <- t(mv5)
ras5
plot(ras5, asp=1)

# Raster List (stack)
Rasters.group <- stack(ras1, ras2, ras3, ras4,  ras5)
plot(Rasters.group)
graphics.off()

#Translating ratsers in images
im.ras1 <- as.im.RasterLayer(ras1); unitname(im.ras1) <- c('km'); summary(im.ras1)
im.ras2 <- as.im.RasterLayer(ras2); unitname(im.ras2) <- c('km'); summary(im.ras2)
im.ras3 <- as.im.RasterLayer(ras3); unitname(im.ras3) <- c('km'); summary(im.ras3)
im.ras4 <- as.im.RasterLayer(ras4); unitname(im.ras4) <- c('km'); summary(im.ras4)
im.ras5 <- as.im.RasterLayer(ras5); unitname(im.ras5) <- c('km'); summary(im.ras5)

#standization of covariates
norm.im.ras1 <- (im.ras1- mean(im.ras1))/sd(im.ras1) ; summary(norm.im.ras1)
norm.im.ras2 <- (im.ras2- mean(im.ras2))/sd(im.ras2) ; summary(norm.im.ras2)
norm.im.ras3 <- (im.ras3- mean(im.ras3))/sd(im.ras3) ; summary(norm.im.ras3)
norm.im.ras4 <- (im.ras4- mean(im.ras4))/sd(im.ras4) ; summary(norm.im.ras4)
norm.im.ras5 <- (im.ras5- mean(im.ras5))/sd(im.ras5) ; summary(norm.im.ras5)


# Plots of images
par(mfrow=c(2, 3))
image.plot(list(x=im.ras1$xcol, y=im.ras1$yrow, z=t(im.ras1$v)), main= "Solar radiation sim", asp=1)
image.plot(list(x=im.ras2$xcol, y=im.ras2$yrow, z=t(im.ras2$v)), main= "Precipitation sim", asp=1)
image.plot(list(x=im.ras3$xcol, y=im.ras3$yrow, z=t(im.ras3$v)), main= "Temperature sim", asp=1)
image.plot(list(x=im.ras4$xcol, y=im.ras4$yrow, z=t(norm.im.ras4$v)), main= "Wind speed sim", asp=1)
image.plot(list(x=im.ras5$xcol, y=im.ras5$yrow, z=t(norm.im.ras5$v)), main= "Water vapour pressure sim", asp=1)

list.covariates <- list(Solar.radiation.sim=norm.im.ras1 , Precipitation.sim=norm.im.ras2,
                        Temperature.sim=norm.im.ras3, Wind.speed.sim=norm.im.ras4,
                        Water.vapour.pressure.sim=norm.im.ras5)

# Compute the log lambda
alpha0 = 0
Beta <- c(1.3, 0.1, 0.1, 0.12, 0.3)
log.lambda <-alpha0 + Beta[1]*norm.im.ras1 + Beta[2]*norm.im.ras2 + Beta[3]*norm.im.ras3 + Beta[4]*norm.im.ras4^2 + Beta[5]*sin(pi*norm.im.ras5)
summary(log.lambda)

par(mfrow=c(1,2))
image.plot(list(x=log.lambda$xcol, y=log.lambda$yrow, z=t(log.lambda$v)), main= " log-intensity(z) (image)", asp=1)
## Plot the Global environment as Raster
Virtual.species.domaine <- raster(log.lambda)
plot(Virtual.species.domaine, main='log-intensity(z) (Raster)')
graphics.off()


################# Model part ###########
#LGCP simulation
lgcp.sim <- lgcp.sim <- rLGCP("matern", mu=log.lambda, var=0.5, scale=0.05, nu=1)
summary(lgcp.sim)
lam=attr(lgcp.sim, "Lambda")
log.lam=log(lam)
summary(log.lam)

image.plot(list(x=log.lambda$xcol, y=log.lambda$yrow, z=t(log.lambda$v)), main= "Extensive Sampling in defined subregion", asp=1)
points(lgcp.sim$x, lgcp.sim$y, cex=0.8, pch=20, lwd=1)
samp.x=lgcp.sim$x[lgcp.sim$x>=3 & lgcp.sim$x<=9]
sampy=lgcp.sim$y[lgcp.sim$y>=2 & lgcp.sim$y<=9]
samp.y=sample(sampy, length(samp.x))
points(samp.x, samp.y, cex=0.8, pch=20, lwd=1, col="white")

Win.sample <- owin(xrange = c(min(samp.x), max(samp.x)), 
                   yrange = c(min(samp.y), max(samp.y)), 
                   poly = NULL,
                   unitname = "km")
ppp.sam=ppp(samp.x, samp.y, window = Win.sample)


#'  *test the complete spacial rendomness*
q.testp <- quadrat.test(ppp.sam, method = "Chisq")
print(q.testp)

#???' Envelope
Sample.Envelope= envelope(ppp.sam, nsim = 30)
plot(Sample.Envelope)



Solar.rad <- norm.im.ras1[Win.sample, drop=F]
Preci <- norm.im.ras2[Win.sample, drop=F]
Temp <- norm.im.ras3[Win.sample, drop=F]
Wind.s <- norm.im.ras4[Win.sample, drop=F]
Water.va<- norm.im.ras5[Win.sample, drop=F]


# Plots of images
par(mfrow=c(2, 3))
plot(ppp.sam, cex=0.8, pch=20)
image.plot(list(x=Solar.rad$xcol, y=Solar.rad$yrow, z=t(Solar.rad$v)), main= "Solar radiation sim", asp=1)
image.plot(list(x=Preci$xcol, y=Preci$yrow, z=t(Preci$v)), main= "Precipitation sim", asp=1)
image.plot(list(x=Temp$xcol, y=Temp$yrow, z=t(Temp$v)), main= "Temperature sim", asp=1)
image.plot(list(x=Wind.s$xcol, y=Wind.s$yrow, z=t(Wind.s$v)), main= "Wind speed sim", asp=1)
image.plot(list(x=Water.va$xcol, y=Water.va$yrow, z=t(Water.va$v)), main= "Water vapour pressure sim", asp=1)


# Model formula
covar.sample.list <- list(Solar.radiation.sim=Solar.rad , Precipitation.sim=exp(Preci), 
                          Temperature.sim=Temp,
                          Wind.speed.sim=Wind.s^2, Water.vapour.pressure.sim=sin(pi*Water.va))

form <- as.formula(paste("~", paste(names(covar.sample.list), collapse = "+")))

#'Fit the *inhomogeneous poisson process model*
#'#mpl
par(mfrow=c(2,3))

Qu <- pixelquad(ppp.sam, Win.sample)

fit.ipp1 <- ppm(Qu, form, covariates = covar.sample.list,
                interaction=NULL, method = 'mpl',eps=c(1,1))
params.fit.ipp1 <- as.vector(coef(fit.ipp1))
IPP.MLP=envelope(fit.ipp1, nsim = 30)
par(mfrow=c(1,2))
plot(IPP.MLP)
graphics.off()

predict.ip <- predict(fit.ipp1,window=Win, covariates = list.covariates, eps=c(1,1))
summary(predict.ip)
Prediction_IPP.MPL=log(predict.ip)
summary(Prediction_IPP.MPL)
plot(Prediction_IPP.MPL, main="IPP.MPL Prediction")


#VBlogi
fit.ipp2 <- ppm(Qu, form, covariates = covar.sample.list, eps=c(1,1),
                interaction = NULL, method = 'VBlogi')
params.fit.ipp2 <- as.vector(coef(fit.ipp2))
IPP.VBlogi=envelope(fit.ipp2, nsim = 30)
par(mfrow=c(1,2))
plot(IPP.VBlogi)
graphics.off()

predict.ip2 <- predict.ppm(fit.ipp2,window=Win,type = "intensity", covariates = list.covariates, eps=c(1,1))
summary(predict.ip2)
par(mfrow=c(1,2))
Prediction_IPP.VBlogi=log(predict.ip2)
summary(Prediction_IPP.VBlogi)
plot(Prediction_IPP.VBlogi, main="IPP.VBlogi Prediction")

#' #' Fit the *gibbs process model*
#' #mpl

fit.gibbs1 <- ppm(Qu, form, covariates = covar.sample.list,
                  interaction= AreaInter(1), method = "mpl")
params.fit.gibbs1 <- as.vector(coef(fit.gibbs1)[-7])
GPP.MPL=envelope(fit.gibbs1, nsim = 30)
par(mfrow=c(1,2))
plot(GPP.MPL)
graphics.off()

predict.gb <- predict.ppm(fit.gibbs1, window=Win, covariates = list.covariates, type="intensity", eps=c(1,1))
summary(predict.gb)
plot(log.lam)
Predictive_GBB.MPL=log(predict.gb)
summary(Predictive_GBB.MPL)
plot(Predictive_GBB.MPL, main="GPP.MPL Prediction")


#' #VBlogi
fit.gibbs2 <- ppm(Qu, form, covariates = covar.sample.list,
                  interaction=AreaInter(1), method = "VBlogi")
params.fit.gibbs2 <- as.vector(coef(fit.gibbs2)[-7])
GPP.VBlogi=envelope(fit.gibbs2, nsim = 30)
par(mfrow=c(1,2))
plot(GPP.VBlogi)
graphics.off()

predict.gb2 <- predict(fit.gibbs2, window=Win,type="intensity", covariates = list.covariates, eps=c(1,1))
summary(predict.gb2)
Prediction_GPP.VBlogi=log(predict.gb2)
summary(Prediction_GPP.VBlogi)
plot(Prediction_GPP.VBlogi, main="GPP.VBlogi Prediction")

#' #' *Fit the cox process model*
#' #' 
#' #clik2
fit.cox1 <- kppm(Qu, form, covariates = covar.sample.list,
                 clusters = 'LGCP',
                 method = "clik2",  var=0.5, scale=0.05, covmodel=list(model="matern", nu=1), epsilon = c(1,1))
params.fit.cox1 <- as.vector(coef(fit.cox1))
CPP.CLIK=envelope(fit.cox1, nsim = 30)
plot(CPP.CLIK)

predict.cx<- predict(fit.cox1,window=Win,type="intensity", covariates = list.covariates, eps=c(1,1))
summary(predict.cx)
par(mfrow=c(1,2))
prediction_cox.CLIK=log(predict.cx)
summary(prediction_cox.CLIK)
plot(prediction_cox.CLIK, main="cox.CLIK prediction")

#' 
#' #palm
fit.cox2 <- kppm(Qu, form, covariates = covar.sample.list,
                 clusters = 'LGCP',
                 method = "palm",  var=0.5, scale=0.05, covmodel=list(model="matern", nu=1), eps = c(1,1))
params.fit.cox2 <- as.vector(coef(fit.cox2))

CPP.PALM=envelope(fit.cox2, nsim = 30)
par(mfrow=c(1,2))
plot(CPP.PALM)

predict.cx2<- predict(fit.cox2,window=Win,type="intensity", covariates = list.covariates, eps=c(1,1))
summary(predict.cx2)
Prediction_CPP.PALM=log(predict.cx2)
summary(Prediction_CPP.PALM)
plot(Prediction_CPP.PALM, main="CPP.PALM Prediction")


## Computing BIC 
BIC(fit.ipp1)
BIC(fit.ipp2)
BIC(fit.gibbs1)
BIC(fit.gibbs2)
BIC(fit.cox1)
BIC(fit.cox2)

## Computing residuals
residuals.ppm(fit.ipp1)
residuals.ppm(fit.ipp2)
residuals.ppm(fit.gibbs1)
residuals.ppm(fit.gibbs2)
residuals.kppm(fit.cox1)
residuals.kppm(fit.cox2)
