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
#####' @details In this part we're going to use *test of conformity of mean* to
#####' *check whether PPMs are able to determine the true coefficients*

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

list.covariates <- list(Solar.radiation.sim=norm.im.ras1, Precipitation.sim=norm.im.ras2,
                        Temperature.sim=norm.im.ras3, Wind.speed.sim=norm.im.ras4,
                        Water.vapour.pressure.sim=norm.im.ras5)

# Compute the log lambda
alpha0 = 0.1
Beta <- c(1.3, 0.1, 0.1, 0.12, 0.3)
log.lambda <-alpha0 + Beta[1]*norm.im.ras1 + Beta[2]*norm.im.ras2 + Beta[3]*norm.im.ras3 + Beta[4]*norm.im.ras4^2 + Beta[5]*sin(pi*norm.im.ras5)
summary(log.lambda)

par(mfrow=c(1,2))
image.plot(list(x=log.lambda$xcol, y=log.lambda$yrow, z=t(log.lambda$v)), main= " log-intensity(z) (image)", asp=1)
## Plot the Global environment as Raster
Virtual.species.domaine <- raster(log.lambda)
plot(Virtual.species.domaine, main='log-intensity(z) (Raster)')
graphics.off()

###### Simulation  part
Simulated <- function(seed){
  set.seed(seed)
  lgcp.sim <- lgcp.sim <- rLGCP("matern", mu=log.lambda, var=0.5, scale=0.05, nu=1)
  summary(lgcp.sim)
  
  n.list.covariates <- list(Solar.radiation.sim=norm.im.ras1 , Precipitation.sim=norm.im.ras2,
                            Temperature.sim=norm.im.ras3, Wind.speed.sim=norm.im.ras4^2,
                            Water.vapour.pressure.sim=sin(pi*norm.im.ras5))
  
  
  form <- as.formula(paste("~", paste(names(n.list.covariates), collapse = "+")))
  
  Qu <- pixelquad(lgcp.sim, Win)
  
  #'Fit the *inhomogeneous poisson process model*
  #'#mpl
  #'
  fit.ipp1 <- ppm(Qu, form, covariates = n.list.covariates, 
                  interaction=NULL, method = 'mpl',eps=c(1,1))
  params.fit.ipp1 <- as.vector(coef(fit.ipp1))
  
  #VBlogi
  fit.ipp2 <- ppm(Qu, form, covariates = n.list.covariates, 
                  interaction=NULL, method = 'VBlogi',eps=c(1,1))
  params.fit.ipp2 <- as.vector(coef(fit.ipp2))
  
  #' #' Fit the *gibbs process model*
  #' #mpl
  fit.gibbs1 <- ppm.ppp(Qu, form, covariates = n.list.covariates,
                        interaction=AreaInter(1), method = 'mpl', eps=c(1,1))
  params.fit.gibbs1 <- as.vector(coef(fit.gibbs1)[-7])
  #' #VBlogi
  fit.gibbs2 <- ppm(Qu, form, covariates = n.list.covariates,
                    interaction=AreaInter(1), method = 'VBlogi', eps=c(1,1))
  params.fit.gibbs2 <- as.vector(coef(fit.gibbs2)[-7])
  
  #' #' *Fit the cox process model*
  #' #' 
  #' #clik2
  
  fit.cox1 <- kppm(Qu, form, covariates = n.list.covariates,
                   clusters = 'LGCP',
                   method = "clik2",  var=0.5, scale=0.05, covmodel=list(model="matern", nu=1), eps=c(1,1))
  params.fit.cox1 <- as.vector(coef(fit.cox1))
  #' 
  #' 
  #' #palm
  fit.cox2 <- kppm(Qu, form, covariates = n.list.covariates,
                   clusters = 'LGCP',
                   method = "palm",  var=0.5, scale=0.05, covmodel=list(model="matern", nu=1), eps = c(1,1))
  params.fit.cox2 <- as.vector(coef(fit.cox2))

   print(1+0) 
  
  return(list(ipp.mpl = params.fit.ipp1, ipp.VBlogi = params.fit.ipp2,
              gibbs.mpl=params.fit.gibbs1, gibbs.VBlogi = params.fit.gibbs2,
              cox.clik2=params.fit.cox1, cox.palm=params.fit.cox2))
}

all.sim = function(n.sim){
  
  ipp.mpl = ipp.VBlogi=gibbs.mpl = gibbs.VBlogi=cox.clik2 =cox.palm= matrix(numeric(6 * n.sim), ncol = 6)
  sed = numeric(n.sim)
  for(i in 1:n.sim){
    sed[i] = sample(1:5000, 1)
    sims = Simulated(sed[i])
    ipp.mpl[i,] = sims$ipp.mpl
    ipp.VBlogi[i,] = sims$ipp.VBlogi
    gibbs.mpl[i,] = sims$gibbs.mpl
    gibbs.VBlogi[i,] = sims$gibbs.VBlogi
    cox.clik2[i,] = sims$cox.clik2
    cox.palm[i,] = sims$cox.palm
  }
  colnames(ipp.mpl) =  colnames(ipp.VBlogi) =  colnames(gibbs.mpl)=  colnames(gibbs.VBlogi) = colnames(cox.clik2) = colnames(cox.palm) = c("Intercept", 'Sol.rad', 'Rainfall',
                                                                                                                                           'Temp', 'Wind.speed', 'Wat.vap.pres')
  
  ipp.Mpl = data.frame(N.sim = 1:n.sim, ipp.mpl)
  ipp.VBLogi = data.frame(N.sim = 1:n.sim, ipp.VBlogi)
  gibbs.Mpl = data.frame(N.sim = 1:n.sim, gibbs.mpl)
  gibbs.VBLogi = data.frame(N.sim = 1:n.sim, gibbs.VBlogi)
  cox.Clik2 = data.frame(N.sim = 1:n.sim, cox.clik2)
  cox.Palm = data.frame(N.sim=1:n.sim, cox.palm)
  
  return(list(ipp.mpl = ipp.Mpl, ipp.VBlogi = ipp.VBLogi,
              gibbs.mpl = gibbs.Mpl, gibbs.VBlogi = gibbs.VBLogi,
              cox.clik2=cox.Clik2,cox.palm=cox.Palm,
              Seeds = sed))
}

Simul <- all.sim(n.sim = 100)
length(Simul)  
structure(Simul)


SHapiRo.T = function(){
  param1 <- as.data.frame(sapply(Simul, '[', "Intercept"))
  alp0 <- apply(param1[,1:6], 2, shapiro.test) 
  Intercept <- round(as.data.frame(sapply(alp0, '[', 'p.value')), digits = 3)
  colnames(Intercept) <- c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  row.names(Intercept) <- c("P.values:")
  
  param2 <- as.data.frame(sapply(Simul, '[', "Sol.rad"))
  sol <- apply(param2[,1:6], 2, shapiro.test) 
  Sol.rad <- round(as.data.frame(sapply(sol, '[', 'p.value')), digits = 3)
  colnames(Sol.rad) <- c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  row.names(Sol.rad) <- c("P.values:")
  
  param3 <- as.data.frame(sapply(Simul, '[', "Rainfall"))
  rain <- apply(param3[,1:6], 2, shapiro.test) 
  Rainfall <- round(as.data.frame(sapply(rain, '[', 'p.value')), digits = 3)
  colnames(Rainfall) <- c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  row.names(Rainfall) <- c("P.values:")
  
  param4 <- as.data.frame(sapply(Simul, '[', "Temp"))
  Te <- apply(param4[,1:6], 2, shapiro.test) 
  Temp <- round(as.data.frame(sapply(Te, '[', 'p.value')), digits = 3)
  colnames(Temp) <- c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  row.names(Temp) <- c("P.values:")
  
  param5 <- as.data.frame(sapply(Simul, '[', "Wind.speed"))
  Wind <- apply(param5[,1:6], 2, shapiro.test) 
  Wind.speed <- round(as.data.frame(sapply(Wind, '[', 'p.value')), digits = 3)
  colnames(Wind.speed) <- c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  row.names(Wind.speed) <- c("P.values:")
  
  param6 <- as.data.frame(sapply(Simul, '[', "Wat.vap.pres"))
  Wat <- apply(param6[,1:6], 2, shapiro.test) 
  Wat.vap.pres <- round(as.data.frame(sapply(Wat, '[', 'p.value')), digits = 3)
  colnames(Wat.vap.pres) <- c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  row.names(Wat.vap.pres) <- c("P.values:")
  
  P_values <- list(Intercept=Intercept, Solar.rad=Sol.rad, Precipitation=Rainfall, 
                   Temp=Temp,Wind.speed=Wind.speed, Wat.vap.pres=Wat.vap.pres)
  
  return(P_values)
}
SHapiRO <- SHapiRo.T()
SHapiRO


TTest = function(){
  param1 <- as.data.frame(sapply(Simul, '[', "Intercept"))
  # Here the true values of mu (alpha0) is 0
  tab= sapply(param1, "[")
  alp0 <- apply(param1[,1:6], 2, function(x) t.test(tab, mu=alpha0))
  Interceptp <- as.data.frame(sapply(alp0, '[', 'p.value'))
  Interceptc <- as.data.frame(sapply(alp0, '[', 'conf.int'))
  colnames(Interceptp)= colnames(Interceptc)= c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  Intercept=round(rbind(Interceptp, Interceptc), digits = 3)
  row.names(Intercept) <- c("P.values:", "CI.lower.bound:", "CI.upper.bound:")
  
  param2 <- as.data.frame(sapply(Simul, '[', "Sol.rad"))
  tab2 = sapply(param2, "[")
  sol <- apply(param2[,1:6], 2, function(x) t.test(tab2, mu=Beta[1])) 
  Sol.radp <- as.data.frame(sapply(sol, '[', 'p.value'))
  Sol.radc <- as.data.frame(sapply(sol, '[', 'conf.int'))
  colnames(Sol.radp)= colnames(Sol.radc)= c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  Sol.rad = round(rbind(Sol.radp, Sol.radc), digits = 3)
  row.names(Sol.rad) = c("P.values:", "CI.lower.bound:", "CI.upper.bound:")
  
  param3 <- as.data.frame(sapply(Simul, '[', "Rainfall"))
  tab3= sapply(param3, "[")
  Rain <- apply(param2[,1:6], 2, function(x) t.test(tab3, mu=Beta[2])) 
  Rainp <- as.data.frame(sapply(Rain, '[', 'p.value'))
  Rainc <- as.data.frame(sapply(Rain, '[', 'conf.int'))
  colnames(Rainp)= colnames(Rainc)= c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  Rainfall = round(rbind(Rainp, Rainc), digits = 3)
  row.names(Rainfall) = c("P.values:", "CI.lower.bound:", "CI.upper.bound:")
  
  param4 <- as.data.frame(sapply(Simul, '[', "Temp"))
  tab4= sapply(param4, "[")
  Te <- apply(param4[,1:6], 2, function(x) t.test(tab4, mu=Beta[3])) 
  Tempp <- as.data.frame(sapply(Te, '[', 'p.value'))
  Tempc <- as.data.frame(sapply(Te, '[', 'conf.int'))
  colnames(Tempp) =colnames(Tempc)= c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  Temp = round(rbind(Tempp, Tempc), digits = 3)
  row.names(Temp) = c("P.values:", "CI.lower.bound:", "CI.upper.bound:")
  
  param5 <- as.data.frame(sapply(Simul, '[', "Wind.speed"))
  tab5= sapply(param5, "[")
  Wind <- apply(param5[,1:6], 2, function(x) t.test(tab5, mu=Beta[4])) 
  Wind.speedp <- as.data.frame(sapply(Wind, '[', 'p.value'))
  Wind.speedc <- as.data.frame(sapply(Wind, '[', 'conf.int'))
  colnames(Wind.speedp)= colnames(Wind.speedc) = c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  Wind.speed= round(rbind(Wind.speedp, Wind.speedc), digits = 3)
  row.names(Wind.speed) = c("P.values:", "CI.lower.bound:", "CI.upper.bound:")
  
  param6 <- as.data.frame(sapply(Simul, '[', "Wat.vap.pres"))
  tab6= sapply(param6, "[")
  Wat <- apply(param6[,1:6], 2, function(x) t.test(tab6, mu=Beta[5])) 
  Wat.vap.presp <- as.data.frame(sapply(Wat, '[', 'p.value'))
  Wat.vap.presc <- as.data.frame(sapply(Wat, '[', 'conf.int'))
  colnames(Wat.vap.presp)= colnames(Wat.vap.presc)= c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  Wat.vap.pres = round(rbind(Wat.vap.presp, Wat.vap.presc), digits = 3)
  row.names(Wat.vap.pres) = c("P.values:", "CI.lower.bound:", "CI.upper.bound:")
  
  T.t <- list(Intercept=Intercept, Solar.rad=Sol.rad, Rainfall=Rainfall, 
              Temp=Temp,Wind.speed=Wind.speed, Wat.vap.pres=Wat.vap.pres)
  
  return(T.t)
}
T.Test <- TTest()
T.Test


Wilk.Test = function(){
  param1 <- as.data.frame(sapply(Simul, '[', "Intercept"))
  # Here the true values of mu (alpha0) is 0
  tab= sapply(param1, "[")
  alp0 <- apply(param1[,1:6], 2, function(x) wilcox.test(tab, mu=alpha0))
  Interceptp <- as.data.frame(sapply(alp0, '[', 'p.value'))
  colnames(Interceptp)= c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  Intercept=round(Interceptp, digits = 3)
  row.names(Intercept) <- c("P.values:")
  
  param2 <- as.data.frame(sapply(Simul, '[', "Sol.rad"))
  tab2 = sapply(param2, "[")
  sol <- apply(param2[,1:6], 2, function(x) wilcox.test(tab2, mu=Beta[1])) 
  Sol.radp <- as.data.frame(sapply(sol, '[', 'p.value'))
  colnames(Sol.radp)=  c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  Sol.rad = round(Sol.radp, digits = 3)
  row.names(Sol.rad) = c("P.values:")
  
  param3 <- as.data.frame(sapply(Simul, '[', "Rainfall"))
  tab3= sapply(param3, "[")
  Rain <- apply(param2[,1:6], 2, function(x) wilcox.test(tab3, mu=Beta[2])) 
  Rainp <- as.data.frame(sapply(Rain, '[', 'p.value'))
  colnames(Rainp)= c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  Rainfall = round(Rainp, digits = 3)
  row.names(Rainfall) = c("P.values:")
  
  param4 <- as.data.frame(sapply(Simul, '[', "Temp"))
  tab4= sapply(param4, "[")
  Te <- apply(param4[,1:6], 2, function(x) wilcox.test(tab4, mu=Beta[3])) 
  Tempp <- as.data.frame(sapply(Te, '[', 'p.value'))
  colnames(Tempp) =c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  Temp = round(Tempp, digits = 3)
  row.names(Temp) = c("P.values:")
  
  param5 <- as.data.frame(sapply(Simul, '[', "Wind.speed"))
  tab5= sapply(param5, "[")
  Wind <- apply(param5[,1:6], 2, function(x) wilcox.test(tab5, mu=Beta[4])) 
  Wind.speedp <- as.data.frame(sapply(Wind, '[', 'p.value'))
  colnames(Wind.speedp) = c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  Wind.speed= round(Wind.speedp, digits = 3)
  row.names(Wind.speed) = c("P.values:")
  
  param6 <- as.data.frame(sapply(Simul, '[', "Wat.vap.pres"))
  tab6= sapply(param6, "[")
  Wat <- apply(param6[,1:6], 2, function(x) wilcox.test(tab6, mu=Beta[5])) 
  Wat.vap.presp <- as.data.frame(sapply(Wat, '[', 'p.value'))
  colnames(Wat.vap.presp)= c("ipp.mpl",'ipp.VBlogi',"gibbs.mpl",'gibbs.VBlogi', "cox.clik2", 'cox.palm')
  Wat.vap.pres = round(Wat.vap.presp, digits = 3)
  row.names(Wat.vap.pres) = c("P.values:")
  
  W.t <- list(Intercept=Intercept, Solar.rad=Sol.rad, Rainfall=Rainfall, 
              Temp=Temp,Wind.speed=Wind.speed, Wat.vap.pres=Wat.vap.pres)
  
  return(W.t)
}
Wilkoxon.Test <- Wilk.Test()
Wilkoxon.Test
as.data.frame(Wilkoxon.Test)