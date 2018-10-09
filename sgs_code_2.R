#ANPP data for github

#input datafiles
anpp.sgs<-read.csv(file.choose(),header=TRUE) #anpp dataset
vwc<-read.csv(file.choose(),header=TRUE) #soil moisture data
n<-read.csv(file.choose(),header=TRUE) #soil nitrogen data

#Averaging sub plot values to produce single value per plot
anpp.ag<-aggregate(x.100 ~ Plot + mm + percentile + vwc.20cm + subset + Block,mean,data=anpp.sgs)
grass.ag<-aggregate(gras.100 ~ Plot + mm + percentile + vwc.20cm + subset + Block,mean,data=anpp.sgs)
forb.ag<-aggregate(forb.100 ~ Plot + mm + percentile + vwc.20cm + subset + Block,mean,data=anpp.sgs)
vwc.ag<-aggregate(calibrated~  mm,mean,data=vwc)
plot(gras.100~mm,data=grass.ag)
plot(x.100~mm,data=anpp.ag)
#moisture and N linear models
lm.vwc<-lm(calibrated ~  mm,data=vwc.ag)
summary(lm.vwc)
n.lm<-lm(total.N~mm,data=n.sgs)
summary(n.lm)

#after ifentifying influential outlier observation
#precipitation
anpp.ag2<-anpp.ag[-4,]
grass.ag2<-grass.ag[-22,]
anpp.precip<-aggregate(x.100~mm,mean,data=anpp.ag2)


#soil moisture
#none found

#identifying potential outliers
#leverage
lev =hat(model.matrix(linear.anpp.sgs)) 
plot(lev)

#assessing residuals
growth.res = resid(linear.anpp.sgs) 
plot(anpp.sgs$ANPP.g,growth.res,
     ylab="Residuals", xlab="predictor", 
     main="Y Residuals")
abline(0,0)

qqnorm(linear.anpp.sgs$res)
qqline(linear.anpp.sgs$res)

#outliers for precip
library(car)
linear.anpp.sgs<-lm(x.100~mm,data=anpp.ag)
outlierTest(linear.anpp.sgs) #remove one observation

linear.grass.sgs<-lm(gras.100~mm,data=grass.ag)
outlierTest(linear.grass.sgs) #one outlier indentified
plot(gras.100~mm,data=grass.ag2)

#outliers for moisture
head(anpp.ag)
#total ANPP
linear.anpp.sgs.moisture<-lm(x.100~vwc.20cm ,data=anpp.ag)
outlierTest(linear.anpp.sgs.moisture) #no outliers in response to moisture

#grass 
linear.grass.sgs.vwc<-lm(gras.100~vwc.20cm,data=grass.ag)
outlierTest(linear.grass.sgs.vwc) #no outliers
#testing block effects
library(nlme)
anpp.block<-lme(x.100 ~ mm + vwc.20cm + Block,random=~1|Plot,data=anpp.ag2)
summary(anpp.block) #no block effect

#asymm, calculated with all data so comparisons are derived from the same dataset.
head(anpp.ag)
anpp.ag.asymm<-aggregate(x.100  ~ subset,mean,data=anpp.ag) #grass and forb
grass.ag.asymm<-aggregate(gras.100  ~ subset,mean,data=grass.ag) #
forb.ag.asymm<-aggregate(forb.100  ~ subset,mean,data=forb.ag)

#ANPP grass and forb
((101.47 - 70.49) - (70.49 - 51.97))/70.49 #.18, positive asymmetry
#grass biomass
((87.84 - 65.075) - (65.075 - 50.82))/65.075 # 0.13., positive asymmetry
#forb biomass
((13.63 - 5.41) - (5.41 - 0.86))/5.41 #0.68, positive asymmetry

#% increase
(0.18-0.13)/0.18 #%27

#linearmodel with precip
anpp.sgs.linear<-lm(x.100~mm,data=anpp.ag2)
summary(anpp.sgs.linear) #r-square = 0.55

#nonlinear least squares model
library(nls2)
library(aod)
library(nlstools)
library(lattice)
library(boot)
head(anpp)

# 2017 ANPP response to precipitation
plot(x.100~mm,data=anpp.ag2)
abline(lm(x.100~mm,data=anpp.ag2))
#find starting parameters
lm(logit(x.100/300)~mm,anpp.ag2)
#run model: logistic growth since the data appears to conform to a logsitic relationship
anpp.sgs.nonlinear<- nls2(x.100 ~ theta1/(1 + exp(-(theta2 + theta3*mm))),
                          start=list(theta1 = 300, theta2 =   -1.815647   , theta3 = 0.004101),data=anpp.ag2, trace=TRUE)
summary(anpp.sgs.nonlinear)
coef(anpp.2015.model)

#forggplot trendline
?range
r <- range(anpp.ag2$mm)
xNew <- seq(r[1],r[2])
yNew <- predict(anpp.sgs.nonlinear,list(mm = xNew))
predictanpp.sgs<-data.frame(xNew,yNew)
predictanpp.sgs
summary(predictanpp)
lines(xNew,yNew)
#comparing linear forsus nonlinear model:

#AIC weights
AIC(anpp.sgs.nonlinear,anpp.sgs.linear)
x <- c(375.4388, 374.4790)         # stores two AIC values in a vector
delta <- x - min(x)          # AIC differences
L <- exp(-0.5 * delta)      #relative likelihoods 
w <- L/sum(L) #AIC weights or posterier probabilities of the models
w #AIC weights



#ANPP grass and forb responses to vwc 20 cm

#response to soil moisture
anpp.sgs.vwc.linear<-lm(x.100~vwc.20cm,data=anpp.ag)
summary(anpp.sgs.vwc.linear) #r-square is 0.6
plot(x.100~vwc.20cm,data=anpp.ag)
abline(lm(x.100~vwc.20cm,data=anpp.ag))

#find starting parameters
lm(logit(x.100/300)~vwc.20cm,anpp.ag)

#run model: logistic growth since the data appears to conform to a logsitic relationship
anpp.sgs.vwc.nonlinear.model<- nls2(x.100~ theta1/(1 + exp(-(theta2 + theta3*vwc.20cm))),
                                    start=list(theta1 = 300, theta2 = -3.9641   , theta3 = 0.1604  ),data=anpp.ag, trace=TRUE)

#forggplot trendline
?range
r <- range(anpp.ag$vwc.20cm)
xNew <- seq(r[1],r[2])
yNew <- predict(anpp.sgs.vwc.nonlinear.model,list(vwc.20cm = xNew))
predictanpp.sgs.vwc<-data.frame(xNew,yNew)
lines(xNew,yNew)
#comparing linear forsus nonlinear model:

#AIC weights
AIC(anpp.sgs.vwc.nonlinear.model,anpp.sgs.vwc.linear)
x <- c(386.8048, 385.9754)         # stores two AIC values in a vector
delta <- x - min(x)          # AIC differences
L <- exp(-0.5 * delta)      #relative likelihoods 
w <- L/sum(L) #AIC weights or posterier probabilities of the models
w

#looking at grass responses to precipitation

#linearmodel with precip
grass.sgs.linear<-lm(gras.100~mm,data=grass.ag2)
summary(grass.sgs.linear) #r-square = 0.47

#nonlinear least squares model

# 2017 ANPP grass response to precipitation
plot(gras.100~mm,data=grass.ag2)
abline(lm(gras.100~mm,data=grass.ag2))

#find starting parmters
lm(logit(gras.100/300)~mm,grass.ag2)
#run model: logistic growth since the data appears to conform to a logsitic relationship
grass.sgs.nonlinear<- nls2(gras.100 ~ theta1/(1 + exp(-(theta2 + theta3*mm))),
                           start=list(theta1 = 300, theta2 =   -1.779922   , theta3 =  0.003081     ),data=grass.ag2, trace=TRUE)


#forggplot trendline
?range
r <- range(grass.ag2$mm)
xNew <- seq(r[1],r[2])
yNew <- predict(grass.sgs.nonlinear,list(mm = xNew))
predictanpp.grass.sgs<-data.frame(xNew,yNew)
predictanpp
summary(predictanpp)
lines(xNew,yNew)

#comparing linear forsus nonlinear model:

#AIC weights
AIC(grass.sgs.nonlinear,grass.sgs.linear)
x <- c(360.9314, 359.0931)         # stores two AIC values in a vector
delta <- x - min(x)          # AIC differences
L <- exp(-0.5 * delta)      #relative likelihoods 
w <- L/sum(L) #AIC weights or posterier probabilities of the models
w

#ANPP grass responses to vwc 20 cm

#response to soil moisture
anpp.grass.sgs.vwc.linear<-lm(gras.100~vwc.20cm,data=grass.ag)
summary(anpp.grass.sgs.vwc.linear) #r-square is 0.53
plot(gras.100~vwc.20cm,data=grass.ag)
abline(lm(gras.100~vwc.20cm,data=grass.ag))

#find starting parameters
lm(logit(gras.100/300)~vwc.20cm,grass.ag)

#run model: logistic growth since the data appears to conform to a logsitic relationship
anpp.grass.sgs.vwc.nonlinear.model<- nls2(gras.100~ theta1/(1 + exp(-(theta2 + theta3*vwc.20cm))),
                                          start=list(theta1 = 300, theta2 = -3.5022, theta3 = 0.1271  ),data=grass.ag, trace=TRUE)

#forggplot trendline
?range
r <- range(anpp.ag$vwc.20cm)
xNew <- seq(r[1],r[2])
yNew <- predict(anpp.grass.sgs.vwc.nonlinear.model,list(vwc.20cm = xNew))
predictanpp.sgs.grass.vwc<-data.frame(xNew,yNew)
lines(xNew,yNew)
#comparing linear forsus nonlinear model:

#AIC weights
AIC(anpp.grass.sgs.vwc.nonlinear.model,anpp.grass.sgs.vwc.linear)
x <- c(374.4406, 372.8795)         # stores two AIC values in a vector
delta <- x - min(x)          # AIC differences
L <- exp(-0.5 * delta)      #relative likelihoods 
w <- L/sum(L) #AIC weights or posterier probabilities of the models
w

#effect sizes using all data like asymmetry
?cohen.d
library(effsize)
#ANPP grass + forb
mean.anpp.sgs<-subset(anpp.ag,subset=="nominal",na.rm=TRUE)
wettest.anpp.sgs<-subset(anpp.ag,subset=="extreme.wet",na.rm=TRUE)
driest.anpp.sgs<-subset(anpp.ag,subset=="extreme.dry",na.rm=TRUE)

#vector
mean.anpp.sgs.vec<-as.vector(mean.anpp.sgs$x.100)
wettest.anpp.sgs.vec<-as.vector(wettest.anpp.sgs$x.100)
driest.anpp.sgs.vec<-as.vector(driest.anpp.sgs$x.100)

#ANPP grass alone
mean.grass.sgs<-subset(grass.ag,subset=="nominal",na.rm=TRUE)
wettest.grass.sgs<-subset(grass.ag,subset=="extreme.wet",na.rm=TRUE)
driest.grass.sgs<-subset(grass.ag,subset=="extreme.dry",na.rm=TRUE)

#vector
mean.grass.sgs.vec<-as.vector(mean.grass.sgs$gras.100)
wettest.grass.sgs.vec<-as.vector(wettest.grass.sgs$gras.100)
driest.grass.sgs.vec<-as.vector(driest.grass.sgs$gras.100)

cohen.d(driest.grass.sgs.vec, mean.grass.sgs.vec,pooled=TRUE,paired=FALSE,
        na.rm=TRUE, hedges.correction=FALSE,
        conf.level=0.95,noncentral=FALSE)

#species composition analysis

comp.sgs<-read.csv(file.choose(),header=TRUE)
head(comp.sgs)


library(reshape2)
library(BiodiversityR)
library(vegan)
library(MASS)
library(plyr)

#similarity percentages, divergence per species and total;
simp<-(comp.sgs[,7:31])
sim<-with(comp.sgs,simper(simp,group)) 
summary(sim,ordered=TRUE,digits=2) 
lapply(sim,FUN=function(x){x$overall}) 

#permanova
#bray curtis distance analysis
extreme<-subset(comp.sgs,subgroup=="extreme")
adonis(extreme[,7:31]~extreme$group,data=extreme,method="bray") #communities differed by year

#SOIL RESPIRATION

#input datafiles
resp<-read.csv(file.choose(),header=TRUE) #soil respiration datafile
head(resp)

#assess raw data for normality, outliers and influential data points
linear.resp.vwc<-lm(flux~calibrated,data=resp)
linear.resp.mm<-lm(flux~mm,data=resp)
summary(linear.resp)
plot(flux~calibrated,data=resp)
#residuals
qqnorm(linear.resp$res) #implies a resonable degree of normality, except at the tails, indicates nonlinearity
qqline(linear.resp$res)

#leverage
lev =hat(model.matrix(linear.resp)) #suggests certain points have high leveage, though this isn't higher than 0.1
plot(lev)

#outlier assessment
library(car)
outlierTest(linear.resp.vwc) #one observation as clear outlier
resp.2<-resp[-210,] #new dataset without the outler
outlierTest(linear.resp.mm) #none for precipitation

# Influence Plot
#cooks distance
cutoff <- 4/((nrow(resp)-length(linear.resp$coefficients)-2))
plot(linear.resp, which=4, cook.levels=cutoff)

#no indiciation of singificant outliers skewing the regression
head(resp)

#lme
library(nlme)
library(car)
plot(flux~temp.C,data=resp.2)
plot(temp.C~vwc.20.cm,data=resp.2)
#testing block effects

resp.lme<-lme(flux~mm + calibrated + Block, random=~1|Plot,data=resp,na.action = na.exclude)
summary(resp.lm) #no block effect


#nonlinear least squares model
library(nls2)
library(aod)
library(nlstools)
library(lattice)
library(boot)

#soil moisture
plot(flux~calibrated,data=resp.2)
abline(lm(flux~calibrated,data=resp.2))
args(nls)
?nls
#soil moisture models
resp.lm.moisture<-lm(flux~calibrated,data=resp.2)
summary(resp.lm.moisture) #r-sauare = 0.53

lm(logit(flux/10)~calibrated,resp.2)
#run model: logistic growth since the data appears to conform to a logsitic relationship
resp.nonlinearmodel.moisture<- nls2(flux~ theta1/(1 + exp(-(theta2 + theta3*calibrated))),
                                    start=list(theta1 = 10, theta2 = -3.7021  , theta3 = 0.1367 ),data=resp.2, trace=TRUE)
summary(resp.nonlinearmodel.precip)

#model fits
#AIC weights
AIC(resp.nonlinearmodel.moisture,resp.lm.moisture)
x <- c(588.9023, 594.4231)         # stores two AIC values in a vector
delta <- x - min(x)          # AIC differences
L <- exp(-0.5 * delta)       
w <- L/sum(L) 
w #evidence strongest for nonlinear model, for soil moisture
0.94049802/0.05950198

#forggplot trendline
r <- range(resp.2$calibrated,na.rm=TRUE)
resp
xNew <- seq(r[1],r[2])
yNew <- predict(resp.nonlinearmodel.moisture,list(calibrated = xNew))
predictresp.vwc<-data.frame(xNew,yNew)
predictresp.vwc
summary(predictresp)
lines(xNew,yNew)

#precipitation models
#get means for basic linear model
resp.mm<-aggregate(flux~mm,mean,data=resp)
summary(resp.mm)
resp.lm.mm<-lm(flux~mm,data=resp)
summary(resp.lm.mm)
#for model
resp.lm.mm<-lm(flux~mm,data=resp)
summary(resp.lm.mm) #r=square = 0.28
plot(flux~mm,data=resp.mm)
abline(lm(flux~mm,data=resp.2)) 

#find starting parmters
lm(logit(flux/10)~mm,resp)
#run model: logistic growth since the data appears to conform to a logsitic relationship
resp.nonlinear.model.mm<- nls2(flux~ theta1/(1 + exp(-(theta2 + theta3*mm))),
                               start=list(theta1 = 10, theta2 =    -2.116397    , theta3 = 0.005367  ),data=resp, trace=TRUE)

AIC(resp.lm.mm,resp.nonlinear.model.mm)
x <- c(703.0607, 703.7447)         # stores two AIC values in a vector
delta <- x - min(x)          # AIC differences
L <- exp(-0.5 * delta)       
w <- L/sum(L) 
w

#forggplot trendline
r <- range(resp$mm,na.rm=TRUE)
xNew <- seq(r[1],r[2])
yNew <- predict(resp.nonlinear.model.mm,list(mm = xNew))
predictresp.mm<-data.frame(xNew,yNew)
predictresp
summary(predictresp)
lines(xNew,yNew)

#effect sizes

library(effsize)
mean.resp<-subset(resp,subsample=="nominal",na.rm=TRUE)
wettest.resp<-subset(resp,subsample=="extreme.wet",na.rm=TRUE)
driest.resp<-subset(resp,subsample=="extreme.dry",na.rm=TRUE)

#vector
mean.resp.sgs.vec<-as.vector(mean.resp$flux)
wettest.resp.sgs.vec<-as.vector(wettest.resp$flux)
driest.resp.sgs.vec<-as.vector(driest.resp$flux)

cohen.d(treatment.vector, control.vector,pooled=TRUE,paired=FALSE,
        na.rm=TRUE, hedges.correction=FALSE,
        conf.level=0.95,noncentral=FALSE)

#asymmetry signifigance tests

anpp.ag.asymm<-aggregate(x.100  ~ subset,mean,data=anpp.ag) #grass and forb
grass.ag.asymm<-aggregate(gras.100  ~ subset,mean,data=grass.ag) #
forb.ag.asymm<-aggregate(forb.100  ~ subset,mean,data=forb.ag)

#ANPP grass and forb
((101.47 - 70.49) - (70.49 - 51.97))/70.49 #.18, positive asymmetry

#anova
g.f.asym<- function(x) {
  
  diff.g.f <- 70.49 - x
  
  return(diff.g.f)
}
anpp.asym.all<-aggregate(x.100  ~  Plot + subset,g.f.asym,data=anpp.ag)
anpp.asym.all.wet.extremes<-subset(anpp.asym.all,subset==c("extreme.wet"))
anpp.asym.all.wet.extremes$anpp<-abs(anpp.asym.all.wet.extremes$x.100)
par(mfrow=c(2,2)) # init 4 charts in 1 panel
plot(anpp.asym.all.wet.extremes$anpp)
shapiro.test(anpp.asym.all.wet.extremes$anpp) #assumption of normality fair
anpp.asym.all.dry.extremes<-subset(anpp.asym.all,subset==c("extreme.dry"))
anpp.asym.all.dry.extremes$anpp<-abs(anpp.asym.all.dry.extremes$x.100)
shapiro.test(anpp.asym.all.dry.extremes$anpp) #assumption of normality fair 
merge.anpp.extremes<-merge(anpp.asym.all.dry.extremes,anpp.asym.all.wet.extremes,by=c("subset","anpp"),all=TRUE)
anpp.all.anova<-lm(anpp~subset,data=merge.anpp.extremes)
anova(anpp.all.anova)
par(mfrow=c(2,2)) # init 4 charts in 1 panel
plot(anpp.all.anova)

#anova for differences from median
anpp.asym.all.nominal<-subset(anpp.ag,subset==c("nominal"))
shapiro.test(log(anpp.asym.all.nominal$x.100)) #need to log transform to meet assumption
anpp.asym.all.wet<-subset(anpp.ag,subset==c("extreme.wet"))
shapiro.test(log(anpp.asym.all.wet$x.100)) #normality assumptions met
merge.anpp.nominal.wet<-merge(anpp.asym.all.wet,anpp.asym.all.nominal,by=c("subset","x.100"),all=TRUE)
anova.wet<-lm(log(x.100) ~subset,data=merge.anpp.nominal.wet)
anova(anova.wet)
#dry
anpp.asym.all.dry<-subset(anpp.ag,subset==c("extreme.dry"))
shapiro.test(log(anpp.asym.all.dry$x.100)) #normality met
merge.anpp.nominal.dry<-merge(anpp.asym.all.dry,anpp.asym.all.nominal,by=c("subset","x.100"),all=TRUE)
anova.dry<-lm(log(x.100) ~subset,data=merge.anpp.nominal.dry)
anova(anova.dry)

#grass biomass
((87.84 - 65.075) - (65.075 - 50.82))/65.075 # 0.13., positive asymmetry

#anova
grass.asym<- function(x) {
  
  diff.grass <- 65.075 - x
  
  return(diff.grass)
}
anpp.asym.grass<-aggregate(gras.100  ~  Plot + subset,grass.asym,data=grass.ag)
anpp.asym.grass.wet.extremes<-subset(anpp.asym.grass,subset==c("extreme.wet"))
anpp.asym.grass.wet.extremes$anpp<-abs(anpp.asym.grass.wet.extremes$gras.100)
shapiro.test(anpp.asym.grass.wet.extremes$anpp) #assumption of normality fair
anpp.asym.grass.dry.extremes<-subset(anpp.asym.grass,subset==c("extreme.dry"))
anpp.asym.grass.dry.extremes$anpp<-abs(anpp.asym.grass.dry.extremes$gras.100)
shapiro.test(anpp.asym.grass.dry.extremes$anpp) #assumption of normality fair 
merge.grass.extremes<-merge(anpp.asym.grass.dry.extremes,anpp.asym.grass.wet.extremes,by=c("subset","anpp"),all=TRUE)
anpp.grass.anova<-lm(anpp~subset,data=merge.grass.extremes)
anova(anpp.grass.anova)
par(mfrow=c(2,2)) # init 4 charts in 1 panel
plot(anpp.grass.anova)

#anova for differences from median
anpp.asym.grass.nominal<-subset(grass.ag,subset==c("nominal"))
shapiro.test(log(anpp.asym.grass.nominal$gras.100)) #need to log transform to meet assumption
anpp.asym.grass.wet<-subset(grass.ag,subset==c("extreme.wet"))
shapiro.test(log(anpp.asym.grass.wet$gras.100)) #lognormal
merge.anpp.grass.nominal.wet<-merge(anpp.asym.grass.wet,anpp.asym.grass.nominal,by=c("subset","gras.100"),all=TRUE)
anova.grass.wet<-lm(log(gras.100) ~subset,data=merge.anpp.grass.nominal.wet)
anova(anova.grass.wet)
#dry
anpp.asym.grass.nominal<-subset(grass.ag,subset==c("nominal"))
shapiro.test(log(anpp.asym.grass.nominal$gras.100)) #need to log transform to meet assumption
anpp.asym.grass.dry<-subset(grass.ag,subset==c("extreme.dry"))
shapiro.test(log(anpp.asym.grass.dry$gras.100)) #lognormal
merge.anpp.grass.nominal.dry<-merge(anpp.asym.grass.dry,anpp.asym.grass.nominal,by=c("subset","gras.100"),all=TRUE)
anova.grass.dry<-lm(log(gras.100) ~subset,data=merge.anpp.grass.nominal.dry)
anova(anova.grass.dry)

#forb biomass
((13.63 - 5.41) - (5.41 - 0.86))/5.41 #0.68, positive asymmetry

#anova
forb.asym<- function(x) {
  
  diff.forb <- 5.51 - x
  
  return(diff.forb)
}
anpp.asym.forb<-aggregate(forb.100  ~  Plot + subset,forb.asym,data=forb.ag)
anpp.asym.forb.wet.extremes<-subset(anpp.asym.forb,subset==c("extreme.wet"))
anpp.asym.forb.wet.extremes$anpp<-abs(anpp.asym.forb.wet.extremes$forb.100)
shapiro.test(log(anpp.asym.forb.wet.extremes$anpp)) #assumption of normality fair
anpp.asym.forb.dry.extremes<-subset(anpp.asym.forb,subset==c("extreme.dry"))
anpp.asym.forb.dry.extremes$anpp<-abs(anpp.asym.forb.dry.extremes$forb.100)
shapiro.test(anpp.asym.forb.dry.extremes$anpp) #assumptions not met even after log transforming
merge.forb.extremes<-merge(anpp.asym.forb.dry.extremes,anpp.asym.forb.wet.extremes,by=c("subset","anpp"),all=TRUE)
anpp.forb.anova<-lm(log(anpp)~subset,data=merge.forb.extremes) #n.s.
anova(anpp.forb.anova) ##n.s.
#because dry extreme wer enot normal
kruskal.test(anpp~subset, data = merge.forb.extremes)
par(mfrow=c(2,2)) # init 4 charts in 1 panel
plot(anpp.forb.anova)

#anova for differences from median
anpp.asym.forb.nominal<-subset(forb.ag,subset==c("nominal"))
shapiro.test(sqrt(anpp.asym.forb.nominal$forb.100)) #cant make normal
anpp.asym.forb.wet<-subset(forb.ag,subset==c("extreme.wet"))
shapiro.test(log(anpp.asym.forb.wet$forb.100)) #lognormal, but median not normal
merge.anpp.forb.nominal.wet<-merge(anpp.asym.forb.wet,anpp.asym.forb.nominal,by=c("subset","forb.100"),all=TRUE)

#kruskall-wallis test due to lack of normality
kruskal.test(forb.100 ~subset, data = merge.anpp.forb.nominal.wet)
hist(sqrt(anpp.asym.forb.nominal$forb.100))

#dry
anpp.asym.forb.nominal<-subset(forb.ag,subset==c("nominal"))
shapiro.test(log(anpp.asym.forb.nominal$forb.100)) #cant make normal
anpp.asym.forb.dry<-subset(forb.ag,subset==c("extreme.dry"))
shapiro.test(log(anpp.asym.forb.dry$forb.100)) #cant maike normal
merge.anpp.forb.nominal.dry<-merge(anpp.asym.forb.dry,anpp.asym.forb.nominal,by=c("subset","forb.100"),all=TRUE)

#kruskall-wallis test instead
kruskal.test(forb.100 ~subset, data = merge.anpp.forb.nominal.dry)



#% increase
(0.18-0.13)/0.13 #%27
