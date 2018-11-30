# To get a better understanding of the concepts explained below I suggest reading the following articles 
# for the following reasons:
#   
# 1. Here we introduce the landslide intensity concept using R-INLA for the first time.
# 
# Lombardo L, Opitz T, Huser R. (2018) Point process-based modeling of multiple debris flow landslides 
# using INLA: an application to the 2009 Messina disaster. Stochastic Environmental Research and Risk 
# Assessment, 32, pp 2179-2198.  
# 
# 2. Here we explain in details how the code works and what are the output one can generate.
# 
# Lombardo L, Opitz T, Huser R. (2018) Numerical recipes for landslide spatial prediction with R-INLA: 
# a step-by- step tutorial}. Chapter accepted for publication in the book 
# "Spatial modeling in GIS and R for Earth and Environmental Sciences", Elsevier.
# 
# 3. Here we give a better example of how powerful the Latent Spatial Effect can be in retrieving the 
# spatial signal of a trigger.
# 
# Lombardo L, Bakka H, Tanyas H, van Westen CJ, Mai PM, Huser R (2018) Geostatistical modeling to capture 
# seismic-shaking patterns from earthquake-induced landslides. Submitted to JGR - Earth Surface.

#########################################################
#################### NOW LET'S START ####################
#########################################################

# Please uncomment and install these packages before running the code
# install.packages("INLA")
install.packages("fields")
install.packages("MASS")
install.packages("maptools")
install.packages("spdep")

## Preparatory part using GIS data to create a matrix where all the adjacent Slope Units are reported

## Set your local directory... This is my local path, please change it according to your local directories
setwd("G:/Lavoro/Haakon/Dataset4Courses/VectorData/")

## Generate adjacency graph of slope units in file file adjgraph.txt

library(maptools)
SU=readShapeSpatial("SU.shp")
View(SU)

library(spdep)
nb2INLA("G:/Lavoro/Haakon/Dataset4Courses/Rdata/adjgraph.txt",poly2nb(SU,queen=F,row.names=SU$SU_ID))

## Let's explore what an adjacency matrix is
adjacency = read.delim("G:/Lavoro/Haakon/Dataset4Courses/Rdata/adjgraph.txt")
View(adjacency)

library(INLA)
library(fields)
library(MASS)

## Set your local directory
setwd("G:/Lavoro/Haakon/Dataset4Courses/Rdata/")

## Load the data 
mydata=read.delim("datacourse.xyz", header = TRUE, sep = "\t")
mydata[mydata==-9999 | mydata==-99999]=NA
mydata = na.omit(mydata) ## This operation ensures that any missing values coming from raster in GIs will be removed. 
## mydata has no missing values to begin with but I usually add this to be on the safe side. 

## Rescale the covariates in the same unitless scale. We opt here for a mean zero, unit variance rescaling. 
data_scaled=mydata
vars2scale=c("AvgPrecipitation","AvgTempDiff","Dist2Faults","Dist2GeoBoundaries","Elevation","MILushan","PlanCur","ProfCur","RSP")
data_scaled[,vars2scale]=apply(data_scaled[,vars2scale],2,scale)
plot(mydata$Dist2Faults,data_scaled$Dist2Faults)

## We merge together the rescaled covariates (originally continuous) and the categorical ones.
covar.inla=data_scaled[,c("LandslideCounts","AvgPrecipitation","AvgTempDiff","Dist2Faults","Dist2GeoBoundaries","Elevation","MILushan","PlanCur","ProfCur","RSP",
                          "Lithology","SU")]
## We add one column for the intercept and prepare 10 classes to assess nonlinear effects of the Microseismic Intensity
covar.inla=cbind(intercept=1,covar.inla)
covar.inla$MILev=inla.group(covar.inla$MILushan,n=10, method = "cut")
table(covar.inla$MILev)
## We create a set of hyperparameters 
hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)))
hyper.iid = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)))
#hyper.fix = list(theta1 = list(initial = log(100), fixed = TRUE))
my.init = NULL

## The target of our point process framework is to produce landslide intensity models. Thus, we need to express the number of landslides and their relation 
## with the pixel we adopt.
y.count=mydata$LandslideCounts
n.pixels=nrow(mydata)
## area.pixel=rep(100^2,n.pixels) if want the units to be landlides per squared meters then we opt for this line otherwise
## area.pixel=rep(0.1^2,n.pixels) ## this sets the units to be landslides per squared kilometer
area.pixel=rep(1^2,n.pixels) ## this will indicate landslides per pixel instead 
offset=area.pixel
n.landslides=sum(y.count)
avg.global=n.landslides/(n.pixels*area.pixel)

## We create a matrix containing all the linear effects.
Xm = as.data.frame(covar.inla[ , c("AvgPrecipitation","AvgTempDiff","Dist2Faults","Dist2GeoBoundaries","Elevation", "PlanCur","ProfCur","RSP")])
Xm = as.matrix(Xm) 

## Here we merge the information of counts, the intercept to be passed to INLA, Lithology, the 10 Macroseismic Intensity classes and the 
## Slope Units that will be used afterwards to calculate the Laternt Spatial Effect
df = data.frame(y=y.count, intercept = 1, Lithology = covar.inla$Lithology, MILev = covar.inla$MILev, SlopeUnitID = covar.inla$SU) # off.hc = off.hc,

## We define the formula we want to pass to INLA. This includes the fixed (or linear) effects, the intercept, Lithology and 
## 10 Microseismic Intensity classes
formula.Lushan.MIonly = y~ -1 + intercept + Xm +
  f(Lithology,model="iid",hyper=hyper.iid, constr=T) +
  f(MILev,model="rw1",hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal = 1E-4)
  
## We finally run INLA
Lushan.MIonly=inla(formula.Lushan.MIonly, family="poisson",data=c(as.list(df), list(Xm=Xm)),
                           control.fixed=list(prec=1),
                           E=offset,num.threads=2,
                           control.inla = list(int.strategy='eb'),
                           control.mode=list(restart=T, theta=my.init),
                           control.predictor=list(link = 1, compute = TRUE),
                           verbose = TRUE)

## In this example, we will use the same model as before, but we will remove the Microseismic Intensity and use the Latent Spatial 
## Effect instead
formula.Lushan.LSEonly = y~ -1 + intercept + Xm +
  f(Lithology,model="iid",hyper=hyper.iid, constr=T)+
  f(SlopeUnitID, model="besag",graph="adjgraph.txt",hyper=hyper.rw, constr=T,scale.model = TRUE, diagonal = 1E-4)

Lushan.LSEonly=inla(formula.Lushan.LSEonly,family="poisson",data=c(as.list(df), list(Xm=Xm)),
                          control.fixed=list(prec=1),
                          E=offset,num.threads=3,
                          control.inla = list(int.strategy='eb'),
                          control.mode=list(restart=T, theta=my.init),
                          control.predictor=list(link = 1, compute = TRUE),
                          verbose = TRUE)

##############################################
#### Now we will explore the INLA results ####
##############################################

## First we will check the model with the Microintensity

summary(Lushan.MIonly) ## here we can get a quick idea of the information about the model (computational time and structure)
str(Lushan.MIonly,1)
Lushan.MIonly$summary.fixed ## Information on covariates treated linearly
Lushan.MIonly$summary.random ## Information on covariate treated nonlinearly but as pure categorical 
Lushan.MIonly$summary.random$MILev ## Information on covariate treated nonlinearly but with internal dependence 
fitted.intensity.MIonly = Lushan.MIonly$summary.fitted.values ## This is where INLA stores the fitted values
quilt.plot(mydata$X, mydata$Y, fitted.intensity.MIonly$mean)

head(fitted.intensity.MIonly) ## Here we can check the information stored by INLA related to intensities
plot(Lushan.MIonly) ## This is a basic plot command in INLA but please note that some of these plot will not be applicable 
## to this specific experiment. 

## Similarly, we can see the result of the model with the Latent Spatial Effect

str(Lushan.LSEonly$summary.random,1)
Lushan.LSEonly$summary.fixed
Lushan.LSEonly$summary.random$Lithology
head(Lushan.LSEonly$summary.random$SlopeUnitID)
fitted.intensity.LSEonly = Lushan.LSEonly$summary.fitted.values
head(fitted.intensity.LSEonly)

## This is a particularly important step because we can assess the contribution of every covariate to the final model
## For example, we can check which covariate has a significant contribution as follows.
## Please note the difference with respect to the previous steps where all covariates where listed with the $summary.fixed command
extractSignificantEffects=function(results.df){
  print(results.df[sign(results.df$`0.025quant`)==sign(results.df$`0.975quant`),])
}
extractSignificantEffects(Lushan.LSEonly$summary.fixed)
extractSignificantEffects(Lushan.LSEonly$summary.random$Lithology)

## Another way to check the covariate effects is to plot their values and relative credible intervals. 
## If we refer to one of the fixed effects, We can do it as follows: 
plot(Lushan.MIonly$summary.fixed$mean[2],ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
text(1,0.6,"Average Precipitation", font = 2, cex = 1.4)
points(Lushan.MIonly$summary.fixed$`0.025quant`[2],col="blue",pch=18, cex = 2)
points(Lushan.MIonly$summary.fixed$`0.975quant`[2],col="blue",pch=18, cex = 2)
segments(1,Lushan.MIonly$summary.fixed$`0.025quant`[2],1,Lushan.MIonly$summary.fixed$`0.975quant`[2],lwd=2, lty = 3, col = "blue")
abline(h=0,lty=2,lwd=2,col="gray50")


## However, if we want a comprehensive view of all fixed effect, we can plot them in a subpaneled structure as follows:
pdf("testFig.pdf", width = 14, height = 7)
par(mfrow=c(2,4))
## Panel 1
plot(Lushan.MIonly$summary.fixed$mean[2], xlab="", xaxt="n", ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
text(1,0.6,"Average Precipitation", font = 2, cex = 1.4)
points(Lushan.MIonly$summary.fixed$`0.025quant`[2],col="blue",pch=18, cex = 2)
points(Lushan.MIonly$summary.fixed$`0.975quant`[2],col="blue",pch=18, cex = 2)
segments(1,Lushan.MIonly$summary.fixed$`0.025quant`[2],1,Lushan.MIonly$summary.fixed$`0.975quant`[2],lwd=2, lty = 3, col = "blue")
abline(h=0,lty=2,lwd=2,col="gray50")
## Panel 2
plot(Lushan.MIonly$summary.fixed$mean[3], xlab="", xaxt="n",  ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
text(1,0.6,"Average Temperature Difference", font = 2, cex = 1.4)
points(Lushan.MIonly$summary.fixed$`0.025quant`[3],col="blue",pch=18, cex = 2)
points(Lushan.MIonly$summary.fixed$`0.975quant`[3],col="blue",pch=18, cex = 2)
segments(1,Lushan.MIonly$summary.fixed$`0.025quant`[3],1,Lushan.MIonly$summary.fixed$`0.975quant`[3],lwd=2, lty = 3, col = "blue")
abline(h=0,lty=2,lwd=2,col="gray50")
## Panel 3
plot(Lushan.MIonly$summary.fixed$mean[4], xlab="", xaxt="n", ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
text(1,0.6,"Distance to Faults", font = 2, cex = 1.4)
points(Lushan.MIonly$summary.fixed$`0.025quant`[4],col="blue",pch=18, cex = 2)
points(Lushan.MIonly$summary.fixed$`0.975quant`[4],col="blue",pch=18, cex = 2)
segments(1,Lushan.MIonly$summary.fixed$`0.025quant`[4],1,Lushan.MIonly$summary.fixed$`0.975quant`[4],lwd=2, lty = 3, col = "blue")
abline(h=0,lty=2,lwd=2,col="gray50")
## Panel 4
plot(Lushan.MIonly$summary.fixed$mean[5], xlab="", xaxt="n", ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
text(1,0.6,"Distance to GeoBoundaries", font = 2, cex = 1.4)
points(Lushan.MIonly$summary.fixed$`0.025quant`[5],col="blue",pch=18, cex = 2)
points(Lushan.MIonly$summary.fixed$`0.975quant`[5],col="blue",pch=18, cex = 2)
segments(1,Lushan.MIonly$summary.fixed$`0.025quant`[5],1,Lushan.MIonly$summary.fixed$`0.975quant`[5],lwd=2, lty = 3, col = "blue")
abline(h=0,lty=2,lwd=2,col="gray50")
## Panel 5
plot(Lushan.MIonly$summary.fixed$mean[6], xlab="", xaxt="n", ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
text(1,0.6,"Elevation", font = 2, cex = 1.4)
points(Lushan.MIonly$summary.fixed$`0.025quant`[6],col="blue",pch=18, cex = 2)
points(Lushan.MIonly$summary.fixed$`0.975quant`[6],col="blue",pch=18, cex = 2)
segments(1,Lushan.MIonly$summary.fixed$`0.025quant`[6],1,Lushan.MIonly$summary.fixed$`0.975quant`[6],lwd=2, lty = 3, col = "blue")
abline(h=0,lty=2,lwd=2,col="gray50")
## Panel 6
plot(Lushan.MIonly$summary.fixed$mean[7], xlab="", xaxt="n", ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
text(1,0.6,"Planar Curvature", font = 2, cex = 1.4)
points(Lushan.MIonly$summary.fixed$`0.025quant`[7],col="blue",pch=18, cex = 2)
points(Lushan.MIonly$summary.fixed$`0.975quant`[7],col="blue",pch=18, cex = 2)
segments(1,Lushan.MIonly$summary.fixed$`0.025quant`[7],1,Lushan.MIonly$summary.fixed$`0.975quant`[7],lwd=2, lty = 3, col = "blue")
abline(h=0,lty=2,lwd=2,col="gray50")
## Panel 7
plot(Lushan.MIonly$summary.fixed$mean[8], xlab="", xaxt="n", ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
text(1,0.6,"Profile Curvature", font = 2, cex = 1.4)
points(Lushan.MIonly$summary.fixed$`0.025quant`[8],col="blue",pch=18, cex = 2)
points(Lushan.MIonly$summary.fixed$`0.975quant`[8],col="blue",pch=18, cex = 2)
segments(1,Lushan.MIonly$summary.fixed$`0.025quant`[8],1,Lushan.MIonly$summary.fixed$`0.975quant`[8],lwd=2, lty = 3, col = "blue")
abline(h=0,lty=2,lwd=2,col="gray50")
## Panel 8
plot(Lushan.MIonly$summary.fixed$mean[9], xlab="", xaxt="n", ylim = c(-0.6,0.6), ylab="Linear predictor",pch=19, col="black", cex.lab=1.5, cex = 2)
text(1,0.6,"Relative Slope Position", font = 2, cex = 1.4)
points(Lushan.MIonly$summary.fixed$`0.025quant`[9],col="blue",pch=18, cex = 2)
points(Lushan.MIonly$summary.fixed$`0.975quant`[9],col="blue",pch=18, cex = 2)
segments(1,Lushan.MIonly$summary.fixed$`0.025quant`[9],1,Lushan.MIonly$summary.fixed$`0.975quant`[9],lwd=2, lty = 3, col = "blue")
abline(h=0,lty=2,lwd=2,col="gray50")
dev.off()


## With a similar concept in mind, we can asses the effect of non linear covariates, by plotting each regression coefficient
## First of all, we need to load a file where each Lithotype name recorded against the numeric ID we used in the calculations
Lushan.MIonly$summary.random$Lithology ## Let's see directly what I mean
LitoNames = read.delim("Lithologies.txt", sep = "\t", header = T)
View(LitoNames)
plot(1:5+0.22,Lushan.MIonly$summary.random$Lithology$mean,xlab="Lithology",ylab="Linear predictor",lwd=2,ylim=c(-2.7,2.7),pch=19,xaxt="n", col="black", cex.lab=1.5, cex = 2)
axis(1,at=1:5+0.22,labels=c("D","P","Pt","Q","T"), cex.axis = 1.5)
points(1:5+0.22,Lushan.MIonly$summary.random$Lithology$"0.025quant",lwd=2,col="blue",pch=18, cex = 2)
points(1:5+0.22,Lushan.MIonly$summary.random$Lithology$"0.975quant",lwd=2,col="blue",pch=18, cex = 2)
segments(1:5+0.22,Lushan.MIonly$summary.random$Lithology$"0.025quant",1:5+0.22,Lushan.MIonly$summary.random$Lithology$"0.975quant",lwd=2, lty = 3, col = "blue")
abline(h=0,lty=2,lwd=2,col="gray50")
text(1:5+0.15,y=Lushan.MIonly$summary.random$Lithology$mean,as.character(LitoNames$Class),srt=90, font = 2, cex = 1.75)

## Differently from a native categorical covariate, random effect (or nonlinear) with internal dependence should be plotted as follows:
# truehist(mydata$MILushan, xlab = "Microseismic Intensity") ## First, let's understand what I meant before,
covar.inla$MILev=inla.group(covar.inla$MILushan,n=10, method = "cut")  ## In line 46, we cut the MI using 10 classes
## However, the classes are cut from the rescaled version of the Microseismic Intensity
## To get the values of the ten classes in the original MI scale, we can do as follows:
MI.cut.original = inla.group(mydata$MILushan,n=10, method = "cut") ## This command will reclassify in 10 classes, Then,
MI.cut.original.val = as.data.frame(table(MI.cut.original)) ## This command will show the 10 classes and how many pixels fell in each one.
View(MI.cut.original.val) ## Let's check them
table(covar.inla$MILev)
xvalues = MI.cut.original.val$MI.cut.original ## We here extract the column where the 10 MI values are contained
View(xvalues) ## Let's check
xlabels = as.character(xvalues) 
xnumbers = as.numeric(xlabels) ## For internal reasons R handles the values in a funny way, so we transform them back to numbers
## from the labels we extracted. If we would have passed the xvalues directly the plot would have been very different.

## Now we can finally plot the MI random effect
pdf("G:/RINLAcourseSheffield/SheffieldCourse_LL/Figures/testMILev.pdf",width = 5, height = 5)
plot(xnumbers, Lushan.MIonly$summary.random$MILev$mean,type="l",
     lwd=3,xaxt="n",xlab="Microseismic Intensity",ylab="Linear predictor",ylim=c(-0.5,0.5),col="black", cex.lab=1.5,
     cex = 1.5, cex.axis = 1.5)
axis(1,at=xnumbers,
     labels=xlabels, cex.axis = 1.5)
lines(xnumbers, Lushan.MIonly$summary.random$MILev$"0.025quant",
      lwd=3,col="blue",lty=3)
lines(xnumbers, Lushan.MIonly$summary.random$MILev$"0.975quant",
      lwd=3,col="blue",lty=3)
abline(h=0,lty=2,lwd=3,col="gray50")
dev.off()
## You may have noticed that the xvalues are not equally spaces. This is due to the way inla.group() operates
## What it does is to initially create equally spaced classes. However, because each class contains its own distribution of values
## inla.group uses the median and return the coefficient for that specific value, making the plot irregular. 
## An update of the inla.group function is on its way. The new one will use the mean instead of the median, which in turn will
## make the x-axis of plot exactly equally spaced.  

## This part will conculde the analyses of the covariate effects with the exception of the Latent Spatial Effect
## Visualize the Latent Spatial Effect requires an additional step, so let's do it.
setwd("G:/Lavoro/Haakon/Dataset4Courses/Rdata/")
LSE = Lushan.LSEonly$summary.random$SlopeUnitID
View(LSE) ## as you can see for each slope unit we computed a LSE value. Thus, we can nicely visualize it in GIS as follows:
write.table(LSE, "LSE4course.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

## PLOT THE WRONG WAY THE LSE

## Before concluding this part, let's save the results. You can do it as follows:
save(Lushan.MIonly, file = "LushanMIonlyRes.Rdata")
save(Lushan.LSEonly, file = "LushanLSEonlyRes.Rdata")
