## Set your local directory
setwd("G:/Lavoro/Haakon/Dataset4Courses/Rdata/")

## Load original data
mydata=read.delim("datacourse.xyz", header = TRUE, sep = "\t")
mydata[mydata==-9999 | mydata==-99999]=NA
mydata = na.omit(mydata)

## Load the results
load("LushanMIonlyRes.Rdata")
load("LushanLSEonlyRes.Rdata")

## Let's extract the landslide intensities for the two models

fitted.intensity.MIonly = Lushan.MIonly$summary.fitted.values ## This is where INLA stores the fitted values
head(fitted.intensity.MIonly) ## As you may have noticed, INLA stores 6 parameters describing the intensity distribution for each pixel
fitted.intensity.LSEonly = Lushan.LSEonly$summary.fitted.values
head(fitted.intensity.LSEonly)

## Plotting fitted values, we will do this operation in GIS after preparing the data here.
## To prepare the data for SAGA GIS, we juct need to extract the coordinates from the original matrix 
## and paste the mean intensity columns fot the two models right after. We do it as follows:  
intensity.data4GIS = cbind(mydata$X,mydata$Y,fitted.intensity.MIonly$mean,fitted.intensity.LSEonly$mean)
head(intensity.data4GIS) ## Let's see what we just did. We lost the headers in the process, so we can fix it by:
colnames(intensity.data4GIS) <- c("X","Y","MeanIntensityMIonly","MeanIntensityLSEonly")
head(intensity.data4GIS) ## Let's check
write.table(intensity.data4GIS,"intensity4maps.txt", sep = "\t", col.names = TRUE, row.names = FALSE) 
## We will now operate in GIS to create landslide intensity maps at pixel scale
## However, as mentioned during the seminar in Day1, one of the best properties of a Point Process is that intensities 
## are not dependent on the mapping unit we chose. Thus, we can immediately project from pixel to any other 
## spatial unit. We will demonstrate this, creating intensity maps for Slope Units, summing all values at pixels 
## contained in each Slope Unit.

# install.packages("data.table")
require(data.table)

SU.table = cbind(mydata$LandslideCounts, fitted.intensity.MIonly$mean, fitted.intensity.LSEonly$mean, mydata$SU)   
head(SU.table)
colnames(SU.table) <- c("ObservedCounts","IntensityMIonly","IntensityLSEonly","SU.ID")
SU.count = data.table(SU.table, key = "SU.ID")
SU.count.sum = SU.count[,list(ObservedCounts=sum(ObservedCounts),IntensityMIonly=sum(IntensityMIonly),
                              IntensityLSEonly=sum(IntensityLSEonly)), by="SU.ID"]
head(SU.count.sum)
write.table(SU.count.sum,"intensity4SU.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

## The landslide intensity is a recent concept in the literature, thus we can explore it according to some metrics 
## we have proposed and also according to more traditional metrics borrowed from landslide susceptibility studies.
## Do we all know what landslide susceptibility means? If yes, let's move forward. If not, 5-10 mins summary.

susceptibility.MIonly = 1- exp(- fitted.intensity.MIonly$mean) ## This relation converts intensities into classic susceptibilities
## Because a Poisson distribution is hierarchically above a Bernoulli, we can always extract the probability of at 
## least one lanslide from the Poisson, thus moving back to the classic susceptibilty concept. 
susceptibility.LSEonly = 1- exp(- fitted.intensity.LSEonly$mean)
head(susceptibility.MIonly)
head(susceptibility.LSEonly)

binarydata = mydata$LandslideCounts ## As we will use some binary metrics, we also need to convert the landslide counts
## into landslide presence/absence data. We do this as follows:
binarydata[binarydata>0] = 1

## Once the two transformations are done, we can use the Receiver Characteristic Operating Curves to assess the quality of
## our two models
install.packages("pROC")
require(pROC)

ROCtest.MIonly = roc(binarydata~susceptibility.MIonly)
AUCtest.MIonly = data.frame((ROCtest.MIonly[[9]])[1])
View(AUCtest.MIonly)

ROCtest.LSEonly = roc(binarydata~susceptibility.LSEonly)
AUCtest.LSEonly = data.frame((ROCtest.LSEonly[[9]])[1])
View(AUCtest.LSEonly)

## Let's now take a look at the result from the perspective of count data. 
## It is worth mentioning that count data at the pixel level do not differ significantly from a binary situation. 
## For this reason, we suggest using ROC curves at the pixel level and the following count metrics for any other mapping unit.

## Because we will use Slope Units as the only example here, we need to compute how many landslide count fall in a given slope
## unit, both for the original landslides and the estimated landslides. Let's recall what we did between lines 35 and 44.
## We can proceed in the same way:

# install.packages("data.table")
require(data.table)

SU.table = cbind(mydata$LandslideCounts, fitted.intensity.MIonly$mean, fitted.intensity.LSEonly$mean, mydata$SU)   
head(SU.table)
colnames(SU.table) <- c("ObservedCounts","IntensityMIonly","IntensityLSEonly","SU.ID")
SU.count = data.table(SU.table, key = "SU.ID")
SU.count.sum = SU.count[,list(ObservedCounts=sum(ObservedCounts),IntensityMIonly=sum(IntensityMIonly),
                              IntensityLSEonly=sum(IntensityLSEonly)), by="SU.ID"]
head(SU.count.sum)

## This operation ensured the aggregation of the pixel counts at the slope unit level.
## Once we have the new aggregated data, we can easily plot it showing the relation between observed landlsides and estimated ones

plot(SU.count.sum$ObservedCounts,SU.count.sum$IntensityMIonly, main = "Slope Units", xlab="Observed Landslide Counts", ylab = "Estimated Landslide Counts",  pch = 19, xlim = c(0,50), ylim = c(0,50), col="black", cex.lab=1.5, cex = 1.5, cex.axis = 1.5)
points(SU.count.sum$ObservedCounts, SU.count.sum$IntensityLSEonly, pch = 19, col = "red", cex.lab=1.5, cex = 1.5, cex.axis = 1.5)
lines(seq(0,50, length=1000),qpois(0.025,lambda=seq(0,50,length=1000)), col = "lightgrey")
lines(seq(0,50, length=1000),qpois(0.975,lambda=seq(0,50,length=1000)), col = "lightgrey")

## There are several ways to compute numerical metrics for count data. However, the appropiate choice may depend on the 
## specific dataset. For instance, for aggregated count data where the majority of values is greater than 5, we could use a 
## standard R squared. However, being the choice conditioned by the dataset you will be working on, in this course, we opted 
## to leave the comparison between observed and estimated counts only via a visual representation.

## Now, lets see how to handle cross-validation procedures with INLA
## In many cases, while working in R, one can calibrate the model, store the corresponding object and use the function predict
## to validate the model over an unknown dataset. In INLA, cross-validation is even simpler!
## We simply have to assign Not-A-Number (NA, in R) to any row of the dependent variable that we require INLA to predict.  
## Let's see one single CV example where we want INLA to use 90% of the data for calibration and 10% of the data, selected at 
## random, to validate.
## Let's reprocess the data as we did for the fitting during the morning. From line 101 to 131 it's the same code as before.  

require(INLA)

data_scaled=mydata
vars2scale=c("AvgPrecipitation","AvgTempDiff","Dist2Faults","Dist2GeoBoundaries","Elevation","MILushan","PlanCur","ProfCur","RSP")
data_scaled[,vars2scale]=apply(data_scaled[,vars2scale],2,scale)

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
hyper.fix = list(theta1 = list(initial = log(100), fixed = TRUE))
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

## However, here we will assign NA to a random subset to perform CV in INLA

n = n.pixels
m = floor(n/10)
set.seed(1)

ind1 = sample(c(1:n),m)
y.count[ind1] = NA

## Here we merge the information of counts, the intercept to be passed to INLA, Lithology, the 10 Macroseismic Intensity classes and the 
## Slope Units that will be used afterwards to calculate the Laternt Spatial Effect
df = data.frame(y=y.count, intercept = 1, Lithology = covar.inla$Lithology, MILev = covar.inla$MILev, SlopeUnitID = covar.inla$SU) # off.hc = off.hc,

## We define the formula we want to pass to INLA. This includes the fixed (or linear) effects, the intercept, Lithology and 
## 10 Microseismic Intensity classes
formula.Lushan.MIonly = y~ -1 + intercept + Xm +
  f(Lithology,model="iid",hyper=hyper.iid, constr=T) +
  f(MILev,model="rw1",hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal = 1E-4)

## We finally run INLA
Lushan.MIonly.CV=inla(formula.Lushan.MIonly, family="poisson",data=c(as.list(df), list(Xm=Xm)),
                   control.fixed=list(prec=1),
                   E=offset,num.threads=2,
                   control.inla = list(int.strategy='eb'),
                   control.mode=list(restart=T, theta=my.init),
                   control.predictor=list(link = 1, compute = TRUE),
                   verbose = TRUE)

## Let's now extract the predicted values at rows were we assigned NA. As you will see, now they are populated with intensity values.

intensity.MI.only.CV = as.data.frame(Lushan.MIonly.CV$summary.fitted.values$mean[ind1])
count.data.CV = as.data.frame(mydata$LandslideCounts[ind1])
View(intensity.MI.only.CV)
## From here we can reuse the code above to measure model performance on a CV subset.
## Do you want to do it or would you like to talk about your own data and figure out how spatial statistics can help you
## to answer your questions. If you want to check CV result, let's write the code together. Do not worry, we just need to 
## copy & paste few lines from above but it's a good feedback for us already to assess if you all understood some key steps.

## If you want to discuss your own data instead, raise your hand and come forward.


