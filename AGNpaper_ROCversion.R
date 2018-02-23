library(data.table)
library(survival)

###### type 1
#
#tempWriteFile <- paste("~/R/GlCoSy/dataForSubmissions/attd2017/admissionDataDT_T1DM.csv",sep=""); reportingDF<-read.csv(tempWriteFile); reportingDF<-data.table(reportingDF); diabetesType="Type 1 Diabetes. "
#

###### type 1 - hba1c revised code
#
# tempWriteFile <- paste("../GlCoSy/source/admissionDataDT_T1DM_hba1cRevision.csv",sep=""); reportingDF<-read.csv(tempWriteFile); reportingDF<-data.table(reportingDF); diabetesType="Type 1 Diabetes. "
#
#####################################################################################################
###### type 2
#
# 
tempWriteFile <- paste("~/R/GlCoSy/source/admissionDataDT_T2DM.csv",sep=""); reportingDF<-read.csv(tempWriteFile); reportingDF<-data.table(reportingDF); diabetesType="Type 2 Diabetes. "
#

# type 2 - hba1c revised code
#
# tempWriteFile <- paste("../GlCoSy/source/admissionDataDT_T2DM_hba1cRevision.csv",sep=""); reportingDF<-read.csv(tempWriteFile); reportingDF<-data.table(reportingDF); diabetesType="Type 2 Diabetes. "

#####################################################################################################
###### type 2 with drug data
#
# tempWriteFile <- paste("../GlCoSy/source/admissionDataDT_T2DM_withDrugs.csv",sep=""); reportingDF<-read.csv(tempWriteFile); reportingDF<-data.table(reportingDF); diabetesType="Type 2 Diabetes. "; reportingDF$nCBGperAdmission<-reportingDF$nCBGperAdmission.x; reportingDF$yyyy<-reportingDF$yyyy.x; reportingDF$age<-reportingDF$age.x; reportingDF$countHypo3<-reportingDF$countHypo3.x; reportingDF$admissionDurationDays<-reportingDF$admissionDurationDays.x; reportingDF$IQR<-reportingDF$IQR.x
#


reportingDF$eAGyyyyDiff<-reportingDF$eAG-reportingDF$yyyy
# reportingDF$eAGyyyyDiff_inFrame<-reportingDF$eAG_inFrame-reportingDF$yyyy
reportingDF$AGN2<-reportingDF$eAG-reportingDF$medianFirst2CBGs
reportingDF$AGN3<-reportingDF$eAG-reportingDF$medianFirst3CBGs
#
reportingDF$AGN4<-reportingDF$eAG-reportingDF$medianFirst4CBGs


reportingDF$hypo<-ifelse(reportingDF$ID_ADMISSIONhypoEpisodes4.60>0,1,0)

#####################################################################################################
## apply conditions for all analyses
# plotReportingDF<-subset(reportingDF,nCBGperAdmission>=2)   # for admisisons: remove single CBG admissions
plotReportingDF<-subset(reportingDF,admissionDurationDays>=1)
plotReportingDF<-subset(plotReportingDF,nHbA1cValuesInFrame>0)  # for HbA1c perior: remove those without a CBG in 15 month window

plotReportingDF$meanCBG <- mean(as.numeric(unlist(regmatches(plotReportingDF$cbgList, gregexpr("-?\\d+\\.\\d+", plotReportingDF$cbgList)))))
plotReportingDF$sdCBG <- sd(as.numeric(unlist(regmatches(plotReportingDF$cbgList, gregexpr("-?\\d+\\.\\d+", plotReportingDF$cbgList)))))

plotReportingDF<-data.table(plotReportingDF)



plotReportingDF$diabetesDurationYears <- (plotReportingDF$dateplustime1 - plotReportingDF$diagnosisDateUnix) / (60*60*24*365.25)


#plotfilename <- paste("./_Plots.pdf",sep="")
#pdf(plotfilename, width=16, height=9)

## IQR
## most recent HbA1c in range (15 months)
boxplot(plotReportingDF$IQR ~ cut(plotReportingDF$lastHbA1cInFrame,breaks=seq(30,200,10)),las=3,varwidth=T,ylim=c(0,10),plot=T,main="IQR vs last measured HbA1c (x axis)")


plot(plotReportingDF$lastHbA1cInFrame,plotReportingDF$IQR)
cor.test(plotReportingDF$lastHbA1cInFrame,plotReportingDF$IQR)
abline(lm(plotReportingDF$IQR ~ plotReportingDF$lastHbA1cInFrame),col="red")


attdAbstractIQR<-boxplot(plotReportingDF$IQR ~ cut(plotReportingDF$eAGyyyyDiff,breaks=seq(-30,30,2)),las=3,varwidth=T,ylim=c(0,10),plot=T,main="IQR vs AGN ATTD abstract 1 (x axis)")

attdAbstractIQRdecile<-boxplot(plotReportingDF$IQR ~cut(plotReportingDF$eAGyyyyDiff,breaks=quantile(plotReportingDF$eAGyyyyDiff, prob = seq(0, 1, length = 11), type = 5)),plot=T,main="IQR vs AGN ATTD abstract 1 (x axis)",las=3)

cut(plotReportingDF$eAGyyyyDiff,breaks=quantile(plotReportingDF$eAGyyyyDiff, prob = seq(0, 1, length = 11), type = 5))

# for paper - convert U shape to a distance from 0 measure
# gives linear association
distanceFromZero_IQR <- boxplot(plotReportingDF$IQR ~ cut(sqrt(plotReportingDF$eAGyyyyDiff^2),breaks=seq(0,30,1)),las=3,varwidth=T,ylim=c(0,15),plot=T,main="IQR vs AGN distance from 0 (x axis)")


plot(sqrt(plotReportingDF$eAGyyyyDiff^2),plotReportingDF$IQR, pch = 16, cex = 2, col = rgb(0,0,0,0.05, maxColorValue = 1))
# plot(sqrt(plotReportingDF$eAGyyyyDiff^2),plotReportingDF$IQR, pch = 16, cex = 2, col = ifelse(plotReportingDF$eAGyyyyDiff > 0, rgb(1,0,0,0.05, maxColorValue = 1), rgb(0,0,0,0.05, maxColorValue = 1)))

cor.test(sqrt(plotReportingDF$eAGyyyyDiff^2),plotReportingDF$IQR)
abline(lm(plotReportingDF$IQR ~ sqrt(plotReportingDF$eAGyyyyDiff^2)),col="red")

# adjust for nCBG
n = 6 # set minimum number of CBGs per admission to calculate the IQR
moreThan_n_CBG <- plotReportingDF[nCBGperAdmission > n]

# flag first admission for each ID in dataset
moreThan_n_CBG[, c("flagFirstAdmission") := ifelse(admissionNumberFlag == min(admissionNumberFlag), 1, 0),by=.(ID)]

# optional single admission per ID
moreThan_n_CBG <- moreThan_n_CBG[flagFirstAdmission == 1]


# calculate QCD
moreThan_n_CBG$QCD = moreThan_n_CBG$IQR / moreThan_n_CBG$medianGlu

# correlation with last HbA1c in this subset
boxplot(moreThan_n_CBG$IQR ~ cut(moreThan_n_CBG$lastHbA1cInFrame,breaks=seq(30,200,10)),las=3,varwidth=T,ylim=c(0,10),plot=T,main="IQR vs last measured HbA1c (x axis)")

plot(moreThan_n_CBG$lastHbA1cInFrame,moreThan_n_CBG$IQR)
cor.test(moreThan_n_CBG$lastHbA1cInFrame,moreThan_n_CBG$IQR)
abline(lm(moreThan_n_CBG$IQR ~ moreThan_n_CBG$lastHbA1cInFrame),col="red")

  # using QCD rather than IQR
  # correlation with last HbA1c in this subset
  boxplot(moreThan_n_CBG$QCD ~ cut(moreThan_n_CBG$lastHbA1cInFrame,breaks=seq(30,200,10)),las=3,varwidth=T,ylim=c(0,2),plot=T,main="IQR vs last measured HbA1c (x axis)")
  
  plot(moreThan_n_CBG$lastHbA1cInFrame,moreThan_n_CBG$QCD)
  cor.test(moreThan_n_CBG$lastHbA1cInFrame,moreThan_n_CBG$QCD)
  abline(lm(moreThan_n_CBG$QCD ~ moreThan_n_CBG$lastHbA1cInFrame),col="red")

# correlation with first CBG in this subset
boxplot(moreThan_n_CBG$IQR ~ cut(moreThan_n_CBG$yyyy1,breaks=seq(1,28,1)),las=3,varwidth=T,plot=T,main="IQR vs last measured HbA1c (x axis)")

plot(moreThan_n_CBG$yyyy1,moreThan_n_CBG$IQR)
cor.test(moreThan_n_CBG$yyyy1,moreThan_n_CBG$IQR)
abline(lm(moreThan_n_CBG$IQR ~ moreThan_n_CBG$yyyy1),col="red")

    # using QCD rather than IQR
    # correlation with first CBG in this subset
    boxplot(moreThan_n_CBG$QCD ~ cut(moreThan_n_CBG$yyyy1,breaks=seq(1,28,1)),las=3,varwidth=T,plot=T,main="IQR vs yyyy1 (x axis)")
    
    plot(moreThan_n_CBG$yyyy1,moreThan_n_CBG$QCD)
    cor.test(moreThan_n_CBG$yyyy1,moreThan_n_CBG$QCD)
    abline(lm(moreThan_n_CBG$QCD ~ moreThan_n_CBG$yyyy1),col="red")

# correlation of AGN with adjusted HbA1c
#adjusted_IQR <-  moreThan_n_CBG$IQR / moreThan_n_CBG$nCBGperAdmission

distanceFromZero_IQR <- boxplot(moreThan_n_CBG$IQR ~ cut(sqrt(moreThan_n_CBG$eAGyyyyDiff^2),breaks=seq(0,22,1)),las=3,varwidth=T,ylim=c(0,15),plot=T,main="IQR vs AGN distance from 0 (x axis)")

plot(sqrt(moreThan_n_CBG$eAGyyyyDiff^2), moreThan_n_CBG$IQR, pch = 16, cex = 2, col = rgb(0,0,0,0.05, maxColorValue = 1), xlim = c(0, 22))
# plot(sqrt(plotReportingDF$eAGyyyyDiff^2),plotReportingDF$IQR, pch = 16, cex = 2, col = ifelse(plotReportingDF$eAGyyyyDiff > 0, rgb(1,0,0,0.05, maxColorValue = 1), rgb(0,0,0,0.05, maxColorValue = 1)))

cor.test(sqrt(moreThan_n_CBG$eAGyyyyDiff^2), moreThan_n_CBG$IQR)
abline(lm(moreThan_n_CBG$IQR ~ sqrt(moreThan_n_CBG$eAGyyyyDiff^2)),col="red")


# using QCD rather than IQR
    distanceFromZero_IQR <- boxplot(moreThan_n_CBG$QCD ~ cut(sqrt(moreThan_n_CBG$eAGyyyyDiff^2),breaks=seq(0,22,1)),las=3,varwidth=T,ylim=c(0,1),plot=T,main="IQR vs AGN distance from 0 (x axis)")
    
    plot(sqrt(moreThan_n_CBG$eAGyyyyDiff^2), moreThan_n_CBG$QCD, pch = 16, cex = 2, col = rgb(0,0,0,0.05, maxColorValue = 1), xlim = c(0, 22))
    # plot(sqrt(plotReportingDF$eAGyyyyDiff^2),plotReportingDF$IQR, pch = 16, cex = 2, col = ifelse(plotReportingDF$eAGyyyyDiff > 0, rgb(1,0,0,0.05, maxColorValue = 1), rgb(0,0,0,0.05, maxColorValue = 1)))
    
    cor.test(sqrt(moreThan_n_CBG$eAGyyyyDiff^2), moreThan_n_CBG$QCD)
    abline(lm(moreThan_n_CBG$QCD ~ sqrt(moreThan_n_CBG$eAGyyyyDiff^2)),col="red")


## ROC for hypo and mortality
# adjust for nCBG
ROCframe <- moreThan_n_CBG
# ensure that all have at least 1y follow up
ROCframe <- ROCframe[dateplustime1 < (max(ROCframe$dateplustime1) - (60*60*24*365.25))]
# calculate hypo rate
ROCframe$hypoRate <- ROCframe$ID_ADMISSIONhypoEpisodes4.60 / ROCframe$admissionDurationDays
# flag top quartile hypo rate
ROCframe$hypoRate_topQ <- ifelse(ROCframe$hypoRate > (quantile(ROCframe$hypoRate)[4]), 1, 0)
ROCframe$hypoRate_topQ <- ifelse(ROCframe$hypoRate > 1, 1, 0)

# split into train / test
library(caTools)
ROCframe$spl = sample.split(ROCframe$ID, SplitRatio=0.7)

# 
# # mortality
# ROCframe$dead <- ifelse(ROCframe$deathDateUnix>0,1,0)
# ROCframe$inHospitalDeath <- ifelse(ROCframe$dead==1 & 
#                                      (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration) > 0) &
#                                      (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration) < (5*24*60*60)),1,0)
# ROCframe$deadWithin_1Y <- ifelse(ROCframe$dead==1 &
#                                    (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration) > 0) &
#                                    (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration) < (365.25*24*60*60)),1,0)


# ROCframe <- subset(ROCframe, yyyy1 > 3.9)
 
library("ROCR")    
roc_function <- function(outcomeParameter, testParameter, inputFrame) {
  
  mylogit <- glm(outcomeParameter ~ (testParameter + age), data = inputFrame, family = "binomial")
  
  summary(mylogit)     
  prob=predict(mylogit,type=c("response"))
  
  pred <- prediction(prob, outcomeParameter)    
  perf <- performance(pred, measure = "tpr", x.measure = "fpr")     
  plot(perf, col=rainbow(7), main="ROC curve Admissions", xlab="Specificity", 
       ylab="Sensitivity")    
  abline(0, 1) #add a 45 degree line
  
  auc = performance(pred, "auc")
  print(auc)
}

roc_function(ROCframe$anyHypo3, sqrt(ROCframe$eAGyyyyDiff^2), ROCframe)
roc_function(ROCframe$anyHypo3, ROCframe$yyyy1, ROCframe)
roc_function(ROCframe$anyHypo3, ROCframe$eAG, ROCframe)

roc_function(ROCframe$anyHypo4, sqrt(ROCframe$eAGyyyyDiff^2), ROCframe)
roc_function(ROCframe$anyHypo4, ROCframe$yyyy1, ROCframe)
roc_function(ROCframe$anyHypo4, ROCframe$eAG, ROCframe)

roc_function(ROCframe$hypoRate_topQ, sqrt(ROCframe$eAGyyyyDiff^2), ROCframe)
roc_function(ROCframe$hypoRate_topQ, ROCframe$yyyy1, ROCframe)
roc_function(ROCframe$hypoRate_topQ, ROCframe$lastHbA1cInFrame, ROCframe)

roc_function(ROCframe$inHospitalDeath, sqrt(ROCframe$eAGyyyyDiff^2), ROCframe)
roc_function(ROCframe$inHospitalDeath, ROCframe$yyyy1, ROCframe)
roc_function(ROCframe$inHospitalDeath, ROCframe$eAG, ROCframe)

roc_function(ROCframe$deadWithin_1Y, sqrt(ROCframe$eAGyyyyDiff^2), ROCframe)
roc_function(ROCframe$deadWithin_1Y, ROCframe$yyyy1, ROCframe)
roc_function(ROCframe$deadWithin_1Y, ROCframe$lastHbA1cInFrame, ROCframe)


#
hypoRateFrame <- plotReportingDF[admissionDurationDays > 1]
hypoRateFrame <- hypoRateFrame[yyyy1 > 3]

hypoRateFrame$hypoRate <- hypoRateFrame$countHypo3 / hypoRateFrame$admissionDurationDays

point5ToOneHypo <- ifelse(hypoRateFrame$hypoRate > 0.5 & hypoRateFrame$hypoRate <=1, 1, 0)
oneHypo <- ifelse(hypoRateFrame$hypoRate > 0 & hypoRateFrame$hypoRate <=1, 1, 0)
topQuart <- ifelse(hypoRateFrame$hypoRate > (quantile(hypoRateFrame$hypoRate)[4]), 1, 0)


nHypo <- topQuart

roc_function(nHypo, sqrt(hypoRateFrame$eAGyyyyDiff^2), hypoRateFrame)
roc_function(nHypo, hypoRateFrame$yyyy1, hypoRateFrame)
roc_function(nHypo, hypoRateFrame$eAG, hypoRateFrame)




