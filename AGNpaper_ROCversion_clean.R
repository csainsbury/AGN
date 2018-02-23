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

plotReportingDF<-data.table(plotReportingDF)
plotReportingDF$uniqueAdmissionID <- seq(1, nrow(plotReportingDF), 1)
plotReportingDF[, c("meanCBG") := mean(as.numeric(unlist(regmatches(cbgList, gregexpr("-?\\d+\\.\\d+", cbgList))))) ,by=.(uniqueAdmissionID)]
plotReportingDF[, c("sdCBG") := sd(as.numeric(unlist(regmatches(cbgList, gregexpr("-?\\d+\\.\\d+", cbgList))))) ,by=.(uniqueAdmissionID)]




plotReportingDF$diabetesDurationYears <- (plotReportingDF$dateplustime1 - plotReportingDF$diagnosisDateUnix) / (60*60*24*365.25)

##
## correlation plots

# adjust for nCBG
n = 6 # set minimum number of CBGs per admission to calculate the IQR
moreThan_n_CBG <- plotReportingDF[nCBGperAdmission > n]

# flag first admission for each ID in dataset
moreThan_n_CBG[, c("flagFirstAdmission") := ifelse(admissionNumberFlag == min(admissionNumberFlag), 1, 0),by=.(ID)]

# optional single admission per ID
moreThan_n_CBG <- moreThan_n_CBG[flagFirstAdmission == 1]


# calculate CV
moreThan_n_CBG$CV = moreThan_n_CBG$sdCBG / moreThan_n_CBG$meanCBG

# using CV rather than IQR
# correlation with last HbA1c in this subset
boxplot(moreThan_n_CBG$CV ~ cut(moreThan_n_CBG$lastHbA1cInFrame,breaks=seq(30,200,10)),las=3,varwidth=T,ylim=c(0,0.5),plot=T,main="CV vs last measured HbA1c (x axis)")

plot(moreThan_n_CBG$lastHbA1cInFrame,moreThan_n_CBG$CV)
cor.test(moreThan_n_CBG$lastHbA1cInFrame,moreThan_n_CBG$CV)
abline(lm(moreThan_n_CBG$CV ~ moreThan_n_CBG$lastHbA1cInFrame),col="red")

# using CV rather than IQR
# correlation with first CBG in this subset
boxplot(moreThan_n_CBG$CV ~ cut(moreThan_n_CBG$yyyy1,breaks=seq(1,28,1)),las=3,varwidth=T,plot=T,main="IQR vs yyyy1 (x axis)")

plot(moreThan_n_CBG$yyyy1,moreThan_n_CBG$CV)
cor.test(moreThan_n_CBG$yyyy1,moreThan_n_CBG$CV)
abline(lm(moreThan_n_CBG$CV ~ moreThan_n_CBG$yyyy1),col="red")

# correlation of AGN with adjusted HbA1c
# using CV rather than IQR
distanceFromZero_IQR <- boxplot(moreThan_n_CBG$CV ~ cut(sqrt(moreThan_n_CBG$eAGyyyyDiff^2),breaks=seq(0,22,1)),las=3,varwidth=T,ylim=c(0,1),plot=T,main="CV vs AGN distance from 0 (x axis)")

plot(sqrt(moreThan_n_CBG$eAGyyyyDiff^2), moreThan_n_CBG$CV, pch = 16, cex = 2, col = rgb(0,0,0,0.05, maxColorValue = 1), xlim = c(0, 22))
# plot(sqrt(plotReportingDF$eAGyyyyDiff^2),plotReportingDF$IQR, pch = 16, cex = 2, col = ifelse(plotReportingDF$eAGyyyyDiff > 0, rgb(1,0,0,0.05, maxColorValue = 1), rgb(0,0,0,0.05, maxColorValue = 1)))

cor.test(sqrt(moreThan_n_CBG$eAGyyyyDiff^2), moreThan_n_CBG$CV)
abline(lm(moreThan_n_CBG$CV ~ sqrt(moreThan_n_CBG$eAGyyyyDiff^2)),col="red")


## ROC for hypo and mortality
# adjust for nCBG
ROCframe <- moreThan_n_CBG
# ensure that all have at least 1y follow up
ROCframe <- ROCframe[dateplustime1 < (max(ROCframe$dateplustime1) - (60*60*24*365.25))]
# calculate hypo rate
ROCframe$hypoRate <- ROCframe$ID_ADMISSIONhypoEpisodes4.60 / ROCframe$admissionDurationDays
# flag top quartile hypo rate
ROCframe$hypoRate_topQ <- ifelse(ROCframe$hypoRate > (quantile(ROCframe$hypoRate)[4]), 1, 0)

ROCframe$anyHypo3 <- ifelse(ROCframe$minGlu < 3, 1, 0)
ROCframe$anyHypo4 <- ifelse(ROCframe$minGlu < 4, 1, 0)

# 
# # mortality
ROCframe$dead <- ifelse(ROCframe$deathDateUnix>0,1,0)
ROCframe$inHospitalDeath <- ifelse(ROCframe$dead==1 & 
                                      (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration) > 0) &
                                      (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration) < (5*24*60*60)),1,0)
 ROCframe$deadWithin_1Y <- ifelse(ROCframe$dead==1 &
                                    (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration) > 0) &
                                    (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration) < (365.25*24*60*60)),1,0)


# ROCframe <- subset(ROCframe, yyyy1 > 3.9)

library("ROCR")
library(caTools)
 set.seed(123)
 
 
 
roc_function <- function(outcomeParameterTrain, testParameterTrain, inputFrameTrain, outcomeParameterTest, testParameterTest, inputFrameTest) {
  
  trainFrame <- data.frame(outcomeParameterTrain, testParameterTrain, inputFrameTrain$age, inputFrameTrain$medianGlu); colnames(trainFrame) <- c("outcomeParameter", "testParameter", "age", "meanCBG")
  testFrame <- data.frame(outcomeParameterTest, testParameterTest, inputFrameTest$age, inputFrameTest$medianGlu); colnames(testFrame) <- c("outcomeParameter", "testParameter", "age", "meanCBG")
  
  classifier <- glm(outcomeParameter ~ ., data = trainFrame, family = "binomial")
  
  summary(classifier)
  
  prob_pred = predict(classifier, type = 'response', newdata = testFrame)
  
  pred <- prediction(prob_pred, testFrame$outcomeParameter)
  roc.perf <- performance(pred, measure = 'tpr', x.measure = 'fpr')
  plot(roc.perf, col=rainbow(7), main="ROC curve Admissions", xlab="Specificity", ylab="Sensitivity")   
  
  auc.perf = performance(pred, measure = "auc")
  auc.perf@y.values
  
   
  abline(0, 1) #add a 45 degree line
  
  auc = performance(pred, "auc")
  print(auc)
}


ROCframe$spl = sample.split(ROCframe$ID, SplitRatio=0.7)

ROCframe_train <- ROCframe[spl == TRUE]
ROCframe_test <- ROCframe[spl == FALSE]



roc_function(ROCframe_train$anyHypo3, sqrt(ROCframe_train$eAGyyyyDiff^2), ROCframe_train, ROCframe_test$anyHypo3, sqrt(ROCframe_test$eAGyyyyDiff^2), ROCframe_test)
roc_function(ROCframe_train$anyHypo3, ROCframe_train$yyyy1, ROCframe_train, ROCframe_test$anyHypo3, ROCframe_test$yyyy1, ROCframe_test)
roc_function(ROCframe_train$anyHypo3, ROCframe_train$eAG, ROCframe_train, ROCframe_test$anyHypo3, ROCframe_test$eAG, ROCframe_test)

roc_function(ROCframe_train$anyHypo4, sqrt(ROCframe_train$eAGyyyyDiff^2), ROCframe_train, ROCframe_test$anyHypo4, sqrt(ROCframe_test$eAGyyyyDiff^2), ROCframe_test)
roc_function(ROCframe_train$anyHypo4, ROCframe_train$yyyy1, ROCframe_train, ROCframe_test$anyHypo4, ROCframe_test$yyyy1, ROCframe_test)
roc_function(ROCframe_train$anyHypo4, ROCframe_train$eAG, ROCframe_train, ROCframe_test$anyHypo4, ROCframe_test$eAG, ROCframe_test)

roc_function(ROCframe_train$hypoRate_topQ, sqrt(ROCframe_train$eAGyyyyDiff^2), ROCframe_train, ROCframe_test$hypoRate_topQ, sqrt(ROCframe_test$eAGyyyyDiff^2), ROCframe_test)
roc_function(ROCframe_train$hypoRate_topQ, ROCframe_train$yyyy1, ROCframe_train, ROCframe_test$hypoRate_topQ, ROCframe_test$yyyy1, ROCframe_test)
roc_function(ROCframe_train$hypoRate_topQ, ROCframe_train$eAG, ROCframe_train, ROCframe_test$hypoRate_topQ, ROCframe_test$eAG, ROCframe_test)


roc_function(ROCframe_train$inHospitalDeath, sqrt(ROCframe_train$eAGyyyyDiff^2), ROCframe_train, ROCframe_test$inHospitalDeath, sqrt(ROCframe_test$eAGyyyyDiff^2), ROCframe_test)
roc_function(ROCframe_train$inHospitalDeath, ROCframe_train$yyyy1, ROCframe_train, ROCframe_test$inHospitalDeath, ROCframe_test$yyyy1, ROCframe_test)
roc_function(ROCframe_train$inHospitalDeath, ROCframe_train$eAG, ROCframe_train, ROCframe_test$inHospitalDeath, ROCframe_test$eAG, ROCframe_test)

roc_function(ROCframe_train$deadWithin_1Y, sqrt(ROCframe_train$eAGyyyyDiff^2), ROCframe_train, ROCframe_test$deadWithin_1Y, sqrt(ROCframe_test$eAGyyyyDiff^2), ROCframe_test)
roc_function(ROCframe_train$deadWithin_1Y, ROCframe_train$yyyy1, ROCframe_train, ROCframe_test$deadWithin_1Y, ROCframe_test$yyyy1, ROCframe_test)
roc_function(ROCframe_train$deadWithin_1Y, ROCframe_train$eAG, ROCframe_train, ROCframe_test$deadWithin_1Y, ROCframe_test$eAG, ROCframe_test)


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




