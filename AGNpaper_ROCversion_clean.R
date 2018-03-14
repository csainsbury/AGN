library(data.table)
library(survival)

tableSummary <- function(inputFrame) {
  
  # n
  print(nrow(inputFrame))

  # age
  print("age")
  print(summary(inputFrame$age))
  
  # sex
  print("sex")
  inputFrame$sexDigit<-ifelse(nchar(inputFrame$charID==9),as.numeric(substr(inputFrame$charID,8,8)),as.numeric(substr(inputFrame$charID,9,9)))
  inputFrame$binarySexDigit <- ifelse(inputFrame$sexDigit %% 2 == 0, 0, 1)
  print(sum(inputFrame$binarySexDigit))
  
  # last HbA1c
  print("last HbA1c")
  print(summary(inputFrame$lastHbA1cInFrame))
  
  # first CBG
  print("first CBG")
  print(summary(inputFrame$yyyy1))
  
  # n CBG per admission
  print("n CBG per admission")
  print(summary(inputFrame$nCBGperAdmission))
  
  # admission Duration
  print("admission duration days")
  print(summary(inputFrame$admissionDurationDays))
  
  # mean / sd CBG
  print("mean CBG")
  print(summary(inputFrame$meanCBG))
  
  print("sd CBG")
  print(summary(inputFrame$sdCBG))
  
  # CV
  print("CV")
  print(summary(inputFrame$CV))
  
  
}

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
plotReportingDF<-subset(reportingDF,admissionDurationDays>=1)
plotReportingDF<-subset(plotReportingDF,nHbA1cValuesInFrame>0)  # for HbA1c perior: remove those without a CBG in 15 month window
# plotReportingDF<-subset(plotReportingDF,nCBGperAdmission>5)   # for admisisons: remove single CBG admissions


plotReportingDF<-data.table(plotReportingDF)
plotReportingDF$uniqueAdmissionID <- seq(1, nrow(plotReportingDF), 1)
plotReportingDF[, c("meanCBG") := mean(as.numeric(unlist(regmatches(cbgList, gregexpr("-?\\d+\\.\\d+", cbgList))))) ,by=.(uniqueAdmissionID)]
plotReportingDF[, c("sdCBG") := sd(as.numeric(unlist(regmatches(cbgList, gregexpr("-?\\d+\\.\\d+", cbgList))))) ,by=.(uniqueAdmissionID)]

plotReportingDF$meanCBG[is.na(plotReportingDF$meanCBG)] <- 1000
plotReportingDF <- plotReportingDF[meanCBG<1000]
plotReportingDF$sdCBG[is.na(plotReportingDF$sdCBG)] <- 1000
plotReportingDF <- plotReportingDF[sdCBG<1000]


plotReportingDF$diabetesDurationYears <- (plotReportingDF$dateplustime1 - plotReportingDF$diagnosisDateUnix) / (60*60*24*365.25)

##
## correlation plots

# adjust for nCBG
#n = 6 # set minimum number of CBGs per admission to calculate the IQR
#moreThan_n_CBG <- plotReportingDF[nCBGperAdmission > n]
moreThan_n_CBG <- plotReportingDF

# flag first admission for each ID in dataset
moreThan_n_CBG[, c("flagFirstAdmission") := ifelse(admissionNumberFlag == min(admissionNumberFlag), 1, 0),by=.(ID)]

# optional single admission per ID
moreThan_n_CBG <- moreThan_n_CBG[flagFirstAdmission == 1]


# calculate CV
moreThan_n_CBG$CV = moreThan_n_CBG$sdCBG / moreThan_n_CBG$meanCBG

tableSummary(moreThan_n_CBG)
moreThan_n_CBG$adjustedAGN <- sqrt(moreThan_n_CBG$eAGyyyyDiff^2)

tableSummary(moreThan_n_CBG[adjustedAGN > quantile(adjustedAGN)[3]])
tableSummary(moreThan_n_CBG[adjustedAGN <= quantile(adjustedAGN)[3]])

wilcox.test(moreThan_n_CBG[adjustedAGN > quantile(adjustedAGN)[3]]$age, moreThan_n_CBG[adjustedAGN <= quantile(adjustedAGN)[3]]$age)

wilcox.test(moreThan_n_CBG[adjustedAGN > quantile(adjustedAGN)[3]]$admissionDurationDays, moreThan_n_CBG[adjustedAGN <= quantile(adjustedAGN)[3]]$admissionDurationDays)



# using CV rather than IQR
# correlation with first CBG in this subset
yyyy1_distanceFromMin <- moreThan_n_CBG$yyyy1 - 5
adjustedFirstCBG <- sqrt((yyyy1_distanceFromMin)^2)
boxplot(moreThan_n_CBG$CV ~ cut(adjustedFirstCBG,breaks=seq(1,28,1)),las=3,varwidth=T,plot=T,main="IQR vs yyyy1 (x axis)")

plot(adjustedFirstCBG,moreThan_n_CBG$CV, pch = 16, cex = 1, col = rgb(0,0,0,1, maxColorValue = 1))
cor.test(adjustedFirstCBG,moreThan_n_CBG$CV, method = "spearman")
abline(lm(moreThan_n_CBG$CV ~ adjustedFirstCBG),col="red", lwd = 3)

# using CV rather than IQR
# correlation with last HbA1c in this subset
boxplot(moreThan_n_CBG$CV ~ cut(moreThan_n_CBG$lastHbA1cInFrame,breaks=seq(30,200,10)),las=3,varwidth=T,ylim=c(0,0.5),plot=T,main="CV vs last measured HbA1c (x axis)")

plot(moreThan_n_CBG$lastHbA1cInFrame,moreThan_n_CBG$CV, pch = 16, cex = 1, col = rgb(0,0,0,1, maxColorValue = 1))
cor.test(moreThan_n_CBG$lastHbA1cInFrame,moreThan_n_CBG$CV, method = "spearman")
abline(lm(moreThan_n_CBG$CV ~ moreThan_n_CBG$lastHbA1cInFrame),col="red", lwd = 3)


# correlation of AGN with adjusted HbA1c
# using CV rather than IQR
distanceFromZero_IQR <- boxplot(moreThan_n_CBG$CV ~ cut(sqrt(moreThan_n_CBG$eAGyyyyDiff^2),breaks=seq(0,22,1)),las=3,varwidth=T,ylim=c(0,1),plot=T,main="CV vs AGN distance from 0 (x axis)")

plot(sqrt(moreThan_n_CBG$eAGyyyyDiff^2), moreThan_n_CBG$CV, pch = 16, cex = 1, col = rgb(0,0,0,1, maxColorValue = 1), xlim = c(0, 22))
# plot(sqrt(plotReportingDF$eAGyyyyDiff^2),plotReportingDF$IQR, pch = 16, cex = 2, col = ifelse(plotReportingDF$eAGyyyyDiff > 0, rgb(1,0,0,0.05, maxColorValue = 1), rgb(0,0,0,0.05, maxColorValue = 1)))

cor.test(sqrt(moreThan_n_CBG$eAGyyyyDiff^2), moreThan_n_CBG$CV, method = "spearman")
abline(lm(moreThan_n_CBG$CV ~ sqrt(moreThan_n_CBG$eAGyyyyDiff^2)),col="red", lwd = 3)

# tests accounting for nCBG
if(!require("ppcor")){
  install.packages("ppcor", repos='http://cran.us.r-project.org')
  library(ppcor)
}

pcor.test(adjustedFirstCBG, moreThan_n_CBG$CV, moreThan_n_CBG$nCBGperAdmission, method = "spearman")
pcor.test(moreThan_n_CBG$lastHbA1cInFrame, moreThan_n_CBG$CV, moreThan_n_CBG$nCBGperAdmission, method = "spearman")
pcor.test(sqrt(moreThan_n_CBG$eAGyyyyDiff^2), moreThan_n_CBG$CV, moreThan_n_CBG$nCBGperAdmission, method = "spearman")



## ROC for hypo and mortality
# adjust for nCBG
ROCframe <- moreThan_n_CBG
# ensure that all have at least 1y follow up
ROCframe <- ROCframe[dateplustime1 < (max(ROCframe$dateplustime1) - (60*60*24*365.25))]
# calculate hypo rate
ROCframe$hypoRate <- ROCframe$ID_ADMISSIONhypoEpisodes4.60 / ROCframe$admissionDurationDays
# flag top quartile hypo rate
ROCframe$hypoRate_topQ <- ifelse(ROCframe$hypoRate > (quantile(ROCframe$hypoRate)[4]), 1, 0)

min(as.numeric(unlist(regmatches(ROCframe$cbgList[1], gregexpr("-?\\d+\\.\\d+", ROCframe$cbgList[1])))))

# ROCframe[, c("anyHypo3") := ifelse(min(as.numeric(unlist(regmatches(cbgList, gregexpr("-?\\d+\\.\\d+", cbgList))))) < 3, 1, 0) ,by=.(uniqueAdmissionID)]
# ROCframe[, c("anyHypo4") := ifelse(min(as.numeric(unlist(regmatches(cbgList, gregexpr("-?\\d+\\.\\d+", cbgList))))) < 4, 1, 0) ,by=.(uniqueAdmissionID)]

# 
# # mortality
ROCframe$dead <- ifelse(ROCframe$deathDateUnix>0,1,0)
ROCframe$inHospitalDeath <- ifelse(ROCframe$dead==1 & 
                                      (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration) > 0) &
                                      (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration) < (5*24*60*60)),1,0)

ROCframe$deadWithin_1Y <- ifelse(ROCframe$dead==1 &
                                    (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration) > 0) &
                                    (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration) < (365.25*24*60*60)),1,0)

ROCframe$deadWithin_3Y <- ifelse(ROCframe$dead==1 &
                                   (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration) > 0) &
                                   (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration) < (3 * 365.25*24*60*60)),1,0)

ROCframe$deadWithin_100days <- ifelse(ROCframe$dead==1 &
                                   (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration) > 0) &
                                   (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration) < (100*24*60*60)),1,0)


# ROCframe <- subset(ROCframe, yyyy1 > 3.9)

 '''
library("ROCR")
library(caTools)
 set.seed(123)
 
roc_function <- function(outcomeParameterTrain, testParameterTrain, inputFrameTrain, outcomeParameterTest, testParameterTest, inputFrameTest) {
  
  trainFrame <- data.frame(outcomeParameterTrain, testParameterTrain, inputFrameTrain$age, inputFrameTrain$meanCBG, inputFrameTrain$sdCBG); colnames(trainFrame) <- c("outcomeParameter", "testParameter", "age", "meanCBG", "sdCBG")
  testFrame <- data.frame(outcomeParameterTest, testParameterTest, inputFrameTest$age, inputFrameTest$meanCBG, inputFrameTest$sdCBG); colnames(testFrame) <- c("outcomeParameter", "testParameter", "age", "meanCBG", "sdCBG")
  
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
'''
## 
simpleSurvivalPlot<-function(inputFrame, eventFlag, timeToEvent, testParam, testParamThresh, endDateUnix,ylimMin) {
  # inputFrame <- ROCframe
  # endDateUnix <- quantile(inputFrame$admissionDurationDays)[4]

  SurvivalData<-inputFrame
  
  #shortCensorPeriodStartDay  <- DaySeconds
  #shortCensorPeriodEndDay    <- DaySeconds*10000
  
  lastDOD<-quantile(inputFrame$admissionDurationDays)[4]
  
  SurvivalData$dateOfDischarge <- 0
  SurvivalData$timeToDeath <- ifelse(eventFlag == 1, timeToEvent, lastDOD + 1)

  SurvivalData$timeToDeathInterval <- SurvivalData$timeToDeath 
    
  #SurvivalData$timeToDeathInterval[is.na(SurvivalData$timeToDeathInterval)]<-0
  #SurvivalData<-subset(SurvivalData,timeToDeathInterval>0)
  
  # SurvivalData$timeToDeathInterval<-SurvivalData$timeToDeathInterval/(60*60*24*365.25)
  
  SurvivalData$shortDeathEvent <- eventFlag
  # SurvivalData$shortDeathEvent <- ifelse(SurvivalData$isDead==1 & SurvivalData$timeToDeath>=(shortCensorPeriodStartDay) & SurvivalData$timeToDeath<(shortCensorPeriodEndDay),1,0)	
  
  #  SurvivalData$sexDigit<-ifelse(nchar(SurvivalData$charID==9),as.numeric(substr(SurvivalData$charID,8,8)),as.numeric(substr(SurvivalData$charID,9,9)))
  # SurvivalData$sexNumber<-ifelse(SurvivalData$sexDigit%%2==0,1,0)
  #  SurvivalData$sex<-factor(1*(SurvivalData$sexNumber <1),levels=0:1,labels=c("F","M"))
  
  
  mfitAge50<-survfit(Surv(timeToDeathInterval, shortDeathEvent) ~ (testParam > testParamThresh), data = SurvivalData)
  shortPlotTitle <- paste("Mortality, time ",round(shortCensorPeriodStartDay)/DaySeconds," to ",round(max(SurvivalData$timeToDeathInterval))/DaySeconds," days\n n= ",nrow(SurvivalData),", threshold: ",quantile(SurvivalData$hba1cIQRinRange)[3],sep="")
  plot(mfitAge50,mark.time=T,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),xlim=c(shortCensorPeriodStartDay,round(max(SurvivalData$timeToDeathInterval))),lwd=5,ylim=c(ylimMin,1))
  
  # mfitAge50.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ age_atSampleTime+medianHbA1cInRange+nValsPerIDinRange+(hba1cIQRinRange>=quantile(SurvivalData$hba1cIQRinRange)[3]), data = SurvivalData)
  
  mfitAge50.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ (testParam > testParamThresh), data = SurvivalData)
  pVal <- summary(mfitAge50.coxph)$coef[,5]; HR <- round(exp(coef(mfitAge50.coxph)),2)
  legendText <- paste("p = ",pVal," | HR = ",HR,sep="")
  summarySurvfit <- summary(mfitAge50); legendNames <- row.names(summarySurvfit$table)
  legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("topright",legendText,cex=0.6)
  
  print(mfitAge50.coxph)
  
}

timeToHypoFunction  <- function(cbgList, cbgTimeList, hypoThreshold) {
  # cbgList = as.numeric(unlist(regmatches(ROCframe$cbgList[1], gregexpr("-?\\d+\\.\\d+", ROCframe$cbgList[1]))))
  # cbgTimeList = as.numeric(unlist(regmatches(ROCframe$cbgTimeList[1], gregexpr("-?\\d+", ROCframe$cbgTimeList[1]))))
  numericCBGlist <- as.numeric(unlist(regmatches(cbgList, gregexpr("-?\\d+\\.\\d+", cbgList))))
  numericCBGTimeList <- as.numeric(unlist(regmatches(cbgTimeList, gregexpr("-?\\d+", cbgTimeList))))
  
  diffTime <- (numericCBGTimeList - numericCBGTimeList[1]) / (60*60*24)
  firstCBGbelowThresh <- diffTime[min(which(numericCBGlist < hypoThreshold))]
  
  return(firstCBGbelowThresh)

}

##

survivalFunction_hypo <- function(survivalFrame, hypoThreshold, maxAdmissionDuration_forAnalysis, ylimMin) {

# survivalFrame = ROCframe; hypoThreshold = 4; maxAdmissionDuration_forAnalysis = 10; ylimMin = 0

survivalFrame = survivalFrame[yyyy1 > hypoThreshold]

# survivalFrame[, c("minGlu_duringAdmission") := min(as.numeric(unlist(regmatches(cbgList, gregexpr("-?\\d+\\.\\d+", cbgList))))) ,by=.(uniqueAdmissionID)]

survivalFrame[, c("anyHypo_thresh") := ifelse(min(as.numeric(unlist(regmatches(cbgList, gregexpr("-?\\d+\\.\\d+", cbgList))))) < hypoThreshold, 1, 0) ,by=.(uniqueAdmissionID)]

survivalFrame[, c("timeToFirstHypo_thresh") := ifelse(anyHypo_thresh == 1, timeToHypoFunction(cbgList, cbgTimeList, hypoThreshold), 0) ,by=.(uniqueAdmissionID)]

AGNdistanceFrom0 = sqrt(survivalFrame$eAGyyyyDiff ^ 2)

survivalFrame$timeToFirstHypo_thresh <- ifelse(survivalFrame$timeToFirstHypo_thresh == 0, max(survivalFrame$admissionDurationDays), survivalFrame$timeToFirstHypo_thresh)

# survivalFrame$timeToFirstHypo3 <- ifelse(survivalFrame$timeToFirstHypo3 > 10, 10, survivalFrame$timeToFirstHypo3)
survivalFrame$anyHypo_thresh <- ifelse(survivalFrame$timeToFirstHypo_thresh > maxAdmissionDuration_forAnalysis, 0, survivalFrame$anyHypo_thresh)

y_bmt <- Surv(survivalFrame$timeToFirstHypo_thresh, survivalFrame$anyHypo_thresh)

fit1_bmt <- survfit(y_bmt ~ AGNdistanceFrom0 > quantile(AGNdistanceFrom0)[3])
#summary(fit1_bmt)

plot(fit1_bmt, xlim = c(0, maxAdmissionDuration_forAnalysis), ylim = c(ylimMin, 1), col = c(1:2), lwd = 3, main = paste("n = ", nrow(survivalFrame), ". hypo threshold = ", hypoThreshold, sep = ""))

fit1.coxph <- coxph(y_bmt ~ (AGNdistanceFrom0 > quantile(AGNdistanceFrom0)[3]) + survivalFrame$age)
summary(fit1.coxph)

# test proportional hazards assumption
test.ph <- cox.zph(fit1.coxph)
print(test.ph)
#plot(test.ph)

}

survivalFunction_hypo_timesplit <- function(survivalFrame, hypoThreshold, maxAdmissionDuration_forAnalysis, ylimMin) {
  
  
  survivalFrame = survivalFrame[yyyy1 > hypoThreshold]
  
  survivalFrame[, c("anyHypo_thresh") := ifelse(min(as.numeric(unlist(regmatches(cbgList, gregexpr("-?\\d+\\.\\d+", cbgList))))) < hypoThreshold, 1, 0) ,by=.(uniqueAdmissionID)]
  
  survivalFrame[, c("timeToFirstHypo_thresh") := ifelse(anyHypo_thresh == 1, timeToHypoFunction(cbgList, cbgTimeList, hypoThreshold), 0) ,by=.(uniqueAdmissionID)]
  
  AGNdistanceFrom0 = sqrt(survivalFrame$eAGyyyyDiff ^ 2)
  
  survivalFrame$timeToFirstHypo_thresh <- ifelse(survivalFrame$timeToFirstHypo_thresh == 0, max(survivalFrame$admissionDurationDays), survivalFrame$timeToFirstHypo_thresh)
  
  survivalFrame = data.frame(survivalFrame$timeToFirstHypo_thresh, survivalFrame$anyHypo_thresh, survivalFrame$age, AGNdistanceFrom0)
  colnames(survivalFrame) <- c("time", "eventMarker", "age", "AGN")
  
  # set marker for analysis
  survivalFrame$highAGN <- ifelse(survivalFrame$AGN > quantile(survivalFrame$AGN)[3], 1, 0)
  # survivalFrame$eventMarker <- factor(survivalFrame$eventMarker)
  
  # survivalFrame$timeToFirstHypo3 <- ifelse(survivalFrame$timeToFirstHypo3 > 10, 10, survivalFrame$timeToFirstHypo3)
  survivalFrame$eventMarker <- ifelse(survivalFrame$time > maxAdmissionDuration_forAnalysis, 0, survivalFrame$eventMarker)
  
  y_bmt <- Surv(survivalFrame$time, survivalFrame$eventMarker)
  
  fit1_bmt <- survfit(y_bmt ~ survivalFrame$highAGN)
  #summary(fit1_bmt)
  
  plot(fit1_bmt, xlim = c(0, maxAdmissionDuration_forAnalysis), ylim = c(ylimMin, 1), col = c(1:2))
  
  fit1.coxph <- coxph(y_bmt ~ (survivalFrame$highAGN) + survivalFrame$age)
  print(summary(fit1.coxph))
  
  # test proportional hazards assumption
  test.ph <- cox.zph(fit1.coxph)
  print(test.ph)
  
  library(Greg)
  library(dplyr)
  
  spl_surv <-
    survivalFrame %>% 
    timeSplitter(by = 1,
                 event_var = "eventMarker",
                 event_start_status = 0,
                 time_var = "time",
                 time_related_vars = c("age"))
  
  interval_model <-
    update(fit1.coxph, 
           Surv(Start_time, Stop_time, eventMarker == 1) ~ age + highAGN,
           data = spl_surv)
            
  summary(interval_model)
  
  time_int_model <- 
    update(interval_model,
           .~.+age:Start_time)
  print(summary(time_int_model))
  
  time.ph <- cox.zph(time_int_model)
  print(time.ph)
  
}
 
survivalFunction_hypo(ROCframe, 4, 10, 0.8)
survivalFunction_hypo(ROCframe, 3, 10, 0.9)

# call survival function with timesplit
survivalFunction_hypo_timesplit(ROCframe, 4, 10, 0.8)
survivalFunction_hypo_timesplit(ROCframe, 3, 10, 0.9)

##

survivalFunction_mort_postDischargeSurvival <- function(survivalFrame, deathEventMarker, rightCensorPointDays, maxAdmissionDuration_forAnalysis, ylimMin) {
  
  # survivalFrame = surviveToDischarge; deathEventMarker = surviveToDischarge$deadWithin_3Y; rightCensorPointDays = 365.25 * 3; maxAdmissionDuration_forAnalysis = 365.25*3
  
  
  timeToDeath = ifelse(deathEventMarker == 1, (survivalFrame$deathDateUnix - (survivalFrame$dateplustime1 + survivalFrame$admissionDuration)) / (60*60*24), 0)
  
  AGNdistanceFrom0 = sqrt(survivalFrame$eAGyyyyDiff ^ 2)
  
  survivalFrame$timeToDeath <- ifelse(timeToDeath == 0, rightCensorPointDays, timeToDeath)
  deathEventMarker <- ifelse(survivalFrame$timeToDeath > maxAdmissionDuration_forAnalysis, 0, deathEventMarker)
  
  
  y_bmt <- Surv(survivalFrame$timeToDeath, deathEventMarker)
  
  fit1_bmt <- survfit(y_bmt ~ AGNdistanceFrom0 > quantile(AGNdistanceFrom0)[3])
  #summary(fit1_bmt)
  
  plot(fit1_bmt, xlim = c(0, maxAdmissionDuration_forAnalysis), ylim = c(ylimMin, 1), col = c(1:2), lwd = 3, main = paste("n = ", nrow(survivalFrame), ". n events = ", sum(fit1_bmt$n.event), sep = ""))
  
  fit1.coxph <- coxph(y_bmt ~ (AGNdistanceFrom0 > quantile(AGNdistanceFrom0)[3]) + survivalFrame$age)
  print(summary(fit1.coxph))
  
  test.ph <- cox.zph(fit1.coxph)
  print(test.ph)
  
}

survivalFunction_mort <- function(survivalFrame, deathEventMarker, rightCensorPointDays, maxAdmissionDuration_forAnalysis, ylimMin) {
  
  # survivalFrame = ROCframe
  
  timeToDeath = ifelse(deathEventMarker == 1, (survivalFrame$deathDateUnix - survivalFrame$dateplustime1) / (60*60*24), 0)
  
  AGNdistanceFrom0 = sqrt(survivalFrame$eAGyyyyDiff ^ 2)
  
  survivalFrame$timeToDeath <- ifelse(timeToDeath == 0, rightCensorPointDays, timeToDeath)
  deathEventMarker <- ifelse(survivalFrame$timeToDeath > maxAdmissionDuration_forAnalysis, 0, deathEventMarker)
  
  
  y_bmt <- Surv(survivalFrame$timeToDeath, deathEventMarker)
  
  fit1_bmt <- survfit(y_bmt ~ AGNdistanceFrom0 > quantile(AGNdistanceFrom0)[3])
  #summary(fit1_bmt)
  
  plot(fit1_bmt, xlim = c(0, maxAdmissionDuration_forAnalysis), ylim = c(ylimMin, 1), col = c(1:2), lwd = 3, main = paste("n = ", nrow(survivalFrame), ". max admission time for analysis = ", maxAdmissionDuration_forAnalysis, ". n events = ", sum(fit1_bmt$n.event), sep = ""))
  
  fit1.coxph <- coxph(y_bmt ~ (AGNdistanceFrom0 > quantile(AGNdistanceFrom0)[3]) + survivalFrame$age)
  print(summary(fit1.coxph))
  
  # test proportional hazards assumption
  test.ph <- cox.zph(fit1.coxph)
  print(test.ph)
  #plot(test.ph)
  
}

survivalFunction_mort(ROCframe, ROCframe$inHospitalDeath, max(ROCframe$admissionDurationDays), 30, 0.98)

# remove those who die in hospital
surviveToDischarge <- ROCframe[deathDateUnix > 0 & (ROCframe$deathDateUnix - (ROCframe$dateplustime1 + ROCframe$admissionDuration)) > (5 * (60*60*24)) | ROCframe$deathDateUnix == 0]

survivalFunction_mort_postDischargeSurvival(surviveToDischarge, surviveToDischarge$deadWithin_3Y, 365.25* 3, 365.25 * 3, 0.8)
survivalFunction_mort_postDischargeSurvival(surviveToDischarge, surviveToDischarge$deadWithin_1Y, 365.25, 365.25, 0.9)
survivalFunction_mort_postDischargeSurvival(surviveToDischarge, surviveToDischarge$deadWithin_100days, 100, 100, 0.95)



