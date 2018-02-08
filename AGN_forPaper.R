library(data.table)
library(survival)
#source("./_T1DM_hba1c_admissionAndMortality.R")
#####################################################################################################
#####################################################################################################
## read in diabetes / hba1c prediction data, analyse and plot
#####################################################################################################
#####################################################################################################

###### type 1
#
#
tempWriteFile <- paste("~/R/GlCoSy/dataForSubmissions/attd2017/admissionDataDT_T1DM.csv",sep=""); reportingDF<-read.csv(tempWriteFile); reportingDF<-data.table(reportingDF); diabetesType="Type 1 Diabetes. "
#

###### type 1 - hba1c revised code
#
# tempWriteFile <- paste("../GlCoSy/source/admissionDataDT_T1DM_hba1cRevision.csv",sep=""); reportingDF<-read.csv(tempWriteFile); reportingDF<-data.table(reportingDF); diabetesType="Type 1 Diabetes. "
#
#####################################################################################################
###### type 2
#
# tempWriteFile <- paste("~/R/GlCoSy/source/admissionDataDT_T2DM.csv",sep=""); reportingDF<-read.csv(tempWriteFile); reportingDF<-data.table(reportingDF); diabetesType="Type 2 Diabetes. "
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
plotReportingDF<-subset(reportingDF,nCBGperAdmission>=2)   # for admisisons: remove single CBG admissions
plotReportingDF<-subset(plotReportingDF,nHbA1cValuesInFrame>0)  # for HbA1c perior: remove those without a CBG in 15 month window
plotReportingDF<-data.table(plotReportingDF)

plotReportingDF$diabetesDurationYears <- (plotReportingDF$dateplustime1 - plotReportingDF$diagnosisDateUnix) / (60*60*24*365.25)


plotfilename <- paste("./_Plots.pdf",sep="")
pdf(plotfilename, width=16, height=9)

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
  n = 16 # set minimum number of CBGs per admission to calculate the IQR
  moreThan_n_CBG <- plotReportingDF[nCBGperAdmission > n]
  
  # demonstrate linear association between IQR and log(nCBG)
  plot(log(moreThan_n_CBG$nCBGperAdmission), moreThan_n_CBG$IQR)
  boxplot(moreThan_n_CBG$IQR ~ cut(log(moreThan_n_CBG$nCBGperAdmission), breaks= 40), varwidth= T)
  
  # correlation with last HbA1c in this subset
  boxplot(moreThan_n_CBG$IQR ~ cut(moreThan_n_CBG$lastHbA1cInFrame,breaks=seq(30,200,10)),las=3,varwidth=T,ylim=c(0,10),plot=T,main="IQR vs last measured HbA1c (x axis)")
  
  plot(moreThan_n_CBG$lastHbA1cInFrame,moreThan_n_CBG$IQR)
  cor.test(moreThan_n_CBG$lastHbA1cInFrame,moreThan_n_CBG$IQR)
  abline(lm(moreThan_n_CBG$IQR ~ moreThan_n_CBG$lastHbA1cInFrame),col="red")
  
    # correlation with first CBG in this subset
    boxplot(moreThan_n_CBG$IQR ~ cut(moreThan_n_CBG$yyyy1,breaks=seq(1,28,1)),las=3,varwidth=T,plot=T,main="IQR vs last measured HbA1c (x axis)")
    
    plot(moreThan_n_CBG$yyyy1,moreThan_n_CBG$IQR)
    cor.test(moreThan_n_CBG$yyyy1,moreThan_n_CBG$IQR)
    abline(lm(moreThan_n_CBG$IQR ~ moreThan_n_CBG$yyyy1),col="red")
  
      # correlation of AGN with adjusted HbA1c
      adjusted_IQR <-  moreThan_n_CBG$IQR / moreThan_n_CBG$nCBGperAdmission
      
      distanceFromZero_IQR <- boxplot(moreThan_n_CBG$IQR ~ cut(sqrt(moreThan_n_CBG$eAGyyyyDiff^2),breaks=seq(0,22,1)),las=3,varwidth=T,ylim=c(0,15),plot=T,main="IQR vs AGN distance from 0 (x axis)")
      
      plot(sqrt(moreThan_n_CBG$eAGyyyyDiff^2), moreThan_n_CBG$IQR, pch = 16, cex = 2, col = rgb(0,0,0,0.05, maxColorValue = 1), xlim = c(0, 22))
      # plot(sqrt(plotReportingDF$eAGyyyyDiff^2),plotReportingDF$IQR, pch = 16, cex = 2, col = ifelse(plotReportingDF$eAGyyyyDiff > 0, rgb(1,0,0,0.05, maxColorValue = 1), rgb(0,0,0,0.05, maxColorValue = 1)))
      
      cor.test(sqrt(moreThan_n_CBG$eAGyyyyDiff^2), moreThan_n_CBG$IQR)
      abline(lm(moreThan_n_CBG$IQR ~ sqrt(moreThan_n_CBG$eAGyyyyDiff^2)),col="red")

##### plotting the correlation of the three options for admissions with variable nCBG
      nn = 100
      
      reportMatrix <- as.data.frame(matrix(0, nrow = nn, ncol = 5))
      colnames(reportMatrix) <- c("min_CBG", "n_IDs", "lastHbA1c", "firstCBG", "AGN")
      
      for (i in seq(1, nn, 1)) {
        # adjust for nCBG
        n = i # set minimum number of CBGs per admission to calculate the IQR
        moreThan_n_CBG <- plotReportingDF[nCBGperAdmission > n]
        
        reportMatrix$min_CBG[i] <- i
        reportMatrix$n_IDs[i] <- nrow(moreThan_n_CBG)
        
        #plot(moreThan_n_CBG$lastHbA1cInFrame,moreThan_n_CBG$IQR)
        hba1cTest <- cor.test(moreThan_n_CBG$lastHbA1cInFrame,moreThan_n_CBG$IQR)
        reportMatrix$lastHbA1c[i] <- hba1cTest$estimate
        
        yyyyTest <- cor.test(moreThan_n_CBG$yyyy1,moreThan_n_CBG$IQR)
        reportMatrix$firstCBG[i] <- yyyyTest$estimate
        
        AGNtest <- cor.test(sqrt(moreThan_n_CBG$eAGyyyyDiff^2), moreThan_n_CBG$IQR)
        reportMatrix$AGN[i] <- AGNtest$estimate
        
      }
      
      plot(reportMatrix$min_CBG, reportMatrix$lastHbA1c, pch = 16, cex = 1, col = "black", ylim = c(0, 0.3))
      lines(reportMatrix$min_CBG, reportMatrix$lastHbA1c, col = "black")
      points(reportMatrix$min_CBG, reportMatrix$firstCBG, pch = 16, cex = 1, col = "green")
      lines(reportMatrix$min_CBG, reportMatrix$firstCBG, col = "green")
      points(reportMatrix$min_CBG, reportMatrix$AGN, pch = 16, cex = 1, col = "red")
      lines(reportMatrix$min_CBG, reportMatrix$AGN, col = "red")
      #####


## type 1 cor - 0.203, type 2 - 0.298

# boxplot(plotReportingDF$IQR ~ cut(log(plotReportingDF$lastHbA1cInFrame),breaks=30),las=3,varwidth=T,ylim=c(0,10),plot=T,main="IQR vs log of last measured HbA1c (x axis)")

# plot(log(subset(plotReportingDF,IQR>0)$lastHbA1cInFrame),subset(plotReportingDF,IQR>0)$IQR,cex=sqrt(plotReportingDF$nCBGperAdmission)/2,col=ifelse(plotReportingDF$nCBGperAdmission<5,"black","black"),main="IQR vs last measured log HbA1c. circle size represents nCBG in admission")

# IQR_hba1cDF<-as.data.frame(matrix(0,nrow=0,ncol=4)); colnames(IQR_hba1cDF)<-c("nCBG","nAdmissions","pval","cor")
# for (ii in seq(2,20,1)) {
#  plot(log(subset(plotReportingDF,nCBGperAdmission==ii)$lastHbA1cInFrame),subset(plotReportingDF,nCBGperAdmission==ii)$IQR)
#  ct<-cor.test(log(subset(plotReportingDF,nCBGperAdmission==ii)$lastHbA1cInFrame),subset(plotReportingDF,nCBGperAdmission==ii)$IQR)

#  output<-data.frame(ii,length(log(subset(plotReportingDF,nCBGperAdmission==ii)$lastHbA1cInFrame)),ct$p.value,ct$estimate)
#  colnames(output)<-c("nCBG","nAdmissions","pval","cor")
#  IQR_hba1cDF<-rbind(IQR_hba1cDF,output)
# }

for (i in seq(2,20,1)) {
  x<-boxplot(subset(plotReportingDF,nCBGperAdmission==i)$IQR ~ cut(subset(plotReportingDF,nCBGperAdmission==i)$lastHbA1cInFrame,breaks=seq(30,200,5)),las=3,varwidth=T,ylim=c(0,10),plot=F)
  if (i==4) {plot(log(c(1:20)),x$stats[3,1:20],ylim=c(0,10),pch=16,cex=3*(sqrt(x$n))/max(sqrt(x$n)),col=i,main="IQR vs log last measured HbA1c. each line represents an admission with n CBGs. point size corresponding to number of values\nthis controls for the issue of IQR increasing with increasing numbers of CBGs measured during an admission")
    lines(log(c(1:20)),x$stats[3,1:20],col=i,lwd=2)
  }
  if (i>4)  {points(log(c(1:20)),x$stats[3,1:20],pch=16,cex=3*(sqrt(x$n))/max(sqrt(x$n)),col=i)
    lines(log(c(1:20)),x$stats[3,1:20],col=i,lwd=2)
  }
  
}

## HbA1c IQR
# boxplot(plotReportingDF$IQR ~ cut(plotReportingDF$IQRHbA1cInFrame,breaks=seq(1,30,1)),las=3,varwidth=T,ylim=c(0,10),plot=T,main="IQR vs IQR of HbA1c measures during 15 months prior to admission")

# plot(log(subset(plotReportingDF,IQR>0)$IQRHbA1cInFrame),subset(plotReportingDF,IQR>0)$IQR,cex=sqrt(plotReportingDF$nCBGperAdmission)/2,col=ifelse(plotReportingDF$nCBGperAdmission<5,"black","black"))

# for (i in seq(2,10,1)) {
#  x<-boxplot(subset(plotReportingDF,nCBGperAdmission==i)$IQR ~ cut(subset(plotReportingDF,nCBGperAdmission==i)$IQRHbA1cInFrame,breaks=seq(1,80,1)),las=3,varwidth=T,ylim=c(0,10),plot=F)

# if (i==4) {plot(log(c(1:20)),x$stats[3,1:20],ylim=c(0,10),pch=16,cex=4*(sqrt(x$n))/max(sqrt(x$n)),col=i,main="IQR vs log median prior HbA1c. each line represents an admission with n CBGs. point size corresponding to number of values\nthis controls for the issue of IQR increasing with increasing numbers of CBGs measured during an admission"); lines(log(c(1:20)),x$stats[3,1:20],col=i,lwd=3)}
# if (i>4)  {points(log(c(1:20)),x$stats[3,1:20],pch=16,cex=4*(sqrt(x$n))/max(sqrt(x$n)),col=i); lines(log(c(1:20)),x$stats[3,1:20],col=i,lwd=3)}

# }

######################################################
## admission duration

## most recent HbA1c in range (15 months)
# boxplot(plotReportingDF$admissionDurationDays ~ cut(plotReportingDF$lastHbA1cInFrame,breaks=seq(30,200,5)),las=3,varwidth=T,ylim=c(0,4),plot=T,main="admission duration vs last measured HbA1c value")

admissionDbox<-boxplot(plotReportingDF$admissionDurationDays ~ cut(plotReportingDF$lastHbA1cInFrame,breaks=seq(30,200,10)),las=3,varwidth=T,ylim=c(0,5),plot=T,main="admission duration vs last measured HbA1c value. increments of 10mmol/mol")
  
admissionDbox<-boxplot(plotReportingDF$admissionDurationDays ~ cut(plotReportingDF$yyyy1,breaks=20),las=3,varwidth=T,ylim=c(0,5),plot=T,main="admission duration vs last measured HbA1c value. increments of 10mmol/mol")

admissionDboxAGN<-boxplot(plotReportingDF$admissionDurationDays ~ cut(plotReportingDF$eAGyyyyDiff,breaks=seq(-30,30,2)),las=3,varwidth=T,ylim=c(0,10),plot=T,main="IQR vs last measured HbA1c (x axis)")



attdAbstractLOS<-boxplot(subset(plotReportingDF,admissionDurationDays>0.2)$admissionDurationDays ~ cut(subset(plotReportingDF,admissionDurationDays>0.2)$eAGyyyyDiff,breaks=seq(-22,22,2)),las=3,varwidth=T,ylim=c(0,8),plot=T,main="LOS vs AGN ATTD abstract 1 (x axis)")

attdAbstractLOSdecile<-boxplot(subset(plotReportingDF,admissionDurationDays>1)$admissionDurationDays ~cut(subset(plotReportingDF,admissionDurationDays>1)$eAGyyyyDiff,breaks=quantile(plotReportingDF$eAGyyyyDiff, prob = seq(0, 1, length = 11), type = 5)),plot=T,main="IQR vs AGN ATTD abstract 1 (x axis)",las=3,ylim=c(0,8))


boxplot(plotReportingDF$admissionDurationDays ~ cut(sqrt(plotReportingDF$eAGyyyyDiff^2),breaks=seq(0,30,1)),las=3,varwidth=T,ylim=c(0,10),plot=T,main="IQR vs last measured HbA1c (x axis)")

######################################################
## hypoEpisodes/day- for type 1 attd abstract
plotReportingDF$hypoEpPerDay<-plotReportingDF$ID_ADMISSIONhypoEpisodes4.60 / plotReportingDF$admissionDurationDays

attdAbstractHYD<-boxplot(subset(plotReportingDF,admissionDurationDays>1)$hypoEpPerDay ~ cut(subset(plotReportingDF,admissionDurationDays>1)$eAGyyyyDiff,breaks=seq(-30,30,2)),las=3,varwidth=T,ylim=c(0,1),plot=T,main="IQR vs AGN ATTD abstract 1 (x axis)")

attdAbstractHYD<-boxplot(subset(plotReportingDF,admissionDurationDays>1 & yyyy>3.9)$hypoEpPerDay ~ cut(subset(plotReportingDF,admissionDurationDays>1 & yyyy>3.9)$eAGyyyyDiff,breaks=10),las=3,varwidth=T,ylim=c(0,1),plot=T,main="IQR vs AGN ATTD abstract 1 (x axis)")

# plot episodes/day
#sequencePlotSet<-subset(plotReportingDF,admissionDurationDays>4)
sequencePlotSet<-plotReportingDF

sequence<-seq(-30,30,2)
epDayMatrix<-as.data.frame(matrix(0,nrow=length(sequence),ncol=4));colnames(epDayMatrix)<-c("AGNlow","AGNhigh","hypEps","admisssonDuration")
epDayMatrix$AGNlow<-sequence; epDayMatrix$AGNhigh[1:nrow(epDayMatrix)-1]<-epDayMatrix$AGNlow[2:nrow(epDayMatrix)]; epDayMatrix$AGNhigh[nrow(epDayMatrix)]<-100
for (i in seq(1,length(sequence),1)) {
  subsetFrame<-subset(sequencePlotSet,eAGyyyyDiff>epDayMatrix$AGNlow[i] & eAGyyyyDiff<=epDayMatrix$AGNhigh[i])
  epDayMatrix$hypEps[i]<-sum(subsetFrame$ID_ADMISSIONhypoEpisodes4.60)
  epDayMatrix$admisssonDuration[i]<-sum(subsetFrame$admissionDurationDays)
}
epDayMatrix$ratio<-epDayMatrix$hypEps/epDayMatrix$admisssonDuration

epDayMatrix_nMoreThan100<-subset(epDayMatrix,hypEps>100)
for (j in seq(1,nrow(epDayMatrix_nMoreThan100),1)) {
  epDayMatrix_nMoreThan100$propTestP[j]<-prop.test(c(epDayMatrix_nMoreThan100$hypEps[j],1082),c(epDayMatrix_nMoreThan100$admisssonDuration[j],7.733240e+03))$p.value
}

plot(epDayMatrix_nMoreThan100$AGNlow,epDayMatrix_nMoreThan100$ratio,main="AGN vs hypo episodes per day for all admissions", pch=16, cex=(sqrt(epDayMatrix_nMoreThan100$hypEps)/10), yaxt="n", xaxt="n", ylab="", xlab="")
lines(epDayMatrix_nMoreThan100$AGNlow,epDayMatrix_nMoreThan100$ratio, lwd=2)
axis(2,cex.axis=2)
axis(1,cex.axis=2)
mtext("episodes hypoglycemia / day", side=2, line=2.5, cex=2)
mtext("AGN (mmol/l)", side=1, line=2.5, cex=2)






#for (i in seq(4,30,1)) {
#  x<-boxplot(subset(plotReportingDF,nCBGperAdmission==i)$admissionDurationDays ~ cut(subset(plotReportingDF,nCBGperAdmission==i)$lastHbA1cInFrame,breaks=seq(30,200,5)),las=3,varwidth=T,ylim=c(0,5),plot=F)

#  if (i==4) {plot(c(1:30),x$stats[3,1:30],ylim=c(0,20),pch=16,cex=4*(sqrt(x$n))/max(sqrt(x$n)),col=i,main="admission duration vs last measured HbA1c. each line represents an admission with n CBGs. point size corresponding to number of values\nthis controls for the issue of admission duration increasing with increasing numbers of CBGs measured during an admission"); lines(c(1:30),x$stats[3,1:30],col=i,lwd=3)}
#  if (i>4) {points(c(1:30),x$stats[3,1:30],ylim=c(0,10),pch=16,cex=4*(sqrt(x$n))/max(sqrt(x$n)),col=i); lines(c(1:30),x$stats[3,1:30],col=i,lwd=3)}

#}

######################################################
## median glucose during admission 

## most recent HbA1c in range (15 months)
boxplot(plotReportingDF$medianGlu ~ cut(plotReportingDF$lastHbA1cInFrame,breaks=seq(30,200,10)),las=3,varwidth=T,ylim=c(5,20),plot=T,main="admission median glucose vs most recent HbA1c value")

plot(plotReportingDF$lastHbA1cInFrame,plotReportingDF$medianGlu)
cor.test(plotReportingDF$lastHbA1cInFrame,plotReportingDF$medianGlu)
abline(lm(plotReportingDF$medianGlu ~ plotReportingDF$lastHbA1cInFrame),col="red")
## cor for t1 0.25, for t2 0.46

boxplot(plotReportingDF$medianGlu ~ cut(plotReportingDF$eAGyyyyDiff,breaks=seq(-30,30,2)),las=3,varwidth=T,ylim=c(0,20),plot=T,main="IQR vs last measured HbA1c (x axis)")

boxplot(plotReportingDF$medianGlu ~ cut(sqrt(plotReportingDF$eAGyyyyDiff^2),breaks=seq(0,30,1)),las=3,varwidth=T,ylim=c(0,20),plot=T,main="IQR vs last measured HbA1c (x axis)")



for (i in seq(1,30,1)) {
  x<-boxplot(subset(plotReportingDF,nCBGperAdmission==i)$medianGlu ~ cut(subset(plotReportingDF,nCBGperAdmission==i)$lastHbA1cInFrame,breaks=seq(30,200,5)),las=3,varwidth=T,ylim=c(0,5),plot=F)
  if (i==1) {plot(c(1:20),x$stats[3,1:20],ylim=c(0,20),pch=16,cex=4*(sqrt(x$n))/max(sqrt(x$n)),col=i,main="admission median glucose vs most recent HbA1c. each line represents an admission with n CBGs. point size corresponding to number of values\nthis controls for the issue of admission duration increasing with increasing numbers of CBGs measured during an admission"); lines(c(1:20),x$stats[3,1:20],col=i,lwd=3)}
  if (i>1)  {points(c(1:20),x$stats[3,1:20],pch=16,cex=4*(sqrt(x$n))/max(sqrt(x$n)),col=i); lines(c(1:20),x$stats[3,1:20],col=i,lwd=3)}
  
}

######################################################
## minimum glucose during admission 

## most recent HbA1c in range (15 months)
boxplot(plotReportingDF$minGlu ~ cut(plotReportingDF$lastHbA1cInFrame,breaks=seq(30,200,10)),las=3,varwidth=T,ylim=c(1,10),plot=T,main="admission minimum glucose vs most recent HbA1c value")

boxplot(plotReportingDF$minGlu ~ cut(plotReportingDF$yyyy,breaks=seq(1,28,1)),las=3,varwidth=T,ylim=c(1,10),plot=T,main="admission minimum glucose vs first CBG value")

boxplot(plotReportingDF$minGlu ~ cut(plotReportingDF$eAGyyyyDiff,breaks=seq(-22,22,2)),las=3,varwidth=T,ylim=c(1,9),plot=T,main="minGlu vs AGN (x axis)")

boxplot(plotReportingDF$minGlu ~ cut(sqrt(plotReportingDF$eAGyyyyDiff^2),breaks=seq(0,30,1)),las=3,varwidth=T,ylim=c(0,8),plot=T,main="minGlu vs distance from AGN==0 (x axis)")

plot(sqrt(plotReportingDF$eAGyyyyDiff^2), plotReportingDF$minGlu)
fit <- lm(plotReportingDF$minGlu ~ sqrt(plotReportingDF$eAGyyyyDiff^2))
abline(fit, col="red")




## plot days with hypo per total days per hba1c, and straightforward proportion of admissions with hypoglycaemia
## plot survival

######################################################
## maximum glucose during admission 

## most recent HbA1c in range (15 months)
boxplot(plotReportingDF$maxGlu ~ cut(plotReportingDF$lastHbA1cInFrame,breaks=seq(30,200,10)),las=3,varwidth=T,ylim=c(10,28),plot=T,main="admission maximum glucose vs most recent HbA1c value")

boxplot(plotReportingDF$maxGlu ~ cut(plotReportingDF$eAGyyyyDiff,breaks=seq(-30,30,2)),las=3,varwidth=T,ylim=c(5,28),plot=T,main="IQR vs last measured HbA1c (x axis)")


######################################################
## propertion of admissions with hypoglycaemia

numberOfDivision<-11

for (hy in seq(3,4,1)) {
  
  hypoThresh<-hy
  
  reportHypoPropDF<-as.data.frame(matrix(0,nrow=numberOfDivision,ncol=6))
  colnames(reportHypoPropDF)<-c("incrementLowerHbA1c","incrementUpperHbA1c","n","hypoN","hypoProp")
  reportHypoPropDF$incrementLowerHbA1c<-quantile(plotReportingDF$lastHbA1cInFrame,prob = seq(0, 1, length = numberOfDivision))
  reportHypoPropDF$incrementUpperHbA1c[1:nrow(reportHypoPropDF)-1]<-reportHypoPropDF$incrementLowerHbA1c[2:nrow(reportHypoPropDF)]
  reportHypoPropDF$incrementUpperHbA1c[nrow(reportHypoPropDF)]<-1000
  
  for (h in seq(1,nrow(reportHypoPropDF),1)) {
    testSub<-subset(plotReportingDF,lastHbA1cInFrame>=reportHypoPropDF$incrementLowerHbA1c[h] & lastHbA1cInFrame<reportHypoPropDF$incrementUpperHbA1c[h])
    testSub$hypoT<-ifelse(testSub$minGlu<hypoThresh,1,0)
    n<-nrow(testSub)
    hypoN<-sum(testSub$hypoT)
    hypoProp<-hypoN/n
    
    reportHypoPropDF$n[h]<-n
    reportHypoPropDF$hypoN[h]<-hypoN
    reportHypoPropDF$hypoProp[h]<-hypoProp
    
  }
  
  plot(reportHypoPropDF$incrementLowerHbA1c[1:nrow(reportHypoPropDF)-1],reportHypoPropDF$hypoProp[1:nrow(reportHypoPropDF)-1],main=paste("proportion of admissions with hypoglycaemia <",hypoThresh,"mmol/L by deciles of last hba1c (x axis)",sep=""),xlab="last measured HbA1c",ylab="proportion of admissions with >=1 hypoglycaemic episodes",pch=16,cex=2)
  
}

dev.off()



plotfilename <- paste("../GlCoSy/plots/hba1c_variability_admissionCharacteristics_type1.pdf",sep="")
pdf(plotfilename, width=16, height=9)


decilesOfParameter<-function(frame,numberOfDivision,hypoThresh,parameter,xlab,ylab) {
  
  # testPlotReportingDF<-subset(frame,inWindow_nVals>2)
  testPlotReportingDF<-frame
  
  reportHypoPropDF<-as.data.frame(matrix(0,nrow=numberOfDivision,ncol=6))
  colnames(reportHypoPropDF)<-c("incrementLowerHbA1c","incrementUpperHbA1c","n","hypoN","hypoProp")
  reportHypoPropDF$incrementLowerHbA1c<-quantile(parameter,prob = seq(0, 1, length = numberOfDivision))
  reportHypoPropDF$incrementUpperHbA1c[1:nrow(reportHypoPropDF)-1]<-reportHypoPropDF$incrementLowerHbA1c[2:nrow(reportHypoPropDF)]
  reportHypoPropDF$incrementUpperHbA1c[nrow(reportHypoPropDF)]<-1000
  
  for (h in seq(1,nrow(reportHypoPropDF),1)) {
    testSub<-subset(testPlotReportingDF,parameter>=reportHypoPropDF$incrementLowerHbA1c[h] & parameter<reportHypoPropDF$incrementUpperHbA1c[h])
    testSub$hypoT<-ifelse(testSub$minGlu<hypoThresh,1,0)
    n<-nrow(testSub)
    hypoN<-sum(testSub$hypoT)
    hypoProp<-hypoN/n
    
    reportHypoPropDF$n[h]<-n
    reportHypoPropDF$hypoN[h]<-hypoN
    reportHypoPropDF$hypoProp[h]<-hypoProp
    
  }
  
  
  plot(reportHypoPropDF$incrementLowerHbA1c[1:nrow(reportHypoPropDF)-1],reportHypoPropDF$hypoProp[1:nrow(reportHypoPropDF)-1],main=paste("proportion of admissions with hypoglycaemia <",hypoThresh,"mmol/L by deciles of x axis parameter",sep=""),xlab=xlab,ylab=ylab,pch=16,cex=2)
  
}

# hypo and hba1c IQR
decilesOfParameter(subset(plotReportingDF,inWindow_nVals>2),11,3,subset(plotReportingDF,inWindow_nVals>2)$inWindow_hbIQR,"IQR of hba1c in 15 month window prior to admission","proportion of admissions with >=1 hypoglycaemic episodes")
decilesOfParameter(subset(plotReportingDF,inWindow_nVals>2),11,4,subset(plotReportingDF,inWindow_nVals>2)$inWindow_hbIQR,"IQR of hba1c in 15 month window prior to admission","proportion of admissions with >=1 hypoglycaemic episodes")

# admission duration and hba1c IQR
boxplot(subset(plotReportingDF,inWindow_nVals>2)$admissionDurationDays ~ cut(subset(plotReportingDF,inWindow_nVals>2)$inWindow_hbIQR,breaks=  quantile(subset(plotReportingDF,inWindow_nVals>2)$inWindow_hbIQR,prob = seq(0, 1, length = 11)) ),las=3,varwidth=T,ylim=c(0,5),plot=T,main="admission duration vs hba1c IQR. 15 month window. deciles")

# admission IQR and hba1c IQR
boxplot(subset(plotReportingDF,inWindow_nVals>2)$IQR ~ cut(subset(plotReportingDF,inWindow_nVals>2)$inWindow_hbIQR,breaks=  quantile(subset(plotReportingDF,inWindow_nVals>2)$inWindow_hbIQR,prob = seq(0, 1, length = 11)) ),las=3,varwidth=T,ylim=c(0,10),plot=T,main="admission IQR vs hba1c IQR. 15 month window. deciles")

# admission IQR and hba1c IQR
boxplot(subset(plotReportingDF,inWindow_nVals>2)$medianGlu ~ cut(subset(plotReportingDF,inWindow_nVals>2)$inWindow_hbIQR,breaks=  quantile(subset(plotReportingDF,inWindow_nVals>2)$inWindow_hbIQR,prob = seq(0, 1, length = 11)) ),las=3,varwidth=T,ylim=c(0,20),plot=T,main="admission median glucose vs hba1c IQR. 15 month window. deciles")



decilesOfParameter(subset(plotReportingDF,inWindow_nVals>2),11,3,subset(plotReportingDF,inWindow_nVals>2)$inWindow_timeIntervalIQR,"time interval IQR of hba1c in 15 month window prior to admission","proportion of admissions with >=1 hypoglycaemic episodes")
decilesOfParameter(subset(plotReportingDF,inWindow_nVals>2),11,4,subset(plotReportingDF,inWindow_nVals>2)$inWindow_timeIntervalIQR,"time interval IQR of hba1c in 15 month window prior to admission","proportion of admissions with >=1 hypoglycaemic episodes")


decilesOfParameter(subset(plotReportingDF,allPrior_hbIQR>2),11,3,subset(plotReportingDF,inWindow_nVals>2)$allPrior_hbIQR,"IQR of all prior hba1c","proportion of admissions with >=1 hypoglycaemic episodes")
decilesOfParameter(subset(plotReportingDF,allPrior_hbIQR>2),11,4,subset(plotReportingDF,inWindow_nVals>2)$allPrior_hbIQR,"IQR of all prior hba1c","proportion of admissions with >=1 hypoglycaemic episodes")

dev.off()


######################################################################################################
## survival plots - come from analysis 3 of _T1DM_hba1c_admissionAndMortality_

simpleSurvivalPlot<-function(inputFrame,parameterToTest,postDischargeStartDay,analysisPlotTitle,quantileDivisions) {
  
  # inputFrame<-survivalPlotReportingDF; parameterToTest<-survivalPlotReportingDF$inWindow_timeIntervalIQR; postDischargeStartDay<-30; analysisPlotTitle<-""; quantileDivisions<-2
  
  SurvivalData<-inputFrame
  #  SurvivalData<-survivalPloReportingDF[eAGyyyyDiff<0]
  
  # SurvivalData<-subset(SurvivalData,diagnosisDateUnix>(-2208988800))
  SurvivalData$diabetesDuration<-SurvivalData$dateplustime1 - SurvivalData$diagnosisDateUnix
  
  DaySeconds<-(60*60*24)
  shortCensorPeriodStartDay  <- DaySeconds*postDischargeStartDay
  shortCensorPeriodEndDay    <- DaySeconds*10000
  
  lastDOD<-max(SurvivalData$deathDateUnix)
  SurvivalData$dateOfDischarge<-SurvivalData$dateplustime1+SurvivalData$admissionDuration
  SurvivalData$deathEvent<-ifelse(SurvivalData$deathDateUnix>0,1,0)
  
  SurvivalData$timeToDeath<-ifelse(SurvivalData$deathEvent==1,(SurvivalData$deathDateUnix-SurvivalData$dateOfDischarge),0)
  #		SurvivalData$timeToDeath<-SurvivalData$timeToDeath/DaySeconds
  SurvivalData$timeToDeathInterval<-ifelse(SurvivalData$deathEvent==0,(lastDOD-SurvivalData$dateOfDischarge),SurvivalData$timeToDeath)
  SurvivalData$timeToDeathInterval[is.na(SurvivalData$timeToDeathInterval)]<-0 # ; SurvivalData<-subset(SurvivalData,timeToDeathInterval>0)
  #		SurvivalData$timeToDeathInterval<-SurvivalData$timeToDeathInterval/DaySeconds
  
  SurvivalData$shortDeathEvent <- SurvivalData$deathEvent
  SurvivalData$shortDeathEvent <- ifelse(SurvivalData$deathEvent==1 & SurvivalData$timeToDeath>=(shortCensorPeriodStartDay) & SurvivalData$timeToDeath<(shortCensorPeriodEndDay),1,0)	
  
  SurvivalData$sexDigit<-ifelse(nchar(SurvivalData$charID==9),as.numeric(substr(SurvivalData$charID,8,8)),as.numeric(substr(SurvivalData$charID,9,9)))
  SurvivalData$sexNumber<-ifelse(SurvivalData$sexDigit%%2==0,1,0)
  SurvivalData$sex<-factor(1*(SurvivalData$sexNumber <1),levels=0:1,labels=c("F","M"))
  
  
  mfitAge50<-survfit(Surv(timeToDeathInterval, shortDeathEvent) ~ (parameterToTest>=quantile(parameterToTest,prob = seq(0, 1, length = (quantileDivisions+1)), type = 5)[quantileDivisions]), data = SurvivalData)
  shortPlotTitle <- paste("Mortality, time ",round(shortCensorPeriodStartDay)/DaySeconds," to ",round(max(SurvivalData$timeToDeathInterval))/DaySeconds,"days\nParameter tested:",analysisPlotTitle,", threshold: ",quantile(parameterToTest)[3],"\nn=",nrow(SurvivalData),". covariables: age, number of hba1c values during period of interest, duration of diabetes",sep="")
  plot(mfitAge50,mark.time=F,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,round(max(SurvivalData$timeToDeathInterval))),lwd=3)
  mfitAge50.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ age+admissionDuration+diabetesDuration+(parameterToTest>=quantile(parameterToTest,prob = seq(0, 1, length = (quantileDivisions+1)), type = 5)[quantileDivisions]), data = SurvivalData)
  pVal <- summary(mfitAge50.coxph)$coef[,5]; HR <- round(exp(coef(mfitAge50.coxph)),2)
  legendText <- paste("p = ",pVal," | HR = ",HR,sep="")
  summarySurvfit <- summary(mfitAge50); legendNames <- row.names(summarySurvfit$table)
  legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("bottomright",legendText,cex=0.6)
  
}


findFirstAdmission<-function(admissionNumberFlag) {
  reportDF_anf<-as.data.frame(matrix(0,nrow=length(admissionNumberFlag),ncol=2))
  colnames(reportDF_anf)<-c("admissionNumberFlag","minFlag")
  
  reportDF_anf$admissionNumberFlag<-admissionNumberFlag
  reportDF_anf$minFlag<-ifelse(reportDF_anf$admissionNumberFlag==min(admissionNumberFlag),1,0)
  
  return(reportDF_anf$minFlag)
}

orderedPloReportingDF<-plotReportingDF[order(plotReportingDF$ID,plotReportingDF$dateplustime1),]
orderedPloReportingDF[, c("flagFirstAdmission") := findFirstAdmission(admissionNumberFlag) , by=.(ID)]

survivalPlotReportingDF<-orderedPloReportingDF[flagFirstAdmission==1]
survivalPlotReportingDF<-subset(survivalPlotReportingDF,diagnosisDateUnix>(-2208988800))

survivalPlotReportingDF$diabetesDuration<-survivalPlotReportingDF$dateplustime1 - survivalPlotReportingDF$diagnosisDateUnix

simpleSurvivalPlot(survivalPlotReportingDF,survivalPlotReportingDF$eAGyyyyDiff,90,"T1DM AGN",2)
simpleSurvivalPlot(survivalPlotReportingDF,sqrt((survivalPlotReportingDF$eAGyyyyDiff - (quantile(survivalPlotReportingDF$eAGyyyyDiff)[3]))^2),60,"T1DM AGN",2)


## ad hoc plot generation for talk
hist(plotReportingDF$eAGyyyyDiff, breaks=seq(-28,28,0.5), yaxt="n", xaxt="n", ylab="", xlab="", main="")
axis(2,cex.axis=2)
axis(1,cex.axis=2)
mtext("Frequency", side=2, line=2.5, cex=2)
mtext("AGN (mmol/l)", side=1, line=2.5, cex=2)

## IQR plots
attdAbstractIQR<-boxplot(plotReportingDF$IQR ~ cut(plotReportingDF$eAGyyyyDiff,breaks=seq(-22,22,2)),las=3,varwidth=T,ylim=c(0,10),plot=T,main="IQR vs AGN ATTD abstract 1 (x axis)", yaxt="n")
axis(2,cex.axis=2)
mtext("IQR (mmol/l)", side=2, line=2.5, cex=2)

    attdAbstractIQR<-boxplot(plotReportingDF$IQR ~ cut(plotReportingDF$yyyy,breaks=seq(0,28,1)),las=3,varwidth=T,ylim=c(0,10),plot=T,main="IQR vs initial CBG", yaxt="n")
    axis(2,cex.axis=2)
    mtext("IQR (mmol/l)", side=2, line=2.5, cex=2)
    
    attdAbstractIQR<-boxplot(plotReportingDF$IQR ~ cut(plotReportingDF$eAG,breaks=seq(0,28,1)),las=3,varwidth=T,ylim=c(0,10),plot=T,main="eAG vs initial CBG", yaxt="n")
    axis(2,cex.axis=2)
    mtext("IQR (mmol/l)", side=2, line=2.5, cex=2)
    
    boxplot(plotReportingDF$IQR ~ cut(sqrt(plotReportingDF$eAGyyyyDiff^2), breaks=seq(0,22,1)), varwidth=T, las=3, ylim=c(0,10), yaxt="n")
    axis(2,cex.axis=2)
    mtext("IQR (mmol/l)", side=2, line=2.5, cex=2)
    
## LOS plots
attdAbstractLOS<-boxplot(subset(plotReportingDF,admissionDurationDays>0.2)$admissionDurationDays ~ cut(subset(plotReportingDF,admissionDurationDays>0.2)$eAGyyyyDiff,breaks=seq(-22,22,4)),las=3,varwidth=T,ylim=c(0,8),plot=T,main="LOS vs AGN ATTD abstract 1 (x axis)", yaxt="n")
axis(2,cex.axis=2)
mtext("LOS (days)", side=2, line=2.5, cex=2)

    attdAbstractLOS<-boxplot(subset(plotReportingDF,admissionDurationDays>0.2)$admissionDurationDays ~ cut(subset(plotReportingDF,admissionDurationDays>0.2)$yyyy,breaks=seq(0,28,1)),las=3,varwidth=T,ylim=c(0,8),plot=T,main="LOS vs AGN ATTD abstract 1 (x axis)", yaxt="n")
    axis(2,cex.axis=2)
      mtext("LOS (days)", side=2, line=2.5, cex=2)
      
      attdAbstractLOS<-boxplot(subset(plotReportingDF,admissionDurationDays>0.2)$admissionDurationDays ~ cut(subset(plotReportingDF,admissionDurationDays>0.2)$eAG,breaks=seq(0,28,1)),las=3,varwidth=T,ylim=c(0,8),plot=T,main="LOS vs AGN ATTD abstract 1 (x axis)", yaxt="n")
      axis(2,cex.axis=2)
      mtext("LOS (days)", side=2, line=2.5, cex=2)
      
      attdAbstractLOS<-boxplot(subset(plotReportingDF,admissionDurationDays>0.2)$admissionDurationDays ~ cut(sqrt(subset(plotReportingDF,admissionDurationDays>0.2)$eAGyyyyDiff^2),breaks=seq(0,28,1)),las=3,varwidth=T,ylim=c(0,8),plot=T,main="LOS vs AGN ATTD abstract 1 (x axis)", yaxt="n")
      axis(2,cex.axis=2)
      mtext("LOS (days)", side=2, line=2.5, cex=2)
      
## minGlu
boxplot(plotReportingDF$minGlu ~ cut(plotReportingDF$eAGyyyyDiff,breaks=seq(-22,22,2)),las=3,varwidth=T,ylim=c(1,9),plot=T,main="minGlu vs AGN (x axis)", yaxt="n")
axis(2,cex.axis=2)
mtext("minimum glucose recorded (mmol/l)", side=2, line=2.5, cex=2)

    attdAbstractIQR<-boxplot(plotReportingDF$minGlu ~ cut(plotReportingDF$yyyy,breaks=seq(0,28,1)),las=3,varwidth=T,ylim=c(0,10),plot=T,main="minGlu vs initial CBG", yaxt="n")
    axis(2,cex.axis=2)
    mtext("minimum glucose (mmol/l)", side=2, line=2.5, cex=2)
    
    attdAbstractIQR<-boxplot(plotReportingDF$minGlu ~ cut(plotReportingDF$eAG,breaks=seq(0,28,1)),las=3,varwidth=T,ylim=c(0,10),plot=T,main="minGlu vs eAG", yaxt="n")
    axis(2,cex.axis=2)
    mtext("minimum glucose (mmol/l)", side=2, line=2.5, cex=2)
      
    
    # mortality
    plotReportingDF$dead <- ifelse(plotReportingDF$deathDateUnix>0,1,0)
    plotReportingDF$inHospitalDeath <- ifelse(plotReportingDF$dead==1 & (plotReportingDF$deathDateUnix - (plotReportingDF$dateplustime1 + plotReportingDF$admissionDuration) < (5*24*60*60)),1,0)
    
    fit <- glm(plotReportingDF$inHospitalDeath ~ sqrt(plotReportingDF$eAGyyyyDiff^2) + plotReportingDF$age)
    summary(fit)
    
    fit <- glm(plotReportingDF$inHospitalDeath ~ plotReportingDF$yyyy + plotReportingDF$age)
    summary(fit)

