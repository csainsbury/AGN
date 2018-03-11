# First we start with loading the dataset
data("melanoma", package = "boot")

# Then we munge it according to ?boot::melanoma
library(dplyr)
library(magrittr)
melanoma %<>% 
  mutate(status = factor(status,
                         levels = 1:3,
                         labels = c("Died from melanoma", 
                                    "Alive", 
                                    "Died from other causes")),
         ulcer = factor(ulcer,
                        levels = 0:1,
                        labels = c("Absent", "Present")),
         time = time/365.25, # All variables should be in the same time unit
         sex = factor(sex,
                      levels = 0:1,
                      labels = c("Female", "Male")))

library(survival)
regular_model <- coxph(Surv(time, status == "Died from melanoma") ~
                         age + sex + year + thickness + ulcer,
                       data = melanoma,
                       x = TRUE, y = TRUE)
summary(regular_model)

spl_melanoma <-
  melanoma %>% 
  timeSplitter(by = .5,
               event_var = "status",
               event_start_status = "Alive",
               time_var = "time",
               time_related_vars = c("age", "year"))

interval_model <-
  update(regular_model, 
         Surv(Start_time, Stop_time, status == "Died from melanoma") ~ .,
         data = spl_melanoma)

summary(interval_model)


## 
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
summary(fit1.coxph)

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
summary(time_int_model)

cox.zph(time_int_model)
