######################################################################################
##                              HW 4:  Survival Analysis 
######################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~Load the necessory Libraries~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(survival)
#install.packages("OIsurv")
library(OIsurv) # to get the simultaneous confidence ban
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~ Load the data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pdat=read.table("http://people.usm.maine.edu/cpeng/STA587/SurvHWdataset.txt",header=T)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

######################################################################################
# Kaplan-Meier Estimation
######################################################################################
fit_km_plain = survfit(Surv(t, status)~1, data = pdat, type = c("kaplan-meier"), conf.type = c("plain"))

# direct output to a file 
sink("fit_km_plain_summary.txt", append=FALSE, split=FALSE)
summary(fit_km_plain)
# return output to the terminal 
sink()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot Kaplan-Meier Estimation and Point wise and simulatenous confidence Interval
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(fit_km_plain, main = "Kaplan-Meier Estimation curve", xlab = "months", ylab = "S(t)", col=c(4,2,2))
surv.obj <- Surv(pdat$t, pdat$status)
my.cb <- confBands(surv.obj, confLevel=0.95, type="hall")
lines(my.cb$time, my.cb$lower, lty=3, type="s", col=6)
lines(my.cb$time, my.cb$upper, lty=3, type="s", col=6)
legend(100, 1, legend=c("K-M survival estimate", "pointwise intervals","confidence bands"), lty=c(1,2,3), col = c(4,2,6))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Analyse KM plots for different conf.type 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fit_km_log = survfit(Surv(t, status)~1, data = pdat, type = c("kaplan-meier"), conf.type = c("log"))
fit_km_log_log = survfit(Surv(t, status)~1, data = pdat, type = c("kaplan-meier"), conf.type = c("log-log"))

plot(fit_km_plain, main = "Kaplan-Meier Estimation curve", xlab = "months", ylab = "S(t)", col=c(4,2,2))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#get log lower and upper CI 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
log_lower = summary(fit_km_log)$lower
log_upper = summary(fit_km_log)$upper
life.time = summary(fit_km_log)$time
lines(life.time,log_lower, type="s", lty=3, col=3)
lines(life.time,log_upper, type="s", lty=3, col=3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#get log-log lower and upper CI 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

log_log_lower = summary(fit_km_log_log)$lower
log_log_upper = summary(fit_km_log_log)$upper
life.time = summary(fit_km_log_log)$time
lines(life.time,log_log_lower, type="s", lty=4, col=9)
lines(life.time,log_log_upper, type="s", lty=4, col=9)

legend(80, 1, legend=c("K-M survival estimate", "plain CI","log CI", "Log-Log CI"), lty=c(1,2,3,4), col = c(4,2,3,9))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mean Survival Time
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# direct output to a file 
sink("mean_survival_hW4.txt", append=FALSE, split=FALSE)
cat("Mean survival time:")
print(fit_km_plain, print.rmean = TRUE )
print(fit_km_log_log, print.rmean = TRUE )
print(fit_km_log, print.rmean = TRUE )

cat("Confiden Interval for mean Survival time :")
lower_Mean = 72.58 - 1.96 * 2.85
upper_Mean = 72.58 + 1.96 * 2.85
cat("Confiden Interval for mean is: [",lower_Mean,",",upper_Mean,"]")

# return output to the terminal 
sink()

########################################################################################
## Kaplan-Meier Estimation of the survival function for the patient stratified by stages
########################################################################################
colnames(pdat)

fit_km_plain_stage = survfit(Surv(t, status)~stage, data = pdat, type = c("kaplan-meier"), conf.type = c("plain"))
sink("KM_stratified.txt", append=FALSE, split=FALSE)
cat("Kaplan-Meier Estimation of the survival function for the patient stratified by stages:\n")
summary(fit_km_plain_stage)
sink()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# KM plots Stratefied by patient stages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(fit_km_plain_stage, col= c(1:4), xlab="time in months", ylab= "Survival Probability")
legend(100, 1, c("stage 1","stage 2","stage 3","stage 4"), col=c(1:4), lwd=0.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Log Rank test
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fit_km_stage.test=survdiff(Surv(t, status)~stage, data = pdat)
sink("LogRankTest.txt", append=FALSE, split=FALSE)
cat("output for fit_km_plain_stage using survdff()\n")
fit_km_stage.test
cat("output for fit_km_plain_stage$var using survdff()\n")
fit_km_stage.test$var
sink()


########################################################################################
## Nelson and Aalen Estimation of the Hazard function 
########################################################################################

NA.H01=cumsum(fit_km_plain$n.event/fit_km_plain$n.risk)
NA_noDup = NA.H01[which(!duplicated(NA.H01))]
# calculate variance 
n.event <- fit_km_plain$n.event
n.risk <- fit_km_plain$n.risk

NA.V01=cumsum(n.event/n.risk^2)
NA_V = NA.V01[which(!duplicated(NA.V01))]

Lower.CI = NA_noDup - (1.96 * sqrt(NA_V))
Upper.CI = NA_noDup + (1.96 * sqrt(NA_V))
hazard_estimation = cbind("interval"= life.time, "hazard rate" = NA_noDup, "variance_H" = NA_V, "lower" = Lower.CI, "upper" = Upper.CI)

sink("hazardFunction.txt", append=FALSE, split=FALSE)
cat("Nelson and Aalen Estimation of the Hazard function \n")
hazard_estimation
sink()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot the Graph for Nelson and Aalen Estimation of the Hazard function along 95% CI
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
idx = which(duplicated(NA.H01))
life.time = c(0, fit_km_plain$time[-idx])
NA_noDup = c(0, NA_noDup)
lower = c(0, Lower.CI)
upper = c(0, Upper.CI)
plot(life.time, NA_noDup , type="s", lty=1, col=1, xlab = "months", ylab = "H(t)")
lines(life.time, lower, type="s", lty=2, col=3)
lines(life.time, upper, type="s", lty=2, col=3)

#####################################################################################
##Simple Bar Plot 
#####################################################################################

counts <- table(pdat$status)
barplot(counts, main="patients distribution by satus", 
        xlab="0 - censored, 1 - death", col="blue", width = c(0.5, 0.5))

counts <- table(pdat$stage)
barplot(counts, main="patients distribution by satus", 
        xlab="stage", col=c(1, 2, 3, 4))

print(fit_km_stage.test, print.rmean = TRUE )


cat("Confiden Interval for mean Survival time :")
lower_Mean = 72.58 - 1.96 * 2.85
upper_Mean = 72.58 + 1.96 * 2.85
cat("Confiden Interval for mean is: [",lower_Mean,",",upper_Mean,"]")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To observe the density plots for the different confidence intervals
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(2,2))

LB_plain <- summary(fit_km_plain)$lower
UB_plain <- summary(fit_km_plain)$upper
LB_log <- summary(fit_km_log)$lower
UB_log <- summary(fit_km_log)$upper
LB_log_log <- summary(fit_km_log_log)$lower
UB_log_log <- summary(fit_km_log_log)$upper

library(lattice)
densityplot(~ LB_plain + UB_plain, auto.key = TRUE,xlab = "lower/upper CI bounds")
densityplot(~ LB_log + UB_log, auto.key = TRUE,xlab = "lower/upper CI bounds")
densityplot(~ LB_log_log + UB_log_log, auto.key = TRUE,xlab = "lower/upper CI bounds")
