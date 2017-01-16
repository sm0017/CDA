
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Final Project - Logistic Regression Analysis
# Smita Sukhhadeve
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
diabetes=  read.csv("E:/CDA/Final_project/diabetes.csv")
cat("Dimension of data", dim(diabetes))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data Preparation 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create bp_status

diabetes$bpStatus  = factor(NA, levels = c('low', 'normal', 'prehigh', 'high'))
BloodPressure <- diabetes$BloodPressure
diabetes$bpStatus[BloodPressure < 60 ] <- 'low'
diabetes$bpStatus[BloodPressure >= 60 & BloodPressure <80 ] <- 'normal'
diabetes$bpStatus[BloodPressure >= 80 & BloodPressure <90 ] <- 'prehigh'
diabetes$bpStatus[BloodPressure >= 90 ] <- 'high'
rm(BloodPressure) 


# Create age_grp

diabetes$ageGrp  = factor(NA, levels = c('age35', 'age50', 'age>50'))
age <- diabetes$Age
diabetes$ageGrp[age >= 20 & age < 35 ] <- 'age35'
diabetes$ageGrp[age >= 35 & age <50] <- 'age50'
diabetes$ageGrp[age >= 50] <- 'age>50'
rm(age)   

# Create bmi_grp
# https://www.nhlbi.nih.gov/health/educational/lose_wt/BMI/bmi_dis.htm

diabetes$bmiGrp  = factor(NA, levels = c('underwt', 'normal', 'overwt', 'obess' ))
bmi <- diabetes$BMI
diabetes$bmiGrp[bmi <18.5 ] <- 'underwt'
diabetes$bmiGrp[bmi >= 18.5 & bmi < 25] <- 'normal'
diabetes$bmiGrp[bmi >= 25 & bmi <30] <- 'overwt'
diabetes$bmiGrp[bmi >= 30] <- 'obess'

#Convert categorical as factor
diabetes$bpStatus <- as.factor(diabetes$bpStatus)
diabetes$bmiGrp <- as.factor(diabetes$bmiGrp)
diabetes$ageGrp <- as.factor(diabetes$ageGrp)  
diabetes$Outcome <- as.factor(diabetes$Outcome)


#### Remove Age and BloodPressures, BMI columns

diabetes$BloodPressure <- NULL
diabetes$Age <- NULL
diabetes$BMI <- NULL

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data exploration 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(diabetes, 3)
str(diabetes)
summary(diabetes[,c(1,2,3,4,5,6,8,7)])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data exploration 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cor.display <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(diabetes[,2:8],
      lower.panel=panel.smooth, upper.panel=cor.display, 
      pch=20, main="Pairwise Scatterplot - Diabetes Data", col="cornflowerblue")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create train and test(For new value prediction)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
smp_size <- floor(0.85 * nrow(diabetes))
## set the seed to make your partition reproductible
set.seed(123)
train_ind <- sample(seq_len(nrow(diabetes)), size = smp_size)
train <- diabetes[train_ind, ]
test <- diabetes[-train_ind, ]
cat("Dimension of test set:", dim(test), "\n")
cat("Dimension of train set:", dim(train))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Backward and Subset Selection
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
full.model = glm(Outcome ~.,data = train, family = binomial("logit"))
null.model = glm(Outcome ~ 1, data = train, family = binomial("logit"))
### Backward Variable Selection Approach: We give full model as input
back.select = step(full.model, data = train, direction = 'backward')
sink("C:/Users/Smita/Desktop/STAbkwdSel.txt")
summary(back.select)
sink()

### Best subset Selection
#install.packages('bestglm', repos='http://cran.us.r-project.org')

library(bestglm)
colnames(train)[6]
colnames(train)[6]<- "y"
subset_selection <- train[, c("Pregnancies", "Glucose", "SkinThickness", "Insulin", "DiabetesPedigreeFunction", "bp_status", 
                              "age_grp" ,"bmi_grp", "y" )]

## Perform
best.subset.model <- bestglm(Xy = subset_selection, family = binomial, IC = "AIC", method = "exhaustive")
summary(best.subset.model$BestModel)
sink("E:/STA/subset.txt")
sink()
colnames(train)[6]<- "Outcome"

colnames(train)

Model1 = glm(Outcome ~ Pregnancies + Glucose + Insulin + DiabetesPedigreeFunction + 
                       bmiGrp, data = train, family = binomial("logit")
             )
summary(Model1)

cat("model 1 Diagnostis\n")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model Diagnotisc Test - Model 1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1. Pearson Goodness of fit test 

cat("Pearson goodness of fit test", 1-pchisq(Model1$deviance, Model1$df.residual), "\n")

#2. Deviance test using Anova
anova(Model1, test = 'Chisq')


anova(Model1)

colnames(train)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model With Interaction term
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

interaction.model = glm(Outcome ~ Pregnancies + Glucose + Insulin + DiabetesPedigreeFunction +
                          bmiGrp + ageGrp + bpStatus + ageGrp:bmiGrp + 
                          ageGrp*bpStatus + Glucose:Insulin + ageGrp:Insulin , 
                        data = train, family = binomial("logit"))

back.model = step(interaction.model, data = train, direction = 'backward')
summary(back.model)

Model2 = glm(formula = Outcome ~ Pregnancies + Glucose + Insulin + DiabetesPedigreeFunction + 
            bmiGrp + ageGrp + Insulin:ageGrp, family = binomial("logit"), 
            data = train)
summary(Model2)

anova(Model2)

library(ResourceSelection)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model Diagnotisc Test - Model 2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1. Pearson Goodness of fit test 

cat("Pearson goodness of fit test", 1-pchisq(Model2$deviance, Model2$df.residual), "\n")

#2. Deviance test using Anova
cat("Deviance test using Anova for Model 2:\n")
anova(Model2, test = 'Chisq')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compare Model1 and Model 2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2. Deviance test using Anova
cat("Comparision between Model 1 and Model 2:\n")
anova(Model1, Model2, test = 'Chisq')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Perform Lack of fit test
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ResourceSelection)
h1=hoslem.test(Model2$y, fitted(Model2), g=10)
cat("result of lack of fit test:",  "\n")
h1


?hoslem.test

result_Hosmer = cbind("Observed" = h1$observed, "Expected"= h1$expected)
result_Hosmer

##coefficient
dat = rbind(t(coef(Model2)), exp(coef(Model2)))
coefficient_exp = t(dat)
colnames(coefficient_exp) = c("reg.coeff", "exp_reg.coeff")
coefficient_exp

colnames(test)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Predict New values : test
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
testVar = c('Pregnancies', 'Glucose', 'Insulin' , 'DiabetesPedigreeFunction', 'bmiGrp', 'ageGrp')
preddat <- predict(Model2, newdata=subset(test, select = testVar), se.fit=TRUE)
testProb <- exp(preddat$fit)/(1+exp(preddat$fit))
testPron.uci <- exp(preddat$fit+1.96*preddat$se.fit)/(1+exp(preddat$fit+1.96*preddat$se.fit))
testPron.lci <- exp(preddat$fit-1.96*preddat$se.fit)/(1+exp(preddat$fit-1.96*preddat$se.fit))

a = sort(test$DiabetesPedigreeFunction)
b = sort(testProb)
b.ui=sort(testPron.uci)
b.li=sort(testPron.lci)

png("E:/plot.png")
plot(a, test$Outcome, type="n", ylim=c(0, 1), 
     ylab="Probability of Occurance of Diabetes", 
     xlab="probabilities")
lines(a, b, col="blue")
lines(a, b.ui, lty=2)
lines(a, b.li, lty=2)
dev.off()

## odds ratios and 95% CI
cat("Odd Ratio with 95% confidence interval")
exp(cbind(OR = coef(Model2), confint(Model2)))

fitted.results <- ifelse(testProb  > 0.5,1,0)
misClasificError <- mean(fitted.results != test$Outcome)
print(paste('Accuracy',1-misClasificError))

install.packages('ROCR', repos='http://cran.us.r-project.org')

library(ROCR)
prob <- predict(Model2, newdata=subset(test, select = testVar), type="response")
pr <- prediction(prob, test$Outcome)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc
