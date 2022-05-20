#rm(list=ls())


install.packages("tidyverse")
library(tidyverse)
# Load packages



setwd("C:/Juv_vib")
dir()
all.summary <- read.csv("animal.summary.csv", header = T)
head(all.summary)

# Tidy data
all.summary <- all.summary %>% filter(id!="dec") %>% droplevels() 
all.summary$class <- ifelse(all.summary$class=="pup","Pup","Juvenile") %>% 
	factor(.,levels=c("Pup","Juvenile"))


###############################################################################
### Figure 1
###############################################################################

pdf("Fig1.pdf")
par(mfrow = c(2,2))


###
### Vibrissae length vs. age class
###

boxplot(vibrissae.length~class,data=all.summary,ylab="Vibrissae length (mm)",xlab="Age class",las=1)


# T-test
mod <- t.test(vibrissae.length~class,data=all.summary,alternative="two.sided")
mod$p.value  # the age classes are significantly different

# Linear regression as an alternative to the t-test above
mod <- lm(vibrissae.length~class,data=all.summary)  
summary(mod)  # coefficients quantify the difference between age classes

# Predictions from linear regression (for reporting in manuscript?)
dat <- data.frame(class=c("Pup","Juvenile"))  # new data for prediction
predict(mod,dat,se.fit=TRUE)[1:2] %>% data.frame() %>%  # mean and standard error for each age class
	mutate(lcl=fit-2*se.fit,ucl=fit+2*se.fit)  # calculate confidence intervals for each age class


#add labels to the upper left corner, you can use this code after each individual figure:
text(x=par("usr")[1],y=par("usr")[4],labels="a.)",adj=c(-0.5,2))


###	
### Vibrissae length vs. age
###

plot(vibrissae.length~age,data=all.summary,ylab="Vibrissae length (mm)", xlab="Age (years)",las=1)


# Linear regression
mod <- lm(vibrissae.length~age,data=all.summary)
summary(mod)  # coefficient for 'age' quantifies the increase in length for each additional year

# Predictions from linear regression for plotting
dat <- data.frame(age=seq(0,4,by=0.1))  # new data for prediction
preds <- predict(mod,dat,se.fit=TRUE)[1:2] %>% data.frame() %>%  # mean and standard error for age
	mutate(lcl=fit-2*se.fit,ucl=fit+2*se.fit)  # calculate confidence intervals

# Add model results to plot
lines(dat$age,preds$fit,lty=1)  # add model fit
polygon(c(dat$age,rev(dat$age)),c(preds$lcl,rev(preds$ucl)),border=NA,col=rgb(0,0,0,0.15))  # add confidence band

text(x=par("usr")[1],y=par("usr")[4],labels="b.)",adj=c(-0.5,2))


###
### Vibrissae length vs. animal.length
###

plot(vibrissae.length~animal.length,data=all.summary,ylab="Vibrissae length (mm)",xlab="Animal length (cm)",las=1)
text(x=par("usr")[1],y=par("usr")[4],labels="c.)",adj=c(-0.5,2))

# Linear regression
mod <- lm(vibrissae.length~animal.length,data=all.summary)
summary(mod)  # coefficient for 'animal.length' quantifies the increase in whisker length for each additional unit of animal length

# Predictions from linear regression for plotting
dat <- data.frame(animal.length=seq(min(all.summary$animal.length,na.rm=TRUE),  # new data for prediction
	max(all.summary$animal.length,na.rm=TRUE),by=1))
preds <- predict(mod,dat,se.fit=TRUE)[1:2] %>% data.frame() %>%  # mean and standard error for animal length
	mutate(lcl=fit-2*se.fit,ucl=fit+2*se.fit)  # calculate confidence intervals

# Add model results to plot
lines(dat$animal.length,preds$fit,lty=1)  # add model fit to plot
polygon(c(dat$animal.length,rev(dat$animal.length)),c(preds$lcl,rev(preds$ucl)),border=NA,col=rgb(0,0,0,0.15))  # add confidence band


###
### Vibrissae length vs. animal mass
###

plot(vibrissae.length~animal.mass,data=all.summary,ylab="Vibrissae length (mm)",xlab="Animal mass (kg)",las=1)
text(x=par("usr")[1],y=par("usr")[4],labels="d.)",adj=c(-0.5,2))

# Linear regression
mod <- lm(vibrissae.length~animal.mass,data=all.summary)
summary(mod)  # coefficient for 'animal.mass' quantifies the increase in whisker length for each additional unit of animal mass

# Predictions from linear regression for plotting
dat <- data.frame(animal.mass=seq(min(all.summary$animal.mass,na.rm=TRUE),  # new data for prediction
	max(all.summary$animal.mass,na.rm=TRUE),by=0.1))
preds <- predict(mod,dat,se.fit=TRUE)[1:2] %>% data.frame() %>%  # mean and standard error for animal length
	mutate(lcl=fit-2*se.fit,ucl=fit+2*se.fit)  # calculate confidence intervals

# Add model results to plot
lines(dat$animal.mass,preds$fit,lty=1)  # add model fit to plot
polygon(c(dat$animal.mass,rev(dat$animal.mass)),c(preds$lcl,rev(preds$ucl)),border=NA,col=rgb(0,0,0,0.15))  # add confidence band


dev.off()

