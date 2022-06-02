###
###
### Fit GAMS to isotope data from paired short-long vibrissae
###
###

rm(list=ls())
setwd("/Users/brian.brost/Documents/sandbox/zeppelin/vibrissae/")

# save.image("/Users/brian.brost/Documents/sandbox/zeppelin/vibrissae/analysis/paired_analysis.Rdata")
# load("/Users/brian.brost/Documents/sandbox/zeppelin/vibrissae/analysis/paired_analysis.Rdata")

ls()


#####################################################################################
### Libraries and subroutines
#####################################################################################

library(mgcv)
library(gratia)
library(rootSolve)
library(numDeriv)
library(tidyverse)
library(furrr)  # for parallel computings
library(rsample)  # for cross validation functions
library(grid)  # for re-arranging figures for the appendix
library(gridExtra)  # for re-arranging figures for the appendix



#####################################################################################
### Import data
#####################################################################################

dat <- read.csv("data/derived/dat.csv")

idx <- which(dat$id=="26.2")  # seal "26.2" is actually seal "26"
dat$id[idx] <- "26"
dat$length_class[idx] <- "short"  # "re-run" is a "short" whisker

# Retain data for individuals with long and short whiskers only
dat <- dat %>% filter(id%in%c("26","T3","T11","T7"))

# Convert to long format and nest data
dat <- dat %>% gather(key="isotope",value="y",d13C,d15N) %>%
  group_by(id,age_n,isotope,length_class) %>% nest() %>%
  mutate(data=map(data,~.x %>% arrange(length_class,segment_percentile))) %>%
  arrange(isotope,age_n,length_class)

# Plot data
dat %>% unnest(data) %>% 
  ggplot(aes(x=segment_percentile,y=y,color=length_class)) +
  facet_grid(isotope~age_n,scales="free_y") +
  geom_line()


#####################################################################################
### Cross-validation to select value for gamma, the extra penalty term in gam(...)
#####################################################################################

get.gam <- function(splits,k,gamma){  # fit GAM and get predictions, effective degrees of freedom, etc.
  
  train <- analysis(splits)  # training data set
  test <- assessment(splits)  # testing data set
  
  # Fit GAM with factor-smooth interaction
  out <- gam(y~s(segment_percentile,k=k),data=train,method="REML",select=FALSE,gamma=gamma)

  # Get predictions for out-of-sample data
  pred <- predict(out,test,se=TRUE) %>% do.call(cbind,.)  # predictions for test data
  list(pred=data.frame(test %>% select(y,segment_percentile),pred),edf=gratia::edf(out))
}


###
### Setup 10-fold cross validation with 10 repeats
###

cv <- dat %>% mutate(cv=map(data,~.x %>% vfold_cv(group=id,v=10,repeats=10)),
                     cv=map(cv,~.x %>% rename("replicate"="id","fold"="id2"))) %>%
  unnest(cv) %>% mutate(k=map_int(splits,nrow),gamma=1,.before=splits) %>%
  group_by(length_class,age_n,isotope) %>% mutate(k=min(k)-1) %>%  # retain minimum values for k
  complete(gamma=seq(1,2,0.1),nesting(id,length_class,age_n,isotope,data,k,splits,replicate,fold))  # complete values for gamma           

cv %>% print(n=1000)  
cv %>% group_by(length_class,age_n,isotope) %>% summarize(k=unique(k))


###
### Perform cross validation 
###

plan("multisession",workers=availableCores()-1)
start <- Sys.time()  # this 10-fold cross validation with 10 repeats takes 12 minutes to complete
cv <- cv %>% mutate(model=future_pmap(list(splits,k,gamma), ~get.gam(..1,..2,..3)))
difftime(Sys.time(),start,units="min") %>% round()
# save.image("/Users/brian.brost/Documents/sandbox/zeppelin/vibrissae/analysis/paired_analysis.Rdata")
cv


###
### Examine results from cross validation
###

# Boxplot of effective degrees of freedom for all replicates
edf <- cv %>% group_by(age_n,isotope,length_class,gamma,replicate) %>% 
  mutate(edf=map(model,~.x$edf)) %>% unnest(edf) %>% 
  filter(smooth=="s(segment_percentile)") %>%
  summarize(edf=mean(edf,na.rm=TRUE))

edf %>% ggplot(aes(x=as.factor(gamma),y=edf)) +
  facet_grid(isotope+length_class~age_n,scales="free_y") +
  geom_boxplot()
  
# Boxplot of mean square error for all replicates
mse <- cv %>% group_by(age_n,isotope,length_class,gamma,replicate) %>% 
  mutate(pred=map(model,~.x$pred)) %>% unnest(pred) %>% 
  summarize(mse=mean((y-fit)^2,na.rm=TRUE))  # mean MSE across cross validation folds

mse %>% ggplot(aes(x=as.factor(gamma),y=mse)) +
  facet_grid(isotope+length_class~age_n,scales="free") +
  geom_boxplot()

# Plot mean MSE across all replicates
mse <- mse %>% group_by(age_n,isotope,length_class,gamma) %>% summarize(mse=mean(mse))

mse %>% group_by(age_n,isotope,length_class) %>% mutate(min=mse==min(mse)) %>%
  ggplot(aes(x=as.factor(gamma),y=mse,color=min)) +
  facet_grid(isotope+length_class~age_n,scales="free") +
  geom_point()



#####################################################################################
### Fit single-level models
#####################################################################################

# Identify best values for gamma
gamma <- mse %>% filter(mse==min(mse)) %>% select(-mse)
gamma

# Get maximum k for factor smooth
k <- cv %>% group_by(age_n,isotope,length_class) %>% summarize(k=unique(k))
k

out <- dat %>% left_join(k) %>% left_join(gamma) %>%  # add max. dimension for basis expansion and gamma
  mutate(model=pmap(list(data,k,gamma),
            ~gam(y~s(segment_percentile,k=..2),data=..1,method="REML",select=FALSE,gamma=..3)))


#####################################################################################
### Examine model fits
#####################################################################################

# Model checking
par(mfrow=c(2,2))
gam.check(out$model[[5]])
k.check(out$model[[6]])

draw(out$model[[8]])
plot.gam(out$model[[4]])
summary(out$model[[3]])


###
### Get predictions from all models
### 

# Get predictions
preds <- out %>% mutate(newdata=map(data,~.x %>% select(segment_percentile)),.after=data) %>%
  mutate(model=map2(model,newdata,~predict(.x,.y)))

# Convert to long format
preds <- preds %>% select(-newdata) %>% unnest(c(data,starts_with("model"))) %>%
  gather("model","fit",starts_with("model"))


###
### Compare modeled estimates
### 

preds %>% ggplot(aes(x=segment_percentile,y=fit,color=length_class)) +
  facet_grid(isotope~age_n,scales="free") +
  geom_point(aes(y=y,x=segment_percentile)) +
  geom_line(show.legend=TRUE)


###
### Examine observed vs. predicted
###

preds %>% ggplot(aes(x=y,y=fit,color=length_class)) +
  facet_grid(isotope~age_n,scales="free") +
  geom_abline(intercept=0,slope=1,color="black",lty=2) +
  geom_point()


###
### Examine residuals
###

preds %>% ggplot(aes(x=segment_percentile,y=y-fit,color=length_class)) +
  facet_grid(isotope~age_n,scales="free") +
  geom_hline(yintercept=0,color="black",lty=2) +
  geom_point() 



#####################################################################################
### Get predictions and 'zeros' from single-level GAMs
#####################################################################################

get_zeros <- function(mod,nd){  # find locations where first derivative is zero

  interval <- c(min(nd$segment_percentile),max(nd$segment_percentile))  # interval over which to find derivatives
  
  f <- function(x){
    predict(mod,newdata=data.frame(segment_percentile=x))
  }
  
  fx <- function(x){
    grad(f,x)
  }
  
  uniroot.all(fx,interval)  # locations where derivative is zero
}


###
### Get model predictions and zeros
### 

# Get individual-level predictions and zeros
preds <- out %>% select(id,age_n,isotope,length_class,data,model) %>%
  mutate(newdata=map(data,~.x %>% summarize(segment_percentile=
        seq(min(segment_percentile),max(segment_percentile),length.out=100))),.after=data) %>%
  mutate(preds=map2(model,newdata,~predict(.x,.y,type="response",se.fit=TRUE)),
         preds=map(preds,~do.call(cbind,.x) %>% as_tibble()),
         preds=map2(newdata,preds,bind_cols)) %>%
  mutate(zeros=map2(model,newdata,~get_zeros(.x,.y)))


#####################################################################################
### Create figures for publication
#####################################################################################

# See ~/figures/figures.r


#####################################################################################
### Save workspace
#####################################################################################

# save.image("/Users/brian.brost/Documents/sandbox/zeppelin/vibrissae/analysis/paired_analysis.Rdata")



