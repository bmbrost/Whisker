###
###
### Fit GAMS to isotope data from long vibrissae
###
###

rm(list=ls())
setwd("/Users/brian.brost/Documents/sandbox/zeppelin/vibrissae/")

# save.image("/Users/brian.brost/Documents/sandbox/zeppelin/vibrissae/analysis/long_analysis.Rdata")
# load("/Users/brian.brost/Documents/sandbox/zeppelin/vibrissae/analysis/long_analysis.Rdata")

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

# Retain long whiskers only and convert to long format
dat <- dat %>% filter(length_class=="long")

# Convert to long format and nest data
dat <- dat %>% gather(key="isotope",value="y",d13C,d15N) %>%
  group_by(age_n,isotope) %>% nest() %>%
  mutate(data=map(data,~.x %>% mutate(id=forcats::as_factor(as.character(id)))),
         data=map(data,~.x %>% arrange(id,segment_percentile))) %>%
  arrange(isotope,age_n)

# Plot data
dat %>%  unnest(data) %>% 
  ggplot(aes(x=segment_percentile,y=y,color=id)) +
  facet_grid(isotope~age_n,scales="free_y") +
  geom_line()


#####################################################################################
### Cross-validation to select value for gamma, the extra penalty term in gam(...)
#####################################################################################

get.k <- function(x,k=10:100,d=1,k_min=15){  # determine maximum size of basis expansion
  
  train <- analysis(x)  # training data
  n <- table(train$id)  # number observations per individual
  m <- length(n)  # number of individuals
  
  # Total df required for model with factor-smooth interaction
  k_pop <- ifelse(round(k/d)<k_min,k_min,round(k/d))
  df <- m+k_pop+k*m
  
  idx <- df<=sum(n)
  list(k_pop=c(k_pop[idx] %>% last()),k_ind=c(k[idx] %>% last()))
}  

get.gam <- function(splits,k_pop,k_ind,gamma){  # fit GAM and get predictions, effective degrees of freedom, etc.
  
  train <- analysis(splits)  # training data set
  test <- assessment(splits)  # testing data set
  
  # Fit GAM with factor-smooth interaction
  out <- gam(y~s(segment_percentile,k=k_pop)+s(segment_percentile,id,bs='fs',m=1,k=k_ind),
                  data=train,method="REML",select=FALSE,gamma=gamma)

  # Get predictions for out-of-sample data
  pred <- predict(out,test,se=TRUE) %>% do.call(cbind,.)  # predictions for test data
  list(pred=data.frame(test %>% select(y,segment_percentile),pred),edf=gratia::edf(out))
}


###
### Setup 10-fold cross validation with 10 repeats
###

cv <- dat %>% mutate(cv=map(data,~.x %>% vfold_cv(group=id,v=10,repeats=10,strata=id))) %>%
  unnest(cv) %>% rename("replicate"="id","fold"="id2") %>%
  mutate(k=case_when(  # identify basis expansion for population- and individual-level smooth
                  age_n==0 ~ map(splits, ~get.k(.x,k=10:100,d=2,k_min=5,fs=TRUE)),  # k for age_n==0
                  age_n>0 ~ map(splits, ~get.k(.x,k=10:100,d=2,k_min=15,fs=TRUE))),  # k for age_n > 0
         k_pop=map_dbl(k,~.x$k_pop),k_ind=map_dbl(k,~.x$k_ind),gamma=1,.before=splits) %>% select(-k) %>%
  group_by(age_n,isotope) %>% mutate(k_pop=min(k_pop),k_ind=min(k_ind)) %>%   # retain minimum values for k
  complete(gamma=seq(1,2,0.1),nesting(age_n,isotope,data,k_pop,k_ind,splits,replicate,fold))  # complete values for gamma

cv %>% print(n=1000)  
cv %>% group_by(age_n,isotope) %>% summarize(k_pop=unique(k_pop),k_ind=unique(k_ind))


###
### Perform cross validation 
###

plan("multisession",workers=availableCores()-1)
start <- Sys.time()  # this 10-fold cross validation with 10 repeats takes 2.6 hours to complete
cv <- cv %>% mutate(model_3=future_pmap(list(splits,k_pop,k_ind,gamma), ~get.gam(..1,..2,..3,..4,fs=TRUE)))
difftime(Sys.time(),start,units="min") %>% round()
# save.image("/Users/brian.brost/Documents/sandbox/zeppelin/vibrissae/analysis/long_analysis.Rdata")
cv


###
### Examine results from cross validation
###

# Boxplot of effective degrees of freedom for all replicates
edf <- cv %>% group_by(age_n,isotope,gamma,replicate) %>% 
  mutate(edf=map(model_3,~.x$edf)) %>% unnest(edf) %>% 
  filter(smooth=="s(segment_percentile)") %>%
  # filter(smooth=="s(segment_percentile,id)") %>%
  summarize(edf=mean(edf,na.rm=TRUE))

edf %>% ggplot(aes(x=as.factor(gamma),y=edf)) +
  facet_grid(isotope~age_n,scales="free_y") +
  geom_boxplot()
  
# Boxplot of mean square error for all replicates
mse <- cv %>% group_by(age_n,isotope,gamma,replicate) %>% 
  mutate(pred=map(model_3,~.x$pred)) %>% unnest(pred) %>% 
  summarize(mse=mean((y-fit)^2,na.rm=TRUE))  # mean MSE across cross validation folds

mse %>% ggplot(aes(x=as.factor(gamma),y=mse)) +
  facet_grid(isotope~age_n,scales="free") +
  geom_boxplot()

# Plot mean MSE across all replicates
mse <- mse %>% group_by(age_n,isotope,gamma) %>% summarize(mse=mean(mse))

mse %>% group_by(age_n,isotope) %>% mutate(min=mse==min(mse)) %>%
  ggplot(aes(x=as.factor(gamma),y=mse,color=min)) +
  facet_grid(isotope~age_n,scales="free") +
  geom_point()



#####################################################################################
### Fit single- and multi-level models
#####################################################################################

# Identify best values for gamma
gamma <- mse %>% filter(mse==min(mse)) %>% select(-mse)
gamma

# Get maximum k for factor smooth
k <- cv %>% group_by(age_n,isotope) %>% summarize(k_pop=unique(k_pop),k_ind=unique(k_ind))
k

out <- dat %>% left_join(k) %>% left_join(gamma) %>%  # add max. dimension for basis expansion and gamma
  mutate(model_1=pmap(list(data,k_ind,gamma),  # constant model
                      ~gam(y~s(segment_percentile,k=..2),
                           data=..1,method="REML",select=FALSE,gamma=..3)),
         model_2=pmap(list(data,k_ind,gamma),  # model with random effect
                      ~gam(y~s(id,bs="re")+s(segment_percentile,k=..2),
                           data=..1,method="REML",select=FALSE,gamma=..3)),
         model_3=pmap(list(data,k_pop,k_ind,gamma),  # factor-smooth interaction
                      ~gam(y~s(segment_percentile,k=..2)+s(segment_percentile,id,bs='fs',m=1,k=..3),
                           data=..1,method="REML",select=FALSE,gamma=..4)))



#####################################################################################
### Examine model fits
#####################################################################################

# mod <- gam(y~s(segment_percentile)+s(segment_percentile,id,bs='fs'),data=dat$data[[5]],method="REML",select=FALSE)
# gam.check(mod)
# k.check(mod)
# 
# dat$data[[5]] %>% select(id,y,segment_percentile) %>% na.omit() %>% mutate(fit=predict(mod)) %>%
#   ggplot(aes(x=segment_percentile,y=y-fit)) +
#   facet_grid(~id) +
#   geom_hline(yintercept=0,color="black",lty=2) +
#   geom_line() 
# 
# dat$data[[5]] %>% select(id,y,segment_percentile) %>% na.omit() %>% mutate(fit=predict(mod)) %>%
#   ggplot(aes(x=fit,y=y-fit)) +
#   facet_grid(~id) +
#   geom_abline(slope=1,yintercept=0,color="black",lty=2) +
#   geom_point() 

# Model checking
par(mfrow=c(2,2))
gam.check(out$model_3[[1]])
k.check(out$model_3[[6]])

draw(out$model_3[[8]])
plot.gam(out$model_3[[4]])
summary(out$model_3[[3]])


###
### Get predictions from all models
### 

# Get predictions
preds <- out %>% mutate(newdata=map(data,~.x %>% select(id,segment_percentile)),.after=data) %>%
  mutate(model_1=map2(model_1,newdata,~predict(.x,.y)),
         model_2=map2(model_2,newdata,~predict(.x,.y)),
         model_3=map2(model_3,newdata,~predict(.x,.y)))

# Convert to long format
preds <- preds %>% select(-newdata) %>% unnest(c(data,starts_with("model"))) %>%
  gather("model","fit",starts_with("model"))


###
### Compare modeled estimates
### 

# Nitrogen
preds %>% filter(age_n==0,isotope=="d15N") %>% 
  ggplot(aes(x=segment_percentile,y=fit,color=model)) +
  facet_wrap(~id) + 
  geom_point(aes(y=y,x=segment_percentile),color="black") +
  geom_line(show.legend=TRUE)

# Carbon
preds %>% filter(age_n==2,isotope=="d13C") %>%
  ggplot(aes(x=segment_percentile,y=fit,color=model)) +
  facet_wrap(~id) + 
  geom_point(aes(y=y,x=segment_percentile),color="black") +
  geom_line(show.legend=TRUE)


###
### Examine observed vs. predicted
###

# Nitrogen
preds %>% filter(age_n==0,isotope=="d15N") %>%
  ggplot(aes(x=y,y=fit)) +
  facet_grid(id~model) +
  geom_abline(intercept=0,slope=1,color="black",lty=2) +
  geom_point()

# Carbon
preds %>% filter(age_n==4,isotope=="d13C") %>%
  ggplot(aes(x=y,y=fit)) +
  facet_grid(id~model) +
  geom_abline(intercept=0,slope=1,color="black",lty=2) +
  geom_point()


###
### Examine residuals
###

# Nitrogen
preds %>% filter(age_n==0,isotope=="d15N") %>% 
  ggplot(aes(x=segment_percentile,y=y-fit)) +
  facet_grid(id~model) +
  geom_hline(yintercept=0,color="black",lty=2) +
  geom_point() 

# Carbon
preds %>% filter(age_n==2,isotope=="d13C") %>%
  ggplot(aes(x=segment_percentile,y=y-fit)) +
  facet_grid(id~model) +
  geom_hline(yintercept=0,color="black",lty=2) +
  geom_point() 



#####################################################################################
### Get predictions and 'zeros' from multi-level GAM with factor-smooth interaction
#####################################################################################

get_zeros <- function(mod,nd,exclude=NULL){  # find locations where first derivative is zero
  out <- data.frame(id=as.character(),zeros=numeric())
  id <- unique(nd$id)
  
  for(i in id){  # loop through individual whiskers
    
    tmp <- nd %>% filter(id==i)  # data for individual whisker
    interval <- c(min(tmp$segment_percentile),max(tmp$segment_percentile))  # interval over which to find derivatives

    f <- function(x){
      predict(mod,newdata=data.frame(segment_percentile=x,id=i),exclude=exclude)
    }
    fx <- function(x){
      grad(f,x)
    }
    zeros <- uniroot.all(fx,interval)  # locations where derivative is zero
 
    if(length(zeros)>0){
      out <- rbind(out,data.frame(id=i,zeros=zeros))  
    }
  }  
  out %>% mutate(id=factor(id,levels=levels(nd$id)))
}


###
### Get individual- and population-level predictions and zeros
### 

plan("multisession",workers=availableCores()-1)

# Get individual-level predictions and zeros
preds <- out %>% select(age_n,isotope,data,model_3) %>% rename("model"="model_3") %>%
  mutate(newdata=map(data,~.x %>% group_by(id) %>% summarize(segment_percentile=
        seq(min(segment_percentile),max(segment_percentile),length.out=100))),.after=data) %>%
  mutate(id_preds=map2(model,newdata,~predict(.x,.y,type="response",se.fit=TRUE)),
         id_preds=map(id_preds,~do.call(cbind,.x) %>% as_tibble()),
         id_preds=map2(newdata,id_preds,bind_cols)) %>%
  mutate(id_zeros=future_map2(model,newdata,~get_zeros(.x,.y)))

# Get population-level predictions
preds <- preds %>% mutate(newdata=map(data,~.x %>% ungroup() %>% summarize(id=first(id),
         segment_percentile=seq(min(segment_percentile),max(segment_percentile),length.out=100))),.after=data) %>%
  mutate(pop_preds=map2(model,newdata,~predict(.x,.y,type="response",se.fit=TRUE,exclude=c("s(segment_percentile,id)"))),
         pop_preds=map(pop_preds,~do.call(cbind,.x) %>% as_tibble()),
         pop_preds=map2(newdata,pop_preds,~cbind(.x %>% select(-id),.y))) %>%
  mutate(pop_zeros=future_map2(model,newdata,~get_zeros(.x,.y,exclude=c("s(segment_percentile,id)")) %>% select(-id)))


###
### Write locations where first derivative is zero for secondary analysis
###

# Individual-level zeros
write.csv(preds %>% select(age_n,isotope,id_zeros) %>% unnest(id_zeros),"results/id_zeros.csv",row.names=FALSE)

# Population-level zeros
write.csv(preds %>% select(age_n,isotope,pop_zeros) %>% unnest(pop_zeros),"results/pop_zeros.csv",row.names=FALSE)


#####################################################################################
### Create figures for publication
#####################################################################################

# See ~/figures/figures.r


#####################################################################################
### Save workspace
#####################################################################################

# save.image("/Users/brian.brost/Documents/sandbox/zeppelin/vibrissae/analysis/long_analysis.Rdata")



