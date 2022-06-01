###
###
### Setup vibrissae isotope data for modeling
###
###

rm(list=ls())
setwd("/Users/brian.brost/Documents/sandbox/zeppelin/vibrissae/")


###
### Libraries and subroutines
###

library(mgcv)
library(gratia)
library(rootSolve)
library(numDeriv)
library(tidyverse)


###
### Import data
###

dat <- read.csv("data/raw/all.csv")


###
### Tidy data
###

# Retain relevant variables
dat <- dat %>% select(id,class,age,long.short,segment,seg.length,seg.end.mm,d13C,d15N)

# Rename variables
names(dat) <- c("id","age_class","age","length_class","segement","segment_length","segment_end","d13C","d15N")

# Define segment midpoints and percentiles
dat <- dat %>% group_by(id,length_class) %>% mutate(segment_midpoint=segment_end-0.5*segment_length,  # segment midpoint
         segment_percentile=1-segment_midpoint/max(segment_end),.after=segment_end)   # segment percentile
   
# Define age to nearest year
dat <- dat %>% mutate(age_n=round(age),.after=age)

# Exclude data for some individuals
dat %>% gather(key="isotope",value="x",d13C,d15N) %>% 
  ggplot(aes(x=segment_percentile,y=x,groups=isotope),color=isotope) +
  facet_grid(age_n~id) + 
  geom_line(show.legend=FALSE)

dat <- dat %>% filter(!id%in%c("CU-2","CU-3","10","dec")) %>% droplevels()


###
### Write tidy data
###

write.csv(dat,"data/derived/dat.csv",row.names=FALSE)
