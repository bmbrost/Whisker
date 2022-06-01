###
###
### Setup vibrissae isotope data for modeling
###
###

rm(list=ls())
setwd("/Users/brian.brost/Documents/sandbox/zeppelin/vibrissae/")
# setwd("C:/Juv_vib")


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
# all <- dat

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





# #define seg.midpoint and midpoint.pctle
# head(all)
# all$seg.midpoint = all$seg.end.mm - .5*all$seg.length
# 
# all$midpoint.pctle = 1 - (all$seg.midpoint/all$whisk.len)
# 
# all$age.n = round(all$age)
# 
# #Set up all2: all variable from  all with id's not used removed
# #and only a subset of variable
# #exclude pup CU2,CU3, 10, and dec)
# 
# all2 <- subset(all, !(id %in% c("CU-2","CU-3","10","dec")))
# #drop unused levels in id
# all2 = droplevels(all2)
# 
# #Now trim all2 to the variable we need:
# 
# myvars=c("id","class","age","age.n","long.short","segment","d15N",
#          "d13C","seg.length","seg.midpoint","whisk.len","midpoint.pctle")
# 
# all2 <- all2[myvars] 
# 
# # Now for all3, use only long whiskers; this can be changed easily or #all2 can be used if we want to model long and short alll at once.
# all3 = all2[all2$long.short == "long",]
# 
# #check that we have done what we think 
# names(all3); dim(all3)
# unique(all3$id)
# summary(all3)
# 
# #subset small vibrissae and calculate means for d13C and d15N
# all4 = all2[all2$long.short == "short",]
# d13C_mean <- aggregate(d13C ~ id, data = dat %>% filter(length_class=="short"), mean)
# d13C_mean
# d15N_mean <- aggregate(d15N ~ id, data = all4, mean)
# d15N_mean
# 
# dat %>% gather(key="isotope",value="x",d13C,d15N) %>% group_by(id,length_class,isotope) %>%
#   summarize(mean=mean(x)) %>% arrange(isotope,id,length_class) %>% print(n=nrow(.))
