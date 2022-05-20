library(mgcv)
library(gratia)
library(rootSolve)
library(numDeriv)
library(tidyverse)

#Script for setting up the vibrissae data for gam modelling
# hopefully only run once

setwd("C:/Juv_vib")
allorig  = read.table("all.csv",header = T, sep ="," )
#make a copy of the original ; leave as is 
all = allorig

#define seg.midpoint and midpoint.pctle
all$seg.midpoint = all$seg.end.mm - .5*all$seg.len
all$midpoint.pctle = 1 - (all$seg.midpoint/all$whisk.len)
all$age.n = round(all$age)

#Set up all2: all variable from  all with id's not used removed
#and only a subset of variable
#exclude pup CU2,CU3, 10, and dec)

all2 <- subset(all, !(id %in% c("CU-2","CU-3","10","dec")))
#drop unused levels in id
all2 = droplevels(all2)

#Now trim all2 to the variable we need:

myvars=c("id","class","age","age.n","long.short","segment","d15N",
         "d13C","seg.length","seg.midpoint","whisk.len","midpoint.pctle")

all2 <- all2[myvars] 

# Now for all3, use only long whiskers; this can be changed easily or #all2 can be used if we want to model long and short alll at once.
all3 = all2[all2$long.short == "long",]

#check that we have done what we think 
names(all3); dim(all3)
unique(all3$id)
summary(all3)

#subset small vibrissae and calculate means for d13C and d15N
all4 = all2[all2$long.short == "short",]
d13C_mean <- aggregate(d13C ~ id, data = all4, mean)
d13C_mean
d15N_mean <- aggregate(d15N ~ id, data = all4, mean)
d15N_mean
