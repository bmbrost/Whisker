library(lubridate)
library(tidyverse)
library(lattice)

setwd("C:/Juv_vib")
dir()
pups2 <- read.csv("pups2.csv", header = T)

# Create multipanel pup plots with reruns of vibrissae 9 and 26. 
unique(c(as.character(pups2$id)))
xyplot(d15N~whisk.pctle|id,data=pups2,type="b",pch=19,cex=0.5) 

xyplot(d15N~whisk.pctle|id,type="b",pch=19,cex=0.5,
       data=pups2 %>% filter(id!="CU-2"&id!="CU-3"&id!="10"),  # exclude pups
       panel=function(x,y,...){
         panel.abline(v=0.58,lty=3)  # add vertical line at 0.55
         panel.xyplot(x,y,...)
       })



#####
xyplot(dC/dN~total.length,data=tmp,group=pup,type="b",pch=19,cex=0.5,
       panel=function(x,y,...){
         panel.xyplot(x,y,...)
         # panel.loess(x,y,...,col=2,span=0.1)
       })

xyplot(dC/dN~total.length|pup,data=tmp,type="b",pch=19,cex=0.5,
       panel=function(x,y,...){
         panel.xyplot(x,y,...)
         # panel.loess(x,y,...,col=2,span=0.1)
       })

## plot dN
xyplot(dN~total.length|pup,data=tmp,type="b",pch=19,cex=0.5,
       panel=function(x,y,...){
         panel.xyplot(x,y,...)
         # panel.loess(x,y,...,col=2,span=0.1)
       })

