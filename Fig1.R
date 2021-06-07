setwd("C:/Juv_vib")
dir()
all.summary <- read.csv("animal.summary2.csv", header = T)
head(all.summary)

#Animal length and vibrissae length comparison
summary(lm(vibrissae.length~animal.length,data =
             all.summary))

#Exclude DEC
all.summary <- subset(all.summary, !(id %in% c("dec")))
#drop unused levels in id
all.summary$id = droplevels(all.summary$id)

#make plots to 
pdf("Fig1.pdf")
par(mfrow = c(2,2))
boxplot(vibrissae.length~class, data = all.summary)
plot(vibrissae.length~age, data = all.summary)
plot(vibrissae.length~animal.length, data = all.summary)
plot(vibrissae.length~animal.mass, data = all.summary)
abline(-58,1.56)
dev.off()

