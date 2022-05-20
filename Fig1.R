#Code for comparing animal length, mass and age to whisker length.

setwd("C:/Juv_vib")
dir()
all.summary <- read.csv("animal.summary.csv", header = T)
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

#t-test to see if whisker lengths are longer for juveniles than pups
t.test(vibrissae.length ~ class, data = all.summary)


# fit linear regressions and reportresults on the graphs (y=mx+b),p_value,adjusted R2,etc.
# for Animal length and vibrissae length comparison
animalLength.mod <- lm(animal.length ~ vibrissae.length, data=all.summary)
print(animalLength.mod)
summary(animalLength.mod)
# for Animal age and vibrissae length comparison
animalage.mod <- lm(age ~ vibrissae.length, data=all.summary)
print(animalage.mod)
summary(animalage.mod)
# for Animal mass and vibrissae length comparison
animalmass.mod <- lm(animal.mass ~ vibrissae.length, data=all.summary)
print(animalmass.mod)
summary(animalmass.mod)
