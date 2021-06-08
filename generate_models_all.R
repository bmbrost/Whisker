library(tidyverse)
source("helper.R")
source('data.setup.R')

#################################################
#model 1: constant 
#################################################
form.list.1 = list(
  N = d15N ~ s(midpoint.pctle),
  C = d13C ~ s(midpoint.pctle)
)
pup.models.1 = make.models(form.list=form.list.1)
age2.models.1 = make.models(whichclass = "juvenile", whichage.n = 2, form.list=form.list.1)
age3.models.1 = make.models(whichclass = "juvenile", whichage.n = 3, form.list=form.list.1)
age4.models.1 = make.models(whichclass = "juvenile", whichage.n = 4, form.list=form.list.1)

#model 2: individual random intercepts
all3$id <- factor(all3$id)
form.list.2 = list(
  N = d15N ~ s(id, bs='re') + s(midpoint.pctle),
  C = d13C ~ s(id, bs='re') + s(midpoint.pctle)
)
pup.models.2 = make.models(form.list=form.list.2)
age2.models.2 = make.models(whichclass = "juvenile", whichage.n = 2, form.list=form.list.2)
age3.models.2 = make.models(whichclass = "juvenile", whichage.n = 3, form.list=form.list.2)
age4.models.2 = make.models(whichclass = "juvenile", whichage.n = 4, form.list=form.list.2)

#model 3: individual random smooths
form.list.3 = list(
  N = d15N ~ s(midpoint.pctle) + s(midpoint.pctle, id, bs='fs'),
  C = d13C ~ s(midpoint.pctle) + s(midpoint.pctle, id, bs='fs')
)
pup.models.3 = make.models(form.list=form.list.3)
age2.models.3 = make.models(whichclass = "juvenile", whichage.n = 2, form.list=form.list.3)
age3.models.3 = make.models(whichclass = "juvenile", whichage.n = 3, form.list=form.list.3)
age4.models.3 = make.models(whichclass = "juvenile", whichage.n = 4, form.list=form.list.3)
