#fitlen is a function to estimate the length of the whisker
#for a given id at a given percentile (pctle)

fitlen = 
  function(id = "7", pctle = .5461) {
    set = all3[all3$id == id,]
    #set = all3[c("whisk.len","seg.midpoint")]
    foo = lm(I(whisk.len-seg.midpoint)~midpoint.pctle, data = set)
    predict(foo, newdata = data.frame(midpoint.pctle =  pctle))
  }

#+++++++++++++++++++++++++++++++++
# calculate derivates

cal.derivative.roots.mp = 
  function( model, id = "7"){
    #model is a gam model
    #new.df is a data frame with all variables model requires
    #    
    require(mgcv)
    require(numDeriv)
    require(rootSolve)
    #want x = percentiles and y prdicted values from same model in a 
    #newdata data frame which must include all the variables in the model
    #this one is for a model that only has percentile and a fixed id as the
    #derivative does not change with id
    f0=function(x){predict(model, newdata = data.frame (midpoint.pctle = x, id = id))}
    
    fprime = function(x){grad(f0, x)}
    #0 solve where derivative = 0
    uniroot.all(fprime, c(0,1))
  }

#+++++++++++++++++++++++++++++++++++++++++++


make.models = function(dataset = all3, whichclass = "pup", whichage.n= 0, form.list, whichlength.class = "long", ...){
  #+++++++++++++++++++++
  #this function generates carbon and nitrogen models, predictions, and #derivaatives  for given class and age
  #at present length class assumed to be "long"
  #+++++++++++++++++++++++++++++++++++++++++++
  require(mgcv)
  require(numDeriv)
  require(rootSolve)
  
  data.used = dataset[dataset$class == whichclass& 
                        dataset$age.n == whichage.n & whichlength.class == "long", ]
  
  # data.used = na.omit(data.used)
  
  mod15 = gam(form.list$N , data = data.used, method="REML")
  mod13 = gam(form.list$C, data = data.used, method="REML")
  
  fit.15 = predict(mod15, type = "terms", se = T , id = data.used$id[1], newdata=data.used, select=TRUE)
  fit.13 = predict(mod13, type = "terms" , se = T, id = data.used$id[1], newdata=data.used, select=TRUE)
  
  data.used$n.smooth = fit.15[[1]] [,'s(midpoint.pctle)']
  data.used$n.smooth.se = fit.15 [[2]] [,'s(midpoint.pctle)']
  data.used$c.smooth = fit.13[[1]] [,'s(midpoint.pctle)']
  data.used$c.smooth.se = fit.13 [[2]] [,'s(midpoint.pctle)'] 
  
  newdf = data.frame(midpoint.pctle = (0:100)/100)
  newdf$id = data.used$id[1]
  
  
  fit.15 = predict(mod15, newdata = newdf,type = "terms", se = T)
  fit.13 = predict(mod13, newdata = newdf, type = "terms" , se = T)
  
  newdf$n.smooth = fit.15[[1]] [,'s(midpoint.pctle)']
  newdf$n.smooth.se = fit.15 [[2]] [,'s(midpoint.pctle)']
  newdf$c.smooth = fit.13 [[1]] [,'s(midpoint.pctle)']
  newdf$c.smooth.se = fit.13[[2]][,'s(midpoint.pctle)']
  newdf = newdf[,-2]
  #names(newdf = c("midpoint.pctle","n.smooth", 
  #    "n.smooth.se","c.smooth". "c.smooth.se"))
  
  newdf$midpoint.mean.whisker.length = round(apply(
    sapply(unique(data.used$id), fitlen,  
           pctle = newdf$midpoint.pctle), 1,mean),2)
  
  mod15.derivatives = cal.derivative.roots.mp(mod15, id = data.used$id[1])
  
  mod13.derivatives = cal.derivative.roots.mp(mod13, id = data.used$id[1])
  
  axis.trans.model = lm(midpoint.mean.whisker.length~midpoint.pctle, data = newdf)
  
  out = list( data.used = data.used,  newdf= newdf , mod15 = mod15, 
              mod13 = mod13, mod15.derivatives = mod15.derivatives, 
              mod13.derivatives= mod13.derivatives, axis.trans.model = axis.trans.model)
}