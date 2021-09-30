###
### Fit GAM with adaptive smooth
###

# Nest data
out <- all2 %>% filter(id%in%c("26","T3","T11","T7")) %>%
  select(d15N,d13C,long.short,id,age.n,midpoint.pctle) %>%
  gather(key="isotope",value="y",d15N,d13C) %>%
  group_by(id,long.short,age.n,isotope) %>%
  nest(data=c(y,midpoint.pctle)) %>%
  arrange(isotope,age.n)

# Fit GAMs
model.wrap <- function(data,age){
  if(age==0){
    mod <- gam(y~s(midpoint.pctle,bs="ad",k=10,m=3),gamma=1,data=data)
  }
  if(age>0){
    mod <- gam(y~s(midpoint.pctle,bs="ad"),gamma=1,data=data)
  }
  mod
}

out <- out %>% mutate(data=map(data,na.omit),  # remove NAs
                      gam=map2(data,age.n,~model.wrap(.x,.y)),  # fit GAMs
                      preds=map(gam,predict),  # get predictions
                      deriv=map(gam,cal.derivative.roots.mp))  # get derivatives


###
### Figure: Overlay of 0-, 2-, 3-, and 4-year-old
###

lmar <- 10  # left margin width
mmar <- 1  # middle margin width
rmar <- 2  # right margin width
tmar <- 8  # top margin width
bmar <- 13  # bottom margin width

pw <- (100-lmar-rmar-mmar*3)/3  # plot width
ph <- (100-tmar-bmar-mmar)/2  # plot height

nf <- layout(mat=matrix(c(0,1,0,2,0,3,0,4,0,
                          13,5,0,6,0,7,0,8,0,
                          0,0,0,0,0,0,0,0,0,
                          14,9,0,10,0,11,0,12,0,
                          0,15,0,15,0,15,0,15,0),5,9,byrow=TRUE),
             widths=c(lmar,pw,mmar,pw,mmar,pw,mmar,pw,rmar),
             heights=c(tmar,ph,mmar,ph,bmar))

layout.show(nf)
par(mar=c(0,0,0,0),xpd=NA)

plot(1,1,axes=FALSE,pch="",ylab="",xlab="")
text(1,1,"Pup",adj=c(0.5,0.5),cex=1.25)
plot(1,1,axes=FALSE,pch="",ylab="",xlab="")
text(1,1,"2-year-old",adj=c(0.5,0.5),cex=1.25)
plot(1,1,axes=FALSE,pch="",ylab="",xlab="")
text(1,1,"3-year-old",adj=c(0.5,0.5),cex=1.25)
plot(1,1,axes=FALSE,pch="",ylab="",xlab="")
text(1,1,"4-year-old",adj=c(0.5,0.5),cex=1.25)

for(i in c(1,3,5,7,9,11,13,15)){
  
  if(i<8) ylims <- range(do.call(rbind,out$data[1:8])$y)
  if(i>8) ylims <- range(do.call(rbind,out$data[9:16])$y)
  
  plot(out$data[[i]]$midpoint.pctle,out$preds[[i]],type="l",  # long whisker fit
       ylab="",xlab="",col=2,xlim=c(-0.05,1.05),ylim=ylims,
       xaxt=ifelse(i %in% c(1,3,5,7),'n','s'),yaxt="n")
  lines(out$data[[i+1]]$midpoint.pctle,out$preds[[i+1]],col=1)  # short whisker fit
  
  points(out$data[[i]]$midpoint.pctle,out$data[[i]]$y,  # long whisker data
         pch=19,col=rgb(1,0,0,0.25),cex=0.75)
  points(out$data[[i]]$midpoint.pctle,out$data[[i]]$y,
         pch=1,col=rgb(1,0,0,0.25),cex=0.75)
  
  points(out$data[[i+1]]$midpoint.pctle,out$data[[i+1]]$y,  # short whisker data
         pch=19,col=rgb(0,0,0,0.25),cex=0.75)
  points(out$data[[i+1]]$midpoint.pctle,out$data[[i+1]]$y,
         pch=1,col=rgb(0,0,0,0.25),cex=0.75)
  
  if(i==1){
    mt <- seq(round(ylims[1]),round(ylims[2]),0.5)
    axis(side=2,at=mt,labels=sprintf("%.1f",mt),las=1)
  }
  if(i==9){
    mt <- seq(round(ylims[1]),round(ylims[2]),2)
    axis(side=2,at=mt,labels=sprintf("%.1f",mt),las=1)
  }
}

plot(1,1,axes=FALSE,pch="",ylab="",xlab="")
text(1,1,expression(paste(delta^13,C)),adj=c(0.5,-1.15),srt=90,cex=1.25)
plot(1,1,axes=FALSE,pch="",ylab="",xlab="")
text(1,1,expression(paste(delta^13,N)),adj=c(0.5,-1.15),srt=90,cex=1.25)
plot(1,1,axes=FALSE,pch="",ylab="",xlab="")
text(1,1,"Percentile",adj=c(0.5,1.5),las=0,cex=1.25)

quartz.save("overlay.pdf",type="pdf")