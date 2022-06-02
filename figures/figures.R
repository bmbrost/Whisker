###
###
### Create figures for publication
###
###


rm(list=ls())
setwd("/Users/brian.brost/Documents/sandbox/zeppelin/vibrissae/")


#####################################################################################
### Figure 1
#####################################################################################

# Add code here



#####################################################################################
### Long vibrissae
#####################################################################################

rm(list=ls())
load("/Users/brian.brost/Documents/sandbox/zeppelin/vibrissae/analysis/long_analysis.Rdata")
ls()

# Update factor labels for panel strips
preds <- preds %>% mutate(xage_n=factor(age_n,levels=c("0","2","3","4"),
                                       labels=c("Age 0","Age 2","Age 3","Age 4")),
                         xisotope=factor(isotope,levels=c("d13C","d15N"),
                                        labels=c(expression(paste(delta^13,"C",sep="")),expression(paste(delta^15,"N",sep="")))))

                                        
###
### Figure for main body of manuscript
###

pdf("figures/isotopes_by_age.pdf",width=10,height=7)
preds %>% unnest(id_preds) %>% 
  ggplot(aes(x=segment_percentile,y=fit,group=id)) +
  facet_grid(xisotope~xage_n,scales="free_y",labeller=labeller(.rows=label_parsed),switch="y") +
  geom_line(show.legend = FALSE,color="gray10",size=0.1) +
  geom_vline(data=preds %>% unnest(pop_zeros) %>% mutate(id=NA),aes(xintercept=zeros),color="red",lty=3,size=0.5) +
  geom_line(data=preds %>% unnest(pop_preds) %>% mutate(id=NA),aes(x=segment_percentile,y=fit),color="red",size=0.5) +
  geom_ribbon(data=preds %>% unnest(pop_preds) %>% mutate(id=NA),
              aes(ymin=fit-2*se.fit,ymax=fit+2*se.fit,x=segment_percentile),fill="red",alpha=0.15) +
  labs(x="Segment percentile",y=NULL) +
  theme(panel.background = element_rect(fill="white",colour=NA),  # recreating theme_bw()
        panel.border=element_rect(fill=NA,colour="grey20"),
        panel.grid=element_line(colour="grey92"), 
        panel.grid.minor=element_line(size=rel(0.5)), 
        legend.key=element_rect(fill="white",colour=NA),complete=TRUE,
        strip.background.x=element_rect(fill="grey85",colour="grey20"),  # x-axis strips only
        strip.background.y=element_blank(),strip.placement="outside",  # updates to theme_bw()...
        axis.text.x.bottom=element_text(size=rel(1),angle=0),
        # axis.title.x.bottom=element_text(vjust=-1),
        panel.spacing.x=unit(0.15,units="inches"))
dev.off()



###
### Pup figures for Appendix
###

preds

# Pups: d13C figure
i <- 1  # age_n==0 and isotope=="d13C"
my.ylab <- expression(paste(delta^13,"C",sep=""))

# Pups: d15N figure
i <- 5  # age_n==0 and isotope=="d15N"
my.ylab <- expression(paste(delta^15,"N",sep=""))


p <- preds$id_preds[[i]] %>% ggplot(aes(x=segment_percentile,y=fit)) +
  facet_wrap(~id) +
  geom_vline(data=preds$id_zeros[[i]],aes(xintercept=zeros),color="black",lty=3,size=0.5) +
  geom_vline(data=preds$pop_zeros[[i]],aes(xintercept=zeros),color="red",lty=3,size=0.5) +
  geom_point(data=preds$data[[i]],aes(x=segment_percentile,y=y),color=rgb(0,0,0,0.35),pch=19,size=0.75) +
  geom_point(data=preds$data[[i]],aes(x=segment_percentile,y=y),color="gray50",pch=1,size=0.75) +
  geom_line(show.legend = FALSE,color="gray10",size=0.1) +
  geom_ribbon(aes(ymin=fit-2*se.fit,ymax=fit+2*se.fit,x=segment_percentile),fill="black",alpha=0.15)+
  geom_line(data=preds$pop_preds[[i]],aes(x=segment_percentile,y=fit),color="red",size=0.5) +
  geom_ribbon(data=preds$pop_preds[[i]],
              aes(ymin=fit-2*se.fit,ymax=fit+2*se.fit,x=segment_percentile),fill="red",alpha=0.15)+
  labs(x="Segment percentile",y=my.ylab) +
  xlim(0,1) +
  theme(panel.background = element_rect(fill="white",colour=NA),  # recreating theme_bw()
        panel.border=element_rect(fill=NA,colour="grey20"),
        panel.grid=element_line(colour="grey92"), 
        panel.grid.minor=element_line(size=rel(0.5)), 
        legend.key=element_rect(fill="white",colour=NA),complete=TRUE,
        strip.background.x=element_rect(fill="grey85",colour="grey20"),  # x-axis strips only
        strip.background.y=element_blank(),strip.placement="outside",  # updates to theme_bw()...
        axis.text.x.bottom=element_text(size=rel(1),angle=0),
        # axis.title.x.bottom=element_text(vjust=-1),
        panel.spacing.x=unit(0.15,units="inches"))

# Save d13C figure
pdf("figures/appendix_d13C_age_0.pdf",width=10,height=7)
p 
dev.off()

# Save d15N figure
pdf("figures/appendix_d15N_age_0.pdf",width=10,height=7)
p
dev.off()

rm(p)  # object size is big


###
### Non-pup figures for Appendix
###

# Non-pups: d13C figures
idx <- 2:4  # age_n==2,3,4 and isotope=="d13C"
my.ylab <- expression(paste(delta^13,"C",sep=""))
my.breaks <- -20:-15

# Non-pups: d15N figures
idx <- 6:8  # age_n==2,3,4 and isotope=="d15N"
my.ylab <- expression(paste(delta^15,"N",sep=""))
my.breaks <- seq(10,22,by=2)

for(i in 1:3){  # loop through ages
  j <- idx[i]  # row in preds
  
  if(preds$age_n[j]!=4){
    tmp <- preds[j,] %>% unnest(id_preds) %>% droplevels()  # ages 2 and 3
  }
  if(preds$age_n[j]==4){  # add extra panel to age 4 figure so it lines up with figures for ages 2 and 3
    tmp <- preds[j,] %>% unnest(id_preds) %>% droplevels() %>% mutate(id=factor(id,levels=c("A",levels(id))))  # age 4
  }

  p <- tmp %>%  ggplot(aes(x=segment_percentile,y=fit)) + 
    facet_grid(xage_n~id,drop=FALSE) +
    geom_vline(data=preds$id_zeros[[j]],aes(xintercept=zeros),color="black",lty=3,size=0.5) +
    geom_vline(data=preds$pop_zeros[[j]],aes(xintercept=zeros),color="red",lty=3,size=0.5) +
    geom_point(data=preds$data[[j]],aes(x=segment_percentile,y=y),color=rgb(0,0,0,0.35),pch=19,size=0.75) +
    geom_point(data=preds$data[[j]],aes(x=segment_percentile,y=y),color="gray50",pch=1,size=0.75) +
    geom_line(show.legend = FALSE,color="gray10",size=0.1) +
    geom_ribbon(aes(ymin=fit-2*se.fit,ymax=fit+2*se.fit,x=segment_percentile),fill="black",alpha=0.15)+
    geom_line(data=preds$pop_preds[[j]],aes(x=segment_percentile,y=fit),color="red",size=0.5) +
    geom_ribbon(data=preds$pop_preds[[j]],
                aes(ymin=fit-2*se.fit,ymax=fit+2*se.fit,x=segment_percentile),fill="red",alpha=0.15)+
    scale_y_continuous(breaks=my.breaks) +
    labs(x="Segment percentile",y=my.ylab) +
    xlim(0,1) +
    theme_bw()
  
  assign(paste0("p",i),p)
}

# Remove left panel in age 4 plots
g <- ggplotGrob(p3)
rm_grobs <- g$layout$name %in% c("panel-1-1", "strip-t-1","axis-t-1","axis-b-1")  # identify grobs that should be removed
g$grobs[rm_grobs] <- NULL  # remove grobs
g$layout <- g$layout[!rm_grobs, ]  # remove grobs
g$layout[g$layout$name=="ylab-l","l"] <- 7  # move y-axis label closer to panel
g$layout[g$layout$name=="axis-l-1","l"] <- 6  # move y-axis tick marks closer to panel
grid.newpage()  # check new plot
grid.draw(g)

# Save d13C figure
pdf("figures/appendix_d13C_age_2_3_4.pdf",width=10,height=7)
grid.arrange(p1,p2,g,nrow=3)
dev.off()

# Save d15N figure
pdf("figures/appendix_d15N_age_2_3_4.pdf",width=10,height=7)
grid.arrange(p1,p2,g,nrow=3)
dev.off()

rm(p,p1,p2,p3)  # object sizes are big




#####################################################################################
### Paired vibrissae
#####################################################################################

rm(list=ls())
# load("/Users/brian.brost/Documents/sandbox/zeppelin/vibrissae/analysis/paired_analysis.Rdata")
ls()

# Update factor labels for panel strips
preds <- preds %>% mutate(xage_n=factor(age_n,levels=c("0","2","3","4"),
                                       labels=c("Age 0","Age 2","Age 3","Age 4")),
    xisotope=factor(isotope,levels=c("d13C","d15N"),
      labels=c(expression(paste(delta^13,"C",sep="")),expression(paste(delta^15,"N",sep="")))))


###
### Figure for main body of manuscript
###

scaleFUN <- function(x) sprintf("%.1f", x)

# With zeros: red = long, blue = short
pdf("figures/paired_whiskers_zeros.pdf",width=10,height=7)
preds %>% unnest(preds) %>% 
  ggplot(aes(x=segment_percentile,y=fit,groups=length_class,color=length_class)) +
  facet_grid(xisotope~xage_n,scales="free_y",labeller=labeller(.rows=label_parsed),switch="y") +
  geom_vline(data=preds %>% unnest(zeros),aes(xintercept=zeros,color=length_class),lty=3,size=0.25,show.legend=FALSE) +
  geom_ribbon(data=preds %>% unnest(preds),aes(ymin=fit-2*se.fit,ymax=fit+2*se.fit,x=segment_percentile,
                                               fill=length_class,linetype=NA),alpha=0.25,show.legend=FALSE) +
  geom_line(show.legend=FALSE,size=0.5) +
  labs(x="Segment percentile",y=NULL) +
  theme(panel.background = element_rect(fill="white",colour=NA),  # recreating theme_bw()
        panel.border=element_rect(fill=NA,colour="grey20"),
        panel.grid=element_line(colour="grey92"), 
        panel.grid.minor=element_line(size=rel(0.5)), 
        legend.key=element_rect(fill="white",colour=NA),complete=TRUE,
        strip.background.x=element_rect(fill="grey85",colour="grey20"),  # x-axis strips only
        strip.background.y=element_blank(),strip.placement="outside",  # updates to theme_bw()...
        axis.text.x.bottom=element_text(size=rel(1),angle=0),
        panel.spacing.x=unit(0.15,units="inches")) +
  scale_y_continuous(labels=scaleFUN)
dev.off()

# Without zeros: red = long, blue = short
pdf("figures/paired_whiskers.pdf",width=10,height=7)
preds %>% unnest(preds) %>% 
  ggplot(aes(x=segment_percentile,y=fit,groups=length_class,color=length_class)) +
  facet_grid(xisotope~xage_n,scales="free_y",labeller=labeller(.rows=label_parsed),switch="y") +
  geom_ribbon(data=preds %>% unnest(preds),aes(ymin=fit-2*se.fit,ymax=fit+2*se.fit,x=segment_percentile,
                                               fill=length_class,linetype=NA),alpha=0.25,show.legend=FALSE) +
  geom_line(show.legend=FALSE,size=0.5) +
  labs(x="Segment percentile",y=NULL) +
  theme(panel.background = element_rect(fill="white",colour=NA),  # recreating theme_bw()
        panel.border=element_rect(fill=NA,colour="grey20"),
        panel.grid=element_line(colour="grey92"), 
        panel.grid.minor=element_line(size=rel(0.5)), 
        legend.key=element_rect(fill="white",colour=NA),complete=TRUE,
        strip.background.x=element_rect(fill="grey85",colour="grey20"),  # x-axis strips only
        strip.background.y=element_blank(),strip.placement="outside",  # updates to theme_bw()...
        axis.text.x.bottom=element_text(size=rel(1),angle=0),
        panel.spacing.x=unit(0.15,units="inches")) +
  scale_y_continuous(labels=scaleFUN)
dev.off()




