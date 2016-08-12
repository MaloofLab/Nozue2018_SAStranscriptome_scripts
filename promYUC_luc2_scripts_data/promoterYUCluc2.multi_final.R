############
##
## promoterYUCluc2.multi_final.R
##    finalized by Kazu
#########
setwd("/Volumes/Data6/data_JM4/promYUC_luc2/scripts_data/")
library(ggplot2);library(reshape2)
source("function.promYUCluc2.R") # run functions used in this scrip
####################################################
# sun + shade trace in one plat version (061716)
####################################################
# start from this line
all<-read.csv("all.csv") # 061716
all<-all[,-1]

## sample number (070714)
ftable(all[all$ZT<7,c("treatment","promoter")],row.vars="treatment",col.vars="promoter")
# promoter YUC5 YUC8 YUC9
# treatment                        
# shade                11   16   17
# sun                  14   11   16

# how to have different scale for different promoter? (092315)
all$TimeAfterTreatment<-all$ZT - 7 # time after treatment (hr)

# for ratio (revert log2(differences) to differeneces)
all$value2<-(2^(all$value) -1)*(10^2) # percent signal increased (071516)

## replace "sun" into "high \nR/FR", "shade into "low \nR/FR"
all$treatment<-as.character(all$treatment)
all$treatment<-gsub("sun","high \nR/FR",all$treatment)
all$treatment<-gsub("shade","low \nR/FR",all$treatment)
all$treatment<-factor(all$treatment,levels=c("high \nR/FR","low \nR/FR"))
# # for ratio (all$value2) for v4_ratio.png (061716)
# day<-data.frame(x1=rep(c(0,24,48,0,7,24,48),3),
#                 x2=rep(c(16,40,57,7,16,40,57),3),
#                 y1=rep(c(0.9997,0.9997-0.00004,0.9999,0.9999-0.000028,0.9999-0.00007,0.9999-0.00007-0.000036),c(3,4,3,4,3,4)),
#                 y2=rep(c(0.9997+0.000037,0.9997-0.00004+0.000037,0.9999+0.000025,0.9999-0.000003,0.9999-0.00007+0.000032,0.9999-0.00007-0.000036+0.000032),c(3,4,3,4,3,4)),
#                 treatment=rep(c("high \nR/FR","high \nR/FR","high \nR/FR","low \nR/FR","low \nR/FR","low \nR/FR","low \nR/FR"),3),
#                 promoter=rep(c("YUC5","YUC8","YUC9"),each=7),
#                 color=rep(c("high \nR/FR","high \nR/FR","high \nR/FR","high \nR/FR","low \nR/FR","low \nR/FR","low \nR/FR"),3))
# night<-data.frame(x1=rep(c(16,40),6),
#                   x2=rep(c(24,48),6),
#                   y1=rep(c(0.9997,0.9997-0.00004,0.9999,0.9999-0.000028,0.9999-0.00007,0.9999-0.00007-0.000036),each=2),
#                   y2=rep(c(0.9997+0.000037,0.9997-0.00004+0.000037,0.9999+0.000025,0.9999-0.000003,0.9999-0.00007+0.000032,0.9999-0.00007-0.000036+0.000032),each=2),
#                   treatment=rep(c("high \nR/FR","high \nR/FR","low \nR/FR","low \nR/FR"),3),
#                   promoter=rep(c("YUC5","YUC8","YUC9"),each=4)                  
#                   )
# for percent increase version (v5)
# for ratio (all$value2) for v4_ratio.png (061716)
day<-data.frame(x1=rep(c(0,24,48,0,7,24,48),3),
                x2=rep(c(16,40,57,7,16,40,57),3),
                y1=rep((c(0.9997,0.9997-0.00004,0.9999,0.9999-0.000028,0.9999-0.00007,0.9999-0.00007-0.000036)-1)*10^2,c(3,4,3,4,3,4)),
                y2=rep((c(0.9997+0.000037,0.9997-0.00004+0.000037,0.9999+0.000025,0.9999-0.000003,0.9999-0.00007+0.000032,0.9999-0.00007-0.000036+0.000032)-1)*10^2,c(3,4,3,4,3,4)),
                treatment=rep(c("high \nR/FR","high \nR/FR","high \nR/FR","low \nR/FR","low \nR/FR","low \nR/FR","low \nR/FR"),3),
                promoter=rep(c("YUC5","YUC8","YUC9"),each=7),
                color=rep(c("high \nR/FR","high \nR/FR","high \nR/FR","high \nR/FR","low \nR/FR","low \nR/FR","low \nR/FR"),3))
night<-data.frame(x1=rep(c(16,40),6),
                  x2=rep(c(24,48),6),
                  y1=rep((c(0.9997,0.9997-0.00004,0.9999,0.9999-0.000028,0.9999-0.00007,0.9999-0.00007-0.000036)-1)*10^2,each=2),
                  y2=rep((c(0.9997+0.000037,0.9997-0.00004+0.000037,0.9999+0.000025,0.9999-0.000003,0.9999-0.00007+0.000032,0.9999-0.00007-0.000036+0.000032)-1)*10^2,each=2),
                  treatment=rep(c("high \nR/FR","high \nR/FR","low \nR/FR","low \nR/FR"),3),
                  promoter=rep(c("YUC5","YUC8","YUC9"),each=4)                  
)

day[,1:2]<-day[,1:2]-7
night[,1:2]<-night[,1:2]-7
# 


## YUC5
YUC5<-ggplot(all[all$promoter=="YUC5",], aes(x=TimeAfterTreatment,y=value2,color=treatment))  + geom_smooth(stat="smooth",method="loess",span = 0.2,aes(fill=treatment)) #+ geom_point()
YUC5<-YUC5 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                   panel.background = element_blank(), axis.line = element_line(colour = "black"))
# drawing night period
# night period (for bottom bar)
YUC5 <- YUC5 + geom_rect(data=night[night$promoter=="YUC5",],mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2),fill="black", inherit.aes = FALSE)# draw dark/light cycle by white and gray rectangle    see http://stackoverflow.com/questions/4733182/how-to-highlight-time-ranges-on-a-plot
# drawing day period (high or low)
YUC5 <-YUC5 + geom_rect(data=day[day$promoter=="YUC5",],mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=color), inherit.aes = FALSE) 
# change bar color
YUC5 <- YUC5 + scale_fill_manual( values = c("high \nR/FR" = "red","low \nR/FR" = "darkred")) +
  scale_colour_manual( values = c("high \nR/FR" = "red","low \nR/FR" = "darkred"))

YUC5<-YUC5 + theme(legend.position = "none", axis.text.y=element_text(size=10)) + labs(x="Time after treatment (h)",y="Percent increase of signal to -1h")
# change range of y
#YUC5<-YUC5 + scale_y_continuous(limits=c(0.9995,1.0003))
# write promoter name
YUC5<-YUC5 + annotate("text",label="YUC5",x=-2,y=(1.00017-1)*10^2)
YUC5
## YUC8
YUC8<-ggplot(all[all$promoter=="YUC8",], aes(x=TimeAfterTreatment,y=value2,color=treatment))  + geom_smooth(stat="smooth",method="loess",span = 0.2,aes(fill=treatment)) #+ geom_point()
YUC8<-YUC8 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                   panel.background = element_blank(), axis.line = element_line(colour = "black"))
# drawing night period
# night period (for bottom bar)
YUC8 <- YUC8 + geom_rect(data=night[night$promoter=="YUC8",],mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2),fill="black", inherit.aes = FALSE)# draw dark/light cycle by white and gray rectangle    see http://stackoverflow.com/questions/4733182/how-to-highlight-time-ranges-on-a-plot
# drawing day period (high or low)
YUC8 <-YUC8 + geom_rect(data=day[day$promoter=="YUC8",],mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=color), inherit.aes = FALSE) 
# change bar color
YUC8 <- YUC8 + scale_fill_manual( values = c("high \nR/FR" = "red","low \nR/FR" = "darkred")) +
  scale_colour_manual( values = c("high \nR/FR" = "red","low \nR/FR" = "darkred"))

YUC8<-YUC8 + theme(legend.position = "none", axis.text.y=element_text(size=10)) + labs(x="Time after treatment (h)",y="Percent increase of signal to -1h")
# write promoter name
YUC8<-YUC8 + annotate("text",label="YUC8",x=-2,y=(1.00021-1)*10^2)
YUC8

## YUC9
YUC9<-ggplot(all[all$promoter=="YUC9",], aes(x=TimeAfterTreatment,y=value2,color=treatment))  + geom_smooth(stat="smooth",method="loess",span = 0.2,aes(fill=treatment)) #+ geom_point()
YUC9<-YUC9 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                   panel.background = element_blank(), axis.line = element_line(colour = "black"))
# drawing night period
# night period (for bottom bar)
YUC9 <- YUC9 + geom_rect(data=night[night$promoter=="YUC9",],mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2),fill="black", inherit.aes = FALSE)# draw dark/light cycle by white and gray rectangle    see http://stackoverflow.com/questions/4733182/how-to-highlight-time-ranges-on-a-plot
# drawing day period (high or low)
YUC9 <-YUC9 + geom_rect(data=day[day$promoter=="YUC9",],mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=color), inherit.aes = FALSE) 
# change bar color
YUC9 <- YUC9 + scale_fill_manual( values = c("high \nR/FR" = "red","low \nR/FR" = "darkred")) +
  scale_colour_manual( values = c("high \nR/FR" = "red","low \nR/FR" = "darkred"))

YUC9<-YUC9 + theme(axis.text.y=element_text(size=10),legend.text=element_text(size=10),legend.title=element_blank()) + labs(x="Time after treatment (h)",y="Percent increase of signal to -1h")
# write promoter name
YUC9<-YUC9 + annotate("text",label="YUC9",x=-2,y=(1.00023-1)*10^2)
YUC9
# cowplot
library(cowplot)
plot.YUC589<-plot_grid(YUC5,YUC8,YUC9,ncol=3,labels=c("A","B","C"),
                       scale=0.9,vjust=0,rel_widths=c(1,1,1.35))
# save_plot("Fig1_promYUC589_LUC2v4_ratio.png", plot.YUC589,
#           ncol = 3, # we're saving a grid plot of 2 columns
#           nrow = 1, # and 2 rows
#           # each individual subplot should have an aspect ratio of 1.3
#           base_aspect_ratio = 0.8)
save_plot("Fig1_promYUC589_LUC2v5_ratio.png", plot.YUC589,
          ncol = 3, # we're saving a grid plot of 2 columns
          nrow = 1, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 0.8)

