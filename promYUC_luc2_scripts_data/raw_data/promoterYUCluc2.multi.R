############
##
## promoterYUCluc2.multi.R
##
#########
# get all file information (ls -lt > exp4files.txt)
# keep original files (compress them)
# separate set1 and set2 (if necessary, rotate/transform images to correct wrong plate posiion/orientaiton) and delete blank images
#
# measurment note; radias is "50".
# set scale
# have an extra image (usually the last image) for labeling samples
# record plant name (promoter, line number, sun/shade)
# have output with two ways ("measure" (save as luc2.data) and "multi" (save as luc2.data.multi))
# for exp3, 4, 5, right LED chamber was used as "sun" and left chamber was used as "shade".
# for exp6 I switched chamber condition. (see spectram)
# 
# (051916) Finding optimal "sapn" in geom_smooth() to visualize peak in shade. (J's suggestion)
# (061816) Plotting both sun and shade in one graph.
# (071516) reorganize folders and make repository.
setwd("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/")
library(ggplot2);library(reshape2)
source("function.promYUCluc2.R") # run functions used in this scrip
setwd("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/")

# exp3
exp3files<-read.table("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp3/exp3files.txt",skip=2)

first_day_ZT0<-"Jun/11/06:00"

## exp3 set1
plant.name<-c(
paste("shade_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1"),sep=""),
paste("shade_plate2_YUC8_line2_rep1",sep=""),
paste("shade_plate3_",c("YUC9_line3_rep1","YUC5_line3_rep1","YUC8_line3_rep2"),sep=""),
paste("sun_plate1_",c("YUC9_line2_rep1","YUC5_line2_rep1","YUC8_line2_rep1"),sep=""),
paste("sun_plate2_",c("YUC9_line1_rep1","YUC5_line1_rep1"),sep="")
)
luc2.data<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp3/exp3set1data.csv") 
luc2.data.multi<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp3/exp3set1data.multi.csv") 
exp3set1.luc2.data.multi.log2<-format.luc2data.multi.exp3set1special.log2(plant.name,luc2.data,luc2.data.multi,files=exp3files,first_day_ZT0=first_day_ZT0) # under construction (070214)
#exp3set1.luc2.data.multi.log2<-temp.multi.dif.melt

head(exp3set1.luc2.data.multi.log2)
exp3set1.luc2.data.multi.log2$value
exp3set1.luc2.data.multi.log2<-exp3set1.luc2.data.multi.log2[!exp3set1.luc2.data.multi.log2$value=="-Inf",]
exp3set1.luc2.data.multi.log2$variable<-paste("exp3_set1_",exp3set1.luc2.data.multi.log2$variable,sep="")

p<-ggplot(exp3set1.luc2.data.multi.log2, aes(x=ZT,y=value,color=line)) + geom_point()
p<-p + facet_wrap(treatment~promoter)
p # 



## exp3 set2
plant.name<-c(
  paste("shade_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1"),sep=""),
  paste("shade_plate2_",c("YUC9_line2_rep1","YUC5_line2_rep1","YUC8_line2_rep1"),sep=""),
  paste("shade_plate3_",c("YUC9_line3_rep1","YUC5_line3_rep1"),sep=""),
  paste("sun_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1"),sep=""),
  paste("sun_plate2_",c("YUC9_line2_rep1","YUC8_line2_rep1"),sep=""),
  paste("sun_plate3_",c("YUC9_line3_rep1","YUC5_line3_rep1"),sep="")
  )
luc2.data<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp3/exp3set2data.csv") 
luc2.data.multi<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp3/exp3set2data.multi.csv") # need to scale!
#exp3set2.luc2.data.multi<-format.luc2data.multi.exp3set2special(plant.name,luc2.data,luc2.data.multi,files=exp3files,first_day_ZT0=first_day_ZT0) # under construction (070214)
exp3set2.luc2.data.multi.log2<-format.luc2data.multi.exp3set2special.log2(plant.name,luc2.data,luc2.data.multi,files=exp3files,first_day_ZT0=first_day_ZT0) # under construction (070314)
exp3set2.luc2.data.multi.log2$variable<-paste("exp3_set2_",exp3set2.luc2.data.multi.log2$variable,sep="")

head(exp3set2.luc2.data.multi.log2)
# remove fake "0"
exp3set2.luc2.data.multi.log2<-exp3set2.luc2.data.multi.log2[!(exp3set2.luc2.data.multi.log2$value==0&exp3set2.luc2.data.multi.log2$ZT>7),]
p<-ggplot(exp3set2.luc2.data.multi.log2, aes(x=ZT,y=value,color=line)) + geom_point()
p<-p + facet_wrap(treatment~promoter)
p # 
# 
# exp4
# exp4 set1 (add "rep" info to avoid same column name)
plant.name<-c(
paste("shade_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1"),sep=""),
paste("shade_plate2_YUC8_line2_rep1",sep=""),
paste("shade_plate3_YUC5_line3_rep1",sep=""),
paste("sun_plate1_",c("YUC9_line1_rep1","YUC5_line1_rep1"),sep=""),
paste("sun_plate3_",c("YUC9_line2_rep1","YUC8_line2_rep1","YUC5_line2_rep1"),sep=""),
paste("sun_plate4_YUC5_line3_rep1",sep="")
)
luc2.data<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp4/exp4set1data.csv") 
luc2.data.multi<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp4/exp4set1data.multi.csv") # need to scale!
#luc2.data.multi<-luc2.data.multi[grep("lab",luc2.data.multi$Label,invert=TRUE),]
exp4files<-read.table("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp4/exp4files.txt",skip=2)
first_day_ZT0<-"Jun/18/06:00"
#exp4set1.luc2.data.multi<-format.luc2data.multi(plant.name,luc2.data,luc2.data.multi,files=exp4files,first_day_ZT0=first_day_ZT0) # under construction (070214)
exp4set1.luc2.data.multi.log2<-format.luc2data.multi.log2(plant.name,luc2.data,luc2.data.multi,files=exp4files,first_day_ZT0=first_day_ZT0) # under construction (070214)
exp4set1.luc2.data.multi.log2$variable<-paste("exp4_set1_",exp4set1.luc2.data.multi.log2$variable,sep="")

head(exp4set1.luc2.data.multi.log2)
library(ggplot2)
p<-ggplot(exp4set1.luc2.data.multi.log2, aes(x=ZT,y=value,color=line)) + geom_point()
p<-p + facet_wrap(treatment~promoter)
p # 


# exp4 set2
plant.name<-c(
paste("shade_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1","YUC5_line1_rep1"),sep=""),
paste("shade_plate2_",c("YUC9_line2_rep1","YUC8_line2_rep1"),sep=""),
paste("shade_plate3_",c("YUC9_line3_rep1","YUC8_line2_rep1","YUC5_line3_rep1"),sep=""),
paste("sun_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1","YUC5_line1_rep1"),sep=""),
paste("sun_plate2_",c("YUC9_line2_rep1","YUC8_line2_rep1","YUC5_line2_rep1"),sep=""),
"sun_plate3_YUC9_line3_rep1"
)
luc2.data<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp4/exp4set2data.csv") 
luc2.data.multi<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp4/exp4set2data.multi.csv") # need to scale!
#luc2.data.multi<-luc2.data.multi[grep("lab",luc2.data.multi$Label,invert=TRUE),]
exp4files<-read.table("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp4/exp4files.txt",skip=2)
first_day_ZT0<-"Jun/18/06:00"
#exp4set2.luc2.data.multi<-format.luc2data.multi(plant.name,luc2.data,luc2.data.multi,files=exp4files,first_day_ZT0=first_day_ZT0) # under construction (070214)
exp4set2.luc2.data.multi.log2<-format.luc2data.multi.log2(plant.name,luc2.data,luc2.data.multi,files=exp4files,first_day_ZT0=first_day_ZT0) # under construction (070214)
exp4set2.luc2.data.multi.log2$variable<-paste("exp4_set2_",exp4set2.luc2.data.multi.log2$variable,sep="")

head(exp4set2.luc2.data.multi.log2)
library(ggplot2)
p<-ggplot(exp4set2.luc2.data.multi.log2, aes(x=ZT,y=value,color=line)) + geom_point()
p<-p + facet_wrap(treatment~promoter)
p # 
# summary
p<-ggplot(rbind(exp4set1.luc2.data.multi.log2,exp4set2.luc2.data.multi.log2), aes(x=ZT,y=value,color=line)) + geom_point()
p<-p + facet_wrap(treatment~promoter)
p # 


# exp5
# ls -lt > exp5set2files.txt # done in set2
exp5files<-read.table("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp5/exp5files.txt", skip=2) # exp5label.tif is for recording lables of each plants
#files<-files[!files$V9=="exp5label.tif",]  # remove "exp5label.tif"
first_day_ZT0<-"Jun/25/06:00"

# exp5 set1
plant.name<-c(
paste("shade_plate1_",c("YUC9_line1_rep1","YUC9_line1_rep2","YUC8_line1_rep1","YUC8_line1_rep1","YUC5_line1_rep1"),sep=""),
paste("shade_plate2_",c("YUC9_line3_rep1","YUC9_line3_rep2","YUC5_line3_rep1"),sep=""),
paste("sun_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1","YUC8_line1_rep2"),sep=""),
paste("sun_plate3_","YUC5_line3_rep1",sep="")
)
luc2.data<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp5/exp5set1data.csv") 
luc2.data.multi<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp5/exp5set1data.multi.csv") # need to scale!
#luc2.data.multi<-luc2.data.multi[grep("lab",luc2.data.multi$Label,invert=TRUE),]
#exp5set1.luc2.data.multi<-format.luc2data.multi(plant.name,luc2.data,luc2.data.multi,files=exp5files,first_day_ZT0=first_day_ZT0) # under construction (070214)
# log2 ratio
exp5set1.luc2.data.multi.log2<-format.luc2data.multi.log2(plant.name,luc2.data,luc2.data.multi,files=exp5files,first_day_ZT0=first_day_ZT0) # under construction (070214)
exp5set1.luc2.data.multi.log2<-exp5set1.luc2.data.multi.log2[exp5set1.luc2.data.multi.log2$ZT< 48,]
exp5set1.luc2.data.multi.log2$variable<-paste("exp5_set1_",exp5set1.luc2.data.multi.log2$variable,sep="")
head(exp5set1.luc2.data.multi.log2)
p<-ggplot(exp5set1.luc2.data.multi.log2, aes(x=ZT,y=value,color=line)) + geom_point()
p<-p + facet_wrap(treatment~promoter)
p # 
# exp5 set2
# enter experimental information
# input genotype and rep by looking "exp5lable.tif" by mannually. Transgenic "line" has to be consistent throughout experiments.
plant.name<-c(
paste("shade_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1","YUC5_line1_rep1","YUC8_line1_rep2","YUC5_line1_rep2"),sep=""),
paste("shade_plate2_",c("YUC9_line2_rep1","YUC9_line2_rep2","YUC5_line2_rep1","YUC5_line2_rep2"),sep=""),
paste("shade_plate3_",c("YUC9_line3_rep1","YUC9_line3_rep2","YUC5_line3_rep1"),sep=""),
paste("sun_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1","YUC5_line1_rep1","YUC5_line1_rep2"),sep=""),
paste("sun_plate2_",c("YUC8_line2_rep1","YUC5_line2_rep1","YUC5_line2_rep2"),sep=""),
paste("sun_plate3_",c("YUC9_line3_rep1","YUC5_line3_rep1","YUC5_line3_rep2"),sep="")) 
#
luc2.data<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp5/exp5set2data.csv") 
luc2.data.multi<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp5/exp5set2data.multi.csv") # need to scale!

#luc2.data.multi<-luc2.data.multi[grep("lab",luc2.data.multi$Label,invert=TRUE),]
exp5set2.luc2.data.multi<-format.luc2data.multi(plant.name,luc2.data,luc2.data.multi,files=exp5files,first_day_ZT0=first_day_ZT0) # under construction (070214)
exp5set2.luc2.data.multi.log2<-format.luc2data.multi.log2(plant.name,luc2.data,luc2.data.multi,files=exp5files,first_day_ZT0=first_day_ZT0) # under construction (070214)

head(exp5set2.luc2.data.multi.log2)

## eliminate day3 (I did not turn on light ...)
exp5set2.luc2.data.multi.log2<-exp5set2.luc2.data.multi.log2[exp5set2.luc2.data.multi.log2$ZT< 48,]
exp5set2.luc2.data.multi.log2$variable<-paste("exp5_set2_",exp5set2.luc2.data.multi.log2$variable,sep="")

p<-ggplot(exp5set2.luc2.data.multi.log2, aes(x=ZT,y=value,color=line)) + geom_point() 
p<-p + facet_wrap(treatment~promoter)
p # 

# exp6 
exp6files<-read.table("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp6/exp6files.txt", skip=2) # exp5label.tif is for recording lables of each plants
first_day_ZT0<-"July/4/06:00"

## exp6 set1
plant.name<-c(
paste("shade_plate1_",c("YUC5_line1_rep1","YUC5_line1_rep2"),sep=""),
paste("shade_plate2_",c("YUC8_line2_rep1","YUC8_line2_rep2","YUC5_line2_rep1","YUC5_line2_rep2"),sep=""),
paste("shade_plate3_",c("YUC9_line3_rep1","YUC5_line3_rep1","YUC5_line3_rep2"),sep=""),
paste("sun_plate1_",c("YUC9_line1_rep1","YUC5_line1_rep1"),sep=""),
paste("sun_plate2_",c("YUC9_line2_rep1","YUC8_line2_rep1","YUC8_line2_rep2","YUC5_line2_rep1"),sep=""),
paste("sun_plate3_",c("YUC9_line3_rep1","YUC9_line3_rep2","YUC8_line3_rep1","YUC5_line3_rep1","YUC5_line3_rep2"),sep="")
)
luc2.data<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp6/exp6set1data.csv") 
luc2.data.multi<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp6/exp6set1data.multi.csv") 
exp6set1.luc2.data.multi.log2<-format.luc2data.multi.log2(plant.name,luc2.data,luc2.data.multi,files=exp6files,first_day_ZT0=first_day_ZT0)
exp6set1.luc2.data.multi.log2$variable<-paste("exp6_set1_",exp6set1.luc2.data.multi.log2$variable,sep="")
p<-ggplot(exp6set1.luc2.data.multi.log2, aes(x=ZT,y=value,color=line)) + geom_point() 
p<-p + facet_wrap(treatment~promoter)
p # 
## exp6 set2
plant.name<-c(
paste("shade_plate1_",c("YUC5_line1_rep1","YUC5_line1_rep2"),sep=""),
paste("shade_plate2_",c("YUC9_line2_rep1","YUC9_line2_rep2","YUC8_line2_rep1","YUC5_line2_rep1","YUC5_line2_rep2"),sep=""),
paste("shade_plate3_",c("YUC9_line3_rep1","YUC9_line3_rep2","YUC5_line3_rep1","YUC5_line3_rep2"),sep=""),
paste("sun_plate1_",c("YUC5_line1_rep1","YUC5_line1_rep2"),sep=""),
paste("sun_plate2_",c("YUC5_line2_rep1","YUC5_line2_rep2"),sep=""),
paste("sun_plate3_",c("YUC9_line3_rep1","YUC9_line3_rep2","YUC5_line3_rep1","YUC5_line3_rep2"),sep="")
)
luc2.data<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp6/exp6set2data.csv") 
luc2.data.multi<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp6/exp6set2data.multi.csv") 
exp6set2.luc2.data.multi.log2<-format.luc2data.multi.log2(plant.name,luc2.data,luc2.data.multi,files=exp6files,first_day_ZT0=first_day_ZT0)
exp6set2.luc2.data.multi.log2$variable<-paste("exp6_set2_",exp6set2.luc2.data.multi.log2$variable,sep="")
p<-ggplot(exp6set2.luc2.data.multi.log2, aes(x=ZT,y=value,color=line)) + geom_point() 
p<-p + facet_wrap(treatment~promoter)
p # 

# summary (log2 fold change version)
all<-rbind(exp3set1.luc2.data.multi.log2,exp3set2.luc2.data.multi.log2,exp4set1.luc2.data.multi.log2,exp4set2.luc2.data.multi.log2,exp5set1.luc2.data.multi.log2,exp5set2.luc2.data.multi.log2,exp6set1.luc2.data.multi.log2,exp6set2.luc2.data.multi.log2)
#all.YUC8and9<-all[!all$promoter=="YUC5",]
all$variable
day<-data.frame(x1=rep(c(0,24,48,0,7,24,48),2),x2=rep(c(16,40,57,7,16,40,57),2),y1=rep(min(all$value),14),y2=rep(max(all$value),14),treatment=rep(c("sun","sun","sun","shade","shade","shade","shade"),2),promoter=rep(c("YUC8","YUC9"),each=7),color=rep(c("high","high","high","high","low","low","low"),2))
night<-data.frame(x1=c(16,40),x2=c(24,48),y1=min(all$value,all$value),y2=max(all$value,all$value))

# for scale_y_continuous(limits=c(-0.001,0.0013))
day<-data.frame(x1=rep(c(0,24,48,0,7,24,48),2),x2=rep(c(16,40,57,7,16,40,57),2),y1=rep(-0.001,14),y2=rep(-0.0009,14),treatment=rep(c("sun","sun","sun","shade","shade","shade","shade"),2),promoter=rep(c("YUC8","YUC9"),each=7),color=rep(c("high","high","high","high","low","low","low"),2))
night<-data.frame(x1=c(16,40),x2=c(24,48),y1=c(-0.001,-0.001),y2=c(-0.0009,-0.0009))

# and adding YUC5
# for scale_y_continuous(limits=c(-0.001,0.0013))
day<-data.frame(x1=rep(c(0,24,48,0,7,24,48),3),x2=rep(c(16,40,57,7,16,40,57),3),y1=rep(-0.001,7*3),y2=rep(-0.0009,7*3),treatment=rep(c("sun","sun","sun","shade","shade","shade","shade"),3),promoter=rep(c("YUC5","YUC8","YUC9"),each=7),color=rep(c("high","high","high","high","low","low","low"),3))
# 

write.csv(all,file="all.csv")

# start from this line
all<-read.csv("all.csv") # 061716
all<-all[,-1]
# drawing graph
p<-ggplot(all, aes(x=ZT,y=value))  # #
p <- p + labs(x="ZT (h)",y="log2 fold change of LUC activities (relative to ZT7)",title="Time-course promoter activities of auxin biosynthesis genes upon shade treatment")
# p<- p + theme(panel.background = element_rect(fill = "white"),axis.ticks = element_blank(), axis.text.x = element_text( angle = 330, vjust=1,hjust = 0, colour = "grey50"))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black")) # see http://felixfan.github.io/rstudy/2013/11/27/ggplot2-remove-grid-background-margin/
p <- p + scale_y_continuous(limits=c(-0.001,0.0013)) # unable to draw geom_rect...

# drawing night period
#p <- p + geom_rect(data=night,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2),fill="grey20", alpha=0.5, inherit.aes = FALSE)# draw dark/light cycle by white and gray rectangle    see http://stackoverflow.com/questions/4733182/how-to-highlight-time-ranges-on-a-plot
# night period (for bottom bar)
p <- p + geom_rect(data=night,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2),fill="black", inherit.aes = FALSE)# draw dark/light cycle by white and gray rectangle    see http://stackoverflow.com/questions/4733182/how-to-highlight-time-ranges-on-a-plot

# drawing day period (sun or shade)
p <- p + geom_rect(data=day,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=color), inherit.aes = FALSE) 
p <- p + facet_grid(treatment~promoter)
#p <- p + facet_wrap(treatment~promoter,scale="free")

# cf http://stackoverflow.com/questions/10097615/how-can-i-have-two-different-scale-fill-manual-active-in-a-ggplot-command
# cf http://stackoverflow.com/questions/19424944/how-to-set-background-color-of-ggplot-and-facet-grid-in-r

p <- p + scale_fill_manual( values = c("high" = "red","low" = "darkred"))
#p <- p + geom_point(aes(color=variable),alpha = 0.3)  #+ geom_smooth(method="ML")
# spline implimentation in ggplot2 http://www.lgbmi.com/2011/10/using-smooth-spline-in-stat-scale-in-ggplot2/
smooth.spline2 <- function(formula, data, ...) {
mat <- model.frame(formula, data)
smooth.spline(mat[, 2], mat[, 1])
}
predictdf.smooth.spline <- function(model, xseq, se, level) {
pred <- predict(model, xseq)
data.frame(x = xseq, y = pred$y)
}

# p <- p + geom_path(aes(color=variable))   + geom_smooth(method="smooth.spline2",se=F)
## another metohd http://f.briatte.org/teaching/ida/092_smoothing.html
# Load packages.
packages <- c("changepoint", "downloader", "ggplot2", "MASS", "reshape", "splines", 
    "XML")
packages <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x)
        library(x, character.only = TRUE)
    }
})

# p <- p + geom_path(aes(color=variable),linetype="dashed")   +geom_smooth(method = "rlm", se = FALSE, formula = y ~ ns(x, 8),size=2)
p <- p + geom_path(aes(alpha=variable))   +geom_smooth(method = "rlm", se = FALSE, formula = y ~ ns(x, 8),size=2,alpha=0.3)

p <- p + theme(legend.position="none")
p # 
#setwd("/Volumes/Data6/data_JM4/promYUC_luc2")
#ggsave(filename=paste("test2",Sys.time(),".pdf",sep=""),width=8,height=8) 

### dif could be better to be as log fold changes (work on functions) (see above)#
# better to elimiate YUC5 (not clear, 070314)

# variable needs to have exp name for geom_path()
## sample number (070714)
 ftable(all[all$ZT<7,c("treatment","promoter")],row.vars="treatment",col.vars="promoter")
          # promoter YUC5 YUC8 YUC9
# treatment                        
# shade                11   16   17
# sun                  14   11   16

#
#p<-ggplot(all, aes(x=ZT,y=value,color=line)) + geom_point() + geom_smooth()
#p<-ggplot(all, aes(x=ZT,y=value,color=line))  + geom_smooth(method="lm",se=F) #+ geom_point()
#p<-ggplot(all, aes(x=ZT,y=value))  + geom_smooth(method="lm",se=F) #+ geom_point() linier
p<-ggplot(all, aes(x=ZT,y=value))  + geom_smooth(stat="smooth",method="loess") #+ geom_point()

p<-p + facet_wrap(treatment~promoter,scales="free_y")
p<-p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
             panel.background = element_blank(), axis.line = element_line(colour = "black"))
night<-data.frame(ZT=c(16,40),value=c(-2e-04,2e-04))
#p<-p + geom_rect(data=night,aes(xmin=ZT,xmax=ZT+8,ymin=-0.00025,ymax=0.00025,alpha=0.1))
#p <- p + scale_y_continuous(limits=c(-5000,10000))
p
#ggsave(filename="test.pdf",width=11,height=20) 

# how to have different scale for different promoter? (092315)
all$TimeAfterTreatment<-all$ZT - 7 # time after treatment (hr)

# for ratio (revert log2(differences) to differeneces)
all$value2<-2^(all$value)

## plot separately and combined with cowplot?
# log2 differences (all$value)
# day<-data.frame(x1=rep(c(0,24,48,0,7,24,48),3),x2=rep(c(16,40,57,7,16,40,57),3),y1=rep(c(-0.00028,-0.0001,-0.0001),each=7),y2=rep(c(-0.00025,-7.8e-05,-7.8e-05),each=7),treatment=rep(c("sun","sun","sun","shade","shade","shade","shade"),3),promoter=rep(c("YUC5","YUC8","YUC9"),each=7),color=rep(c("high","high","high","high","low","low","low"),3))
# night<-data.frame(x1=rep(c(16,40),3),x2=rep(c(24,48),3),y1=rep(c(-0.00028,-0.0001,-0.0001),each=2),y2=rep(c(-0.00025,-7.8e-05,-7.8e-05),each=2),promoter=rep(c("YUC5","YUC8","YUC9"),each=2))

# for ratio (all$value2) for v3_ratio.png (092415)
day<-data.frame(x1=rep(c(0,24,48,0,7,24,48),3),x2=rep(c(16,40,57,7,16,40,57),3),y1=rep(c(0.9997,0.9999,0.9999),each=7),y2=rep(c(0.9997+0.00004,0.9999+0.000025,0.9999+0.000025),each=7),treatment=rep(c("sun","sun","sun","shade","shade","shade","shade"),3),promoter=rep(c("YUC5","YUC8","YUC9"),each=7),color=rep(c("high","high","high","high","low","low","low"),3))
night<-data.frame(x1=rep(c(16,40),3),x2=rep(c(24,48),3),y1=rep(c(0.9997,0.9999,0.9999),each=2),y2=rep(c(0.9997+0.00004,0.9999+0.000025,0.9999+0.000025),each=2),promoter=rep(c("YUC5","YUC8","YUC9"),each=2))


day[,1:2]<-day[,1:2]-7
night[,1:2]<-night[,1:2]-7
####
YUC5<-ggplot(all[all$promoter=="YUC5",], aes(x=TimeAfterTreatment,y=value2))  + geom_smooth(stat="smooth",method="loess",span = 0.2) #+ geom_point()
YUC5<-YUC5 + facet_grid(treatment~.)
YUC5<-YUC5 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
             panel.background = element_blank(), axis.line = element_line(colour = "black"))
# drawing night period
# night period (for bottom bar)
YUC5 <- YUC5 + geom_rect(data=night[night$promoter=="YUC5",],mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2),fill="black", inherit.aes = FALSE)# draw dark/light cycle by white and gray rectangle    see http://stackoverflow.com/questions/4733182/how-to-highlight-time-ranges-on-a-plot
# drawing day period (sun or shade)
YUC5 <-YUC5 + geom_rect(data=day[day$promoter=="YUC5",],mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=color), inherit.aes = FALSE) 
# change bar color
YUC5 <- YUC5 + scale_fill_manual( values = c("high" = "red","low" = "darkred"))
YUC5<-YUC5 + theme(legend.position = "none", axis.text.y=element_text(size=7)) + labs(x="Time after treatment (h)",y="Relative signal to -1h")
YUC5
## YUC8
YUC8<-ggplot(all[all$promoter=="YUC8",], aes(x=TimeAfterTreatment,y=value2))  + geom_smooth(stat="smooth",method="loess",span = 0.2) #+ geom_point()

YUC8<-YUC8 + facet_grid(treatment~.)
YUC8<-YUC8 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                   panel.background = element_blank(), axis.line = element_line(colour = "black"))
# drawing night period
# night period (for bottom bar)
YUC8 <- YUC8 + geom_rect(data=night[night$promoter=="YUC8",],mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2),fill="black", inherit.aes = FALSE)# draw dark/light cycle by white and gray rectangle    see http://stackoverflow.com/questions/4733182/how-to-highlight-time-ranges-on-a-plot
# drawing day period (sun or shade)
YUC8 <-YUC8 + geom_rect(data=day[day$promoter=="YUC8",],mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=color), inherit.aes = FALSE) 
# change bar color
YUC8 <- YUC8 + scale_fill_manual( values = c("high" = "red","low" = "darkred"))
#p <- p + scale_y_continuous(limits=c(-5000,10000))
YUC8<-YUC8 + theme(legend.position = "none", axis.text.y=element_text(size=7),axis.title.y=element_blank()) +labs(x="Time after treatment (h)",y="Relative signal to -1h")
YUC8
## YUC9
YUC9<-ggplot(all[all$promoter=="YUC9",], aes(x=TimeAfterTreatment,y=value2))  + geom_smooth(stat="smooth",method="loess",span = 0.2) #+ geom_point()
YUC9<-YUC9 + facet_grid(treatment~.)
YUC9<-YUC9 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                   panel.background = element_blank(), axis.line = element_line(colour = "black"))
# drawing night period
# night period (for bottom bar)
YUC9 <- YUC9 + geom_rect(data=night[night$promoter=="YUC9",],mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2),fill="black", inherit.aes = FALSE)# draw dark/light cycle by white and gray rectangle    see http://stackoverflow.com/questions/4733182/how-to-highlight-time-ranges-on-a-plot
# drawing day period (sun or shade)
YUC9 <-YUC9 + geom_rect(data=day[day$promoter=="YUC9",],mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=color), inherit.aes = FALSE) 
# change bar color
YUC9 <- YUC9 + scale_fill_manual( values = c("high" = "red","low" = "darkred"))
YUC9<-YUC9 + theme(legend.position = "none", axis.text.y=element_text(size=7),axis.title.y=element_blank()) + labs(x="Time after treatment (h)",y="Relative signal to -1h")
#
# cowplot
library(cowplot)
plot.YUC589<-plot_grid(YUC5,YUC8,YUC9,ncol=3,labels=c("A","B","C"))
save_plot("promYUC589_LUC2v3_ratio.png", plot.YUC589,
          ncol = 3, # we're saving a grid plot of 2 columns
          nrow = 1, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 0.8)
# save_plot("promYUC589_LUC2v2.png", plot.YUC589,
#           ncol = 3, # we're saving a grid plot of 2 columns
#           nrow = 1, # and 2 rows
#           # each individual subplot should have an aspect ratio of 1.3
#           base_aspect_ratio = 0.8)


# how many plants measured?
all.unique<-all[!duplicated(all$variable),]
tapply(all.unique$variable,list(all.unique$treatment,all.unique$promoter),length)
#       YUC5 YUC8 YUC9
# shade   25   16   23
# sun     25   13   19
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
all$value2<-2^(all$value)


## replace "sun" into "high R/FR", "shade into "low R/FR"
all$treatment<-as.character(all$treatment)
all$treatment<-gsub("sun","high R/FR",all$treatment)
all$treatment<-gsub("shade","low R/FR",all$treatment)
all$treatment<-factor(all$treatment,levels=c("high R/FR","low R/FR"))
# for ratio (all$value2) for v4_ratio.png (061716)
day<-data.frame(x1=rep(c(0,24,48,0,7,24,48),3),
                x2=rep(c(16,40,57,7,16,40,57),3),
                y1=rep(c(0.9997,0.9997-0.00004,0.9999,0.9999-0.000028,0.9999-0.00007,0.9999-0.00007-0.000036),c(3,4,3,4,3,4)),
                y2=rep(c(0.9997+0.000037,0.9997-0.00004+0.000037,0.9999+0.000025,0.9999-0.000003,0.9999-0.00007+0.000032,0.9999-0.00007-0.000036+0.000032),c(3,4,3,4,3,4)),
                treatment=rep(c("high R/FR","high R/FR","high R/FR","low R/FR","low R/FR","low R/FR","low R/FR"),3),
                promoter=rep(c("YUC5","YUC8","YUC9"),each=7),
                color=rep(c("high R/FR","high R/FR","high R/FR","high R/FR","low R/FR","low R/FR","low R/FR"),3))
night<-data.frame(x1=rep(c(16,40),6),
                  x2=rep(c(24,48),6),
                  y1=rep(c(0.9997,0.9997-0.00004,0.9999,0.9999-0.000028,0.9999-0.00007,0.9999-0.00007-0.000036),each=2),
                  y2=rep(c(0.9997+0.000037,0.9997-0.00004+0.000037,0.9999+0.000025,0.9999-0.000003,0.9999-0.00007+0.000032,0.9999-0.00007-0.000036+0.000032),each=2),
                  treatment=rep(c("high R/FR","high R/FR","low R/FR","low R/FR"),3),
                  promoter=rep(c("YUC5","YUC8","YUC9"),each=4)                  
                  )
day[,1:2]<-day[,1:2]-7
night[,1:2]<-night[,1:2]-7
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
YUC5 <- YUC5 + scale_fill_manual( values = c("high R/FR" = "red","low R/FR" = "darkred")) +
  scale_colour_manual( values = c("high R/FR" = "red","low R/FR" = "darkred"))

YUC5<-YUC5 + theme(legend.position = "none", axis.text.y=element_text(size=7),axis.title.y=element_blank()) + labs(x="Time after treatment (h)",y="Relative signal to -1h")
# change range of y
#YUC5<-YUC5 + scale_y_continuous(limits=c(0.9995,1.0003))
# write promoter name
YUC5<-YUC5 + annotate("text",label="YUC5",x=-3,y=1.00017)
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
YUC8 <- YUC8 + scale_fill_manual( values = c("high R/FR" = "red","low R/FR" = "darkred")) +
  scale_colour_manual( values = c("high R/FR" = "red","low R/FR" = "darkred"))

YUC8<-YUC8 + theme(legend.position = "none", axis.text.y=element_text(size=7),axis.title.y=element_blank()) + labs(x="Time after treatment (h)",y="Relative signal to -1h")
# write promoter name
YUC8<-YUC8 + annotate("text",label="YUC8",x=-3,y=1.00021)
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
YUC9 <- YUC9 + scale_fill_manual( values = c("high R/FR" = "red","low R/FR" = "darkred")) +
  scale_colour_manual( values = c("high R/FR" = "red","low R/FR" = "darkred"))

YUC9<-YUC9 + theme(axis.text.y=element_text(size=7),axis.title.y=element_blank(),legend.text=element_text(size=7),legend.title=element_blank()) + labs(x="Time after treatment (h)",y="Relative signal to -1h")
# write promoter name
YUC9<-YUC9 + annotate("text",label="YUC9",x=-3,y=1.00023)
YUC9
# see http://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
library(ggplot2); library(gridExtra)
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 
legend <- g_legend(my.category.heatmap.graph[[i]]) 


# cowplot
library(cowplot)
plot.YUC589<-plot_grid(YUC5,YUC8,YUC9,ncol=3,labels=c("A","B","C"),
                       scale=0.9,vjust=0,rel_widths=c(1,1,1.38))
save_plot("Fig1_promYUC589_LUC2v4_ratio.png", plot.YUC589,
          ncol = 3, # we're saving a grid plot of 2 columns
          nrow = 1, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 0.8)

##### background subtraction (071216) ######
# start from this line
## exp3 set1
exp3files<-read.table("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp3/exp3files.txt",skip=2)

first_day_ZT0<-"Jun/11/06:00"
plant.name<-c(
  paste("shade_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1"),sep=""),
  paste("shade_plate2_YUC8_line2_rep1",sep=""),
  paste("shade_plate3_",c("YUC9_line3_rep1","YUC5_line3_rep1","YUC8_line3_rep2"),sep=""),
  paste("sun_plate1_",c("YUC9_line2_rep1","YUC5_line2_rep1","YUC8_line2_rep1"),sep=""),
  paste("sun_plate2_",c("YUC9_line1_rep1","YUC5_line1_rep1"),sep="")
)
luc2.data<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp3/exp3set1data.csv") 
luc2.data.multi<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp3/exp3set1data.multi.csv") # need to scale!

exp3set1.luc2.data.multi.log2<-format.luc2data.multi.exp3set1special.log2.BS(plant.name,luc2.data,luc2.data.multi,files=exp3files,first_day_ZT0=first_day_ZT0) # under construction (070214)
head(exp3set1.luc2.data.multi.log2)
exp3set1.luc2.data.multi.log2$value
exp3set1.luc2.data.multi.log2<-exp3set1.luc2.data.multi.log2[!exp3set1.luc2.data.multi.log2$value=="-Inf",]
exp3set1.luc2.data.multi.log2$variable<-paste("exp3_set1_",exp3set1.luc2.data.multi.log2$variable,sep="")

p<-ggplot(exp3set1.luc2.data.multi.log2[exp3set1.luc2.data.multi.log2$value<3,], aes(x=ZT,y=value,color=line)) + geom_point()
p<-p + facet_wrap(treatment~promoter)
p # 
## exp3 set2
plant.name<-c(
  paste("shade_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1"),sep=""),
  paste("shade_plate2_",c("YUC9_line2_rep1","YUC5_line2_rep1","YUC8_line2_rep1"),sep=""),
  paste("shade_plate3_",c("YUC9_line3_rep1","YUC5_line3_rep1"),sep=""),
  paste("sun_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1"),sep=""),
  paste("sun_plate2_",c("YUC9_line2_rep1","YUC8_line2_rep1"),sep=""),
  paste("sun_plate3_",c("YUC9_line3_rep1","YUC5_line3_rep1"),sep=""))
luc2.data<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp3/exp3set2data.csv") 
luc2.data.multi<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp3/exp3set2data.multi.csv") # need to scale!
exp3set2.luc2.data.multi.log2<-format.luc2data.multi.exp3set2special.log2.BS(plant.name,luc2.data,luc2.data.multi,files=exp3files,first_day_ZT0=first_day_ZT0) # under construction (070314)
exp3set2.luc2.data.multi.log2$variable<-paste("exp3_set2_",exp3set2.luc2.data.multi.log2$variable,sep="")

head(exp3set2.luc2.data.multi.log2)
# remove fake "0"
exp3set2.luc2.data.multi.log2<-exp3set2.luc2.data.multi.log2[!(exp3set2.luc2.data.multi.log2$value==0&exp3set2.luc2.data.multi.log2$ZT>7),]
p<-ggplot(exp3set2.luc2.data.multi.log2[exp3set2.luc2.data.multi.log2$value<4,], aes(x=ZT,y=value,color=line)) + geom_point()
p<-p + facet_wrap(treatment~promoter)
p # 
#
# exp4
# exp4 set1 (add "rep" info to avoid same column name)
plant.name<-c(
  paste("shade_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1"),sep=""),
  paste("shade_plate2_YUC8_line2_rep1",sep=""),
  paste("shade_plate3_YUC5_line3_rep1",sep=""),
  paste("sun_plate1_",c("YUC9_line1_rep1","YUC5_line1_rep1"),sep=""),
  paste("sun_plate3_",c("YUC9_line2_rep1","YUC8_line2_rep1","YUC5_line2_rep1"),sep=""),
  paste("sun_plate4_YUC5_line3_rep1",sep="")
)
luc2.data<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp4/exp4set1data.csv") 
luc2.data.multi<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp4/exp4set1data.multi.csv") # need to scale!
exp4files<-read.table("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp4/exp4files.txt",skip=2)
first_day_ZT0<-"Jun/18/06:00"
#exp4set1.luc2.data.multi<-format.luc2data.multi(plant.name,luc2.data,luc2.data.multi,files=exp4files,first_day_ZT0=first_day_ZT0) # under construction (070214)
exp4set1.luc2.data.multi.log2<-format.luc2data.multi.log2.BS(plant.name,luc2.data,luc2.data.multi,files=exp4files,first_day_ZT0=first_day_ZT0) # under construction (070214)
exp4set1.luc2.data.multi.log2$variable<-paste("exp4_set1_",exp4set1.luc2.data.multi.log2$variable,sep="")

head(exp4set1.luc2.data.multi.log2)
library(ggplot2)
p<-ggplot(exp4set1.luc2.data.multi.log2, aes(x=ZT,y=value,color=line)) + geom_point()
p<-p + facet_wrap(treatment~promoter)
p # 

# exp4 set2
plant.name<-c(
  paste("shade_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1","YUC5_line1_rep1"),sep=""),
  paste("shade_plate2_",c("YUC9_line2_rep1","YUC8_line2_rep1"),sep=""),
  paste("shade_plate3_",c("YUC9_line3_rep1","YUC8_line2_rep1","YUC5_line3_rep1"),sep=""),
  paste("sun_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1","YUC5_line1_rep1"),sep=""),
  paste("sun_plate2_",c("YUC9_line2_rep1","YUC8_line2_rep1","YUC5_line2_rep1"),sep=""),
  "sun_plate3_YUC9_line3_rep1"
)
luc2.data<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp4/exp4set2data.csv") 
luc2.data.multi<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp4/exp4set2data.multi.csv") # need to scale!
#luc2.data.multi<-luc2.data.multi[grep("lab",luc2.data.multi$Label,invert=TRUE),]
exp4files<-read.table("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp4/exp4files.txt",skip=2)
first_day_ZT0<-"Jun/18/06:00"
#exp4set2.luc2.data.multi<-format.luc2data.multi(plant.name,luc2.data,luc2.data.multi,files=exp4files,first_day_ZT0=first_day_ZT0) # under construction (070214)
exp4set2.luc2.data.multi.log2<-format.luc2data.multi.log2.BS(plant.name,luc2.data,luc2.data.multi,files=exp4files,first_day_ZT0=first_day_ZT0) # under construction (070214)
exp4set2.luc2.data.multi.log2$variable<-paste("exp4_set2_",exp4set2.luc2.data.multi.log2$variable,sep="")

head(exp4set2.luc2.data.multi.log2)
library(ggplot2)
p<-ggplot(exp4set2.luc2.data.multi.log2, aes(x=ZT,y=value,color=line)) + geom_point()
p<-p + facet_wrap(treatment~promoter)
p # 
# summary
p<-ggplot(rbind(exp4set1.luc2.data.multi.log2,exp4set2.luc2.data.multi.log2), aes(x=ZT,y=value,color=line)) + geom_point()
p<-p + facet_wrap(treatment~promoter)
p # 

# exp5
# ls -lt > exp5set2files.txt # done in set2
exp5files<-read.table("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp5/exp5files.txt", skip=2) # exp5label.tif is for recording lables of each plants
#files<-files[!files$V9=="exp5label.tif",]  # remove "exp5label.tif"
first_day_ZT0<-"Jun/25/06:00"

# exp5 set1
plant.name<-c(
  paste("shade_plate1_",c("YUC9_line1_rep1","YUC9_line1_rep2","YUC8_line1_rep1","YUC8_line1_rep1","YUC5_line1_rep1"),sep=""),
  paste("shade_plate2_",c("YUC9_line3_rep1","YUC9_line3_rep2","YUC5_line3_rep1"),sep=""),
  paste("sun_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1","YUC8_line1_rep2"),sep=""),
  paste("sun_plate3_","YUC5_line3_rep1",sep="")
)
luc2.data<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp5/exp5set1data.csv") 
luc2.data.multi<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp5/exp5set1data.multi.csv") # need to scale!
#luc2.data.multi<-luc2.data.multi[grep("lab",luc2.data.multi$Label,invert=TRUE),]
#exp5set1.luc2.data.multi<-format.luc2data.multi(plant.name,luc2.data,luc2.data.multi,files=exp5files,first_day_ZT0=first_day_ZT0) # under construction (070214)
# log2 ratio
exp5set1.luc2.data.multi.log2<-format.luc2data.multi.log2.BS(plant.name,luc2.data,luc2.data.multi,files=exp5files,first_day_ZT0=first_day_ZT0) # under construction (070214)
exp5set1.luc2.data.multi.log2<-exp5set1.luc2.data.multi.log2[exp5set1.luc2.data.multi.log2$ZT< 48,]
exp5set1.luc2.data.multi.log2$variable<-paste("exp5_set1_",exp5set1.luc2.data.multi.log2$variable,sep="")

head(exp5set1.luc2.data.multi.log2)
library(ggplot2)
p<-ggplot(exp5set1.luc2.data.multi.log2, aes(x=ZT,y=value,color=line)) + geom_point()
p<-p + facet_wrap(treatment~promoter)
p # 



# exp5 set2
# enter experimental information
# input genotype and rep by looking "exp5lable.tif" by mannually. Transgenic "line" has to be consistent throughout experiments.
plant.name<-c(
  paste("shade_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1","YUC5_line1_rep1","YUC8_line1_rep2","YUC5_line1_rep2"),sep=""),
  paste("shade_plate2_",c("YUC9_line2_rep1","YUC9_line2_rep2","YUC5_line2_rep1","YUC5_line2_rep2"),sep=""),
  paste("shade_plate3_",c("YUC9_line3_rep1","YUC9_line3_rep2","YUC5_line3_rep1"),sep=""),
  paste("sun_plate1_",c("YUC9_line1_rep1","YUC8_line1_rep1","YUC5_line1_rep1","YUC5_line1_rep2"),sep=""),
  paste("sun_plate2_",c("YUC8_line2_rep1","YUC5_line2_rep1","YUC5_line2_rep2"),sep=""),
  paste("sun_plate3_",c("YUC9_line3_rep1","YUC5_line3_rep1","YUC5_line3_rep2"),sep="")) 
#
luc2.data<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp5/exp5set2data.csv") 
luc2.data.multi<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp5/exp5set2data.multi.csv") # need to scale!

#luc2.data.multi<-luc2.data.multi[grep("lab",luc2.data.multi$Label,invert=TRUE),]
# exp5set2.luc2.data.multi<-format.luc2data.multi(plant.name,luc2.data,luc2.data.multi,files=exp5files,first_day_ZT0=first_day_ZT0) # under construction (070214)
exp5set2.luc2.data.multi.log2<-format.luc2data.multi.log2.BS(plant.name,luc2.data,luc2.data.multi,files=exp5files,first_day_ZT0=first_day_ZT0) # under construction (070214)

head(exp5set2.luc2.data.multi.log2)

## eliminate day3 (I did not turn on light ...)
exp5set2.luc2.data.multi.log2<-exp5set2.luc2.data.multi.log2[exp5set2.luc2.data.multi.log2$ZT< 48,]
exp5set2.luc2.data.multi.log2$variable<-paste("exp5_set2_",exp5set2.luc2.data.multi.log2$variable,sep="")

p<-ggplot(exp5set2.luc2.data.multi.log2, aes(x=ZT,y=value,color=line)) + geom_point() 
p<-p + facet_wrap(treatment~promoter)
p # 
# exp6 
exp6files<-read.table("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp6/exp6files.txt", skip=2) # exp5label.tif is for recording lables of each plants
first_day_ZT0<-"July/4/06:00"

## exp6 set1
plant.name<-c(
  paste("shade_plate1_",c("YUC5_line1_rep1","YUC5_line1_rep2"),sep=""),
  paste("shade_plate2_",c("YUC8_line2_rep1","YUC8_line2_rep2","YUC5_line2_rep1","YUC5_line2_rep2"),sep=""),
  paste("shade_plate3_",c("YUC9_line3_rep1","YUC5_line3_rep1","YUC5_line3_rep2"),sep=""),
  paste("sun_plate1_",c("YUC9_line1_rep1","YUC5_line1_rep1"),sep=""),
  paste("sun_plate2_",c("YUC9_line2_rep1","YUC8_line2_rep1","YUC8_line2_rep2","YUC5_line2_rep1"),sep=""),
  paste("sun_plate3_",c("YUC9_line3_rep1","YUC9_line3_rep2","YUC8_line3_rep1","YUC5_line3_rep1","YUC5_line3_rep2"),sep="")
)
luc2.data<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp6/exp6set1data.csv") 
luc2.data.multi<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp6/exp6set1data.multi.csv") 
exp6set1.luc2.data.multi.log2<-format.luc2data.multi.log2.BS(plant.name,luc2.data,luc2.data.multi,files=exp6files,first_day_ZT0=first_day_ZT0)
exp6set1.luc2.data.multi.log2$variable<-paste("exp6_set1_",exp6set1.luc2.data.multi.log2$variable,sep="")

p<-ggplot(exp6set1.luc2.data.multi.log2, aes(x=ZT,y=value,color=line)) + geom_point() 
p<-p + facet_wrap(treatment~promoter)
p # 
## exp6 set2
plant.name<-c(
  paste("shade_plate1_",c("YUC5_line1_rep1","YUC5_line1_rep2"),sep=""),
  paste("shade_plate2_",c("YUC9_line2_rep1","YUC9_line2_rep2","YUC8_line2_rep1","YUC5_line2_rep1","YUC5_line2_rep2"),sep=""),
  paste("shade_plate3_",c("YUC9_line3_rep1","YUC9_line3_rep2","YUC5_line3_rep1","YUC5_line3_rep2"),sep=""),
  paste("sun_plate1_",c("YUC5_line1_rep1","YUC5_line1_rep2"),sep=""),
  paste("sun_plate2_",c("YUC5_line2_rep1","YUC5_line2_rep2"),sep=""),
  paste("sun_plate3_",c("YUC9_line3_rep1","YUC9_line3_rep2","YUC5_line3_rep1","YUC5_line3_rep2"),sep="")
)
luc2.data<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp6/exp6set2data.csv") 
luc2.data.multi<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp6/exp6set2data.multi.csv") 
exp6set2.luc2.data.multi.log2<-format.luc2data.multi.log2.BS(plant.name,luc2.data,luc2.data.multi,files=exp6files,first_day_ZT0=first_day_ZT0)
exp6set2.luc2.data.multi.log2$variable<-paste("exp6_set2_",exp6set2.luc2.data.multi.log2$variable,sep="")
p<-ggplot(exp6set2.luc2.data.multi.log2, aes(x=ZT,y=value,color=line)) + geom_point() 
p<-p + facet_wrap(treatment~promoter)
p # 

# summary (log2 fold change version)
all<-rbind(exp3set1.luc2.data.multi.log2,exp3set2.luc2.data.multi.log2,exp4set1.luc2.data.multi.log2,exp4set2.luc2.data.multi.log2,exp5set1.luc2.data.multi.log2,exp5set2.luc2.data.multi.log2,exp6set1.luc2.data.multi.log2,exp6set2.luc2.data.multi.log2)
#all.YUC8and9<-all[!all$promoter=="YUC5",]
all$variable
setwd("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/")
write.csv(all, file="all.background.subtraction.csv")
p<-ggplot(all[all$value<3,], aes(x=ZT,y=value,color=variable)) + geom_point() 
p<-p + facet_wrap(treatment~promoter) + theme(legend.position="none")
p # 
#ggsave(file="promYUC589_LUC2_ratio_background_subtraction_each_trace.png")

# background subtraction version
####################################################
# sun + shade trace in one plat version (071516)
####################################################
# start from this line
setwd("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/")
all<-read.csv("all.background.subtraction.csv") 
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
all$value2<-2^(all$value)


## replace "sun" into "high R/FR", "shade into "low R/FR"
all$treatment<-as.character(all$treatment)
all$treatment<-gsub("sun","high R/FR",all$treatment)
all$treatment<-gsub("shade","low R/FR",all$treatment)
all$treatment<-factor(all$treatment,levels=c("high R/FR","low R/FR"))
# for ratio (all$value2) for v4_ratio.png (061716)
day<-data.frame(x1=rep(c(0,24,48,0,7,24,48),3),
                x2=rep(c(16,40,57,7,16,40,57),3),
                y1=rep(c(0.9997,0.9997-0.00004,0.9999,0.9999-0.000028,0.9999-0.00007,0.9999-0.00007-0.000036),c(3,4,3,4,3,4)),
                y2=rep(c(0.9997+0.000037,0.9997-0.00004+0.000037,0.9999+0.000025,0.9999-0.000003,0.9999-0.00007+0.000032,0.9999-0.00007-0.000036+0.000032),c(3,4,3,4,3,4)),
                treatment=rep(c("high R/FR","high R/FR","high R/FR","low R/FR","low R/FR","low R/FR","low R/FR"),3),
                promoter=rep(c("YUC5","YUC8","YUC9"),each=7),
                color=rep(c("high R/FR","high R/FR","high R/FR","high R/FR","low R/FR","low R/FR","low R/FR"),3))
night<-data.frame(x1=rep(c(16,40),6),
                  x2=rep(c(24,48),6),
                  y1=rep(c(0.9997,0.9997-0.00004,0.9999,0.9999-0.000028,0.9999-0.00007,0.9999-0.00007-0.000036),each=2),
                  y2=rep(c(0.9997+0.000037,0.9997-0.00004+0.000037,0.9999+0.000025,0.9999-0.000003,0.9999-0.00007+0.000032,0.9999-0.00007-0.000036+0.000032),each=2),
                  treatment=rep(c("high R/FR","high R/FR","low R/FR","low R/FR"),3),
                  promoter=rep(c("YUC5","YUC8","YUC9"),each=4)                  
)
day[,1:2]<-day[,1:2]-7
night[,1:2]<-night[,1:2]-7
# remove outliers
all<-all[all$value2<10,]
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
YUC5 <- YUC5 + scale_fill_manual( values = c("high R/FR" = "red","low R/FR" = "darkred")) +
  scale_colour_manual( values = c("high R/FR" = "red","low R/FR" = "darkred"))

YUC5<-YUC5 + theme(legend.position = "none", axis.text.y=element_text(size=7)) + labs(x="Time after treatment (h)",y="Relative signal to -1h")
# change range of y
#YUC5<-YUC5 + scale_y_continuous(limits=c(0.9995,1.0003))
# write promoter name
YUC5<-YUC5 + annotate("text",label="YUC5",x=-2,y=1.00017)
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
YUC8 <- YUC8 + scale_fill_manual( values = c("high R/FR" = "red","low R/FR" = "darkred")) +
  scale_colour_manual( values = c("high R/FR" = "red","low R/FR" = "darkred"))

YUC8<-YUC8 + theme(legend.position = "none", axis.text.y=element_text(size=7)) + labs(x="Time after treatment (h)",y="Relative signal to -1h")
# write promoter name
YUC8<-YUC8 + annotate("text",label="YUC8",x=-2,y=1.00021)
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
YUC9 <- YUC9 + scale_fill_manual( values = c("high R/FR" = "red","low R/FR" = "darkred")) +
  scale_colour_manual( values = c("high R/FR" = "red","low R/FR" = "darkred"))

YUC9<-YUC9 + theme(axis.text.y=element_text(size=7),legend.text=element_text(size=7),legend.title=element_blank()) + labs(x="Time after treatment (h)",y="Relative signal to -1h")
# write promoter name
YUC9<-YUC9 + annotate("text",label="YUC9",x=-2,y=1.00023)
YUC9
# cowplot
library(cowplot)
plot.YUC589<-plot_grid(YUC5,YUC8,YUC9,ncol=3,labels=c("A","B","C"),
                       scale=0.9,vjust=0,rel_widths=c(1,1,1.38))
save_plot("Fig1_promYUC589_LUC2_ratio_background_subtraction_v1.png", plot.YUC589,
          ncol = 3, # we're saving a grid plot of 2 columns
          nrow = 1, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 0.8)


