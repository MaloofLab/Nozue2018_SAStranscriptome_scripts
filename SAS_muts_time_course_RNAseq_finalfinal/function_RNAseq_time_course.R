
#######
##
## function_RNAseq_time_course.R
##  
###################################################
library("annotate") #
library(edgeR)
library(ggplot2);library(reshape2);library(grid);library(class);library(MASS);library(plyr)
library(kohonen) # for SOM analysis
library(scales) # for muted
library(WGCNA);library(ShortRead);library(goseq);library(GO.db); library("org.At.tair.db")
options(stringsAsFactors = FALSE) # for WGCNA


# see http://www.bioconductor.org/install/ for installation of these packages 
#library(ggdendro) # for dendrogram
library(lmerTest) # for significant analysis
# this does not work while knitting
#TAIR10_gene_descriptions<-read.csv(file.path(homedir2,"../../Nozue2016_SAStranscriptome_data/input/TAIR10_functional_descriptions.csv") )
# how to do?
#setwd("../")
#TAIR10_gene_descriptions<-read.csv("../Nozue2016_SAStranscriptome_data/input/TAIR10_functional_descriptions.csv")
#TAIR10_gene_descriptions<-read.csv("Nozue2016_SAStranscriptome_data/input/TAIR10_functional_descriptions.csv")
TAIR10_gene_descriptions<-read.csv(file.path(homedir,"..","..","Nozue2016_SAStranscriptome_data","input","TAIR10_functional_descriptions.csv"))

#setwd(homedir2)
TAIR10_gene_descriptions$AGI<-gsub("([[:print:]]+)(.[[:digit:]]+)","\\1",TAIR10_gene_descriptions$"Model_name")
# 
#conversion.table.nam<-data.frame(num=20:27, genotype=c("Blh_1","Jea","Shahdara","Col_0","Cvi_0","Bur_0","Oy_0","Ita_0"))
conversion.table<-data.frame(num=c(1,10,13,14,17,19,3,4,5,6,7,8,9), genotype=c("Col","PAR1_RNAi09","pif3","mida9_4","aos","co_9","hy5","jar1","kat1_2","phyB","pif45","spt_11","yuc2589"))
#all.equal(as.character(conversion.table[c(1,10,13,14,17,19,4,5,6,7,8,9),"genotype"]),conversion.table2[,"genotype"])
expression.pattern.graph<-function(data.cpm,target.genes,samples){# require ggplot2, reshape2 packages, this is only for SAS timecourse data
  temp.data<-as.data.frame(t(data.cpm[rownames(data.cpm) %in% target.genes,]))
  temp.data$file<-rownames(temp.data)
  temp.data.samples<-merge(samples,temp.data,by="file")
  print(temp.data.samples)
  temp.data.samples.melt<-melt(temp.data.samples[,-1],id=names(samples)[-1])
  
  temp.data.samples.melt$time<-factor(temp.data.samples.melt$time,levels=c("1","4","16","25","49"))
  q<-ggplot(temp.data.samples.melt,aes(x=time,y=value,color=batch)) + geom_point(alpha = 0.5) + facet_grid(variable~trt,scale="free") + theme(strip.text.y=element_text(angle=0))
  #  q<-q + scale_y_continuous(trans=log2_trans())
  return(q)
}
#expression.pattern.graph(data_Col.int.cpm,target.genes=result.genes[1:10],samples=samples)
#expression.pattern.graph(samples.nolow.Col3.cpm,target.genes=rownames(topTags(dge.nolow.Col.batch.glm.lrt,n=50)$table),samples=samples.nolow.Col3)

expression.pattern.graph2<-function(data.cpm,target.genes,samples){# require ggplot2, reshape2 packages, this is only for SAS timecourse data
  temp.data<-as.data.frame(t(data.cpm[rownames(data.cpm) %in% target.genes,]))
  temp.data$file<-rownames(temp.data)
  temp.data.samples<-merge(samples,temp.data,by="file")
  print(temp.data.samples)
  temp.data.samples.melt<-melt(temp.data.samples[,-1],id=names(samples)[-1])
  
  temp.data.samples.melt$time<-factor(temp.data.samples.melt$time,levels=c("1","4","16","25","49"))
  q<-ggplot(temp.data.samples.melt,aes(x=time,y=value)) + geom_boxplot(alpha = 0.5) + facet_grid(variable~trt,scale="free") + theme(strip.text.y=element_text(angle=0))
  #  q<-q + scale_y_continuous(trans=log2_trans())
  return(q)
}
# fold change (fixed; 060515)
expression.pattern.graph2b<-function(summary.vst.response.kazu,target.genes,plot.order){# require ggplot2, reshape2 packages, this is only for SAS timecourse data  
  temp.data<-as.data.frame(t(summary.vst.response.kazu[rownames(summary.vst.response.kazu) %in% target.genes,]))
  temp.data$genotype<-gsub("([[:digit:]]+)(_)(1|4|25|49)(hrA)","\\1",rownames(temp.data))
  temp.data$time<-gsub("([[:digit:]]+)(_)(1|4|25|49)(hrA)","\\3",rownames(temp.data)) 
  # convert genotype number into actual name
  conversion.table<-data.frame(num=1:27, genotype=c("Col","AT5G02540_1","hy5","jar1","kat1_2","phyB","pif45","spt_11","yuc2589","PAR1_RNAi09","coi1_16","phyAB","pif3","mida9_4","bsk5_1","sto","aos","argos","co_9","Blh_1","Jea","Shahdara","Col_0","Cvi_0","Bur_0","Oy_0","Ita_0"))
  temp.data$genotype<-as.character(temp.data$genotype)
  for(i in 1:27) {
    conversion.table[i,]
    temp.data[temp.data$genotype==conversion.table[i,"num"],"genotype"]<-rep(as.character(conversion.table[i,"genotype"]),sum(as.integer(temp.data$genotype==conversion.table[i,"num"])))
  }
  #
  print(temp.data)
  temp.data.melt<-melt(temp.data,id=c("time","genotype"))
  temp.data.melt$time<-factor(temp.data.melt$time,levels=c("1","4","25","49"))  
  temp.data.melt$variable<-factor(temp.data.melt$variable,levels=plot.order)
  q<-ggplot(temp.data.melt,aes(x=time,y=value)) + geom_point() + facet_grid(variable~genotype,scale="free") + theme(strip.text.y=element_text(angle=0))
  #  q<-q + scale_y_continuous(trans=log2_trans())
  return(q)
}

expression.pattern.graph3<-function(data.cpm,target.genes,plot.order){# require ggplot2, reshape2 packages, this is only for SAS timecourse data (mean type, eg. "summary4")
  temp.data<-as.data.frame(t(data.cpm[rownames(data.cpm) %in% target.genes,]))
  print(temp.data)
  temp.data$genotype<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\1",rownames(temp.data))
  temp.data$trt<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\2",rownames(temp.data))
  temp.data$time<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\3",rownames(temp.data))
  # convert genotype number into actual name
  conversion.table<-data.frame(num=1:27, genotype=c("Col","AT5G02540_1","hy5","jar1","kat1_2","phyB","pif45","spt_11","yuc2589","PAR1_RNAi09","coi1_16","phyAB","pif3","mida9_4","bsk5_1","sto","aos","argos","co_9","Blh_1","Jea","Shahdara","Col_0","Cvi_0","Bur_0","Oy_0","Ita_0"))
  temp.data$genotype<-as.character(temp.data$genotype)
  for(i in 1:27) {
    conversion.table[i,]
    temp.data[temp.data$genotype==conversion.table[i,"num"],"genotype"]<-rep(as.character(conversion.table[i,"genotype"]),sum(as.integer(temp.data$genotype==conversion.table[i,"num"])))
  }
  #
  temp.data.melt<-melt(temp.data,id=c("trt","time","genotype"))
  print(temp.data.melt)
  #temp.data.samples.melt$time<-factor(temp.data.samples.melt$time,levels=c("1","4","16","25","49"))
  temp.data.melt$time<-factor(temp.data.melt$time,levels=c("1","4","25","49"))
  temp.data.melt$variable<-factor(temp.data.melt$variable,levels=plot.order)
  
  q<-ggplot(temp.data.melt,aes(x=time,y=value,color=trt)) + geom_point(alpha = 0.5) 
  q<-q + facet_grid(variable~genotype,scale="free") + theme(strip.text.y=element_text(angle=0))
  #  q<-q + scale_y_continuous(trans=log2_trans())
  return(q)
}
#
expression.pattern.graph3b<-function(data.mean,data.expression,target.genes,plot.order){# require ggplot2, reshape2 packages, this is only for SAS timecourse data (jitter and mean)
  temp.data<-as.data.frame(t(data.mean[rownames(data.mean) %in% target.genes,]))
  print(temp.data)
  temp.data$genotype<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\1",rownames(temp.data))
  temp.data$trt<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\2",rownames(temp.data))
  temp.data$time<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\3",rownames(temp.data))
  # convert genotype number into actual name
  conversion.table<-data.frame(num=1:27, genotype=c("Col","AT5G02540_1","hy5","jar1","kat1_2","phyB","pif45","spt_11","yuc2589","PAR1_RNAi09","coi1_16","phyAB","pif3","mida9_4","bsk5_1","sto","aos","argos","co_9","Blh_1","Jea","Shahdara","Col_0","Cvi_0","Bur_0","Oy_0","Ita_0"))
  temp.data$genotype<-as.character(temp.data$genotype)
  # for(i in 1:27) {
  #   conversion.table[i,]
  #   temp.data[temp.data$genotype==conversion.table[i,"num"],"genotype"]<-rep(as.character(conversion.table[i,"genotype"]),sum(as.integer(temp.data$genotype==conversion.table[i,"num"])))
  # }
  temp.data$genotype<-vlookup(temp.data$genotype,conversion.table,2)
  #
  temp.data.melt<-melt(temp.data,id=c("trt","time","genotype"))
  print(temp.data.melt)
  #temp.data.samples.melt$time<-factor(temp.data.samples.melt$time,levels=c("1","4","16","25","49"))
  temp.data.melt$time<-factor(temp.data.melt$time,levels=c("1","4","25","49"))
  temp.data.melt$variable<-factor(temp.data.melt$variable,levels=plot.order)
  # for individual library 
  temp.data.idv<-as.data.frame(t(data.expression[rownames(data.expression) %in% target.genes,]))
  temp.data.idv$genotype<-gsub("(C|D|E)([[:digit:]]+)(H|L)(1|4|16|25|49)(A.merged.bam)","\\2",rownames(temp.data.idv))
  temp.data.idv$trt<-gsub("(C|D|E)([[:digit:]]+)(H|L)(1|4|16|25|49)(A.merged.bam)","\\3",rownames(temp.data.idv))
  temp.data.idv$time<-gsub("(C|D|E)([[:digit:]]+)(H|L)(1|4|16|25|49)(A.merged.bam)","\\4",rownames(temp.data.idv))
  temp.data.idv$rep<-gsub("(C|D|E)([[:digit:]]+)(H|L)(1|4|16|25|49)(A.merged.bam)","\\1",rownames(temp.data.idv))
  temp.data.idv$genotype<-vlookup(temp.data.idv$genotype,conversion.table,2)
  temp.data.idv.melt<-melt(temp.data.idv,id=c("trt","time","genotype","rep"))
  
  q<-ggplot(temp.data.melt) + geom_bar(aes(x=time,y=value,color=trt), stat="identity",position=position_dodge(),alpha = 0.5) + geom_jitter(data=temp.data.idv.melt,aes(x=time,y=value,color=trt),position=position_jitterdodge(jitter.width=0.3))
  q<-q + facet_grid(variable~genotype) + theme(strip.text.y=element_text(angle=0))
  #  q<-q + scale_y_continuous(trans=log2_trans())
  return(q)
}




expression.pattern.logFC.graph4<-function(data.cpm,target.genes,logFC){# require ggplot2, reshape2 packages, this is only for SAS timecourse data (mean type, eg. "summary4")
  temp.data<-as.data.frame(t(data.cpm[rownames(data.cpm) %in% target.genes,]))
  print(temp.data)
  temp.data$genotype<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\1",rownames(temp.data))
  temp.data$trt<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\2",rownames(temp.data))
  temp.data$time<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\3",rownames(temp.data))
  
  temp.data.melt<-melt(temp.data[,-1],id=c("trt","time","genotype"))
  temp.data.melt$sample<-paste(temp.data.melt$trt,temp.data.melt$time,temp.data.melt$genotype,temp.data.melt$variable,sep=".")
  names(temp.data.melt)[5]<-"expression"
  # combine expression value and logFC (Col vs mutant) data
  temp.logFC<-as.data.frame(t(logFC[rownames(logFC) %in% target.genes,]))
  temp.logFC$genotype<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\1",rownames(temp.logFC))
  temp.logFC$trt<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\2",rownames(temp.logFC))
  temp.logFC$time<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\3",rownames(temp.logFC))
  
  temp.logFC.melt<-melt(temp.logFC[,-1],id=c("trt","time","genotype"))
  temp.logFC.melt$sample<-paste(temp.logFC.melt$trt,temp.logFC.melt$time,temp.logFC.melt$genotype,temp.logFC.melt$variable,sep=".")
  names(temp.logFC.melt)[5]<-"logFC"
  
  # combine
  temp.combined<-merge(temp.data.melt,temp.logFC.melt[,c("logFC","sample")],by="sample",all=TRUE)
  names(temp.combined)<-gsub("variable","AGI",names(temp.combined))
  #temp.combined[is.na(temp.combined)]
  temp.combined$time<-factor(temp.combined$time,levels=c("1","4","25","49"))
  #  temp.combined$logFC2<-temp.combined$logFC +7  
  temp.combined[is.na(temp.combined)]<-0
  
  
  #    q<-ggplot(temp.combined,aes(x=time,y=expression)) + geom_text(label="*",aes(colour=logFC) )  + scale_color_gradient2(low="green",high="magenta",mid="white",midpoint=0)
  # q<-q + geom_point(alpha = 0.5,aes(shape=factor(trt)))
  # q<-q + facet_grid(AGI~genotype,scale="free") + theme(strip.text.y=element_text(angle=0)) 
  #   #q
  #   ggsave(q,file="test.pdf",width=11,height=8)
  temp.combined[temp.combined$logFC>3,"logFC"]<-3
  temp.combined[temp.combined$logFC< -3,"logFC"]<- -3
  
  q<-ggplot(temp.combined,aes(x=time,y=expression))  + geom_point(alpha = 0.5,aes(shape=factor(trt),colour=logFC) ) + scale_color_gradient2(limit=c(-3,3),low="green",high="magenta",mid="black",midpoint=0)
  q<-q + facet_grid(AGI~genotype,scale="free") + theme(strip.text.y=element_text(angle=0)) 
  #q
  #ggsave(q,file="test2.pdf",width=11,height=8)
  #   #  q<-q + scale_y_continuous(trans=log2_trans())
  return(q)
}


expression.mean<-function(data) { # return rowMean results for given data
  # also needs samples.nolow2$file (use samples.nolow3$file for wo16h data)
  x<-samples.nolow5[samples.nolow5$file %in% names(data),] # s change this line only
  
  print(x)
  str(x)
  x$sample<-as.character(x$sample)
  x$sample<-factor(x$sample)
  str(x)
  filename<-as.character(x[x$sample %in% levels(x$sample)[1],"file"])  
  print(paste("filename is",filename))						
  data.s<-data[,names(data) %in% filename]	
  print(head(data.s))
  if(class(data.s)=="data.frame") {
    data.rowMeans<-data.frame(rowMeans(data.s))
    names(data.rowMeans)<-levels(x$sample)[1]
  }
  else {
    data.rowMeans<-data.frame(data.s)
    names(data.rowMeans)<-levels(x$sample)[1]
  }		
  # more samples
  for(sample in levels(x$sample)[-1]) {
    print("sample is")
    print(sample)
    filename<-as.character(x[x$sample %in% sample,"file"])
    print(filename)
    data.s<-data[,names(data) %in% filename]	
    print(head(data.s))	
    if(class(data.s)=="data.frame") {
      data.s.rowMeans<-data.frame(rowMeans(data.s))
      names(data.s.rowMeans)<-sample
      data.rowMeans<-cbind(data.rowMeans,data.s.rowMeans)
    }
    else {
      data.s.rowMeans<-data.frame(data.s)
      names(data.s.rowMeans)<-sample
      data.rowMeans<-cbind(data.rowMeans,data.s.rowMeans)
    }
  }	
  return(data.rowMeans)
}
# 

expression.pattern.graph5<-function(data.cpm,target.genes,plot.order=NULL){# require ggplot2, reshape2 packages, this is only for SAS timecourse data (mean type, eg. "summary4")
  temp.data<-as.data.frame(t(data.cpm[rownames(data.cpm) %in% target.genes,]))
  print(temp.data)
  # convert into relative value Col, sun 1h is 1 (new to this function)
  temp.data<-as.data.frame(t(t(temp.data)/(t(temp.data)[,1])))
  temp.data<-as.data.frame(t(log2(t(temp.data)/(t(temp.data)[,1]))))
  
  #
  temp.data$genotype<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\1",rownames(temp.data))
    temp.data$trt<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\2",rownames(temp.data))
   temp.data$time<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\3",rownames(temp.data))
    temp.data$trt.time<-with(temp.data,paste(trt,time,sep=""))
    #temp.data$trt.time<-factor(rownames(temp.data),levels=c("1H1","1H4","1H25","1H49","1L1","1L4","1L25","1L49")) # only works for Col
  temp.data$trt.time<-factor(temp.data$trt.time,levels=c("H1","H4","H25","H49","L1","L4","L25","L49")) # only works for Col
  
  # convert genotype number into actual name
  conversion.table<-data.frame(num=1:27, genotype=c("Col","AT5G02540_1","hy5","jar1","kat1_2","phyB","pif45","spt_11","yuc2589","PAR1_RNAi09","coi1_16","phyAB","pif3","mida9_4","bsk5_1","sto","aos","argos","co_9","Blh_1","Jea","Shahdara","Col_0","Cvi_0","Bur_0","Oy_0","Ita_0"))
  temp.data$genotype<-as.character(temp.data$genotype)
  for(i in 1:27) {
    conversion.table[i,]
    temp.data[temp.data$genotype==conversion.table[i,"num"],"genotype"]<-rep(as.character(conversion.table[i,"genotype"]),sum(as.integer(temp.data$genotype==conversion.table[i,"num"])))
  }
  #temp.data.melt<-melt(temp.data,id=c("trt","time","genotype"))
  temp.data.melt<-melt(subset(temp.data,select=-time),id=c("trt.time","genotype","trt"))

  str(temp.data.melt)
  temp.data.melt$value<-as.numeric(temp.data.melt$value) 
  
  print(temp.data.melt)
  summary(temp.data.melt)
  #temp.data.samples.melt$time<-factor(temp.data.samples.melt$time,levels=c("1","4","16","25","49"))
  #temp.data.melt$time<-factor(temp.data.melt$time,levels=c("1","4","25","49"))
  if(is.null(plot.order)==TRUE) {temp.data.melt<-temp.data.melt} else {temp.data.melt$variable<-factor(temp.data.melt$variable,levels=plot.order)}
  # graph
  q<-ggplot(temp.data.melt,aes(x=trt.time,y=value)) + geom_violin(aes(fill=trt))
  q<-q + facet_grid(.~genotype,scale="free") + theme(strip.text.y=element_text(angle=0)) + labs(y="log2 normalized value") 
  #  q<-q + scale_y_continuous(trans=log2_trans())
  return(q)
}
# line version
expression.pattern.graph5b<-function(data.cpm,target.genes,plot.order=NULL){# require ggplot2, reshape2 packages, this is only for SAS timecourse data (mean type, eg. "summary4")
  temp.data<-as.data.frame(t(data.cpm[rownames(data.cpm) %in% target.genes,]))
  print(temp.data)
  # convert into relative value Col, sun 1h is 1 (new to this function)
  temp.data<-as.data.frame(t(t(temp.data)/(t(temp.data)[,1])))
  temp.data<-as.data.frame(t(log2(t(temp.data)/(t(temp.data)[,1]))))
  
  #
  temp.data$genotype<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\1",rownames(temp.data))
  temp.data$trt<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\2",rownames(temp.data))
  temp.data$time<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\3",rownames(temp.data))
  temp.data$trt.time<-with(temp.data,paste(trt,time,sep=""))
  #temp.data$trt.time<-factor(rownames(temp.data),levels=c("1H1","1H4","1H25","1H49","1L1","1L4","1L25","1L49")) # only works for Col
  temp.data$trt.time<-factor(temp.data$trt.time,levels=c("H1","H4","H25","H49","L1","L4","L25","L49")) # only works for Col
  
  # convert genotype number into actual name
  conversion.table<-data.frame(num=1:27, genotype=c("Col","AT5G02540_1","hy5","jar1","kat1_2","phyB","pif45","spt_11","yuc2589","PAR1_RNAi09","coi1_16","phyAB","pif3","mida9_4","bsk5_1","sto","aos","argos","co_9","Blh_1","Jea","Shahdara","Col_0","Cvi_0","Bur_0","Oy_0","Ita_0"))
  temp.data$genotype<-as.character(temp.data$genotype)
  for(i in 1:27) {
    conversion.table[i,]
    temp.data[temp.data$genotype==conversion.table[i,"num"],"genotype"]<-rep(as.character(conversion.table[i,"genotype"]),sum(as.integer(temp.data$genotype==conversion.table[i,"num"])))
  }
  temp.data.melt<-melt(temp.data,id=c("trt","time","genotype","trt.time"))
  #temp.data.melt<-melt(subset(temp.data,select=-time),id=c("trt.time","genotype","trt"))
  
  str(temp.data.melt)
  temp.data.melt$value<-as.numeric(temp.data.melt$value) 
  
  print(temp.data.melt)
  summary(temp.data.melt)
  #temp.data.samples.melt$time<-factor(temp.data.samples.melt$time,levels=c("1","4","16","25","49"))
  #temp.data.melt$time<-factor(temp.data.melt$time,levels=c("1","4","25","49"))
  temp.data.melt$time<-as.numeric(temp.data.melt$time)
  if(is.null(plot.order)==TRUE) {temp.data.melt<-temp.data.melt} else {temp.data.melt$variable<-factor(temp.data.melt$variable,levels=plot.order)}
  # graph
  q<-ggplot(temp.data.melt,aes(x=time,y=value,color=trt,group=variable)) + geom_line(alpha=0.1)
  q<-q + facet_grid(trt~genotype,scale="free") + theme(strip.text.y=element_text(angle=0)) + labs(y="log2 normalized value") 
  #  q<-q + scale_y_continuous(trans=log2_trans())
  return(q)
}
# absolute value
expression.pattern.graph5.abs<-function(data.cpm,target.genes,plot.order=NULL){# require ggplot2, reshape2 packages, this is only for SAS timecourse data (mean type, eg. "summary4")
  temp.data<-as.data.frame(t(data.cpm[rownames(data.cpm) %in% target.genes,]))
  print(temp.data)
  # convert into relative value Col, sun 1h is 1 (new to this function)
  #temp.data<-as.data.frame(t(t(temp.data)/(t(temp.data)[,1])))
  #temp.data<-as.data.frame(t(log2(t(temp.data)/(t(temp.data)[,1]))))
  
  #
  temp.data$genotype<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\1",rownames(temp.data))
  temp.data$trt<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\2",rownames(temp.data))
  temp.data$time<-gsub("([[:digit:]]+)(H|L)(1|4|16|25|49)","\\3",rownames(temp.data))
  temp.data$trt.time<-with(temp.data,paste(trt,time,sep=""))
  #temp.data$trt.time<-factor(rownames(temp.data),levels=c("1H1","1H4","1H25","1H49","1L1","1L4","1L25","1L49")) # only works for Col
  temp.data$trt.time<-factor(temp.data$trt.time,levels=c("H1","H4","H25","H49","L1","L4","L25","L49")) # only works for Col
  
  # convert genotype number into actual name
  conversion.table<-data.frame(num=1:27, genotype=c("Col","AT5G02540_1","hy5","jar1","kat1_2","phyB","pif45","spt_11","yuc2589","PAR1_RNAi09","coi1_16","phyAB","pif3","mida9_4","bsk5_1","sto","aos","argos","co_9","Blh_1","Jea","Shahdara","Col_0","Cvi_0","Bur_0","Oy_0","Ita_0"))
  temp.data$genotype<-as.character(temp.data$genotype)
  for(i in 1:27) {
    conversion.table[i,]
    temp.data[temp.data$genotype==conversion.table[i,"num"],"genotype"]<-rep(as.character(conversion.table[i,"genotype"]),sum(as.integer(temp.data$genotype==conversion.table[i,"num"])))
  }
  temp.data.melt<-melt(temp.data,id=c("trt","time","genotype","trt.time"))
  #temp.data.melt<-melt(subset(temp.data,select=-time),id=c("trt.time","genotype","trt"))
  
  str(temp.data.melt)
  #temp.data.melt$value<-as.numeric(temp.data.melt$value) 
  
  print(temp.data.melt)
  summary(temp.data.melt)
  #temp.data.samples.melt$time<-factor(temp.data.samples.melt$time,levels=c("1","4","16","25","49"))
  temp.data.melt$time<-factor(temp.data.melt$time,levels=c("1","4","25","49"))
  #temp.data.melt$time<-as.numeric(temp.data.melt$time)
  if(is.null(plot.order)==TRUE) {temp.data.melt<-temp.data.melt} else {temp.data.melt$variable<-factor(temp.data.melt$variable,levels=plot.order)}
  # graph
  q<-ggplot(temp.data.melt,aes(x=time,y=value,color=variable,group=variable)) + geom_line(alpha=0.1)
  q<-q + facet_grid(genotype~trt,scale="free") + theme(strip.text.y=element_text(angle=0)) + labs(y="value") 
  #  q<-q + scale_y_continuous(trans=log2_trans())
  return(q)
}
# fold chagne (shade responsiveness using summary.vst.response.kazu
expression.pattern.graph5.FC<-function(data.cpm,target.genes,plot.order=NULL){# require ggplot2, reshape2 packages, this is only for SAS timecourse data (mean type, eg. "summary4")
  temp.data<-as.data.frame(t(data.cpm[rownames(data.cpm) %in% target.genes,]))

  temp.data$genotype<-gsub("([[:digit:]]+)(_)(1|4|25|49)(hrA)","\\1",rownames(temp.data))
  temp.data$time<-gsub("([[:digit:]]+)(_)(1|4|25|49)(hrA)","\\3",rownames(temp.data)) 
    # convert genotype number into actual name
  conversion.table<-data.frame(num=1:27, genotype=c("Col","AT5G02540_1","hy5","jar1","kat1_2","phyB","pif45","spt_11","yuc2589","PAR1_RNAi09","coi1_16","phyAB","pif3","mida9_4","bsk5_1","sto","aos","argos","co_9","Blh_1","Jea","Shahdara","Col_0","Cvi_0","Bur_0","Oy_0","Ita_0"))
  temp.data$genotype<-as.character(temp.data$genotype)
  for(i in 1:27) {
    conversion.table[i,]
    temp.data[temp.data$genotype==conversion.table[i,"num"],"genotype"]<-rep(as.character(conversion.table[i,"genotype"]),sum(as.integer(temp.data$genotype==conversion.table[i,"num"])))
  }
  #
  print(temp.data)

    temp.data.melt<-melt(temp.data,id=c("time","genotype"))

  str(temp.data.melt)
  #temp.data.melt$value<-as.numeric(temp.data.melt$value) 
  
  print(temp.data.melt)
  summary(temp.data.melt)
  #temp.data.samples.melt$time<-factor(temp.data.samples.melt$time,levels=c("1","4","16","25","49"))
  temp.data.melt$time<-factor(temp.data.melt$time,levels=c("1","4","25","49"))
  #temp.data.melt$time<-as.numeric(temp.data.melt$time)
  if(is.null(plot.order)==TRUE) {temp.data.melt<-temp.data.melt} else {temp.data.melt$variable<-factor(temp.data.melt$variable,levels=plot.order)}
  # graph
  q<-ggplot(temp.data.melt,aes(x=time,y=value,color=variable,group=variable)) + geom_line()
  q<-q + facet_grid(.~genotype,scale="free") + theme(strip.text.y=element_text(angle=0)) + labs(y="value") 
  #  q<-q + scale_y_continuous(trans=log2_trans())
  return(q)
}

# ORA with GOseq
# prerequisit
#library(ShortRead);library(goseq);library(GO.db);library("org.At.tair.db");library("annotate")

#TIR10_cdna_rep_model<-readDNAStringSet("../../Nozue2016_SAStranscriptome_data/input/TAIR10_cdna_20110103_representative_gene_model") 
#setwd("../")
TIR10_cdna_rep_model<-readDNAStringSet(file.path(homedir,"..","..","Nozue2016_SAStranscriptome_data/input/TAIR10_cdna_20110103_representative_gene_model") )
#setwd(homedir2)
head(TIR10_cdna_rep_model)
bias<-nchar(TIR10_cdna_rep_model)
names(bias)<-substr(names(TIR10_cdna_rep_model),1,9)
length(bias) 
#  bias.data vector must have the same length as DEgenes vector!
###Read in AtGO
Atgo <- toTable(org.At.tairGO)
#head(Atgo)
BP <- TRUE #only keep BP go TERMS
if (BP) Atgo <- Atgo[Atgo$Ontology=="BP",]
#convert to list
Atgo.list <- tapply(Atgo$go_id,Atgo$gene_id,c)

GOseq.ORA<-function(genelist,padjust=0.05) { # return GO enrichment table, padjus, padjust=0.05 , modified 092817
  TF<-(names(bias) %in% genelist)*1
  names(TF)<-names(bias)
  #print(TF)
  pwf<-nullp(TF,bias.data=bias)
  #print(pwf$DEgenes)
  ###Read in AtGO
  Atgo <- toTable(org.At.tairGO)
  #head(Atgo)
  BP <- TRUE #only keep BP go TERMS
  if (BP) Atgo <- Atgo[Atgo$Ontology=="BP",]
  #convert to list
  Atgo.list <- tapply(Atgo$go_id,Atgo$gene_id,c)
  #
  GO.pval <- goseq(pwf,gene2cat=Atgo.list,use_genes_without_cat=TRUE) # format became different in new goseq version (021111)
  #head(GO.pval) 
  GO.pval$over_represented_padjust<-p.adjust(GO.pval$over_represented_pvalue,method="BH")
  #if(GO.pval$over_represented_padjust[1]>padjust) stop("no enriched GO")
  if(GO.pval$over_represented_padjust[1]>padjust) return("no enriched GO")
  
  else {
    enriched.GO<-GO.pval[GO.pval$over_represented_padjust<padjust,] 
    print("enriched.GO is")
    print(enriched.GO)
    
    ## write Term and Definition 
    for(i in 1:dim(enriched.GO)[1]) {
      enriched.GO$Term[i]<-Term(GOTERM[[enriched.GO[i,"category"]]])
      enriched.GO$Definition[i]<-Definition(GOTERM[[enriched.GO[i,"category"]]])
    }
    return(enriched.GO)
  }
}
# example
#enriched.GO.Col.SOMcluster1<-GOseq.ORA(rownames(data.val3.4.SOM.Col.SOM1.all.barcode.s)[rownames(data.val3.4.SOM.Col.SOM1.all.barcode.s) %in% names(bias)]) 
GOseq.ORA_old<-function(genelist,padjust=0.05) { # return GO enrichment table, padjus, padjust=0.05 
  TF<-(names(bias) %in% genelist)*1
  names(TF)<-names(bias)
  #print(TF)
  pwf<-nullp(TF,bias.data=bias)
  #print(pwf$DEgenes)
  ###Read in AtGO
  Atgo <- toTable(org.At.tairGO)
  #head(Atgo)
  BP <- TRUE #only keep BP go TERMS
  if (BP) Atgo <- Atgo[Atgo$Ontology=="BP",]
  #convert to list
  Atgo.list <- tapply(Atgo$go_id,Atgo$gene_id,c)
  #
  GO.pval <- goseq(pwf,gene2cat=Atgo.list,use_genes_without_cat=TRUE) # format became different in new goseq version (021111)
  #head(GO.pval) 
  GO.pval$over_represented_padjust<-p.adjust(GO.pval$over_represented_pvalue,method="BH")
  #if(GO.pval$over_represented_padjust[1]>padjust) stop("no enriched GO")
  if(GO.pval$over_represented_padjust[1]>padjust) return("no enriched GO")
  
  else {
    enriched.GO<-GO.pval[GO.pval$over_represented_padjust<padjust,] 
    print("enriched.GO is")
    print(enriched.GO)
    
    ## write Term and Definition 
    for(i in 1:dim(enriched.GO)[1]) {
      enriched.GO$Term[i]<-Term(GOTERM[[enriched.GO[i,"category"]]])
      enriched.GO$Definition[i]<-Definition(GOTERM[[enriched.GO[i,"category"]]])
    }
    return(enriched.GO)
  }
}



# revised GOseq.CC.ORA (112816)
GOseq.CC.ORA<-function(genelist,padjust=0.05) { # return GO enrichment table, padjus, padjust=0.05 
  TF<-(names(bias) %in% genelist)*1
  names(TF)<-names(bias)
  #print(TF)
  pwf<-nullp(TF,bias.data=bias)
  #print(pwf$DEgenes)
  ###Read in AtGO
  Atgo <- toTable(org.At.tairGO)
  #head(Atgo)
  CC <- TRUE #only keep BP go TERMS
  if (CC) Atgo <- Atgo[Atgo$Ontology=="CC",]
  #convert to list
  Atgo.list <- tapply(Atgo$go_id,Atgo$gene_id,c)
  #
  #GO.pval <- goseq(pwf,gene2cat=Atgo.list) # format became different in new goseq version (021111)
  GO.pval <- goseq(pwf,gene2cat=Atgo.list,use_genes_without_cat=TRUE) # format became different in new goseq version (021111)
  
  #head(GO.pval) 
  GO.pval$over_represented_padjust<-p.adjust(GO.pval$over_represented_pvalue,method="BH")
  if(GO.pval$over_represented_padjust[1]>padjust) return("no enriched GO")
  
  else {
    enriched.GO<-GO.pval[GO.pval$over_represented_padjust<padjust,] 
    print("enriched.GO is")
    print(enriched.GO)
    
    ## write Term and Definition 
    for(i in 1:dim(enriched.GO)[1]) {
      enriched.GO$Term[i]<-Term(GOTERM[[enriched.GO[i,"category"]]])
      enriched.GO$Definition[i]<-Definition(GOTERM[[enriched.GO[i,"category"]]])
    }
    return(enriched.GO)
  }  
}

genes.in.enriched.category<-function(enrich.result,gene.list,category.table=hormone.responsiveness) { # return data frame with all input data with ... which format is good? VennDiagram input format?  
  enrich.category<-enrich.result[enrich.result$over_represented_pvalue<0.05,"category"]
  temp<-category.table[gene.list] # short category.table for enriched.result # gene.list shas to be character (not factor)  
  #names(temp[!is.na(names(temp))])
  temp<-temp[!is.na(names(temp))] # remove non matched AGI
  test<-data.frame(AGI=names(temp))
  for(i in enrich.category) {
    x<-rep(0,length(names(temp)))
    x[grep(i,temp)]<-1
    test<-cbind(test,x)
    names(test)[dim(test)[2]]<-i
  }
  return(test) 
}

my.category.heatmap5<-function(summary3.response,gt.num,gt2,newdata.Col.wGO.selected,plotname,h=5,w=5*4/3,path,save.plot=T,legend.fill=F) { # gt.num (num), gt2 (gt name), for w/o 16h data, newdata.Col.wGO.selected should have "AGI" and "my.category"
  # voom transformed data has been transformd into log2
  ### for temp 
  # calculate newdata.SOM.temp (needs to fix becasue changing into gt = 1 did not reproduce above)
  #newdata.SOM.temp<-newdata[newdata$gt==9,] # not all genes were selected. Start from summary3.response
  summary3.response.log2.temp<-log2(summary3.response[,gsub("([[:digit:]]+)(_)(1|4|16|25|49)hrA","\\1",names(summary3.response))==gt.num]) # 
  print(head(summary3.response.log2.temp))
  summary3.response.log2.temp$AGI<-rownames(summary3.response.log2.temp)  
  print(head(summary3.response.log2.temp$AGI))
  #### select genes in newdata.Col.wGO.selected and look expression pattern in temp
  print(newdata.Col.wGO.selected)
  selected.genes<-newdata.Col.wGO.selected$AGI #
  print(paste("selected.genes are",selected.genes))
  summary3.response.log2.temp.selected<-summary3.response.log2.temp[summary3.response.log2.temp$AGI %in% selected.genes,]
  # not that way. use Col data
  summary3.response.log2.temp.selected<-merge(summary3.response.log2.temp.selected,newdata.Col.wGO.selected[,c("AGI","my.category")],by="AGI")
  #
  names(summary3.response.log2.temp.selected)[2:5]<-as.numeric(gsub("([[:digit:]]+)(_)(1|4|16|25|49)hrA","\\3",names(summary3.response.log2.temp.selected)[2:5]))
  # remove "-Inf" gene
  #summary3.response.log2.temp.selected.noINF<-summary3.response.log2.temp.selected[apply(summary3.response.log2.temp.selected[,2:5], 1, Compose(is.finite, all)),]
  # Error in Compose(is.finite, all) : Argument is not a function
  # see http://stackoverflow.com/questions/15773189/remove-na-nan-inf-in-a-matrix for repair
  summary3.response.log2.temp.selected.noINF<-summary3.response.log2.temp.selected[apply(summary3.response.log2.temp.selected[,2:5], 1, function(x) all(is.finite(x))),]
  
  print(summary3.response.log2.temp.selected.noINF)
  # melt
  summary3.response.log2.temp.selected.noINF.melt<-melt(summary3.response.log2.temp.selected.noINF,id=c("AGI","my.category"))
  print("summary3.response.log2.temp.selected.noINF.melt is")
  print(head(summary3.response.log2.temp.selected.noINF.melt)) # OK
  summary(summary3.response.log2.temp.selected.noINF.melt)
  summary3.response.log2.temp.selected.noINF.melt$variable<-factor(summary3.response.log2.temp.selected.noINF.melt$variable,levels=c("1","4","25","49"))
  
  summary3.response.log2.temp.selected.noINF.melt$my.category<-factor(summary3.response.log2.temp.selected.noINF.melt$my.category,levels=levels(newdata.Col.wGO.selected$my.category))
  print(head(summary3.response.log2.temp.selected.noINF.melt)) # OK? no.
  
  # mean table (= heat map should reflect this table)
  summary.table<-as.data.frame(tapply(summary3.response.log2.temp.selected.noINF.melt$value,list(summary3.response.log2.temp.selected.noINF.melt$my.category,summary3.response.log2.temp.selected.noINF.melt$variable),mean))
  print(paste("summary table for ",gt2))
  print(summary.table)
  summary.table$my.category<-factor(rownames(summary.table),levels=levels(summary3.response.log2.temp.selected.noINF.melt$my.category))
  
  summary.table.melt<-melt(summary.table,id="my.category")
  #names(summary.table)
  #summary.table
  # additional part in this function
  summary.table.melt$variable<-factor(summary.table.melt$variable,levels=c("49","25","4","1"))
  summary.table.melt$my.category<-factor(summary.table.melt$my.category,levels=rev(levels(summary.table.melt$my.category)))
  
  save(summary.table.melt,file=file.path("..","..","Nozue2016_SAStranscriptome_output","output",paste("FigS2.absolute.summary.table.melt",gt.num,".Rdata",sep="")))
  # drawing heatmap
  #p.selected <- ggplot() + geom_tile(sumamry.table.melt,aes(x=variable,y=my.category,fill=value),size=0.3,colour="black") # does not work (052816)
  p.selected <- ggplot(summary.table.melt,aes(x=variable,y=my.category)) + geom_tile(size=0.3,colour="black",aes(fill=value)) 
  
  library(scales) # for muted
  #p.selected <- p.selected + scale_fill_gradient2(limits=c(-1,1),low=muted("green"), high=muted("magenta")) 
  p.selected <- p.selected + scale_fill_gradient2(limits=c(-1,1),low=muted("green"), high=muted("magenta")) #,guide=guide_legend("none")) 
  p.selected <- p.selected  + 
    theme(axis.text.x=element_text(size=30,angle=90),
          axis.text.y=element_text(size=30),
          axis.title=element_text(size=40),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.title=element_text(size=40),
          axis.line=element_blank())
  if(legend.fill==TRUE) {
    #p.selected <- p.selected  + labs(fill="log2\n fold change",legend.text=element_text(size=20)) 
    p.selected <- p.selected  + labs(fill="") + theme(legend.title=element_text(size=40),legend.text=element_text(size=40),legend.key.height=unit(1, "cm")) 
    
  } else {p.selected <- p.selected + theme(legend.position="none")}
  #,strip.background=element_rect(fill=c("red","blue")))
  #p.selected <- p.selected + labs(x="Time (hr)",y="",title=gt2,fill="log2\n fold change")
  p.selected <- p.selected + labs(x="Time (hr)",y="",title=gt2)
  
  # flip xy
  p.selected<-p.selected + coord_flip() 
  if(save.plot==T)   {
    ggsave(file=paste(gt2,plotname,sep="."),p.selected,height=h,width=w,path=path)
  } 
  else {return(p.selected)}
  #### the end of temp
}
# tsne version
my.category.heatmap6<-function(summary3.response,gt.num,gt2,newdata.Col.wGO.selected,plotname,h=5,w=5*4/3,path,save.plot=T,legend.fill=F,file.name="",category.level="") { # gt.num (num), gt2 (gt name), for w/o 16h data, newdata.Col.wGO.selected should have "AGI" and "my.category"
  # voom transformed data has been transformd into log2
  ### for temp 
  # calculate newdata.SOM.temp (needs to fix becasue changing into gt = 1 did not reproduce above)
  #newdata.SOM.temp<-newdata[newdata$gt==9,] # not all genes were selected. Start from summary3.response
  if(exists("category.level")) {category.level<-category.level} else { category.level<-levels(newdata.Col.wGO.selected$my.category)} # 102617
  
  summary3.response.log2.temp<-log2(summary3.response[,gsub("([[:digit:]]+)(_)(1|4|16|25|49)hrA","\\1",names(summary3.response))==gt.num]) # 
  print(head(summary3.response.log2.temp))
  summary3.response.log2.temp$AGI<-rownames(summary3.response.log2.temp)  
  print(head(summary3.response.log2.temp$AGI))
  #### select genes in newdata.Col.wGO.selected and look expression pattern in temp
  print(newdata.Col.wGO.selected)
  selected.genes<-newdata.Col.wGO.selected$AGI #
  print(paste("selected.genes are",selected.genes))
  summary3.response.log2.temp.selected<-summary3.response.log2.temp[summary3.response.log2.temp$AGI %in% selected.genes,]
  # not that way. use Col data
  summary3.response.log2.temp.selected<-merge(summary3.response.log2.temp.selected,newdata.Col.wGO.selected[,c("AGI","my.category")],by="AGI")
  #
  names(summary3.response.log2.temp.selected)[2:5]<-as.numeric(gsub("([[:digit:]]+)(_)(1|4|16|25|49)hrA","\\3",names(summary3.response.log2.temp.selected)[2:5]))
  # remove "-Inf" gene
  #summary3.response.log2.temp.selected.noINF<-summary3.response.log2.temp.selected[apply(summary3.response.log2.temp.selected[,2:5], 1, Compose(is.finite, all)),]
  # Error in Compose(is.finite, all) : Argument is not a function
  # see http://stackoverflow.com/questions/15773189/remove-na-nan-inf-in-a-matrix for repair
  summary3.response.log2.temp.selected.noINF<-summary3.response.log2.temp.selected[apply(summary3.response.log2.temp.selected[,2:5], 1, function(x) all(is.finite(x))),]
  
  print(summary3.response.log2.temp.selected.noINF)
  # melt
  summary3.response.log2.temp.selected.noINF.melt<-melt(summary3.response.log2.temp.selected.noINF,id=c("AGI","my.category"))
  print("summary3.response.log2.temp.selected.noINF.melt is")
  print(head(summary3.response.log2.temp.selected.noINF.melt)) # OK
  summary(summary3.response.log2.temp.selected.noINF.melt)
  summary3.response.log2.temp.selected.noINF.melt$variable<-factor(summary3.response.log2.temp.selected.noINF.melt$variable,levels=c("1","4","25","49"))
  
  summary3.response.log2.temp.selected.noINF.melt$my.category<-factor(summary3.response.log2.temp.selected.noINF.melt$my.category,levels=levels(newdata.Col.wGO.selected$my.category))
  print(head(summary3.response.log2.temp.selected.noINF.melt)) # OK? no.
  
  # mean table (= heat map should reflect this table)
  summary.table<-as.data.frame(tapply(summary3.response.log2.temp.selected.noINF.melt$value,list(summary3.response.log2.temp.selected.noINF.melt$my.category,summary3.response.log2.temp.selected.noINF.melt$variable),mean))
  print(paste("summary table for ",gt2))
  print(summary.table)
  summary.table$my.category<-factor(rownames(summary.table),levels=levels(summary3.response.log2.temp.selected.noINF.melt$my.category))
  
  summary.table.melt<-melt(summary.table,id="my.category")
  #names(summary.table)
  #summary.table
  # additional part in this function
  summary.table.melt$variable<-factor(summary.table.melt$variable,levels=c("49","25","4","1"))
  #summary.table.melt$my.category<-factor(summary.table.melt$my.category,levels=rev(category.level))
  summary.table.melt$my.category<-factor(summary.table.melt$my.category,levels=category.level)
  
  save(summary.table.melt,file=file.path(path,paste("FigS2.absolute.summary.table.melt",gt.num,file.name,".Rdata",sep=""))) # 102617
  # drawing heatmap
  #p.selected <- ggplot() + geom_tile(sumamry.table.melt,aes(x=variable,y=my.category,fill=value),size=0.3,colour="black") # does not work (052816)
  p.selected <- ggplot(summary.table.melt,aes(x=variable,y=my.category)) + geom_tile(size=0.3,colour="black",aes(fill=value)) 
  
  library(scales) # for muted
  #p.selected <- p.selected + scale_fill_gradient2(limits=c(-1,1),low=muted("green"), high=muted("magenta")) 
  p.selected <- p.selected + scale_fill_gradient2(limits=c(-1,1),low=muted("green"), high=muted("magenta")) #,guide=guide_legend("none")) 
  p.selected <- p.selected  + 
    theme(axis.text.x=element_text(size=30,angle=90,hjust=0.95,vjust=0.2),
          axis.text.y=element_text(size=30),
          axis.title=element_text(size=40),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.title=element_text(size=40),
          axis.line=element_blank())
  if(legend.fill==TRUE) {
    #p.selected <- p.selected  + labs(fill="log2\n fold change",legend.text=element_text(size=20)) 
    p.selected <- p.selected  + labs(fill="") + theme(legend.title=element_text(size=40),legend.text=element_text(size=40),legend.key.height=unit(1, "cm")) 
    
  } else {p.selected <- p.selected + theme(legend.position="none")}
  #,strip.background=element_rect(fill=c("red","blue")))
  #p.selected <- p.selected + labs(x="Time (hr)",y="",title=gt2,fill="log2\n fold change")
  p.selected <- p.selected + labs(x="Time (hr)",y="",title=gt2)
  
  # flip xy
  p.selected<-p.selected + coord_flip() 
  if(save.plot==T)   {
    ggsave(file=paste(gt2,plotname,sep="."),p.selected,height=h,width=w,path=path)
  } 
  else {return(p.selected)}
  #### the end of temp
}

## diff version (061816)
## calculate diff at gene level and calculate mean value of the diff among each category
my.category.diff.heatmap1<-function(summary3.response,gt.num,gt2,newdata.Col.wGO.selected,plotname,h=5,w=5*4/3,path,save.plot=T,legend.fill=F,max.min.value=c(-1,1)) { # gt.num (num), gt2 (gt name), for w/o 16h data, newdata.Col.wGO.selected should have "AGI" and "my.category"
  # voom transformed data has been transformd into log2
  ### for temp 
  # calculate newdata.SOM.temp (needs to fix becasue changing into gt = 1 did not reproduce above)
  #newdata.SOM.temp<-newdata[newdata$gt==9,] # not all genes were selected. Start from summary3.response
  # Col data
  summary3.response.log2.Col<-log2(summary3.response[,gsub("([[:digit:]]+)(_)(1|4|16|25|49)hrA","\\1",names(summary3.response))=="1"]) # 
  print(head(summary3.response.log2.Col))  
  # gt2 data
  summary3.response.log2.mut<-log2(summary3.response[,gsub("([[:digit:]]+)(_)(1|4|16|25|49)hrA","\\1",names(summary3.response))==gt.num]) # 
  print(head(summary3.response.log2.mut))
  # differences (new)
  summary3.response.log2.temp<-summary3.response.log2.mut[,1:4]-summary3.response.log2.Col
  
  summary3.response.log2.temp[,1]<-(summary3.response.log2.temp[,1])*sign(summary3.response.log2.Col[,1])
  summary3.response.log2.temp[,2]<-(summary3.response.log2.temp[,2])*sign(summary3.response.log2.Col[,2])
  summary3.response.log2.temp[,3]<-(summary3.response.log2.temp[,3])*sign(summary3.response.log2.Col[,3])
  summary3.response.log2.temp[,4]<-(summary3.response.log2.temp[,4])*sign(summary3.response.log2.Col[,4])
  print(head(summary3.response.log2.temp))
  summary3.response.log2.temp$AGI<-rownames(summary3.response.log2.Col)    
  #### select genes in newdata.Col.wGO.selected and look expression pattern in temp
  print(newdata.Col.wGO.selected)
  selected.genes<-newdata.Col.wGO.selected$AGI
  print(paste("selected.genes are",selected.genes))
  summary3.response.log2.temp.selected<-summary3.response.log2.temp[summary3.response.log2.temp$AGI %in% selected.genes,]
  # not that way. use Col data
  summary3.response.log2.temp.selected<-merge(summary3.response.log2.temp.selected,newdata.Col.wGO.selected[,c("AGI","my.category")],by="AGI")
  #
  names(summary3.response.log2.temp.selected)[2:5]<-as.numeric(gsub("([[:digit:]]+)(_)(1|4|16|25|49)hrA","\\3",names(summary3.response.log2.temp.selected)[2:5]))
  # remove "-Inf" gene
  #summary3.response.log2.temp.selected.noINF<-summary3.response.log2.temp.selected[apply(summary3.response.log2.temp.selected[,2:5], 1, Compose(is.finite, all)),]
  summary3.response.log2.temp.selected.noINF<-summary3.response.log2.temp.selected[apply(summary3.response.log2.temp.selected[,2:5], 1, function(x) all(is.finite(x))),]
  print(summary3.response.log2.temp.selected.noINF)
  # melt
  summary3.response.log2.temp.selected.noINF.melt<-melt(summary3.response.log2.temp.selected.noINF,id=c("AGI","my.category"))
  print("summary3.response.log2.temp.selected.noINF.melt is")
  print(head(summary3.response.log2.temp.selected.noINF.melt)) # OK
  summary(summary3.response.log2.temp.selected.noINF.melt)
  summary3.response.log2.temp.selected.noINF.melt$variable<-factor(summary3.response.log2.temp.selected.noINF.melt$variable,levels=c("1","4","25","49"))
  
  summary3.response.log2.temp.selected.noINF.melt$my.category<-factor(summary3.response.log2.temp.selected.noINF.melt$my.category,levels=levels(newdata.Col.wGO.selected$my.category))
  print(head(summary3.response.log2.temp.selected.noINF.melt)) # OK? no.
  
  # mean table (= heat map should reflect this table)
  summary.table<-as.data.frame(tapply(summary3.response.log2.temp.selected.noINF.melt$value,list(summary3.response.log2.temp.selected.noINF.melt$my.category,summary3.response.log2.temp.selected.noINF.melt$variable),mean))
  print(paste("summary table for ",gt2))
  print(summary.table)
  summary.table$my.category<-factor(rownames(summary.table),levels=levels(summary3.response.log2.temp.selected.noINF.melt$my.category))
  
  summary.table.melt<-melt(summary.table,id="my.category")
  save(summary.table.melt,file=file.path("..","..","Nozue2016_SAStranscriptome_output","output",paste("Fig5B.diff.summary.table.melt",gt.num,".Rdata",sep="")))
  #summary.table # additional part in this function
  summary.table.melt$variable<-factor(summary.table.melt$variable,levels=c("49","25","4","1"))
  summary.table.melt$my.category<-factor(summary.table.melt$my.category,levels=rev(levels(summary.table.melt$my.category)))
  # plotting data (heatmap)
  #p.selected <- ggplot() + geom_tile(summary.table.melt,aes(x=variable,y=my.category,fill=value),size=0.3,colour="black") # does not work (052816)
  p.selected <- ggplot(summary.table.melt,aes(x=variable,y=my.category)) + geom_tile(size=0.3,colour="black",aes(fill=value)) 
  
  library(scales) # for muted
  #p.selected <- p.selected + scale_fill_gradient2(limits=c(-1,1),low=muted("green"), high=muted("magenta")) 
  p.selected <- p.selected + scale_fill_gradient2(limits=max.min.value,low=muted("green"), high=muted("magenta")) #,guide=guide_legend("none")) 
  p.selected <- p.selected  + 
    theme(axis.text.x=element_text(size=30,angle=90),
          axis.text.y=element_text(size=30),
          axis.title=element_text(size=40),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.title=element_text(size=40),
          axis.line=element_blank())
  if(legend.fill==TRUE) {
    p.selected <- p.selected  + labs(fill="Relative to Col",legend.text=element_text(size=20)) 
  } else {p.selected <- p.selected + theme(legend.position="none")}
  #,strip.background=element_rect(fill=c("red","blue")))
  p.selected <- p.selected + labs(x="Time (hr)",y="",title=gt2)
  
  # flip xy
  p.selected<-p.selected + coord_flip() 
  if(save.plot==T)   {
    ggsave(file=paste(gt2,plotname,sep="."),p.selected,height=h,width=w,path=path)
  } 
  else {return(p.selected)}
  #### the end of temp
} # the end of my.category.diff.heatmap1

# diff ver2 (101117)
my.category.diff.heatmap2<-function(summary3.response.diff,gt.num,gt2,newdata.Col.wGO.selected,plotname,h=5,w=5*4/3,path,save.plot=T,legend.fill=F,max.min.value=c(-1,1)) { # gt.num (num), gt2 (gt name), for w/o 16h data, newdata.Col.wGO.selected should have "AGI" and "my.category"
  # voom transformed data has been transformd into log2
  # melt  
  summary3.response.diff.melt<-melt(summary3.response.diff,id=c("time","my.category"))
  # select genotype by gt.num
  summary3.response.diff.melt.s<-summary3.response.diff.melt[summary3.response.diff.melt$variable==gt.num,]
  
  #summary3.response.diff.melt$my.category<-factor(rownames(summary3.response.diff.melt),levels=levels(summary3.response.diff.melt$my.category))
  summary3.response.diff.melt.s$time<-factor(summary3.response.diff.melt.s$time,levels=c("49","25","4","1"))
  # plotting data (heatmap)
  #p.selected <- ggplot() + geom_tile(summary.table.melt,aes(x=variable,y=my.category,fill=value),size=0.3,colour="black") # does not work (052816)
  p.selected <- ggplot(summary3.response.diff.melt.s,aes(x=time,y=my.category)) + geom_tile(size=0.3,colour="black",aes(fill=value)) 
  
  library(scales) # for muted
  #p.selected <- p.selected + scale_fill_gradient2(limits=c(-1,1),low=muted("green"), high=muted("magenta")) 
  p.selected <- p.selected + scale_fill_gradient2(limits=max.min.value,low=muted("green"), high=muted("magenta")) #,guide=guide_legend("none")) 
  p.selected <- p.selected  + 
    theme(axis.text.x=element_text(size=30,angle=90),
          axis.text.y=element_text(size=30),
          axis.title=element_text(size=40),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.title=element_text(size=40),
          axis.line=element_blank())
  if(legend.fill==TRUE) {
    p.selected <- p.selected  + labs(fill="Relative to Col",legend.text=element_text(size=20)) 
  } else {p.selected <- p.selected + theme(legend.position="none")}
  #,strip.background=element_rect(fill=c("red","blue")))
  p.selected <- p.selected + labs(x="Time (hr)",y="",title=gt2)
  
  # flip xy
  p.selected<-p.selected + coord_flip() 
  if(save.plot==T)   {
    ggsave(file=paste(gt2,plotname,sep="."),p.selected,height=h,width=w,path=path)
  } 
  else {return(p.selected)}
  #### the end of temp
} # the end of my.category.diff.heatmap2


# overlapTable GOseq version (from 082515)
#TIR10_cdna_rep_model<-readDNAStringSet("/Volumes/Data8/NGS_related/Arabidopsis_analysis/reference/TAIR10_cdna_20110103_representative_gene_model") 
#TIR10_cdna_rep_model<-readDNAStringSet("Nozue2016_SAStranscriptome_data/input/TAIR10_cdna_20110103_representative_gene_model") 
head(TIR10_cdna_rep_model)
bias<-nchar(TIR10_cdna_rep_model)
names(bias)<-substr(names(TIR10_cdna_rep_model),1,9)

GOseq.overlapTable<-function (labels1, labels2, na.rm = TRUE, ignore = NULL, levels1 = NULL, 
                              levels2 = NULL,genes,permutation = 2000)  # genes is names of all genes in this analysis (usugally rownames of matrix). This considers cDNA length bias in DE.
{ # default is permutaiton = 2000
  # for GOseq
  bias2<-bias[names(bias) %in% genes]
  #
  labels1 = as.vector(labels1)
  labels2 = as.vector(labels2)
  if (na.rm) {
    keep = !is.na(labels1) & !is.na(labels2)
    labels1 = labels1[keep]
    labels2 = labels2[keep]
  }
  if (is.null(levels1)) {
    levels1 = sort(unique(labels1))
    levels1 = levels1[!levels1 %in% ignore]
  }
  if (is.null(levels2)) {
    levels2 = sort(unique(labels2))
    levels2 = levels2[!levels2 %in% ignore]
  }
  n1 = length(levels1)
  n2 = length(levels2)
  countMat = matrix(0, n1, n2)
  pMat = matrix(0, n1, n2)
  for (m1 in 1:n1) for (m2 in 1:n2) {
    m1Members = (labels1 == levels1[m1])
    m2Members = (labels2 == levels2[m2])
    #pMat[m1, m2] = fisher.test(m1Members, m2Members, alternative = "greater")$p.value   
    # use GOseq instead of fisher's exact test 
    print(m1Members)
    print(class(m1Members))
    TF<-m1Members
    names(TF)<-genes  
    print("length of m1Members is")
    print(length(TF))
    print("length of bias is")
    print(length(bias2))
    pwf<-nullp(TF,bias.data=bias2)
    ### custom category
    temp<-list(x=genes[m2Members])
    print(levels2[m2])
    names(temp)<-levels2[m2]
    print("temp is")
    print(temp)
    GO.pval <- goseq(pwf,gene2cat=temp,use_genes_without_cat=TRUE,method="Sampling",repcnt=permutation) # format became different in new goseq version (021111)
    # enter pvalue
    pMat[m1, m2]<-GO.pval$over_represented_pvalue       
    ####################### the end of GOseq 
    countMat[m1, m2] = sum(labels1 == levels1[m1] & labels2 == 
                             levels2[m2])
  }
  dimnames(pMat) = list(levels1, levels2)
  dimnames(countMat) = list(levels1, levels2)
  pMat[is.na(pMat)] = 1
  list(countTable = countMat, pTable = pMat)
}

# test function
# DE.genotype.batch.genes.SOMclusters<-read.csv("/Volumes/Data8/NGS_related/Arabidopsis_analysis/SAS_muts_time_course_RNAseq2_mycomp_output/figs_tables/DE.genotype.batch.genes.SOMclusters.csv")
# DE.genotype.batch.genes.SOMclusters<-DE.genotype.batch.genes.SOMclusters[,-1]
# test<-DE.genotype.batch.genes.SOMclusters
# # overlap table
# i<-29 # yuc2589.shade
# overlap<-overlapTable(test[,30],test[,i]) # cluster 10 is enriched in yuc2589.shade
# overlap.GOseq<-GOseq.overlapTable(test[,30],test[,i],genes=test$AGI) # yes!



#Vlookup in R (https://gist.github.com/jnmaloof/7367450)
#Version 0.3 November 12, 2013
#Return senesical results if return column is a factor
#Version 0.2 November 11, 2013
#Require first column of table to be numeric if range lookup is being done
#Change defaults to larger=FALSE
#Julin Maloof
vlookup <- function(ref, #the value or values that you want to look for
                    table, #the table where you want to look for it; will look in first column
                    column, #the column that you want the return data to come from,
                    range=FALSE, #if there is not an exact match, return the closest?
                    larger=FALSE) #if doing a range lookup, should the smaller or larger key be used?)
{
  if(!is.numeric(column) & !column %in% colnames(table)) {
    stop(paste("can't find column",column,"in table"))
  }
  if(range) {
    if(!is.numeric(table[,1])) {
      stop(paste("The first column of table must be numeric when using range lookup"))
    }
    table <- table[order(table[,1]),] 
    index <- findInterval(ref,table[,1])
    if(larger) {
      index <- ifelse(ref %in% table[,1],index,index+1)
    }
    output <- table[index,column]
    output[!index <= dim(table)[1]] <- NA
    
  } else {
    output <- table[match(ref,table[,1]),column]
    output[!ref %in% table[,1]] <- NA #not needed?
  }
  dim(output) <- dim(ref)
  output
}
