#######
##
## function_RNAseq_time_course.R
##  
###################################################
library(lme4);library(lmerTest) # for significant analysis
library(edgeR);library(ggplot2);library(reshape2);library(ggdendro)
library(grid)
library(class)
library(MASS)
library(kohonen);library(plyr)
#library(scales) # for muted
library(WGCNA)
library(ShortRead)
library(goseq)
library(GO.db)
library("org.At.tair.db")
#library("annotate") # error in whitney
# see http://www.bioconductor.org/install/ for installation of these packages 

conversion.table<-data.frame(num=1:27, genotype=c("Col","AT5G02540_1","hy5","jar1","kat1_2","phyB","pif45","spt_11","yuc2589","PAR1_RNAi09","coi1_16","phyAB","pif3","mida9_4","bsk5_1","sto","aos","argos","co_9","Blh_1","Jea","Shahdara","Col_0","Cvi_0","Bur_0","Oy_0","Ita_0"))

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
expression.pattern.graph2b<-function(summary.vst.response.kazu,target.genes){# require ggplot2, reshape2 packages, this is only for SAS timecourse data  
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

  q<-ggplot(temp.data.melt,aes(x=time,y=value)) + geom_point() + facet_grid(variable~genotype,scale="free") + theme(strip.text.y=element_text(angle=0))
  #  q<-q + scale_y_continuous(trans=log2_trans())
  return(q)
}

expression.pattern.graph3<-function(data.cpm,target.genes){# require ggplot2, reshape2 packages, this is only for SAS timecourse data (mean type, eg. "summary4")
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
  
  q<-ggplot(temp.data.melt,aes(x=time,y=value,color=trt)) + geom_point(alpha = 0.5) 
  q<-q + facet_grid(variable~genotype,scale="free") + theme(strip.text.y=element_text(angle=0))
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
### barcoding
###### for any genotype 
SOM.clustering3<-function(data,gt="Col",cluster) { # "Col","YuQ","Coi","Jar","Pi3","Pi4","Spt","PhB","Co_" (021114)
  # genotype is numeric (061114) 
  data.SOM.gt<-data[grep(gt,rownames(data)),] # gt subset
  #print(data.SOM.gt)
  data.SOM.gt.SOM1<-data.SOM.gt[data.SOM.gt$ssom3.4.unit.classif==cluster,]
  #print(data.SOM.gt.SOM1)
  temp<-paste("(AT[[:alnum:]]+)(\\.[[:digit:]]+)(\\.",gt,")",sep="")
  #print(paste("temp is",temp))
  #print(gsub(temp,"\\1",rownames(data.SOM.gt.SOM1)))
  data.SOM.gt.SOM1.all<-data[gsub("(AT[[:alnum:]]+)(\\.[[:digit:]]+)(\\.[[:print:]]+)","\\1",rownames(data)) %in% gsub(temp,"\\1",rownames(data.SOM.gt.SOM1)),]
  #print(data.SOM.gt.SOM1.all) 
  data.SOM.gt.SOM1.all.s<-data.frame(row.names=rownames(data.SOM.gt.SOM1.all),ssom3.4.unit.classif=data.SOM.gt.SOM1.all[,"ssom3.4.unit.classif"])
  #print(data.SOM.gt.SOM1.all.s) # does not work
  data.SOM.gt.SOM1.all.s$AGI<-gsub("(AT[[:alnum:]]+)(\\.)([[:digit:]]+)(\\.)([[:print:]]+)","\\1",rownames(data.SOM.gt.SOM1.all.s))
  data.SOM.gt.SOM1.all.s$genotype<-gsub("(AT[[:alnum:]]+)(\\.)([[:digit:]]+)(\\.)([[:print:]]+)","\\5",rownames(data.SOM.gt.SOM1.all.s))
  data.SOM.gt.SOM1.all.s.sorted<-data.SOM.gt.SOM1.all.s[order(data.SOM.gt.SOM1.all.s$AGI),1:3]
  data.SOM.gt.SOM1.all.s.sorted$genotype<-as.factor(data.SOM.gt.SOM1.all.s.sorted$genotype)
  #print(data.SOM.Col.SOM1.all.s.sorted)
  # reshape
  data.SOM.gt.SOM1.all.barcode<-reshape(data.SOM.gt.SOM1.all.s.sorted,idvar="AGI",v.names ="ssom3.4.unit.classif",timevar="genotype",direction="wide")
  names(data.SOM.gt.SOM1.all.barcode)[-1]<-gsub("(ssom3.4.unit.classif\\.)([[:print:]]+)","\\2",names(data.SOM.gt.SOM1.all.barcode)[-1])
  rownames(data.SOM.gt.SOM1.all.barcode)<-data.SOM.gt.SOM1.all.barcode$AGI
  data.SOM.gt.SOM1.all.barcode.s<-data.SOM.gt.SOM1.all.barcode[,-1]
  return(data.SOM.gt.SOM1.all.barcode.s)#
}
###
SOM.clustering4<-function(data,gt=1,cluster) { # "Col","YuQ","Coi","Jar","Pi3","Pi4","Spt","PhB","Co_" (021114)
  # genotype is numeric (061114) 
  data.SOM.gt<-data[grep(paste(".",gt,"$",sep=""),rownames(data),value=T),] # gt subset
  #print(data.SOM.gt)
  data.SOM.gt.SOM1<-data.SOM.gt[data.SOM.gt$ssom3.4.unit.classif==cluster,]
  #print(data.SOM.gt.SOM1)
  temp<-paste("(AT[[:alnum:]]+)(\\.[[:digit:]]+)(\\.",gt,")",sep="")
  #print(paste("temp is",temp))
  #print(gsub(temp,"\\1",rownames(data.SOM.gt.SOM1)))
  data.SOM.gt.SOM1.all<-data[gsub("(AT[[:alnum:]]+)(\\.[[:digit:]]+)(\\.[[:print:]]+)","\\1",rownames(data)) %in% gsub(temp,"\\1",rownames(data.SOM.gt.SOM1)),]
  #print(data.SOM.gt.SOM1.all) 
  data.SOM.gt.SOM1.all.s<-data.frame(row.names=rownames(data.SOM.gt.SOM1.all),ssom3.4.unit.classif=data.SOM.gt.SOM1.all[,"ssom3.4.unit.classif"])
  #print(data.SOM.gt.SOM1.all.s) # does not work
  data.SOM.gt.SOM1.all.s$AGI<-gsub("(AT[[:alnum:]]+)(\\.)([[:digit:]]+)(\\.)([[:print:]]+)","\\1",rownames(data.SOM.gt.SOM1.all.s))
  data.SOM.gt.SOM1.all.s$genotype<-gsub("(AT[[:alnum:]]+)(\\.)([[:digit:]]+)(\\.)([[:print:]]+)","\\5",rownames(data.SOM.gt.SOM1.all.s))
  data.SOM.gt.SOM1.all.s.sorted<-data.SOM.gt.SOM1.all.s[order(data.SOM.gt.SOM1.all.s$AGI),1:3]
  data.SOM.gt.SOM1.all.s.sorted$genotype<-factor(data.SOM.gt.SOM1.all.s.sorted$genotype,levels=levels(data$gt))
  #print(data.SOM.gt.SOM1.all.s.sorted)
  #print(levels(data.SOM.gt.SOM1.all.s.sorted$genotype))
  # reshape
  data.SOM.gt.SOM1.all.barcode<-reshape(data.SOM.gt.SOM1.all.s.sorted,idvar="AGI",v.names ="ssom3.4.unit.classif",timevar="genotype",direction="wide")
  #print("data.SOM.gt.SOM1.all.barcode is")
  #print(data.SOM.gt.SOM1.all.barcode)  
  names(data.SOM.gt.SOM1.all.barcode)[-1]<-gsub("(ssom3.4.unit.classif\\.)([[:print:]]+)","\\2",names(data.SOM.gt.SOM1.all.barcode)[-1])
  rownames(data.SOM.gt.SOM1.all.barcode)<-data.SOM.gt.SOM1.all.barcode$AGI
  data.SOM.gt.SOM1.all.barcode.s<-data.SOM.gt.SOM1.all.barcode[,-1]
  data.SOM.gt.SOM1.all.barcode.s<-data.SOM.gt.SOM1.all.barcode.s[,order(names(data.SOM.gt.SOM1.all.barcode.s))]
  return(data.SOM.gt.SOM1.all.barcode.s)#
}



# ORA with GOseq
# prerequisit
#library(ShortRead);library(goseq);library(GO.db);library("org.At.tair.db");library("annotate")

TIR10_cdna_rep_model<-readDNAStringSet("../input/TAIR10_cdna_20110103_representative_gene_model") 
#TIR10_cdna_rep_model<-readDNAStringSet("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/kazu/Documents/SASmutRNAseq_whitney/TAIR10_cdna_20110103_representative_gene_model") 
#TIR10_cdna_rep_model<-readDNAStringSet("/mydata/kazu_data/SASmutRNAseq2/TAIR10_cdna_20110103_representative_gene_model") 

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

GOseq.ORA<-function(genelist,padjust=0.05) { # return GO enrichment table, padjus, padjust=0.05 
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
  GO.pval <- goseq(pwf,gene2cat=Atgo.list) # format became different in new goseq version (021111)
  #head(GO.pval) 
  GO.pval$over_represented_padjust<-p.adjust(GO.pval$over_represented_pvalue,method="BH")
  if(GO.pval$over_represented_padjust[1]>padjust) stop("no enriched GO")
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






#############
###
### GO heatmap
###
#####################
GO.expression<-function(GO,TERM,expression.data,cluster) { 
  # expression<-expression.data[grep(GO[1],expression.data$AtGO),c("X1h.1","X4h.1","X16h.1","X25h.1","X49h.1")]
  expression<-expression.data[grep(GO[1],expression.data$AtGO),c("1h","4h","16h","25h","49h")]
  #print(expression)
  expression.mean3<-colMeans(expression)
  #print(expression.mean3)
  if(length(GO)>1) {	
    for(i in 2:length(GO)) {
      expression<-expression.data[grep(GO[i],expression.data$AtGO),c("1h","4h","16h","25h","49h")]
      expression.mean<-colMeans(expression)
      expression.mean3<-rbind(expression.mean3,expression.mean)		
    }
    expression.mean3<-as.data.frame(expression.mean3)				
    rownames(expression.mean3)<-GO
    names(expression.mean3)<-c("1h","4h","16h","25h","49h")
    expression.mean3$GO<-paste(GO,TERM)
    expression.mean3$cluster<-cluster
    expression.mean3.melt<-melt(expression.mean3,id.var=c("GO","cluster"))
    return(expression.mean3.melt)			
  }
  else {
    #print(expression.mean3)
    expression.mean4<-data.frame(GO=rep(paste(GO,TERM),5),cluster=rep(cluster,5),variable=c("1h","4h","16h","25h","49h"),value=expression.mean3)
    #print(expression.mean4)		
    return(expression.mean4)	
  }	
}
# for w/o 16h
GO.expression<-function(GO,TERM,expression.data,cluster) { 
  # expression<-expression.data[grep(GO[1],expression.data$AtGO),c("X1h.1","X4h.1","X16h.1","X25h.1","X49h.1")]
  expression<-expression.data[grep(GO[1],expression.data$AtGO),c("1h","4h","25h","49h")]
  #print(expression)
  expression.mean3<-colMeans(expression)
  #print(expression.mean3)
  if(length(GO)>1) {	
    for(i in 2:length(GO)) {
      expression<-expression.data[grep(GO[i],expression.data$AtGO),c("1h","4h","25h","49h")]
      expression.mean<-colMeans(expression)
      expression.mean3<-rbind(expression.mean3,expression.mean)		
    }
    expression.mean3<-as.data.frame(expression.mean3)				
    rownames(expression.mean3)<-GO
    names(expression.mean3)<-c("1h","4h","25h","49h")
    expression.mean3$GO<-paste(GO,TERM)
    expression.mean3$cluster<-cluster
    expression.mean3.melt<-melt(expression.mean3,id.var=c("GO","cluster"))
    return(expression.mean3.melt)			
  }
  else {
    #print(expression.mean3)
    expression.mean4<-data.frame(GO=rep(paste(GO,TERM),4),cluster=rep(cluster,4),variable=c("1h","4h","25h","49h"),value=expression.mean3)
    #print(expression.mean4)		
    return(expression.mean4)	
  }	
}



############################################################################
### What are genes in given cluster and GO term?
###                 input ** output  barcode, AGI, GO term, symbol, description
###                    022514 Kazu 
############################################################################
#library(ShortRead);library(goseq);library(GO.db);library("org.At.tair.db");library("annotate")
### modified version of gene list data frame with .... (under construction and not used)
# annotate2 <- function(results,sample) { #resultsis vector of AGI name. sample is sample name for output data frame
#   #  results.annotated <- result$table #need to extract the gene table from results to get started
#   #  results.annotated$AGI <- substr(rownames(results.annotated),1,9) #get rid of the gene model number
#   results.annotated<-data.frame(AGI=results)
#   
#   #add At GO
#   Atgo <- toTable(org.At.tairGO)
#   head(Atgo)
#   
#   BP <- TRUE #only keep BP go TERMS
#   if (BP) Atgo <- Atgo[Atgo$Ontology=="BP",]
#   #convert to list
#   Atgo.list <- tapply(Atgo$go_id,Atgo$gene_id,c)
#   #collapse list
#   Atgo.terms <- data.frame(AtGO=unlist(lapply(Atgo.list,paste,collapse=";")),AGI=names(Atgo.list))
#   results.annotated <- merge(results.annotated,Atgo.terms,by="AGI",all.x=T,sort=F,incomparables=NA)
#   summary(  results.annotated$AGI %in% results)   # OK
#   
#   #add At gene symbol
#   AtSymbol <- toTable(org.At.tairSYMBOL)
#   AtSymbol <- AtSymbol[!duplicated(AtSymbol$gene_id),]
#   results.annotated <- merge(results.annotated,AtSymbol,by.x="AGI",by.y="gene_id",all.x=T,sort=F)
#   
#   #Add At description
#   AtNames <- toTable(org.At.tairGENENAME)
#   AtNames <- AtNames[!duplicated(AtNames$gene_id),]
#   anyDuplicated(AtNames$gene_id)
#   results.annotated <- merge(results.annotated,AtNames,by.x="AGI",by.y="gene_id",all.x=T,sort=F)
#   head(results.annotated)
#   return(results.annotated)
#   #  write.csv(results.annotated,sample)
# }

#####
## barcoding
############
mut.SOM.clustering.heatmap<-function(barcode,cluster) {
  # data.cor<-cor(t(barcode),method="pearson")
  # data.cor.dist<-as.dist(1-data.cor)  
  barcode2<-paste(
    barcode[,1],
    barcode[,2],
    barcode[,3],
    barcode[,4],
    barcode[,5],
    barcode[,6],
    barcode[,7],
    barcode[,8],
    barcode[,9],
    sep="")
  names(barcode2)<-rownames(barcode)
  print(barcode2)
  d  <- adist(barcode2)
  hc <- hclust(as.dist(d))
  
  pdf(file=paste("barcode.color.clustering",cluster,date(),".pdf",sep=""),width=8,height=8)
  heatmap(as.matrix(barcode),Rowv=as.dendrogram(hc),col=rainbow(12),na.rm = TRUE,scale="none",main=paste("cluster",cluster)) #,Colv=NA
  dev.off() # this is not what I want.
  ### BW
  barcode.bw<-(barcode==barcode$Col[1])*1
  barcode.bw2<-paste(
    barcode.bw[,1],
    barcode.bw[,2],
    barcode.bw[,3],
    barcode.bw[,4],
    barcode.bw[,5],
    barcode.bw[,6],
    barcode.bw[,7],
    barcode.bw[,8],
    barcode.bw[,9],
    sep="")
  names(barcode.bw2)<-rownames(barcode)
  d.bw  <- adist(barcode.bw2)
  hc.bw <- hclust(as.dist(d.bw))
  
  pdf(file=paste("barcode.bw.clustering",cluster,date(),".pdf",sep=""),width=8,height=8)
  heatmap(as.matrix(barcode.bw),Rowv=as.dendrogram(hc.bw),col=colorRampPalette(c("green", "magenta"))(n = 2),na.rm = TRUE,scale="none",main=paste("cluster",cluster)) #,Colv=NA
  dev.off() # this is not what I want.
  
}
### do not use
library(ggplot2)
library(reshape2)
library(ggdendro)

mut.SOM.clustering.heatmap.ggplot<-function(barcode,title) {
  #needs to convert clustering number 1:12 into A:L
  nuc.To.alpha.data<-data.frame(nuc=1:12,alpha=c("A","B","C","D","E","F","G","H","I","J","K","K"))
  head(barcode)
  nuc.To.alpha<-function(number) {nuc.To.alpha.data[grep(number,nuc.To.alpha.data[,1])[1],2]}
  barcode.alpha<-data.frame()
  for(i in 1:dim(barcode)[1]) {
    for(n in 1:dim(barcode)[2]) {
      barcode.alpha[i,n] <- nuc.To.alpha(barcode[i,n])  
    }
  }
  rownames(barcode.alpha)<-rownames(barcode)
  colnames(barcode.alpha)<-colnames(barcode)
  # for clustering AGI
  barcode2<-paste(
    barcode.alpha[,1],
    barcode.alpha[,2],
    barcode.alpha[,3],
    barcode.alpha[,4],
    barcode.alpha[,5],
    barcode.alpha[,6],
    barcode.alpha[,7],
    barcode.alpha[,8],
    barcode.alpha[,9],
    barcode.alpha[,10],
    barcode.alpha[,11],
    barcode.alpha[,12],
    barcode.alpha[,13],
    barcode.alpha[,14],
    barcode.alpha[,15],
    sep="")
  names(barcode2)<-rownames(barcode.alpha)
  #	print(barcode2)
  d.col <-adist(barcode2)
  hc.col<-hclust(as.dist(d.col))
  
  dd.col <- as.dendrogram(hc.col)
  col.ord <- order.dendrogram(dd.col)
  
  # calculate hc for genotype 
  xx.alpha<-barcode.alpha
  barcode.row<-vector()
  
  for(t in 1:15) {		
    temp <- xx.alpha[1,t]
    for(m in 2:dim(xx.alpha)[1]) {
      temp2<- xx.alpha[m,t]
      temp <- paste(temp, temp2,sep="")		
    }
    barcode.row[t]<-temp
  }
  
  d.row <-adist(barcode.row)
  hc.row<-hclust(as.dist(d.row))
  dd.row <- as.dendrogram(hc.row)
  row.ord <- order.dendrogram(dd.row)
  
  #### sort data by col.ord and row.ord
  x <- as.matrix(barcode)
  #xx <- scale(mtcars)[col.ord, row.ord] # not need to scale
  xx <- barcode[col.ord,row.ord]
  print("xx (sorted) is")
  print(xx[seq(dim(xx)[1],1),])
  
  xx_names <- attr(xx, "dimnames")
  df <- as.data.frame(xx)
  #colnames(df) <- xx_names[[2]]
  df$AGI <- rownames(xx)
  df$AGI <- with(df, factor(AGI, levels=AGI, ordered=TRUE))
  
  # print("xx is")
  # print(xx)
  
  mdf <- melt(df, id.vars="AGI")
  
  #Extract dendrogram data and create the plots
  ddata_x <- dendro_data(dd.row)
  ddata_y <- dendro_data(dd.col)
  
  ### Set up a blank theme
  theme_none <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(colour=NA),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank()
    #axis.ticks.length = element_blank()
  )
  ### color coding matched up with boxplot 
  #mdf$value<-factor(mdf$value,levels=c(10,7,6,11,3,9,5,1,4,8,12,2))
  #mdf$value<-factor(mdf$value,levels=c(12, 6,7,3,4,11,2,9,8,10,5,1))
  mdf$value<-factor(mdf$value,levels=c(10, 7,12,1,11,4,2,8,5,9,6,3)) #081814 
  
  ### Create plot components ###    
  # Heatmap
  p1 <- ggplot(mdf, aes(x=variable, y=AGI)) +
    geom_tile(aes(fill=factor(value))) + guides(fill = guide_legend(title = "cluster")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  p1 <- p1 + geom_text(aes(label = value))
  
  #+ scale_fill_gradient2()
  
  
  # Dendrogram 1
  p2 <- ggplot(segment(ddata_x)) + labs(title=title) + 
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
    theme_none + theme(axis.title.x=element_blank())
  
  # Dendrogram 2
  p3 <- ggplot(segment(ddata_y))   +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
    coord_flip() + theme_none 
  #Use grid graphics and some manual alignment to position the three plots on the page
  ### Draw graphic ###
  pdf(file=paste("heatmap.ggplot",date(),".pdf",sep=""),width=11,height=18)
  print("barcode gene num is")
  print(dim(barcode))
  adjfact<-(dim(barcode)[1]/12)*0.022 + 1
  
  library(grid)
  grid.newpage()
  print(p1, vp=viewport(0.84, 0.8, x=0.43, y=0.4))
  print(p2, vp=viewport(0.61, 0.2, x=0.43, y=0.89))
  print(p3, vp=viewport(0.19, 0.785*adjfact, x=0.91, y=0.405))
  dev.off()
}
############# graphng is good, but results is messy for all barcodes (021414)
### needs to work on graph position (021514)

# library(ggplot2)
# library(reshape2)
# library(ggdendro)

# mut.WGCNA.module.clustering.heatmap.ggplot<-function(barcode,title) {
#   # there are new modules in some mutant
#   load("/Volumes/Data8/NGS_related/Arabidopsis_analysis/SAS_muts_time_course_RNAseq2_mycomp_output/SAS.expression.vst.s.kazu.largeCV.modules.Rdata") 
#   SAS.expression.vst.s.kazu.largeCV.modules.melt<-melt(SAS.expression.vst.s.kazu.largeCV.modules,id="AGI")
#   all.module<-unique(SAS.expression.vst.s.kazu.largeCV.modules.melt$value)
#   all.module<-all.module[!all.module=="0"]
#   mod.to.apha.data<-data.frame(module=all.module,alpha=c(LETTERS[1:23],letters[1:22]))
#   #
#   head(barcode)
#   mod.To.alpha<-function(module) {mod.to.apha.data[grep(module,mod.to.apha.data[,1])[1],2]}
#   barcode.alpha<-data.frame()
#   for(i in 1:dim(barcode)[1]) {
#     for(n in 1:dim(barcode)[2]) {
#       barcode.alpha[i,n] <- mod.To.alpha(barcode[i,n])  
#     }
#   }
#   rownames(barcode.alpha)<-rownames(barcode)
#   colnames(barcode.alpha)<-colnames(barcode)
#   print("new barcode is")
#   print(barcode.alpha)
#   
#   # for clustering AGI
#   barcode2<-paste(
#     barcode.alpha[,1],
#     barcode.alpha[,2],
#     barcode.alpha[,3],
#     barcode.alpha[,4],
#     barcode.alpha[,5],
#     barcode.alpha[,6],
#     barcode.alpha[,7],
#     barcode.alpha[,8],
#     barcode.alpha[,9],
#     barcode.alpha[,10],
#     barcode.alpha[,11],
#     barcode.alpha[,12],
#     barcode.alpha[,13],
#     #barcode.alpha[,14],
#     #barcode.alpha[,15],
#     sep="")
#   names(barcode2)<-rownames(barcode.alpha)
#   #  print(barcode2)
#   d.col <-adist(barcode2)
#   hc.col<-hclust(as.dist(d.col))
#   
#   dd.col <- as.dendrogram(hc.col)
#   col.ord <- order.dendrogram(dd.col)
#   
#   # calculate hc for genotype 
#   xx.alpha<-barcode.alpha
#   barcode.row<-vector()
#   
#   #for(t in 1:15) {	
#   for(t in 1:13) {  
#       
#     temp <- xx.alpha[1,t]
#     for(m in 2:dim(xx.alpha)[1]) {
#       temp2<- xx.alpha[m,t]
#       temp <- paste(temp, temp2,sep="")		
#     }
#     barcode.row[t]<-temp
#   }
#   
#   d.row <-adist(barcode.row)
#   hc.row<-hclust(as.dist(d.row))
#   dd.row <- as.dendrogram(hc.row)
#   row.ord <- order.dendrogram(dd.row)
#   
#   #### sort data by col.ord and row.ord
#   x <- as.matrix(barcode)
#   #xx <- scale(mtcars)[col.ord, row.ord] # not need to scale
#   xx <- barcode[col.ord,row.ord]
#   print("xx (sorted) is")
#   print(xx[seq(dim(xx)[1],1),])
#   
#   xx_names <- attr(xx, "dimnames")
#   df <- as.data.frame(xx)
#   #colnames(df) <- xx_names[[2]]
#   df$AGI <- rownames(xx)
#   df$AGI <- with(df, factor(AGI, levels=AGI, ordered=TRUE))
#   
#   # print("xx is")
#   # print(xx)
#   
#   mdf <- melt(df, id.vars="AGI")
#   
#   #Extract dendrogram data and create the plots
#   ddata_x <- dendro_data(dd.row)
#   ddata_y <- dendro_data(dd.col)
#   
#   ### Set up a blank theme
#   theme_none <- theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     axis.title.x = element_text(colour=NA),
#     axis.title.y = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.line = element_blank()
#     #axis.ticks.length = element_blank()
#   )
#   ### color coding matched up with boxplot 
#   #mdf$value<-factor(mdf$value,levels=c(10,7,6,11,3,9,5,1,4,8,12,2))
#   #mdf$value<-factor(mdf$value,levels=c(12, 6,7,3,4,11,2,9,8,10,5,1))
#   #mdf$value<-factor(mdf$value,levels=c(10, 7,12,1,11,4,2,8,5,9,6,3)) #081814 
#   
#   ### Create plot components ###    
#   # Heatmap
#   p1 <- ggplot(mdf, aes(x=variable, y=AGI)) + geom_tile(aes(fill=factor(value))) + guides(fill = guide_legend(title = "cluster")) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
#   p1 <- p1 + geom_text(aes(label = value))
#   
#   #+ scale_fill_gradient2()
#   
#   
#   # Dendrogram 1
#   p2 <- ggplot(segment(ddata_x)) + labs(title=title) + 
#     geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
#     theme_none + theme(axis.title.x=element_blank())
#   
#   # Dendrogram 2
#   p3 <- ggplot(segment(ddata_y))   +
#     geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
#     coord_flip() + theme_none 
#   #Use grid graphics and some manual alignment to position the three plots on the page
#   ### Draw graphic ###
#   pdf(file=paste("heatmap.ggplot",date(),".pdf",sep=""),width=29,height=15)
#   print("barcode gene num is")
#   print(dim(barcode2))
#   adjfact<-(dim(barcode)[1]/12)*0.022 + 1
#   
#   library(grid)
#   grid.newpage()
#   print(p1, vp=viewport(0.84, 0.8, x=0.43, y=0.4))
#   print(p2, vp=viewport(0.61, 0.2, x=0.43, y=0.89))
#   print(p3, vp=viewport(0.19, 0.785*adjfact, x=0.91, y=0.405))
#   dev.off()
# } 
## the end

# mut.WGCNA.module.clustering.heatmap.ggplot2<-function(barcode,misexpression,title) { # text modules only misexpressed in mutants
# #   barcode<-salmon.mut.transition.sss # this is only for testing function
# #   misexpression<-salmon.disc.s
#   # there are new modules in some mutant
#   load("/Volumes/Data8/NGS_related/Arabidopsis_analysis/SAS_muts_time_course_RNAseq2_mycomp_output/SAS.expression.vst.s.kazu.largeCV.modules.Rdata") 
#   SAS.expression.vst.s.kazu.largeCV.modules.melt<-melt(SAS.expression.vst.s.kazu.largeCV.modules,id="AGI")
#   all.module<-unique(SAS.expression.vst.s.kazu.largeCV.modules.melt$value)
#   all.module<-all.module[!all.module=="0"]
#   mod.to.apha.data<-data.frame(module=all.module,alpha=c(LETTERS[1:23],letters[1:22]))
#   #
#   head(barcode)
#   mod.To.alpha<-function(module) {mod.to.apha.data[grep(module,mod.to.apha.data[,1])[1],2]}
#   barcode.alpha<-data.frame()
#   for(i in 1:dim(barcode)[1]) {
#     for(n in 1:dim(barcode)[2]) {
#       barcode.alpha[i,n] <- mod.To.alpha(barcode[i,n])  
#     }
#   }
#   rownames(barcode.alpha)<-rownames(barcode)
#   colnames(barcode.alpha)<-colnames(barcode)
#   print("new barcode is")
#   print(barcode.alpha)
#   
#   # for clustering AGI
#   barcode2<-paste(
#     barcode.alpha[,1],
#     barcode.alpha[,2],
#     barcode.alpha[,3],
#     barcode.alpha[,4],
#     barcode.alpha[,5],
#     barcode.alpha[,6],
#     barcode.alpha[,7],
#     barcode.alpha[,8],
#     barcode.alpha[,9],
#     barcode.alpha[,10],
#     barcode.alpha[,11],
#     barcode.alpha[,12],
#     barcode.alpha[,13],
#     #barcode.alpha[,14],
#     #barcode.alpha[,15],
#     sep="")
#   names(barcode2)<-rownames(barcode.alpha)
#   #  print(barcode2)
#   d.col <-adist(barcode2)
#   hc.col<-hclust(as.dist(d.col))
#   
#   dd.col <- as.dendrogram(hc.col)
#   col.ord <- order.dendrogram(dd.col)
#   
#   # calculate hc for genotype 
#   xx.alpha<-barcode.alpha
#   barcode.row<-vector()
#   
#   #for(t in 1:15) {  
#   for(t in 1:13) {  
#     
#     temp <- xx.alpha[1,t]
#     for(m in 2:dim(xx.alpha)[1]) {
#       temp2<- xx.alpha[m,t]
#       temp <- paste(temp, temp2,sep="")		
#     }
#     barcode.row[t]<-temp
#   }
#   
#   d.row <-adist(barcode.row)
#   hc.row<-hclust(as.dist(d.row))
#   dd.row <- as.dendrogram(hc.row)
#   row.ord <- order.dendrogram(dd.row)
#   
#   #### sort data by col.ord and row.ord
#   x <- as.matrix(barcode)
#   #xx <- scale(mtcars)[col.ord, row.ord] # not need to scale
#   xx <- barcode[col.ord,row.ord]
#   print("xx (sorted) is")
#   print(xx[seq(dim(xx)[1],1),])
#   
#   xx_names <- attr(xx, "dimnames")
#   df <- as.data.frame(xx)
#   #colnames(df) <- xx_names[[2]]
#   df$AGI <- rownames(xx)
#   df$AGI <- with(df, factor(AGI, levels=AGI, ordered=TRUE))
#   
#   # print("xx is")
#   # print(xx)
#   df.melt<-melt(df, id.vars="AGI")
#   names(df.melt)[2:3]<-c("genotype","WGCNAmodule")
#   df.melt$genotype_AGI<-paste(df.melt$genotype,df.melt$AGI,sep="_")
#   # misexpressed genes (sun)
#   salmon.disc.s.sun<-salmon.disc.s[,c("AGI",grep(".sun",names(salmon.disc.s),value=TRUE))]
#   misexpress.sun.melt<-melt(salmon.disc.s.sun,id.vars="AGI")
#   misexpress.sun.melt$genotype<-gsub("([[:print:]]+)(.)(sun|shade)","\\1",misexpress.sun.melt$variable)
#   misexpress.sun.melt$light<-gsub("([[:print:]]+)(.)(sun|shade)","\\3",misexpress.sun.melt$variable)
#   misexpress.sun.melt$genotype_AGI<-paste(misexpress.sun.melt$genotype,misexpress.sun.melt$AGI,sep="_")
#   names(misexpress.sun.melt)[names(misexpress.sun.melt)=="value"]<-"sun"
#   misexpress.sun.melt$sun[misexpress.sun.melt$sun=="1"]<-"H"
#   misexpress.sun.melt$sun[misexpress.sun.melt$sun=="0"]<-""
#   # misexpressed genes (shade)
#   salmon.disc.s.shade<-salmon.disc.s[,c("AGI",grep(".shade",names(salmon.disc.s),value=TRUE))]
#   misexpress.shade.melt<-melt(salmon.disc.s.shade,id.vars="AGI")
#   misexpress.shade.melt$genotype<-gsub("([[:print:]]+)(.)(sun|shade)","\\1",misexpress.shade.melt$variable)
#   misexpress.shade.melt$light<-gsub("([[:print:]]+)(.)(sun|shade)","\\3",misexpress.shade.melt$variable)
#   misexpress.shade.melt$genotype_AGI<-paste(misexpress.shade.melt$genotype,misexpress.shade.melt$AGI,sep="_")
#   names(misexpress.shade.melt)[names(misexpress.shade.melt)=="value"]<-"shade"
#   misexpress.shade.melt$shade[misexpress.shade.melt$shade=="1"]<-"L"
#   misexpress.shade.melt$shade[misexpress.shade.melt$shade=="0"]<-""
#   
#   # merge
#   mdf.misexpression<-merge(df.melt,misexpress.sun.melt[,c("genotype_AGI","sun")],by="genotype_AGI")
#   mdf.misexpression<-merge(mdf.misexpression,misexpress.shade.melt[,c("genotype_AGI","shade")],by="genotype_AGI")
#   ### changed module names did not match to misexpression...
# #   > head(mdf.misexpression)
# #   genotype_AGI       AGI genotype WGCNAmodule sun shade
# #   1 aos_AT1G02920 AT1G02920      aos greenyellow          
# #   2 aos_AT1G09415 AT1G09415      aos      salmon          
# #   3 aos_AT1G14880 AT1G14880      aos      salmon   H     L
# #   4 aos_AT1G19960 AT1G19960      aos     magenta   H     L
# #   5 aos_AT1G21250 AT1G21250      aos      salmon   H      
# #   6 aos_AT1G24147 AT1G24147      aos      salmon          
#   ### better to use method for SOM cluster analysis; treating gene + genotype = new gene (062115)
#   ### due to missing values, this way of WGCNA did not work (062215)
#   
#   #mdf <- melt(df, id.vars="AGI")
#   #mdf.misexpression <- melt(df.misexpression, id.vars="AGI") # this is not what I want. Have two columns for values (0/1 and module)
#   
# mdf.misexpression$WGCNAmodule_HL<-paste(mdf.misexpression$WGCNAmodule,  mdf.misexpression$sun,  mdf.misexpression$shade,sep="_")
# mdf.misexpression$WGCNAmodule_HL[grep("__$",mdf.misexpression$WGCNAmodule_HL)]<-""  
#   
#   #Extract dendrogram data and create the plots
#   ddata_x <- dendro_data(dd.row)
#   ddata_y <- dendro_data(dd.col)
#   
#   ### Set up a blank theme
#   theme_none <- theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     axis.title.x = element_text(colour=NA),
#     axis.title.y = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.line = element_blank()
#     #axis.ticks.length = element_blank()
#   )
#   ### color coding matched up with boxplot 
#   #mdf$value<-factor(mdf$value,levels=c(10,7,6,11,3,9,5,1,4,8,12,2))
#   #mdf$value<-factor(mdf$value,levels=c(12, 6,7,3,4,11,2,9,8,10,5,1))
#   #mdf$value<-factor(mdf$value,levels=c(10, 7,12,1,11,4,2,8,5,9,6,3)) #081814 
#   
#   ### Create plot components ###    
#   # Heatmap
#   p1 <- ggplot(mdf.misexpression, aes(x=genotype, y=AGI)) + geom_tile(aes(fill=factor(WGCNAmodule))) + guides(fill = guide_legend(title = "cluster")) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
#   p1 <- p1 + geom_text(aes(label = WGCNAmodule_HL))
#   p1
#   #+ scale_fill_gradient2()
#   
#   
#   # Dendrogram 1
#   p2 <- ggplot(segment(ddata_x)) + labs(title=title) + 
#     geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
#     theme_none + theme(axis.title.x=element_blank())
#   
#   # Dendrogram 2
#   p3 <- ggplot(segment(ddata_y))   +
#     geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
#     coord_flip() + theme_none 
#   #Use grid graphics and some manual alignment to position the three plots on the page
#   ### Draw graphic ###
#   pdf(file=paste("heatmap.ggplot",date(),".pdf",sep=""),width=29,height=29)
#   print("barcode gene num is")
#   print(dim(barcode2))
#   adjfact<-(dim(barcode)[1]/12)*0.022 + 1
#   
#   library(grid)
#   grid.newpage()
#   print(p1, vp=viewport(0.84, 0.8, x=0.43, y=0.4))
#   print(p2, vp=viewport(0.61, 0.2, x=0.43, y=0.89))
#   print(p3, vp=viewport(0.19, 0.785*adjfact, x=0.91, y=0.405))
#   dev.off()
# } 
# ## the end




library(functional) # for Compose function
my.category.heatmap<-function(summary3.response,gt.num,gt2) { # gt.num (num), gt2 (gt name)
### for temp 
# calculate newdata.SOM.temp (needs to fix becasue changing into gt = 1 did not reproduce above)
#newdata.SOM.temp<-newdata[newdata$gt==9,] # not all genes were selected. Start from summary3.response
summary3.response.log2.temp<-log2(summary3.response[,gsub("([[:digit:]]+)(_)(1|4|16|25|49)hrA","\\1",names(summary3.response))==gt.num]) # 
head(summary3.response.log2.temp)
summary3.response.log2.temp$AGI<-gsub("(AT[[:print:]]+).[[:digit:]]","\\1",rownames(summary3.response.log2.temp))

#newdata.SOM.temp.wGO<-merge(newdata.SOM.temp,Atgo.terms,by="AGI") # I should use newdata.SOM.temp for Col. no need to add AtGO
#### select genes in newdata.Col.wGO.selected and look expression pattern in temp
selected.genes<-newdata.Col.wGO.selected$AGI
summary3.response.log2.temp.selected<-summary3.response.log2.temp[summary3.response.log2.temp$AGI %in% selected.genes,]
# not that way. use Col data
summary3.response.log2.temp.selected<-merge(summary3.response.log2.temp.selected,newdata.Col.wGO.selected[,c("AGI","my.category")],by="AGI")
#
names(summary3.response.log2.temp.selected)[2:6]<-as.numeric(gsub("([[:digit:]]+)(_)(1|4|16|25|49)hrA","\\3",names(summary3.response.log2.temp.selected)[2:6]))
# remove "-Inf" gene # alternative http://stackoverflow.com/questions/15773189/remove-na-nan-inf-in-a-matrix
summary3.response.log2.temp.selected.noINF<-summary3.response.log2.temp.selected[apply(summary3.response.log2.temp.selected[,2:6], 1, Compose(is.finite, all)),]

# melt
summary3.response.log2.temp.selected.noINF.melt<-melt(summary3.response.log2.temp.selected.noINF,id=c("AGI","my.category"))
head(summary3.response.log2.temp.selected.noINF.melt)
summary(summary3.response.log2.temp.selected.noINF.melt)
summary3.response.log2.temp.selected.noINF.melt$variable<-factor(summary3.response.log2.temp.selected.noINF.melt$variable,levels=c("1","4","16","25","49"))
summary3.response.log2.temp.selected.noINF.melt$my.category<-factor(summary3.response.log2.temp.selected.noINF.melt$my.category,levels=c("auxin","ABA","defense","SA","JA","energy","ethylene"))
# mean table (= heat map should reflect this table)
summary.table<-as.data.frame(tapply(summary3.response.log2.temp.selected.noINF.melt$value,list(summary3.response.log2.temp.selected.noINF.melt$my.category,summary3.response.log2.temp.selected.noINF.melt$variable),mean))
summary.table$my.category<-factor(rownames(summary.table),levels=c("auxin","brassinosteroid","cell wall","GA","abiotic stress","defense","JA","SAR","energy"))
sumamry.table.melt<-melt(summary.table,id="my.category")
names(summary.table)
summary.table
print(paste("summary table for ",gt2))
print(summary.table)

p.selected <- ggplot(sumamry.table.melt,aes(x=variable,y=my.category,fill=value)) + geom_tile(size=0.3) 
#p.selected <- p.selected + scale_fill_gradient(limit=c(-0.5,3), high="magenta", low="green") 
#p.selected <- p.selected + scale_fill_gradient2(limits=c(-1,3), mid="gray2",high="magenta", low="green") 
library(scales) # for muted
p.selected <- p.selected + scale_fill_gradient2(limits=c(-2,3),low=muted("green"), high=muted("magenta"),guide = guide_legend(title = "log2 fold change"))
p.selected <- p.selected + theme_bw() + theme(legend.position="right",axis.title.y=element_blank(),axis.text.x=element_text(size=18),axis.ticks = element_blank())#,strip.background=element_rect(fill=c("red","blue")))
p.selected <- p.selected + labs(x="Time (hr)",y="",title=gt2)
#p.all <- p.all + facet_grid(variable~.,space="free",scale="free") # scale="free"
#p.selected <- p.selected + guides(color=guide_legend(title="log2 fold change")) # does not work
#p.selected
ggsave(file=paste("/mydata/kazu_data/SASmutRNAseq2/fix_no_outlier_elimination/fixed_timelabel/my.category.",gt2,".expression.pattern.noscaleFDR0.001.pdf",sep=""),width=5,height=2.5)
#### the end of temp
}

my.category.heatmap2<-function(summary3.response,gt.num,gt2) { # gt.num (num), gt2 (gt name), for w/o 16h data
### for temp 
# calculate newdata.SOM.temp (needs to fix becasue changing into gt = 1 did not reproduce above)
#newdata.SOM.temp<-newdata[newdata$gt==9,] # not all genes were selected. Start from summary3.response
summary3.response.log2.temp<-log2(summary3.response[,gsub("([[:digit:]]+)(_)(1|4|16|25|49)hrA","\\1",names(summary3.response))==gt.num]) # 
head(summary3.response.log2.temp)
summary3.response.log2.temp$AGI<-gsub("(AT[[:print:]]+).[[:digit:]]","\\1",rownames(summary3.response.log2.temp))

#newdata.SOM.temp.wGO<-merge(newdata.SOM.temp,Atgo.terms,by="AGI") # I should use newdata.SOM.temp for Col. no need to add AtGO
#### select genes in newdata.Col.wGO.selected and look expression pattern in temp
selected.genes<-newdata.Col.wGO.selected$AGI
summary3.response.log2.temp.selected<-summary3.response.log2.temp[summary3.response.log2.temp$AGI %in% selected.genes,]
# not that way. use Col data
summary3.response.log2.temp.selected<-merge(summary3.response.log2.temp.selected,newdata.Col.wGO.selected[,c("AGI","my.category")],by="AGI")
#
names(summary3.response.log2.temp.selected)[2:5]<-as.numeric(gsub("([[:digit:]]+)(_)(1|4|16|25|49)hrA","\\3",names(summary3.response.log2.temp.selected)[2:5]))
# remove "-Inf" gene
summary3.response.log2.temp.selected.noINF<-summary3.response.log2.temp.selected[apply(summary3.response.log2.temp.selected[,2:5], 1, Compose(is.finite, all)),]

# melt
summary3.response.log2.temp.selected.noINF.melt<-melt(summary3.response.log2.temp.selected.noINF,id=c("AGI","my.category"))
head(summary3.response.log2.temp.selected.noINF.melt)
summary(summary3.response.log2.temp.selected.noINF.melt)
summary3.response.log2.temp.selected.noINF.melt$variable<-factor(summary3.response.log2.temp.selected.noINF.melt$variable,levels=c("1","4","25","49"))
summary3.response.log2.temp.selected.noINF.melt$my.category<-factor(summary3.response.log2.temp.selected.noINF.melt$my.category,levels=c("energy","SAR","JA","defense","abiotic stress","GA","cell wall","brassinosteroid","auxin"))
# mean table (= heat map should reflect this table)
summary.table<-as.data.frame(tapply(summary3.response.log2.temp.selected.noINF.melt$value,list(summary3.response.log2.temp.selected.noINF.melt$my.category,summary3.response.log2.temp.selected.noINF.melt$variable),mean))
summary.table$my.category<-factor(rownames(summary.table),levels=c("energy","SAR","JA","defense","abiotic stress","GA","cell wall","brassinosteroid","auxin"))
sumamry.table.melt<-melt(summary.table,id="my.category")
names(summary.table)
summary.table
print(paste("summary table for ",gt2))
print(summary.table)

p.selected <- ggplot(sumamry.table.melt,aes(x=variable,y=my.category,fill=value)) + geom_tile(size=0.3) 
#p.selected <- p.selected + scale_fill_gradient(limit=c(-0.5,3), high="magenta", low="green") 
#p.selected <- p.selected + scale_fill_gradient2(limits=c(-1,3), mid="gray2",high="magenta", low="green") 
library(scales) # for muted
p.selected <- p.selected + scale_fill_gradient2(limits=c(-2,3),low=muted("green"), high=muted("magenta"),guide = guide_legend(title = "log2 fold change"))
p.selected <- p.selected + theme_bw() + theme(legend.position="bottom",axis.title.y=element_blank(),axis.text.x=element_text(size=18),axis.ticks = element_blank())#,strip.background=element_rect(fill=c("red","blue")))
p.selected <- p.selected + labs(x="Time (hr)",y="",title=gt2)
#p.all <- p.all + facet_grid(variable~.,space="free",scale="free") # scale="free"
#p.selected <- p.selected + guides(color=guide_legend(title="log2 fold change")) # does not work
#p.selected
ggsave(file=paste("my.category.",gt2,".expression.pattern.noscaleFDR0.01.pdf",sep=""),width=3.5,height=2.5)
#### the end of temp
}

###
my.category.heatmap3<-function(summary3.response,gt.num,gt2,newdata.Col.wGO.selected,plotname) { # gt.num (num), gt2 (gt name), for w/o 16h data
  ### for temp 
  # calculate newdata.SOM.temp (needs to fix becasue changing into gt = 1 did not reproduce above)
  #newdata.SOM.temp<-newdata[newdata$gt==9,] # not all genes were selected. Start from summary3.response
  summary3.response.log2.temp<-log2(summary3.response[,gsub("([[:digit:]]+)(_)(1|4|16|25|49)hrA","\\1",names(summary3.response))==gt.num]) # 
  print(head(summary3.response.log2.temp))
  summary3.response.log2.temp$AGI<-gsub("(AT[[:print:]]+).[[:digit:]]","\\1",rownames(summary3.response.log2.temp))
  print(head(summary3.response.log2.temp$AGI))
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
  summary3.response.log2.temp.selected.noINF<-summary3.response.log2.temp.selected[apply(summary3.response.log2.temp.selected[,2:5], 1, Compose(is.finite, all)),]
  print(summary3.response.log2.temp.selected)
  # melt
  summary3.response.log2.temp.selected.noINF.melt<-melt(summary3.response.log2.temp.selected.noINF,id=c("AGI","my.category"))
  head(summary3.response.log2.temp.selected.noINF.melt)
  summary(summary3.response.log2.temp.selected.noINF.melt)
  summary3.response.log2.temp.selected.noINF.melt$variable<-factor(summary3.response.log2.temp.selected.noINF.melt$variable,levels=c("1","4","25","49"))
  summary3.response.log2.temp.selected.noINF.melt$my.category<-factor(summary3.response.log2.temp.selected.noINF.melt$my.category,levels=levels(newdata.Col.wGO.selected$my.category))
  # mean table (= heat map should reflect this table)
  summary.table<-as.data.frame(tapply(summary3.response.log2.temp.selected.noINF.melt$value,list(summary3.response.log2.temp.selected.noINF.melt$my.category,summary3.response.log2.temp.selected.noINF.melt$variable),mean))
  summary.table$my.category<-factor(rownames(summary.table),levels=levels(newdata.Col.wGO.selected$my.category))
  sumamry.table.melt<-melt(summary.table,id="my.category")
  #names(summary.table)
  #summary.table
  print(paste("summary table for ",gt2))
  print(summary.table)
  
  p.selected <- ggplot(sumamry.table.melt,aes(x=variable,y=my.category,fill=value)) + geom_tile(size=0.3) 
   library(scales) # for muted
   p.selected <- p.selected + scale_fill_gradient2(limits=c(-2,3),low=muted("green"), high=muted("magenta"),guide = guide_legend(title = "log2 fold change"))
   #p.selected <- p.selected + theme_bw() + theme(legend.position="none",axis.title.y=element_blank(),axis.text.x=element_text(size=18),axis.ticks = element_blank())#,strip.background=element_rect(fill=c("red","blue")))
      p.selected <- p.selected + theme_bw() + theme(axis.title.y=element_blank(),axis.text.x=element_text(size=18),axis.ticks = element_blank())#,strip.background=element_rect(fill=c("red","blue")))

   p.selected <- p.selected + labs(x="Time (hr)",y="",title=gt2)
   ggsave(file=paste(plotname,gt2,".expression.pattern.pdf",sep=""),width=5,height=2.5)
  #### the end of temp
  return(summary3.response.log2.temp.selected.noINF)
}
#############
my.category.heatmap4<-function(summary3.response,gt.num,gt2,newdata.Col.wGO.selected,plotname,h=5,w=5*4/3,path,save.plot=T) { # gt.num (num), gt2 (gt name), for w/o 16h data, newdata.Col.wGO.selected should have "AGI" and "my.category"
  # voom transformed data has been transformd into log2
  ### for temp 
  # calculate newdata.SOM.temp (needs to fix becasue changing into gt = 1 did not reproduce above)
  #newdata.SOM.temp<-newdata[newdata$gt==9,] # not all genes were selected. Start from summary3.response
  summary3.response.log2.temp<-log2(summary3.response[,gsub("([[:digit:]]+)(_)(1|4|16|25|49)hrA","\\1",names(summary3.response))==gt.num]) # 
  #summary3.response.log2.temp<-summary3.response[,gsub("([[:digit:]]+)(_)(1|4|16|25|49)hrA","\\1",names(summary3.response))==gt.num] # 121315
  #summary3.response.log2.temp<-summary3.response[,gsub("([[:digit:]]+)(_)(1|4|16|25|49)hrA","\\1",names(summary3.response))==gt.num] #   
  print(head(summary3.response.log2.temp))
  #summary3.response.log2.temp$AGI<-gsub("(AT[[:print:]]+).[[:digit:]]","\\1",rownames(summary3.response.log2.temp))
  summary3.response.log2.temp$AGI<-rownames(summary3.response.log2.temp)
  
  print(head(summary3.response.log2.temp$AGI))
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
  summary3.response.log2.temp.selected.noINF<-summary3.response.log2.temp.selected[apply(summary3.response.log2.temp.selected[,2:5], 1, Compose(is.finite, all)),]
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
  
  sumamry.table.melt<-melt(summary.table,id="my.category")
  #names(summary.table)
  #summary.table
  
  p.selected <- ggplot(sumamry.table.melt,aes(x=variable,y=my.category,fill=value)) + geom_tile(size=0.3,colour="black") 
   library(scales) # for muted
  p.selected <- p.selected + scale_fill_gradient2(limits=c(-1,1),low=muted("green"), high=muted("magenta")) 
  # scale_fill_gradient2 itself create "Non Lab interpolation is deprecated " warning (121315)
  #,guide = guide_legend(title = "log2 fold change"
   #p.selected <- p.selected + theme_bw() + theme(legend.position="none",axis.title.y=element_blank(),axis.text.x=element_text(size=18),axis.ticks = element_blank())#,strip.background=element_rect(fill=c("red","blue")))
  p.selected <- p.selected  + 
    theme(axis.title.y=element_blank(),
          axis.text.x=element_text(size=18),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white"))
  #,strip.background=element_rect(fill=c("red","blue")))
   p.selected <- p.selected + labs(x="Time (hr)",y="",title=gt2,fill="log2 fold change")
  if(save.plot==T)   {
    ggsave(file=paste(gt2,plotname,sep="."),p.selected,height=h,width=w,path=path)
    } 
  else {return(sumamry.table.melt)}
  #### the end of temp
  
}
## new version (052716)
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
  
  sumamry.table.melt<-melt(summary.table,id="my.category")
  #names(summary.table)
  #summary.table
  # additional part in this function
  sumamry.table.melt$variable<-factor(sumamry.table.melt$variable,levels=c("49","25","4","1"))
  sumamry.table.melt$my.category<-factor(sumamry.table.melt$my.category,levels=rev(levels(sumamry.table.melt$my.category)))
  
  #p.selected <- ggplot() + geom_tile(sumamry.table.melt,aes(x=variable,y=my.category,fill=value),size=0.3,colour="black") # does not work (052816)
  p.selected <- ggplot(sumamry.table.melt,aes(x=variable,y=my.category)) + geom_tile(size=0.3,colour="black",aes(fill=value)) 
  
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
    p.selected <- p.selected  + labs(fill="log2\n fold change",legend.text=element_text(size=20)) 
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
my.category.diff.heatmap1<-function(summary3.response,gt.num,gt2,newdata.Col.wGO.selected,plotname,h=5,w=5*4/3,path,save.plot=T,legend.fill=F) { # gt.num (num), gt2 (gt name), for w/o 16h data, newdata.Col.wGO.selected should have "AGI" and "my.category"
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
  
  sumamry.table.melt<-melt(summary.table,id="my.category")
  #names(summary.table)
  #summary.table
  # additional part in this function
  sumamry.table.melt$variable<-factor(sumamry.table.melt$variable,levels=c("49","25","4","1"))
  sumamry.table.melt$my.category<-factor(sumamry.table.melt$my.category,levels=rev(levels(sumamry.table.melt$my.category)))
  
  #p.selected <- ggplot() + geom_tile(sumamry.table.melt,aes(x=variable,y=my.category,fill=value),size=0.3,colour="black") # does not work (052816)
  p.selected <- ggplot(sumamry.table.melt,aes(x=variable,y=my.category)) + geom_tile(size=0.3,colour="black",aes(fill=value)) 
  
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


# overlapTable GOseq version (from 082515)
TIR10_cdna_rep_model<-readDNAStringSet("../input/TAIR10_cdna_20110103_representative_gene_model") 
head(TIR10_cdna_rep_model)
bias<-nchar(TIR10_cdna_rep_model)
names(bias)<-substr(names(TIR10_cdna_rep_model),1,9)

GOseq.overlapTable<-function (labels1, labels2, na.rm = TRUE, ignore = NULL, levels1 = NULL, 
                              levels2 = NULL,genes,permutation = 2000)  # genes is names of all genes in this analysis (usugally rownames of matrix). This considers cDNA length bias in DE.
{
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

# for module network
library(glmnet)
module.lasso<-function(data,module,timepoint,actual.value=F) { # timepoint (character); 
  if(timepoint=="1")  {
    data.subset<-data[data$time==timepoint,]
    data.subset$trt<-factor(data.subset$trt,levels=c("H","L"))
    data.subset$genotype<-as.factor(data.subset$genotype)
    #data.subset<-data.frame(lapply(data.subset,as.character), stringsAsFactors=TRUE)
    print(str(data.subset))
    fmla<-as.formula(paste(module,"~",paste(names(data.subset)[5:27][!names(data.subset)[5:27] %in% module],"*trt",collapse="+",sep=""),"+genotype*trt"))
    x2<-model.matrix(fmla,data=data.subset)[,-1]
    save(x2,file=paste(module,timepoint,"x2.Rdata",sep=""))
    print("model matrix is")
    print(x2)
    temp<-data.frame(actual=data.subset[,module],prediction="")
    grid =10^seq(10,-2, length =100)
    # lasso
    train=sample(1:nrow(x2),nrow(x2)/2)
    lasso.mod=glmnet(x2[train,],temp$actual[train],alpha=1,lambda=grid)
    cv.out=cv.glmnet(x2[train,],temp$actual[train],alpha=1)
    plot(cv.out)
    bestlam=cv.out$lambda.min
    out=glmnet(x2,temp$actual,alpha=1,lambda=grid)
    lasso.coef=predict(out,type="coefficients",s=bestlam)
    #lasso.coef=coef(out,s=bestlam)  # same as previous line. how to get coeff names?
    lasso.coef.vector<-as.vector(lasso.coef)
    names(lasso.coef.vector)<-c("Intercept",colnames(x2))
    if(actual.value==F) {return(lasso.coef.vector)} 
    else {
      temp.data<-data.subset[,c("genotype","trt","time","set",module)]
      temp.data$module<-rep(module,dim(temp.data)[1])
      names(temp.data)[names(temp.data)==module]<-"actual.value"
      temp.data.final<-temp.data[,c("genotype","trt","time","set","module","actual.value")]
      return(temp.data.final)}    
  }
  else { # other timepints (4, 25, 49 h)
    timecourse<-c("1","4","25","49")
    time.match<-match(timepoint,timecourse)
    data$time<-factor(data$time,levels=timecourse)
    data.subset<-data[data$time %in% timecourse[c(time.match-1, time.match)],]
    data.subset$trt<-factor(data.subset$trt,levels=c("H","L"))
    data.subset$genotype<-as.factor(data.subset$genotype)
    # needs to reshape 
    data.subset.reshape<-reshape(data.subset,idvar=c("genotype","trt","set"),timevar="time",direction="wide")
    # clean up data that has both 1h and 4h data
    data.subset.reshape.s<-data.subset.reshape[!rowSums(is.na(data.subset.reshape)*1)==23,]
    # 
    print(str(data.subset.reshape.s))
    fmla<-as.formula(paste(module,".",timepoint,"~",paste(names(data.subset.reshape.s)[4:49][!names(data.subset.reshape.s)[4:49] %in% paste(module,".",timepoint,sep="")],"*trt",collapse="+",sep=""),"+genotype*trt",sep=""))
    x2<-model.matrix(fmla,data=data.subset.reshape.s)[,-1]
    save(x2,file=paste(module,timepoint,"x2.Rdata",sep=""))
    print("model matrix is")
    print(colnames(x2))    
    temp<-data.frame(actual=data.subset.reshape.s[,paste(module,".",timepoint,sep="")],prediction="")
    grid =10^seq(10,-2, length =100)
    # lasso
    train=sample(1:nrow(x2),nrow(x2)/2)
    lasso.mod=glmnet(x2[train,],temp$actual[train],alpha=1,lambda=grid)
    cv.out=cv.glmnet(x2[train,],temp$actual[train],alpha=1)
    plot(cv.out)
    bestlam=cv.out$lambda.min
    out=glmnet(x2,temp$actual,alpha=1,lambda=grid)
    lasso.coef=predict(out,type="coefficients",s=bestlam)
    #lasso.coef=coef(out,s=bestlam)  # same as previous line. how to get coeff names?
    lasso.coef.vector<-as.vector(lasso.coef)
    names(lasso.coef.vector)<-c("Intercept",colnames(x2))
    if(actual.value==F) {return(lasso.coef.vector)} 
    else {
      # reformat data.subset.reshape.s 
      temp.data<-data.subset.reshape.s[,c("genotype","trt","set",paste(module,".",timepoint,sep=""))]
      temp.data$time<-rep(timepoint,dim(temp.data)[1])
      temp.data$module<-rep(module,dim(temp.data)[1])
      names(temp.data)[names(temp.data)==paste(module,".",timepoint,sep="")]<-"actual.value"
      temp.data.final<-temp.data[,c("genotype","trt","time","set","module","actual.value")]
      return(temp.data.final)
    }    
  }
}
# module network for subset data (sun or shade)
module.lasso2<-function(data,module,timepoint,actual.value=F,path="/Volumes/Data8/NGS_related/Arabidopsis_analysis/SAS_muts_time_course_RNAseq2_mycomp_output/") { # timepoint (character); 
  if(timepoint=="1")  {
    data.subset<-data[data$time==timepoint,]
    #data.subset$trt<-factor(data.subset$trt,levels=c("H","L"))
    data.subset$genotype<-as.factor(data.subset$genotype)
    #data.subset<-data.frame(lapply(data.subset,as.character), stringsAsFactors=TRUE)
    print(str(data.subset))
    fmla<-as.formula(paste(module,"~",paste(names(data.subset)[5:27][!names(data.subset)[5:27] %in% module],collapse="+",sep=""),"+genotype")) # 
    x3<-model.matrix(fmla,data=data.subset)[,-1]
    save(x3,file=paste(path,module,timepoint,data.subset$trt[1],"x3.Rdata",sep=""))
    print("model matrix is")
    print(x3)
    temp<-data.frame(actual=data.subset[,module],prediction="")
    grid =10^seq(10,-2, length =100)
    # lasso
    train=sample(1:nrow(x3),nrow(x3)/2)
    lasso.mod=glmnet(x3[train,],temp$actual[train],alpha=1,lambda=grid)
    cv.out=cv.glmnet(x3[train,],temp$actual[train],alpha=1,grouped=F)
    plot(cv.out)
    bestlam=cv.out$lambda.min
    out=glmnet(x3,temp$actual,alpha=1,lambda=grid)
    lasso.coef=predict(out,type="coefficients",s=bestlam)
    lasso.coef=coef(out,s=bestlam)  # same as previous line. how to get coeff names?
    lasso.coef.vector<-as.vector(lasso.coef)
    names(lasso.coef.vector)<-c("Intercept",colnames(x3))
    if(actual.value==F) {return(lasso.coef.vector)} 
    else {
      temp.data<-data.subset[,c("genotype","trt","time","set",module)]
      temp.data$module<-rep(module,dim(temp.data)[1])
      names(temp.data)[names(temp.data)==module]<-"actual.value"
      temp.data.final<-temp.data[,c("genotype","trt","time","set","module","actual.value")]
      return(temp.data.final)}    
  }
  else { # other timepints (4, 25, 49 h)
    timecourse<-c("1","4","25","49")
    time.match<-match(timepoint,timecourse)
    data$time<-factor(data$time,levels=timecourse)
    data.subset<-data[data$time %in% timecourse[c(time.match-1, time.match)],]
    #data.subset$trt<-factor(data.subset$trt,levels=c("H","L"))
    data.subset$genotype<-as.factor(data.subset$genotype)
    # needs to reshape 
    data.subset.reshape<-reshape(data.subset,idvar=c("genotype","trt","set"),timevar="time",direction="wide")
    # clean up data that has both 1h and 4h data
    data.subset.reshape.s<-data.subset.reshape[!rowSums(is.na(data.subset.reshape)*1)==23,]
    # 
    print(str(data.subset.reshape.s))
    fmla<-as.formula(paste(module,".",timepoint,"~",paste(names(data.subset.reshape.s)[4:49][!names(data.subset.reshape.s)[4:49] %in% paste(module,".",timepoint,sep="")],collapse="+",sep=""),"+genotype",sep=""))
    x3<-model.matrix(fmla,data=data.subset.reshape.s)[,-1]
    save(x3,file=paste(path,module,timepoint,data.subset$trt[1],"x3.Rdata",sep=""))
    print("model matrix is")
    print(colnames(x3))    
    temp<-data.frame(actual=data.subset.reshape.s[,paste(module,".",timepoint,sep="")],prediction="")
    grid =10^seq(10,-2, length =100)
    # lasso
    train=sample(1:nrow(x3),nrow(x3)/2)
    lasso.mod=glmnet(x3[train,],temp$actual[train],alpha=1,lambda=grid)
    cv.out=cv.glmnet(x3[train,],temp$actual[train],alpha=1,grouped=F)
    plot(cv.out)
    bestlam=cv.out$lambda.min
    out=glmnet(x3,temp$actual,alpha=1,lambda=grid)
    lasso.coef=predict(out,type="coefficients",s=bestlam)
    #lasso.coef=coef(out,s=bestlam)  # same as previous line. how to get coeff names?
    lasso.coef.vector<-as.vector(lasso.coef)
    names(lasso.coef.vector)<-c("Intercept",colnames(x3))
    if(actual.value==F) {return(lasso.coef.vector)} 
    else {
      # reformat data.subset.reshape.s 
      temp.data<-data.subset.reshape.s[,c("genotype","trt","set",paste(module,".",timepoint,sep=""))]
      temp.data$time<-rep(timepoint,dim(temp.data)[1])
      temp.data$module<-rep(module,dim(temp.data)[1])
      names(temp.data)[names(temp.data)==paste(module,".",timepoint,sep="")]<-"actual.value"
      temp.data.final<-temp.data[,c("genotype","trt","time","set","module","actual.value")]
      return(temp.data.final)
    }    
  }
}
# omit genothpe factors in model
module.lasso3<-function(data,module,timepoint,actual.value=F,model.matrix=F,path="/Volumes/Data8/NGS_related/Arabidopsis_analysis/SAS_muts_time_course_RNAseq2_mycomp_output/") 
  { # timepoint (character); 
  if(timepoint=="1")  {
    data.subset<-data[data$time==timepoint,]
    #data.subset$trt<-factor(data.subset$trt,levels=c("H","L"))
    data.subset$genotype<-as.factor(data.subset$genotype)
    #data.subset<-data.frame(lapply(data.subset,as.character), stringsAsFactors=TRUE)
    print(str(data.subset))
    fmla<-as.formula(paste(module,"~",paste(names(data.subset)[5:27][!names(data.subset)[5:27] %in% module],collapse="+",sep=""))) # 
    x3<-model.matrix(fmla,data=data.subset)[,-1]    
    if(model.matrix==T & actual.value==T) {
      temp.data<-data.subset[,c("genotype","trt","time","set",module)]
      temp.data$module<-rep(module,dim(temp.data)[1])
      names(temp.data)[names(temp.data)==module]<-"actual.value"
      temp.data.final<-temp.data[,c("genotype","trt","time","set","module","actual.value")]    
      save(x3,file=paste(path,module,timepoint,data.subset$trt[1],"x3.Rdata",sep=""))
      return(temp.data.final)
      #stop("only need x3")
        } 
    else {
    print("model matrix is")
    print(x3)
    temp<-data.frame(actual=data.subset[,module],prediction="")
    grid =10^seq(10,-2, length =100)
    # lasso
    train=sample(1:nrow(x3),nrow(x3)/2)
    lasso.mod=glmnet(x3[train,],temp$actual[train],alpha=1,lambda=grid)
    cv.out=cv.glmnet(x3[train,],temp$actual[train],alpha=1,grouped=F)
    plot(cv.out)
    bestlam=cv.out$lambda.min
    out=glmnet(x3,temp$actual,alpha=1,lambda=grid)
    lasso.coef=predict(out,type="coefficients",s=bestlam)
    lasso.coef=coef(out,s=bestlam)  # same as previous line. how to get coeff names?
    lasso.coef.vector<-as.vector(lasso.coef)
    names(lasso.coef.vector)<-c("Intercept",colnames(x3))
    return(lasso.coef.vector)
    } 
  }
  else { # other timepints (4, 25, 49 h)
    timecourse<-c("1","4","25","49")
    time.match<-match(timepoint,timecourse)
    data$time<-factor(data$time,levels=timecourse)
    data.subset<-data[data$time %in% timecourse[c(time.match-1, time.match)],]
    # get reid of genotype information
    
    #data.subset$trt<-factor(data.subset$trt,levels=c("H","L"))
    data.subset$genotype<-as.factor(data.subset$genotype)
    # needs to reshape 
    data.subset.reshape<-reshape(data.subset,idvar=c("genotype","trt","set"),timevar="time",direction="wide")
    # clean up data that has both 1h and 4h data
    data.subset.reshape.s<-data.subset.reshape[!rowSums(is.na(data.subset.reshape)*1)==23,]
    # 
    print(str(data.subset.reshape.s))
    fmla<-as.formula(paste(module,".",timepoint,"~",paste(names(data.subset.reshape.s)[4:49][!names(data.subset.reshape.s)[4:49] %in% paste(module,".",timepoint,sep="")],collapse="+",sep=""),sep=""))
    x3<-model.matrix(fmla,data=data.subset.reshape.s)[,-1]
      if(model.matrix==T & actual.value==T) {
        # reformat data.subset.reshape.s 
        temp.data<-data.subset.reshape.s[,c("genotype","trt","set",paste(module,".",timepoint,sep=""))]
        temp.data$time<-rep(timepoint,dim(temp.data)[1])
        temp.data$module<-rep(module,dim(temp.data)[1])
        names(temp.data)[names(temp.data)==paste(module,".",timepoint,sep="")]<-"actual.value"
        temp.data.final<-temp.data[,c("genotype","trt","time","set","module","actual.value")]
        save(x3,file=paste(path,module,timepoint,data.subset$trt[1],"x3.Rdata",sep=""))
        return(temp.data.final)
      } 
    else { 
      print("model matrix is")
      print(colnames(x3))    
      temp<-data.frame(actual=data.subset.reshape.s[,paste(module,".",timepoint,sep="")],prediction="")
      grid =10^seq(10,-2, length =100)
      # lasso
      train=sample(1:nrow(x3),nrow(x3)/2)
      lasso.mod=glmnet(x3[train,],temp$actual[train],alpha=1,lambda=grid)
      cv.out=cv.glmnet(x3[train,],temp$actual[train],alpha=1,grouped=F)
      plot(cv.out)
      bestlam=cv.out$lambda.min
      out=glmnet(x3,temp$actual,alpha=1,lambda=grid)
      lasso.coef=predict(out,type="coefficients",s=bestlam)
      #lasso.coef=coef(out,s=bestlam)  # same as previous line. how to get coeff names?
      lasso.coef.vector<-as.vector(lasso.coef)
      names(lasso.coef.vector)<-c("Intercept",colnames(x3))
      return(lasso.coef.vector)
    }
  } # else for timepint==1
} # close function

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


