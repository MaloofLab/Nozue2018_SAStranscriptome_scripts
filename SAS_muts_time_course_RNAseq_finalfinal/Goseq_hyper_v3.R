
#### enrichment analysis with GOseq analysis
###### origianlly from Matthew young (GOseq developer)
# ORA with GOseq
# prerequisit
library(ShortRead);library(goseq);library(GO.db);library("org.At.tair.db");library("annotate")

TIR10_cdna_rep_model<-readDNAStringSet("/Volumes/Data8/NGS_related/Arabidopsis_analysis/reference/TAIR10_cdna_20110103_representative_gene_model") 
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
# example
#enriched.GO.Col.SOMcluster1<-GOseq.ORA(rownames(data.val3.4.SOM.Col.SOM1.all.barcode.s)[rownames(data.val3.4.SOM.Col.SOM1.all.barcode.s) %in% names(bias)]) 

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

#######################

## regular GO seq analysis for Arabidopsis
library(ShortRead)
TIR10_cdna_rep_model<-readDNAStringSet("/Volumes/Data8/NGS_related/Arabidopsis_analysis/reference/TAIR10_cdna_20110103_representative_gene_model") 
head(TIR10_cdna_rep_model)
####################
bias<-nchar(TIR10_cdna_rep_model) 
names(bias)<-substr(names(TIR10_cdna_rep_model),1,9)
length(bias) 
length(genes) # different from length of bias...

#  bias.data vector must have the same length as DEgenes vector!
###
gene.name<-substr(names(TIR10_cdna_rep_model),1,9)
bias2<-bias[names(bias) %in% names(genes)] # genes came from differential expression analysis (from Plant Mart 7)
#summary(names(bias) %in% names(genes))
summary(names(bias2) %in% names(genes))
## 
library(goseq)


# (2) extract genes in enriched category
###Read in AtGO
Atgo <- toTable(org.At.tairGO)
#head(Atgo)
BP <- TRUE #only keep BP go TERMS
if (BP) Atgo <- Atgo[Atgo$Ontology=="BP",]
#convert to list
Atgo.list <- tapply(Atgo$go_id,Atgo$gene_id,c)

# adding gene description
TIR10_gene_descriptions<-read.csv("/Volumes/Data8/NGS_related/Arabidopsis_analysis/reference/TAIR10_functional_descriptions.csv") 
head(TIR10_gene_descriptions)
TIR10_gene_descriptions$AGI<-gsub("([[:print:]]+)(.[[:digit:]]+)","\\1",TIR10_gene_descriptions$Model_name)
# making custom made category list that is compatible to GOseq
# integrate all hormoe responsiveness in one data
#making empty dataframe called "dummy.cat2" # needs "bias" object
setwd("/Volumes/Data8/NGS_related/Arabidopsis_analysis")
group_comparison2014v1<-read.csv("group_comparison2014v1.csv") # Nemhauser papers plus SA (see above)

# target categories
categories<-c("IAAup","IAAdown","BLup","BLdown","GAinga1.3UP","GAinga1.3DOWN","BLup","BLdown","MJup","MJdown","CKup","CKdown","SA_UP","SA_DOWN",
                "SARplus","SARminus","GRF1_3target","PIFtarget","MYC234_up","MYC234_down","NPR1_up","NPR1_down","WRKY33up","WRKY33down")
# rename categories
group_comparison2014v1.selected<-group_comparison2014v1[,names(group_comparison2014v1[-c(1:3),]) %in% categories]
names(group_comparison2014v1.selected)<-c("BLup","BLdown","IAAup","IAAdown","MJup","MJdown","CKup","CKdown","GAup","GAdown",
                                          "SAup","SAdown","SARminus","SARplus","GRF1_3target","PIFtarget","MYC234up","MYC234down","NPR1up","NPR1down","WRKY33up","WRKY33down")
#dummy.cat2<-data.frame(BLup=rep("",length(bias)),row.names=names(bias))
dummy.cat2<-as.data.frame(matrix(nrow=length(names(bias)),ncol=length(names(group_comparison2014v1.selected))))
dummy.cat2<-as.data.frame(matrix(nrow=length(names(bias)),ncol=1))
rownames(dummy.cat2)<-names(bias)
for(i in names(group_comparison2014v1.selected)) {# 
  print(i)
  x.locus<-group_comparison2014v1.selected[-c(1:3),names(group_comparison2014v1.selected)==i]
  x.locus<-x.locus[!x.locus==""]  
  dummy.cat2[rownames(dummy.cat2) %in% toupper(x.locus),i]<-rep(i,sum(as.integer(rownames(dummy.cat2) %in% toupper(x.locus))))
}
head(dummy.cat2)
dummy.cat2<-dummy.cat2[,-1]
# clean up <NA>
dummy.cat2[is.na(dummy.cat2)]<-""
### convert into list
temp<-list()
for(i in 1:dim(dummy.cat2)[1]) {
  temp[[i]]<-paste(dummy.cat2[i,])
}
names(temp)<-rownames(dummy.cat2)   
hormone.responsiveness2<-temp
head(hormone.responsiveness2)
table(hormone.responsiveness2[1:10])
setwd("/Volumes/Data8/NGS_related/Arabidopsis_analysis/")
save(hormone.responsiveness2,file="hormone.responsiveness2.Rdata")

# add ABA and ACC (022415)
setwd("/Volumes/Data8/NGS_related/Arabidopsis_analysis")
group_comparison2014v1<-read.csv("group_comparison2014v1.csv") # Nemhauser papers plus SA (see above)

# target categories
categories<-c("IAAup","IAAdown","BLup","BLdown","GAinga1.3UP","GAinga1.3DOWN","BLup","BLdown","MJup","MJdown","ABAup","ABAdown","ACCup","ACCdown","CKup","CKdown","SA_UP","SA_DOWN",
              "SARplus","SARminus","GRF1_3target","PIFtarget","MYC234_up","MYC234_down","NPR1_up","NPR1_down","WRKY33up","WRKY33down")
# rename categories
group_comparison2014v1.selected<-group_comparison2014v1[,names(group_comparison2014v1[-c(1:3),]) %in% categories]
names(group_comparison2014v1.selected)<-c("ABAup","ABAdown","ACCup","ACCdown","BLup","BLdown","IAAup","IAAdown","MJup","MJdown","CKup","CKdown","GAup","GAdown",
                                          "SAup","SAdown","SARminus","SARplus","GRF1_3target","PIFtarget","MYC234up","MYC234down","NPR1up","NPR1down","WRKY33up","WRKY33down")
#dummy.cat2<-data.frame(BLup=rep("",length(bias)),row.names=names(bias))
dummy.cat2<-as.data.frame(matrix(nrow=length(names(bias)),ncol=length(names(group_comparison2014v1.selected))))
dummy.cat2<-as.data.frame(matrix(nrow=length(names(bias)),ncol=1))
rownames(dummy.cat2)<-names(bias)
for(i in names(group_comparison2014v1.selected)) {# 
  print(i)
  x.locus<-group_comparison2014v1.selected[-c(1:3),names(group_comparison2014v1.selected)==i]
  x.locus<-x.locus[!x.locus==""]  
  dummy.cat2[rownames(dummy.cat2) %in% toupper(x.locus),i]<-rep(i,sum(as.integer(rownames(dummy.cat2) %in% toupper(x.locus))))
}
head(dummy.cat2)
dummy.cat2<-dummy.cat2[,-1]
# clean up <NA>
dummy.cat2[is.na(dummy.cat2)]<-""
### convert into list
temp<-list()
for(i in 1:dim(dummy.cat2)[1]) {
  temp[[i]]<-paste(dummy.cat2[i,])
}
names(temp)<-rownames(dummy.cat2)   
hormone.responsiveness3<-temp
head(hormone.responsiveness3)
table(hormone.responsiveness3[1:10])
setwd("/Volumes/Data8/NGS_related/Arabidopsis_analysis/")
save(hormone.responsiveness3,file="hormone.responsiveness3.Rdata")



# loading custom made categories
#load("/Volumes/Data8/NGS_related/Arabidopsis_analysis/hormone.responsiveness.Rdata") #old version. check if new version is identical to this
load("/Volumes/Data8/NGS_related/Arabidopsis_analysis/hormone.responsiveness2.Rdata") #new version. check if new version is identical to this

# version 4
# add SA regulated genes(060116)

setwd("/Volumes/Data8/NGS_related/Arabidopsis_analysis")
group_comparison2016v1<-read.csv("group_comparison2016v1.csv") 

# target categories
categories<-c("IAAup","IAAdown","BLup","BLdown","GAinga1.3UP","GAinga1.3DOWN","BLup","BLdown","MJup","MJdown","ABAup","ABAdown","ACCup","ACCdown","CKup","CKdown","SA_UP","SA_DOWN",
              "SARplus","SARminus","GRF1_3target","PIFtarget","MYC234_up","MYC234_down","NPR1_up","NPR1_down","WRKY33up","WRKY33down","PCC1up","PCC1down")
# rename categories
group_comparison2016v1.selected<-group_comparison2016v1[,names(group_comparison2016v1[-c(1:3),]) %in% categories]
names(group_comparison2016v1.selected)<-c("ABAup","ABAdown","ACCup","ACCdown","BLup","BLdown","IAAup","IAAdown","MJup","MJdown","CKup","CKdown","GAup","GAdown",
                                          "SAup","SAdown","SARminus","SARplus","GRF1_3target","PIFtarget","MYC234up","MYC234down","NPR1up","NPR1down","WRKY33up","WRKY33down","PCC1up","PCC1down")
#dummy.cat2<-data.frame(BLup=rep("",length(bias)),row.names=names(bias))
dummy.cat2<-as.data.frame(matrix(nrow=length(names(bias)),ncol=length(names(group_comparison2016v1.selected))))
dummy.cat2<-as.data.frame(matrix(nrow=length(names(bias)),ncol=1))
rownames(dummy.cat2)<-names(bias) # use "bias" object from function_RNAseq_time_course.R 
for(i in names(group_comparison2016v1.selected)) {# 
  print(i)
  x.locus<-group_comparison2016v1.selected[-c(1:3),names(group_comparison2016v1.selected)==i]
  # clean up x.locus
  x.locus<-x.locus[!x.locus==""]  
  # only single AGI name is extracted (no probe naem, no multiple AGI name, such as "AT1G23900 AT1G23940") 
  dummy.cat2[rownames(dummy.cat2) %in% toupper(x.locus),i]<-rep(i,sum(as.integer(rownames(dummy.cat2) %in% toupper(x.locus))))
}
head(dummy.cat2)
dummy.cat2<-dummy.cat2[,-1]
# clean up <NA>
dummy.cat2[is.na(dummy.cat2)]<-""
# 
### convert into list
temp<-list()
for(i in 1:dim(dummy.cat2)[1]) {
  temp[[i]]<-paste(dummy.cat2[i,])
}
names(temp)<-rownames(dummy.cat2)   
hormone.responsiveness4<-temp
head(hormone.responsiveness4)
table(hormone.responsiveness4[1:10])
attributes(hormone.responsiveness4)$category<-names(group_comparison2016v1.selected)
setwd("/Volumes/Data8/NGS_related/Arabidopsis_analysis/")
save(hormone.responsiveness4,file="hormone.responsiveness4.Rdata")

# version 5
# add more Sight (2015) data: NPR1 regulated gene upon SA treatment (24h), 24h SA responsive genes (Col)
# 
setwd("/Volumes/Data8/NGS_related/Arabidopsis_analysis") 
setwd("/Volumes/Data8new/Data8/NGS_related/Arabidopsis_analysis") # for temporal system

group_comparison.Ver5<-read.csv("group_comparison.Ver5.csv") 

# target categories
categories<-c("IAAup","IAAdown","BLup","BLdown","GAinga1.3UP","GAinga1.3DOWN","BLup","BLdown","MJup","MJdown","ABAup","ABAdown","ACCup","ACCdown","CKup","CKdown","SA_UP","SA_DOWN",
              "SARplus","SARminus","GRF1_3target","PIFtarget","MYC234_up","MYC234_down","NPR1_up","NPR1_down","WRKY33up","WRKY33down","PCC1up","PCC1down","SAup_24h","SAdown_24h","NPR1up_uponSA","NPR1down_uponSA")
# rename categories
group_comparison.Ver5.selected<-group_comparison.Ver5[,names(group_comparison.Ver5[-c(1:3),]) %in% categories]
names(group_comparison.Ver5.selected)<-c("ABAup","ABAdown","ACCup","ACCdown","BLup","BLdown","IAAup","IAAdown","MJup","MJdown","CKup","CKdown","GAup","GAdown",
                                          "SAup_2h","SAdown_2h","SARminus","SARplus","GRF1_3target","PIFtarget","MYC234up","MYC234down","NPR1up","NPR1down","WRKY33up","WRKY33down","PCC1up","PCC1down","SAup_24h","SAdown_24h","NPR1up_uponSA","NPR1down_uponSA")
#dummy.cat2<-data.frame(BLup=rep("",length(bias)),row.names=names(bias))
dummy.cat2<-as.data.frame(matrix(nrow=length(names(bias)),ncol=length(names(group_comparison.Ver5.selected))))
dummy.cat2<-as.data.frame(matrix(nrow=length(names(bias)),ncol=1))
rownames(dummy.cat2)<-names(bias) # use "bias" object from function_RNAseq_time_course.R 
for(i in names(group_comparison.Ver5.selected)) {# 
  print(i)
  x.locus<-group_comparison.Ver5.selected[-c(1:3),names(group_comparison.Ver5.selected)==i]
  # clean up x.locus
  x.locus<-x.locus[!x.locus==""]  
  # only single AGI name is extracted (no probe naem, no multiple AGI name, such as "AT1G23900 AT1G23940") 
  dummy.cat2[rownames(dummy.cat2) %in% toupper(x.locus),i]<-rep(i,sum(as.integer(rownames(dummy.cat2) %in% toupper(x.locus))))
}
head(dummy.cat2)
dummy.cat2<-dummy.cat2[,-1]
# clean up <NA>
dummy.cat2[is.na(dummy.cat2)]<-""
# 
### convert into list
temp<-list()
for(i in 1:dim(dummy.cat2)[1]) {
  temp[[i]]<-paste(dummy.cat2[i,])
}
names(temp)<-rownames(dummy.cat2)   
hormone.responsiveness5<-temp
head(hormone.responsiveness5)
table(hormone.responsiveness5[1:10])
attributes(hormone.responsiveness5)$category<-names(group_comparison.Ver5.selected)
setwd(homedir)
setwd("../../Nozue2016_SAStranscriptome_data/input")
save(hormone.responsiveness5,file="hormone.responsiveness5.Rdata")

