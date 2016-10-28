#setwd("~/Dropbox/Kazu's SAS timecourse data/")

setwd("../../Nozue2016_SAStranscriptome_data/input/")
data <- read.delim("SAS_counts_merged_updated_final.bam.tsv",sep="\t",row.names=2)

data[1:10,1:10]
# eliminated.libraries<-c(
#   "D2H16A.merged.bam", # AT5G02540_2
#   "D4H4A.merged.bam", # jar1-1
#   "D4H16A.merged.bam", # jar1-1
#   "D4H49A.merged.bam", # jar1-1
#   "C6H1A.merged.bam", # phyB
#   "C6H16A.merged.bam", # phyB
#   "C7L4A.merged.bam", # pif4pi5
#   "C7H16A.merged.bam", # pif4pif5
#   "C7H49A.merged.bam", # pif4pif5
#   "C7L49A.merged.bam", # pif4pif5
#   "C11H1A.merged.bam", # coi1_16
#   "C11L1A.merged.bam", # coi1_16
#   "C11H4A.merged.bam", # coi1_16
#   "D11H4A.merged.bam", # coi1_16
#   "C11L4A.merged.bam", # coi1_16
#   "C11H16A.merged.bam", # coi1_16
#   "C11L16A.merged.bam", # coi1_16
#   "C11H25A.merged.bam", # coi1_16
#   "C11L25A.merged.bam", # coi1_16
#   "C11H49A.merged.bam", # coi1_16
#   "D11H49A.merged.bam", # coi1_16
#   "C11L49A.merged.bam", # coi1_16
#   "C12H1A.merged.bam", # phyAB
#   "C12L1A.merged.bam", # phyAB
#   "D12L1A.merged.bam", # phyAB
#   "C12H4A.merged.bam", # phyAB
#   "C12L4A.merged.bam", # phyAB
#   "C12H16A.merged.bam", # phyAB
#   "C12L16A.merged.bam", # phyAB
#   "C12L25A.merged.bam", # phyAB
#   "C12H49A.merged.bam", # phyAB
#   "C12L49A.merged.bam", # phyAB
#   "C13H1A.merged.bam", # pif3_1
#   "C13H16A.merged.bam", # pif3_1
#   "E13H16A.merged.bam", # pif3_1
#   "C13L16A.merged.bam", # pif3_1
#   "E13L16A.merged.bam", # pif3_1
#   "C13H25A.merged.bam", # pif3_1
#   "C13L25A.merged.bam", # pif3_1
#   "C13H49A.merged.bam", # pif3_1
#   "D13H49A.merged.bam", # pif3_1
#   "C13L49A.merged.bam", # pif3_1
#   "D14L4A.merged.bam", # mida9_4
#   "C14H25A.merged.bam", # mida9_4
#   "D14L25A.merged.bam", # mida9_4
#   "C16H1A.merged.bam", # sto
#   "C16L1A.merged.bam", # sto
#   "C16H4A.merged.bam", # sto
#   "C16L4A.merged.bam", # sto
#   "C16H16A.merged.bam", # sto
#   "C16H25A.merged.bam", # sto
#   "C16L25A.merged.bam", # sto
#   "C16H49A.merged.bam", # sto
#   "C16L49A.merged.bam", # sto
#   "D16L49A.merged.bam") # sto
# # eliminate 55 libraries within my genotypes (081514)
# 
# # conversion.table<-data.frame(num=1:27, genotype=c("Col","AT5G02540_1","hy5","jar1",
# #                                                   "kat1_2","phyB","pif45","spt_11","yuc2589",
# #                                                   "PAR1_RNAi09","coi1_16","phyAB","pif3",
# #                                                   "AT5G02760_2","AT5G59010_2","sto","aos","argos",
# #                                                   "co_9","Blh_1","Jea","Shahdara","Col_0","Cvi_0",
# #                                                   "Bur_0","Oy_0","Ita_0"))
# # remove non PLoS Genetics paper (052616; revised 071516)
# librarynames<-colnames(data)
# PLoSGenetics2015mutnats<-!gsub("(C|D|E)([[:digit:]]+)(H|L)(1|4|16|25|49)(A.merged.bam)","\\2",librarynames) %in% c("2","15","18")
# # also remove 16h timepoints (from MDS plot after edgeR)
# No_sixteen_hours<-!(gsub("(C|D|E)([[:digit:]]+)(H|L)(1|4|16|25|49)(A.merged.bam)","\\4",librarynames) =="16")
# data<-data[,PLoSGenetics2015mutnats&No_sixteen_hours]
# 
# #
# data <- data[,!colnames(data) %in% eliminated.libraries]
# 
if(rownames(data)[1]=="*" & colnames(data)[1]=="row.names") data <- data[-1,-1] #get rid of unused first column and unmapped reads

head(data[1:10])

data[is.na(data)] <- 0

# samples <- names(data)
# 
# sample.info <- data.frame(sample=samples,
#                           rep=substr(samples,1,1),
#                           gt=sub("[A-Z]([0-9]{1,2})[A-Z].*","\\1",samples),
#                           trt=sub("[A-Z][0-9]{1,2}([A-Z]).*","\\1",samples),
#                           time=sub("[A-Z][0-9]{1,2}[A-Z]([0-9]{1,2}).*","\\1",samples)
# )
# 
# write.csv(sample.info,"sample.info.csv")
# ftable(sample.info,row.vars=c(2:3),col.vars=c(4:5)) # gt=1 is strange only 16h samples, why?

### new version (070516) ##
load("../../Nozue2016_SAStranscriptome_output/output/samples.nolow5.Rdata")
samples.nolow5$genotype2<-""
for(i in 1:27){
  samples.nolow5[samples.nolow5$genotype==conversion.table[i,1],"genotype2"]<-as.character(conversion.table[i,2])  
} 
# all libraries > eliminate low count libraries > eliminate libraries with wrong genotypes > eliminate outliers based on DE (sun vs shade) > samples.nolow5
## double check if this is correct. (092116)
PLoSGenetics2015mutants<-samples.nolow5[samples.nolow5$genotype2 %in% c("Col","aos","co_9","coi1-16","jar1","kat1_2","hy5","mida9_4","PAR1_RNAi09","phyAB","phyB","pif3","pif45","spt_11","sto","yuc2589"),]
# genotypes not working in edgeR 
errorDE<-c("coi1-16","phyAB","sto")
PLoSGenetics2015mutants.errrDE<-PLoSGenetics2015mutants[!PLoSGenetics2015mutants$genotype2 %in% errorDE,]
save(PLoSGenetics2015mutants.errrDE,file="../../Nozue2016_SAStranscriptome_output/output/PLoSGenetics2015mutants.errrDE.Rdata")
data<-data[,colnames(data) %in% PLoSGenetics2015mutants.errrDE$file]

library(DESeq2)
PLoSGenetics2015mutants.errrDE$batch

dds <- DESeqDataSetFromMatrix(countData = data, colData = PLoSGenetics2015mutants.errrDE, design = ~ batch + genotype*trt*time) # modified (070616)
vsd <- varianceStabilizingTransformation(dds)

vstMat <- assay(vsd)
colnames(vstMat) <- colnames(data)
vstMat[1:10,1:10]
#write.csv(vstMat,file="../input/kazu.SAS.expression.VST.070516.csv")
write.csv(vstMat,file="../input/kazu.SAS.expression.VST.082016.csv")

log2(data[1:10,1:10]+1)
# reality check: highest correlations should be of the sample to itself
cor(vstMat[,1:10],log2(data[,1:10]+1),use="complete.obs") # yes
#write.csv(vstMat,"kazu.SAS.expression.VST.112614.csv")
## prune low expressed genes first?

data.small <- data[rowSums(data > 5) > ncol(data)/2,] #more than 5 counts in more than half the samples

dim(data.small)

dim(data)

dds.small <- DESeqDataSetFromMatrix(countData = data.small, colData = sample.info, design = ~ rep + gt*trt*time)

system.time(vsd.small <- varianceStabilizingTransformation(dds.small))

vstMat.small <- assay(vsd.small)

colnames(vstMat.small) <- colnames(data.small)

#write.csv(vstMat.small,"kazu.SAS.expression.VST.SMALL.112614.csv")
setwd("../output/")
write.csv(vstMat.small,"kazu.SAS.expression.VST.SMALL.052616.csv")

# compare old (112614) vs new (052616)
SAS.expression.vst<-read.csv("../input/kazu.SAS.expression.VST.SMALL.112614.csv")
rownames(SAS.expression.vst)<-SAS.expression.vst[,1]
SAS.expression.vst.s<-SAS.expression.vst[,-1] # should be numeric for cor()
# the same columns
SAS.expression.vst.s<-SAS.expression.vst.s[,colnames(SAS.expression.vst.s) %in% colnames(vstMat.small)]
vstMat.small.s<-vstMat.small[,colnames(vstMat.small) %in% colnames(SAS.expression.vst.s)]
# the same genes
SAS.expression.vst.ss<-SAS.expression.vst.s[rownames(SAS.expression.vst.s) %in% rownames(vstMat.small.s),]
vstMat.small.ss<-vstMat.small.s[rownames(vstMat.small.s) %in% rownames(SAS.expression.vst.s),]

X<-sapply(1:ncol(SAS.expression.vst.ss),function(i) cor(SAS.expression.vst.ss[,i],vstMat.small.ss[,i]))

min(X) # 0.9921697


