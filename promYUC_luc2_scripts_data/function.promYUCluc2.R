#####
##
## functions for promYUC::luc2 exp 
##
#################
format.luc2data<-function(plant.name,luc2.data,files) { 
#names(sample.num)<-c("shade","sun")
# extract genotype and sample number, line name
genotype<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)","\\5",plant.name)
sample.num<-c(length(grep("shade",gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)","\\1",plant.name))),length(grep("sun",gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)","\\1",plant.name))))
names(sample.num)<-c("shade","sun")
line<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)","\\7",plant.name)
# format luc2.data
set.files<-levels(as.factor(gsub("([[:print:]]+)(exp[[:digit:]]+)","\\2",luc2.data$Label))) # and needs to eliminate "set2.tif:exp5label"
set.files<-grep("label",set.files,invert=TRUE,value=TRUE)
luc2.data$treatment<-rep(rep(c("shade","sun"),c(sample.num["shade"],sample.num["sun"])),length(set.files)) # this is wrong?
luc2.data$promoter<-rep(genotype,length(set.files))
luc2.data$line<-rep(line,length(set.files))

# add time point info (under construction)
print(files)
files$date.time<-strptime(paste(files$V6,files$V7,files$V8,sep="/"),"%b/%d/%H:%M")
print(files)

luc2.data$date.time<-strptime("","%b/%d/%H:%M")
luc2.data$exp.file<-gsub("([[:print:]]+)(exp[[:digit:]]+)","\\2",luc2.data$Label)
for(i in 1:dim(luc2.data)[1]) {
luc2.data[i,"date.time"]<-as.character(files[gsub("(exp[[:digit:]]+)(.tif)","\\1",files$V9)==luc2.data$exp.file[i],"date.time"])
}

#luc2.data$ZT<-luc2.data$date.time - sort(luc2.data$date.time)[1] # this is sec. how to change into hr?
#as.difftime(c(luc2.data$date.time,rep(sort(luc2.data$date.time)[1],dim(luc2.data$date.time)[1]), units = "hours",format="%b/%d/%H:%M")
luc2.data$ZT<-as.numeric(luc2.data$date.time - sort(luc2.data$date.time)[1], units = "hours",format="%b/%d/%H:%M")	
# reshape? and substract values in before shade treatment
#luc2.data.reshape<-reshape(***********) # under construction

	return(luc2.data)
}

# functin for "multi" output (under construction)
format.luc2data.multi<-function(plant.name,luc2.data,luc2.data.multi,files,first_day_ZT0) { 
	temp.multi<-luc2.data.multi[,c("Slice",grep("^IntDen",names(luc2.data.multi),value=TRUE))]
	# add file name
	set.files<-unique(gsub("([[:print:]]+)(exp[[:digit:]]+)","\\2",luc2.data$Label)) # keeping file order (by time)
	#print(set.files)
	# remove label file
	temp.multi$setfiles<-set.files
	temp.multi<-temp.multi[grep("lab",temp.multi$setfiles,invert=TRUE),]

	# calculate increase/decrease from begining
#	temp.multi.dif<-as.matrix(temp.multi[,-1])-as.vector(temp.multi[1,-1]) # does not work
	temp.multi.dif<-sweep(data.matrix(temp.multi[,grep("IntDen",names(temp.multi))]),2,data.matrix(temp.multi[1,grep("IntDen",names(temp.multi))]))
	temp.multi.dif<-as.data.frame(temp.multi.dif)
	temp.multi.dif<-cbind(temp.multi.dif,temp.multi$setfiles)
	names(temp.multi.dif)[grep("setfiles",names(temp.multi.dif))]<-"setfiles"
	#print(temp.multi.dif)

	# add sample names
	names(temp.multi.dif)[grep("IntDen",names(temp.multi.dif))]<-plant.name
	temp.multi.dif
	print(temp.multi.dif) # Error: id variables not found in data: setfiles
	# melt
	temp.multi.dif.melt<-melt(temp.multi.dif,id.var="setfiles")

	# treatment
	temp.multi.dif.melt$treatment<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\1",temp.multi.dif.melt$variable)
	
	# promoter
	temp.multi.dif.melt$promoter<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\5",temp.multi.dif.melt$variable)
	
	# line
	temp.multi.dif.melt$line<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\7",temp.multi.dif.melt$variable)	
	
# add time point info (under construction)
#print(files)
files$date.time<-strptime(paste(files$V6,files$V7,files$V8,sep="/"),"%b/%d/%H:%M")
#print(files)

temp.multi.dif.melt$date.time<-strptime("","%b/%d/%H:%M")
for(i in 1:dim(temp.multi.dif.melt)[1]) {
temp.multi.dif.melt[i,"date.time"]<-as.character(files[gsub("(exp[[:digit:]]+)(.tif)","\\1",files$V9)==temp.multi.dif.melt$setfiles[i],"date.time"])
}

temp.multi.dif.melt$ZT<-as.numeric(temp.multi.dif.melt$date.time - strptime(first_day_ZT0,"%b/%d/%H:%M"), units = "hours",format="%b/%d/%H:%M") 
	return(temp.multi.dif.melt)
}

# function exp3set1 special
# functin for "multi" output (under construction)
format.luc2data.multi.exp3set1special<-function(plant.name,luc2.data,luc2.data.multi,files,first_day_ZT0) { 
	temp.multi<-luc2.data.multi[,c("Slice",grep("^IntDen",names(luc2.data.multi),value=TRUE))]
	# add file name
	set.files<-unique(gsub("([[:print:]]+)(exp[[:digit:]]+)","\\2",luc2.data$Label)) # keeping file order (by time)
	#print(set.files)
	# remove label file
	temp.multi$setfiles<-set.files 
	temp.multi<-temp.multi[grep("lab",temp.multi$setfiles,invert=TRUE),]
	write.csv(temp.multi, file="temp.multi.csv") # and then edit files
	#temp.multi.mod<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/exp3/temp.multi2.csv")
	temp.multi.mod<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/raw_data/exp3/temp.multi2.csv")
	
	# calculate increase/decrease from begining
	# temp.multi.dif<-as.matrix(temp.multi[,-1])-as.vector(temp.multi[1,-1]) # does not work
	temp.multi.dif<-sweep(data.matrix(temp.multi.mod[,grep("IntDen",names(temp.multi.mod))]),2,data.matrix(temp.multi.mod[1,grep("IntDen",names(temp.multi.mod))]))	
	temp.multi.dif<-as.data.frame(temp.multi.dif)
	temp.multi.dif<-cbind(temp.multi.dif,temp.multi.mod$setfiles)
	names(temp.multi.dif)[grep("setfiles",names(temp.multi.dif))]<-"setfiles"
	print(temp.multi.dif)

	# add sample names
	names(temp.multi.dif)[grep("IntDen",names(temp.multi.dif))]<-plant.name
	temp.multi.dif
	print(temp.multi.dif) # Error: id variables not found in data: setfiles
	# melt
	temp.multi.dif.melt<-melt(temp.multi.dif,id.var="setfiles")
	temp.multi.dif.melt[is.na(temp.multi.dif.melt$value),"value"]<-0
	# treatment
	temp.multi.dif.melt$treatment<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\1",temp.multi.dif.melt$variable)
	
	# promoter
	temp.multi.dif.melt$promoter<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\5",temp.multi.dif.melt$variable)
	
	# line
	temp.multi.dif.melt$line<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\7",temp.multi.dif.melt$variable)	
	
# add time point info (under construction)
#print(files)
files$date.time<-strptime(paste(files$V6,files$V7,files$V8,sep="/"),"%b/%d/%H:%M")
#print(files)

temp.multi.dif.melt$date.time<-strptime("","%b/%d/%H:%M")
for(i in 1:dim(temp.multi.dif.melt)[1]) {
temp.multi.dif.melt[i,"date.time"]<-as.character(files[gsub("(exp[[:digit:]]+)(.tif)","\\1",files$V9)==temp.multi.dif.melt$setfiles[i],"date.time"])
}

temp.multi.dif.melt$ZT<-as.numeric(temp.multi.dif.melt$date.time - strptime(first_day_ZT0,"%b/%d/%H:%M"), units = "hours",format="%b/%d/%H:%M") 
	return(temp.multi.dif.melt)
}

# function exp3set2 special
# functin for "multi" output (under construction)
format.luc2data.multi.exp3set2special<-function(plant.name,luc2.data,luc2.data.multi,files,first_day_ZT0) { 
	temp.multi<-luc2.data.multi[,c("Slice",grep("^IntDen",names(luc2.data.multi),value=TRUE))]
	# add file name
	set.files<-unique(gsub("([[:print:]]+)(exp[[:digit:]]+)","\\2",luc2.data$Label)) # keeping file order (by time)
	#print(set.files)
	# remove label file
	temp.multi$setfiles<-set.files 
	temp.multi<-temp.multi[grep("lab",temp.multi$setfiles,invert=TRUE),]
	write.csv(temp.multi, file="/Volumes/Data6/data_JM4/promYUC_luc2/exp3/temp.multi.set2.csv") # and then edit files
	temp.multi.mod<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/exp3/temp.multi.set2.2.csv")
	# calculate increase/decrease from begining
	# temp.multi.dif<-as.matrix(temp.multi[,-1])-as.vector(temp.multi[1,-1]) # does not work
	temp.multi.dif<-sweep(data.matrix(temp.multi.mod[,grep("IntDen",names(temp.multi.mod))]),2,data.matrix(temp.multi.mod[1,grep("IntDen",names(temp.multi.mod))]))	
	temp.multi.dif<-as.data.frame(temp.multi.dif)
	temp.multi.dif<-cbind(temp.multi.dif,temp.multi.mod$setfiles)
	names(temp.multi.dif)[grep("setfiles",names(temp.multi.dif))]<-"setfiles"
	print(temp.multi.dif)

	# add sample names
	names(temp.multi.dif)[grep("IntDen",names(temp.multi.dif))]<-plant.name
	temp.multi.dif
	print(temp.multi.dif) # Error: id variables not found in data: setfiles
	# melt
	temp.multi.dif.melt<-melt(temp.multi.dif,id.var="setfiles")
	temp.multi.dif.melt[is.na(temp.multi.dif.melt$value),"value"]<-0
	# treatment
	temp.multi.dif.melt$treatment<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\1",temp.multi.dif.melt$variable)
	
	# promoter
	temp.multi.dif.melt$promoter<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\5",temp.multi.dif.melt$variable)
	
	# line
	temp.multi.dif.melt$line<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\7",temp.multi.dif.melt$variable)	
	
# add time point info (under construction)
#print(files)
files$date.time<-strptime(paste(files$V6,files$V7,files$V8,sep="/"),"%b/%d/%H:%M")
#print(files)

temp.multi.dif.melt$date.time<-strptime("","%b/%d/%H:%M")
for(i in 1:dim(temp.multi.dif.melt)[1]) {
temp.multi.dif.melt[i,"date.time"]<-as.character(files[gsub("(exp[[:digit:]]+)(.tif)","\\1",files$V9)==temp.multi.dif.melt$setfiles[i],"date.time"])
}

temp.multi.dif.melt$ZT<-as.numeric(temp.multi.dif.melt$date.time - strptime(first_day_ZT0,"%b/%d/%H:%M"), units = "hours",format="%b/%d/%H:%M") 
	return(temp.multi.dif.melt)
}
###########################################################
################## for log2 fold change caluculation
######
format.luc2data.multi.log2<-function(plant.name,luc2.data,luc2.data.multi,files,first_day_ZT0) { 
	temp.multi<-luc2.data.multi[,c("Slice",grep("^IntDen",names(luc2.data.multi),value=TRUE))]
	# add file name
	set.files<-unique(gsub("([[:print:]]+)(exp[[:digit:]]+)","\\2",luc2.data$Label)) # keeping file order (by time)
	#print(set.files)
	# remove label file
	temp.multi$setfiles<-set.files
	temp.multi<-temp.multi[grep("lab",temp.multi$setfiles,invert=TRUE),]

	# calculate increase/decrease from begining
#	temp.multi.dif<-as.matrix(temp.multi[,-1])-as.vector(temp.multi[1,-1]) # does not work
#	temp.multi.dif<-sweep(data.matrix(temp.multi[,grep("IntDen",names(temp.multi))]),2,data.matrix(temp.multi[1,grep("IntDen",names(temp.multi))]))
#    log2 transform all data
	temp.multi.log2<-log2(temp.multi[,grep("IntDen",names(temp.multi))])
	temp.multi.dif<-sweep(data.matrix(temp.multi.log2),2,data.matrix(temp.multi.log2[1,])) # difference of log2 data = ratio

	temp.multi.dif<-as.data.frame(temp.multi.dif)
	temp.multi.dif<-cbind(temp.multi.dif,temp.multi$setfiles)
	names(temp.multi.dif)[grep("setfiles",names(temp.multi.dif))]<-"setfiles"
	#print(temp.multi.dif)

	# add sample names
	names(temp.multi.dif)[grep("IntDen",names(temp.multi.dif))]<-plant.name
	temp.multi.dif
	print(temp.multi.dif) # Error: id variables not found in data: setfiles
	# melt
	temp.multi.dif.melt<-melt(temp.multi.dif,id.var="setfiles")

	# treatment
	temp.multi.dif.melt$treatment<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\1",temp.multi.dif.melt$variable)
	
	# promoter
	temp.multi.dif.melt$promoter<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\5",temp.multi.dif.melt$variable)
	
	# line
	temp.multi.dif.melt$line<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\7",temp.multi.dif.melt$variable)	
	
# add time point info (under construction)
#print(files)
files$date.time<-strptime(paste(files$V6,files$V7,files$V8,sep="/"),"%b/%d/%H:%M")
#print(files)

temp.multi.dif.melt$date.time<-strptime("","%b/%d/%H:%M")
for(i in 1:dim(temp.multi.dif.melt)[1]) {
temp.multi.dif.melt[i,"date.time"]<-as.character(files[gsub("(exp[[:digit:]]+)(.tif)","\\1",files$V9)==temp.multi.dif.melt$setfiles[i],"date.time"])
}

temp.multi.dif.melt$ZT<-as.numeric(temp.multi.dif.melt$date.time - strptime(first_day_ZT0,"%b/%d/%H:%M"), units = "hours",format="%b/%d/%H:%M") 
	return(temp.multi.dif.melt)
}
# functin for "multi" output (under construction)
format.luc2data.multi.exp3set1special.log2<-function(plant.name,luc2.data,luc2.data.multi,files,first_day_ZT0) { 
	temp.multi<-luc2.data.multi[,c("Slice",grep("^IntDen",names(luc2.data.multi),value=TRUE))]
	# add file name
	set.files<-unique(gsub("([[:print:]]+)(exp[[:digit:]]+)","\\2",luc2.data$Label)) # keeping file order (by time)
	#print(set.files)
	# remove label file
	temp.multi$setfiles<-set.files 
	temp.multi<-temp.multi[grep("lab",temp.multi$setfiles,invert=TRUE),]
	write.csv(temp.multi, file="temp.multi.csv") # and then edit files
	#temp.multi.mod<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/exp3/temp.multi2.csv")
	temp.multi.mod<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp3/temp.multi2.csv")
	
	# calculate increase/decrease from begining
	# temp.multi.dif<-as.matrix(temp.multi[,-1])-as.vector(temp.multi[1,-1]) # does not work
	temp.multi.log2<-log2(temp.multi.mod[,grep("IntDen",names(temp.multi.mod))])
	temp.multi.dif<-sweep(data.matrix(temp.multi.log2),2,data.matrix(temp.multi.log2[1,]))
	temp.multi.dif<-as.data.frame(temp.multi.dif)
	temp.multi.dif<-cbind(temp.multi.dif,temp.multi.mod$setfiles)
	names(temp.multi.dif)[grep("setfiles",names(temp.multi.dif))]<-"setfiles"
	print(temp.multi.dif)

	# add sample names
	names(temp.multi.dif)[grep("IntDen",names(temp.multi.dif))]<-plant.name
	temp.multi.dif
	print(temp.multi.dif) # Error: id variables not found in data: setfiles
	# melt
	temp.multi.dif.melt<-melt(temp.multi.dif,id.var="setfiles")
	temp.multi.dif.melt[is.na(temp.multi.dif.melt$value),"value"]<-0
	# treatment
	temp.multi.dif.melt$treatment<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\1",temp.multi.dif.melt$variable)
	
	# promoter
	temp.multi.dif.melt$promoter<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\5",temp.multi.dif.melt$variable)
	
	# line
	temp.multi.dif.melt$line<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\7",temp.multi.dif.melt$variable)	
	
# add time point info (under construction)
#print(files)
files$date.time<-strptime(paste(files$V6,files$V7,files$V8,sep="/"),"%b/%d/%H:%M")
#print(files)

temp.multi.dif.melt$date.time<-strptime("","%b/%d/%H:%M")
for(i in 1:dim(temp.multi.dif.melt)[1]) {
temp.multi.dif.melt[i,"date.time"]<-as.character(files[gsub("(exp[[:digit:]]+)(.tif)","\\1",files$V9)==temp.multi.dif.melt$setfiles[i],"date.time"])
}

temp.multi.dif.melt$ZT<-as.numeric(temp.multi.dif.melt$date.time - strptime(first_day_ZT0,"%b/%d/%H:%M"), units = "hours",format="%b/%d/%H:%M") 
	return(temp.multi.dif.melt)
}

# function exp3set2 special
# functin for "multi" output (under construction)
format.luc2data.multi.exp3set2special.log2<-function(plant.name,luc2.data,luc2.data.multi,files,first_day_ZT0) { 
	temp.multi<-luc2.data.multi[,c("Slice",grep("^IntDen",names(luc2.data.multi),value=TRUE))]
	# add file name
	set.files<-unique(gsub("([[:print:]]+)(exp[[:digit:]]+)","\\2",luc2.data$Label)) # keeping file order (by time)
	#print(set.files)
	# remove label file
	temp.multi$setfiles<-set.files 
	temp.multi<-temp.multi[grep("lab",temp.multi$setfiles,invert=TRUE),]
	write.csv(temp.multi, file="/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp3/temp.multi.set2.csv") # and then edit files
	temp.multi.mod<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp3/temp.multi.set2.2.csv")
	# calculate increase/decrease from begining
	temp.multi.log2<-log2(temp.multi.mod[,grep("IntDen",names(temp.multi.mod))])
	temp.multi.dif<-sweep(data.matrix(temp.multi.log2),2,data.matrix(temp.multi.log2[1,]))

	temp.multi.dif<-as.data.frame(temp.multi.dif)
	temp.multi.dif<-cbind(temp.multi.dif,temp.multi.mod$setfiles)
	names(temp.multi.dif)[grep("setfiles",names(temp.multi.dif))]<-"setfiles"
	print(temp.multi.dif)

	# add sample names
	names(temp.multi.dif)[grep("IntDen",names(temp.multi.dif))]<-plant.name
	temp.multi.dif
	print(temp.multi.dif) # Error: id variables not found in data: setfiles
	# melt
	temp.multi.dif.melt<-melt(temp.multi.dif,id.var="setfiles")
	temp.multi.dif.melt[is.na(temp.multi.dif.melt$value),"value"]<-0
	# treatment
	temp.multi.dif.melt$treatment<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\1",temp.multi.dif.melt$variable)
	
	# promoter
	temp.multi.dif.melt$promoter<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\5",temp.multi.dif.melt$variable)
	
	# line
	temp.multi.dif.melt$line<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\7",temp.multi.dif.melt$variable)	
	
# add time point info (under construction)
#print(files)
files$date.time<-strptime(paste(files$V6,files$V7,files$V8,sep="/"),"%b/%d/%H:%M")
#print(files)

temp.multi.dif.melt$date.time<-strptime("","%b/%d/%H:%M")
for(i in 1:dim(temp.multi.dif.melt)[1]) {
temp.multi.dif.melt[i,"date.time"]<-as.character(files[gsub("(exp[[:digit:]]+)(.tif)","\\1",files$V9)==temp.multi.dif.melt$setfiles[i],"date.time"])
}

temp.multi.dif.melt$ZT<-as.numeric(temp.multi.dif.melt$date.time - strptime(first_day_ZT0,"%b/%d/%H:%M"), units = "hours",format="%b/%d/%H:%M") 
	return(temp.multi.dif.melt)
}

## background subtraction (071316)
# functin for "multi" output 
format.luc2data.multi.exp3set1special.log2.BS<-function(plant.name,luc2.data,luc2.data.multi,files,first_day_ZT0) { 
  temp.multi<-luc2.data.multi[,c("Slice",grep("^IntDen",names(luc2.data.multi),value=TRUE))]
  # add file name
  set.files<-unique(gsub("([[:print:]]+)(exp[[:digit:]]+)","\\2",luc2.data$Label)) # keeping file order (by time)
  #print(set.files)
  # remove label file
  temp.multi$setfiles<-set.files 
  temp.multi<-temp.multi[grep("lab",temp.multi$setfiles,invert=TRUE),]
  write.csv(temp.multi, file="temp.multi.csv") # and then edit files
  #temp.multi.mod<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/exp3/temp.multi2.csv")
  #temp.multi.mod<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/raw_data/exp3/temp.multi2.csv")
  
  # background subtraction (BS; 071216)
  temp.multi.mod<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp3/temp.multi3.csv")
  temp.multi.mod2<-temp.multi.mod[,grep("IntDen",names(temp.multi.mod))]-temp.multi.mod[,"background"]
  ## the end of background subtraction  
  # calculate increase/decrease from begining
  temp.multi.log2<-log2(temp.multi.mod2[,grep("IntDen",names(temp.multi.mod2))])
  temp.multi.dif<-sweep(data.matrix(temp.multi.log2),2,data.matrix(temp.multi.log2[1,]))
  temp.multi.dif<-as.data.frame(temp.multi.dif)
  temp.multi.dif<-cbind(temp.multi.dif,temp.multi.mod$setfiles)
  names(temp.multi.dif)[grep("setfiles",names(temp.multi.dif))]<-"setfiles"
  print(temp.multi.dif)
  
  # add sample names
  names(temp.multi.dif)[grep("IntDen",names(temp.multi.dif))]<-plant.name
  temp.multi.dif
  print(temp.multi.dif) # Error: id variables not found in data: setfiles
  # melt
  temp.multi.dif.melt<-melt(temp.multi.dif,id.var="setfiles")
  temp.multi.dif.melt[is.na(temp.multi.dif.melt$value),"value"]<-0
  # treatment
  temp.multi.dif.melt$treatment<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\1",temp.multi.dif.melt$variable)
  
  # promoter
  temp.multi.dif.melt$promoter<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\5",temp.multi.dif.melt$variable)  
  # line
  temp.multi.dif.melt$line<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\7",temp.multi.dif.melt$variable)	
  # add time point info (under construction)
  #print(files)
  files$date.time<-strptime(paste(files$V6,files$V7,files$V8,sep="/"),"%b/%d/%H:%M")
  #print(files)
  temp.multi.dif.melt$date.time<-strptime("","%b/%d/%H:%M")
  for(i in 1:dim(temp.multi.dif.melt)[1]) {
    temp.multi.dif.melt[i,"date.time"]<-as.character(files[gsub("(exp[[:digit:]]+)(.tif)","\\1",files$V9)==temp.multi.dif.melt$setfiles[i],"date.time"])
  }  
  temp.multi.dif.melt$ZT<-as.numeric(temp.multi.dif.melt$date.time - strptime(first_day_ZT0,"%b/%d/%H:%M"), units = "hours",format="%b/%d/%H:%M") 
  return(temp.multi.dif.melt)
}

# function exp3set2 special
# functin for "multi" output (under construction)
format.luc2data.multi.exp3set2special.log2.BS<-function(plant.name,luc2.data,luc2.data.multi,files,first_day_ZT0) { 
  temp.multi<-luc2.data.multi[,c("Slice",grep("^IntDen",names(luc2.data.multi),value=TRUE))]
  # add file name
  set.files<-unique(gsub("([[:print:]]+)(exp[[:digit:]]+)","\\2",luc2.data$Label)) # keeping file order (by time)
  #print(set.files)
  # remove label file
  temp.multi$setfiles<-set.files 
  temp.multi<-temp.multi[grep("lab",temp.multi$setfiles,invert=TRUE),]
  #write.csv(temp.multi, file="/Volumes/Data6/data_JM4/promYUC_luc2/exp3/temp.multi.set2.csv") # and then edit files
  temp.multi.mod<-read.csv("/Volumes/Data6/data_JM4/promYUC_luc2/promYUC_luc2_scripts_data/raw_data/exp3/temp.multi.set2.2.csv")
  # background subtraction (BS; 071516)
  temp.multi.mod2<-temp.multi.mod[,grep("IntDen",names(temp.multi.mod))]-temp.multi.mod[,"background"]
  
  # calculate increase/decrease from begining
  # temp.multi.dif<-as.matrix(temp.multi[,-1])-as.vector(temp.multi[1,-1]) # does not work
  temp.multi.log2<-log2(temp.multi.mod2[,grep("IntDen",names(temp.multi.mod2))])
  temp.multi.dif<-sweep(data.matrix(temp.multi.log2),2,data.matrix(temp.multi.log2[1,]))
  
  temp.multi.dif<-as.data.frame(temp.multi.dif)
  temp.multi.dif<-cbind(temp.multi.dif,temp.multi.mod$setfiles)
  names(temp.multi.dif)[grep("setfiles",names(temp.multi.dif))]<-"setfiles"
  print(temp.multi.dif)
  
  # add sample names
  names(temp.multi.dif)[grep("IntDen",names(temp.multi.dif))]<-plant.name
  temp.multi.dif
  print(temp.multi.dif) # Error: id variables not found in data: setfiles
  # melt
  temp.multi.dif.melt<-melt(temp.multi.dif,id.var="setfiles")
  temp.multi.dif.melt[is.na(temp.multi.dif.melt$value),"value"]<-0
  # treatment
  temp.multi.dif.melt$treatment<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\1",temp.multi.dif.melt$variable)
  
  # promoter
  temp.multi.dif.melt$promoter<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\5",temp.multi.dif.melt$variable)
  
  # line
  temp.multi.dif.melt$line<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\7",temp.multi.dif.melt$variable)	
  
  # add time point info (under construction)
  #print(files)
  files$date.time<-strptime(paste(files$V6,files$V7,files$V8,sep="/"),"%b/%d/%H:%M")
  #print(files)
  
  temp.multi.dif.melt$date.time<-strptime("","%b/%d/%H:%M")
  for(i in 1:dim(temp.multi.dif.melt)[1]) {
    temp.multi.dif.melt[i,"date.time"]<-as.character(files[gsub("(exp[[:digit:]]+)(.tif)","\\1",files$V9)==temp.multi.dif.melt$setfiles[i],"date.time"])
  }
  
  temp.multi.dif.melt$ZT<-as.numeric(temp.multi.dif.melt$date.time - strptime(first_day_ZT0,"%b/%d/%H:%M"), units = "hours",format="%b/%d/%H:%M") 
  return(temp.multi.dif.melt)
}
###########################################################
################## for log2 fold change caluculation (backgroudn subtraction version)
######
format.luc2data.multi.log2.BS<-function(plant.name,luc2.data,luc2.data.multi,files,first_day_ZT0) { 
  #temp.multi<-luc2.data.multi[,c("Slice",grep("^IntDen",names(luc2.data.multi),value=TRUE))]
  temp.multi<-luc2.data.multi[,c("Slice",grep("RawIntDen",names(luc2.data.multi),value=TRUE),"background")]
  
  # add file name
  set.files<-unique(gsub("([[:print:]]+)(exp[[:digit:]]+)","\\2",luc2.data$Label)) # keeping file order (by time)
  #print(set.files)
  # remove label file
  temp.multi$setfiles<-set.files
  temp.multi<-temp.multi[grep("lab",temp.multi$setfiles,invert=TRUE),]
  # background subtraction (BS; 071516)
  temp.multi.mod2<-temp.multi[,grep("RawIntDen",names(temp.multi))]-temp.multi[,"background"]
  
  
  # calculate increase/decrease from begining
  #    log2 transform all data
  temp.multi.log2<-log2(temp.multi[,grep("RawIntDen",names(temp.multi))])
  temp.multi.dif<-sweep(data.matrix(temp.multi.log2),2,data.matrix(temp.multi.log2[1,])) # difference of log2 data = ratio
  
  temp.multi.dif<-as.data.frame(temp.multi.dif)
  temp.multi.dif<-cbind(temp.multi.dif,temp.multi$setfiles)
  names(temp.multi.dif)[grep("setfiles",names(temp.multi.dif))]<-"setfiles"
  #print(temp.multi.dif)
  
  # add sample names
  names(temp.multi.dif)[grep("IntDen",names(temp.multi.dif))]<-plant.name
  temp.multi.dif
  print(temp.multi.dif) # Error: id variables not found in data: setfiles
  # melt
  temp.multi.dif.melt<-melt(temp.multi.dif,id.var="setfiles")
  
  # treatment
  temp.multi.dif.melt$treatment<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\1",temp.multi.dif.melt$variable)
  
  # promoter
  temp.multi.dif.melt$promoter<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\5",temp.multi.dif.melt$variable)
  
  # line
  temp.multi.dif.melt$line<-gsub("(sun|shade)(_)(plate[[:digit:]])(_)(YUC[[:digit:]]+)(_)(line[[:digit:]]+)(_)(rep[[:digit:]])","\\7",temp.multi.dif.melt$variable)	
  
  # add time point info (under construction)
  #print(files)
  files$date.time<-strptime(paste(files$V6,files$V7,files$V8,sep="/"),"%b/%d/%H:%M")
  #print(files)
  
  temp.multi.dif.melt$date.time<-strptime("","%b/%d/%H:%M")
  for(i in 1:dim(temp.multi.dif.melt)[1]) {
    temp.multi.dif.melt[i,"date.time"]<-as.character(files[gsub("(exp[[:digit:]]+)(.tif)","\\1",files$V9)==temp.multi.dif.melt$setfiles[i],"date.time"])
  }
  
  temp.multi.dif.melt$ZT<-as.numeric(temp.multi.dif.melt$date.time - strptime(first_day_ZT0,"%b/%d/%H:%M"), units = "hours",format="%b/%d/%H:%M") 
  return(temp.multi.dif.melt)
}



