# function.R
# 	clean up 022614 Kazu
library(lme4);library(lmerTest);library(ggplot2);library(reshape2)

calc.sem <- function(y) {
  y <- y[complete.cases(y)]
  mean.y <- mean(y)
  sem.y <- sd(y)/sqrt(length(y))
  data.frame(ymax=mean.y+sem.y,ymin=mean.y-sem.y)
}



model2c.lmer<-function(phenotype,data,refAccession,refTreatment,formula,saveplot=F) { # specify level before running this function becasue releveling data and use it in lme did not work although data was releveled (in model1 function).
	# model2("hypocotyl",hyp_SAS_data2.healthy,"Col",lme(hypocotyl~plant*treatment+exp,random=~treatment|rep,data=hyp_SAS_data2.healthy))
	# able to work with refTreatment == "sun" or "shade" (052213)
	# requires library(lme4);library(languageR);library(ggplot2)
	# str(data)
	print(paste("phenotype is ",phenotype,".",sep=""))
	print(paste("reference (accession) is ",refAccession,".",sep=""))
	print(paste("reference (treatment) is ",refTreatment,".",sep=""))

	temp.lmer<-formula
	# print(temp.lme)
	# table<-as.data.frame(summary(temp.lme)$tTable)
	table<-as.data.frame(summary(temp.lmer)@coefs)
    # calculate pvalue (require library(linguageR))
	pvals<-pvals.fnc(templ.lmer,addPlot=TRUE,ndigits=3,nsim=1000) # for testing model2c.lmer nsim=1000, for real calc. set nsim=10000)
	table$pvalue<-pvals # is this correct format?
   #### neeeds to work on "treatmentsun" for refAccession (052013). solved (052213)   
	if(refAccession=="Col") {
		rownames(table)<-sub("\\(Intercept\\)","plantCol",rownames(table)) 
		rownames(table)<-sub("^treatmentshade$","plantCol\\:treatmentshade",rownames(table))
		rownames(table)<-sub("^treatmentsun$","plantCol\\:treatmentsun",rownames(table))
		}
	else if(refAccession=="Ler") {
		rownames(table)<-sub("\\(Intercept\\)","plantLer",rownames(table)) 
		rownames(table)<-sub("^treatmentshade$","plantLer\\:treatmentshade",rownames(table))
		rownames(table)<-sub("^treatmentsun$","plantLer\\:treatmentsun",rownames(table))
		}
	else if(refAccession=="Ws") {
		rownames(table)<-sub("\\(Intercept\\)","plantWs",rownames(table)) 
		rownames(table)<-sub("^treatmentshade$","plantWs\\:treatmentshade",rownames(table))
		rownames(table)<-sub("^treatmentsun$","plantWs\\:treatmentsun",rownames(table))
		}
    print(table)
    # extract extra fixed effects (ohter than light, treatment)
    fixed_effects<-unique(sub("([[:print:]]+)\\:([[:print:]]+)","\\2",grep("\\:",rownames(table),value=T)))
    extra_fixed_effects<-grep(pattern="treatment",x=fixed_effects,invert=T,value=T)
    extra_fixed_effects<-c(extra_fixed_effects,unique(sub("([[:alpha:]]+)([[:digit:]]+)","\\1",rownames(table)[grep("plant",rownames(table),invert=T)[grep("plant",rownames(table),invert=T) %in% grep("treatment",rownames(table),invert=T)]])))

    print("fixed effects other than treatment and plant are:")
    print(extra_fixed_effects)
       # refTreatment is "sun" (needs to became more general expression to cope with various formula)
   if(refTreatment=="sun"){
	table.sun<-table[-grep(pattern="treatmentshade",x=rownames(table)),]
		if(length(extra_fixed_effects)==1) {
			table.sun<-table.sun[grep(pattern=extra_fixed_effects,x=rownames(table.sun),invert=T),]
			}else if(length(extra_fixed_effects)==2) {
				table.sun<-table.sun[grep(pattern=extra_fixed_effects[1],x=rownames(table.sun),invert=T),]
				table.sun<-table.sun[grep(pattern=extra_fixed_effects[2],x=rownames(table.sun),invert=T),]				
			}else if(length(extra_fixed_effects)==3) {
				table.sun<-table.sun[grep(pattern=extra_fixed_effects[1],x=rownames(table.sun),invert=T),]
				table.sun<-table.sun[grep(pattern=extra_fixed_effects[2],x=rownames(table.sun),invert=T),]				
				table.sun<-table.sun[grep(pattern=extra_fixed_effects[3],x=rownames(table.sun),invert=T),]				
			}else if(length(extra_fixed_effects)==4) {
				table.sun<-table.sun[grep(pattern=extra_fixed_effects[1],x=rownames(table.sun),invert=T),]
				table.sun<-table.sun[grep(pattern=extra_fixed_effects[2],x=rownames(table.sun),invert=T),]				
				table.sun<-table.sun[grep(pattern=extra_fixed_effects[3],x=rownames(table.sun),invert=T),]				
				table.sun<-table.sun[grep(pattern=extra_fixed_effects[4],x=rownames(table.sun),invert=T),]
			}
	# table.sun<-table.sun[grep(pattern="plant",x=rownames(table.sun)),]
	print("table.sun is")
	print(table.sun)
	table.sun$plant<-gsub("plant","",rownames(table.sun))
	table.shade<-table[grep(pattern="treatmentshade",x=rownames(table)),]
	table.sun$plant<-gsub("plant","",rownames(table.sun))
	#
	table.shade$plant<-gsub("plant","",rownames(table.shade))
	table.shade$plant[1]<-refAccession 
	table.shade$plant<-gsub(":treatmentshade","",table.shade$plant)
	print("table.shade is")
	print(table.shade)	
	}
	else if(refTreatment =="shade") {
		
	table.sun<-table[grep(pattern="treatmentsun",x=rownames(table)),]
		# print(gsub("^treatmentsun","plantCol",rownames(table.sun)))

	# table.sun<-table.sun[grep(pattern="plant",x=rownames(table.sun)),]
	# table.sun<-table.sun[grep(pattern=extra_fixed_effects,x=rownames(table.sun),invert=T),]

	table.sun$plant<-gsub("plant","",rownames(table.sun))	
	table.sun$plant<-gsub(":treatmentsun","",table.sun$plant)
	print("table.sun is")
	print(table.sun)
	table.shade<-table[grep(pattern="treatmentsun",x=rownames(table),invert=T),]
		# print("table.shade is (before grep(pattern=extra_fixed_effects,x=rownames(table.shade),invert=T))")
			# print(table.shade)
		if(length(extra_fixed_effects)==1) {
			table.shade<-table.shade[grep(pattern=extra_fixed_effects,x=rownames(table.shade),invert=T),]
			
		# print("table.shade is (after grep(pattern=extra_fixed_effects,x=rownames(table.shade),invert=T))"	
		# print(table.shade)
        # print("grep result is")
        # print(grep(pattern=extra_fixed_effects,x=rownames(table.shade),invert=T))
			}else if(length(extra_fixed_effects)==2) {
				table.shade<-table.shade[grep(pattern=extra_fixed_effects[1],x=rownames(table.shade),invert=T),]
				table.shade<-table.shade[grep(pattern=extra_fixed_effects[2],x=rownames(table.shade),invert=T),]
			}else if(length(extra_fixed_effects)==3) {
				table.shade<-table.shade[grep(pattern=extra_fixed_effects[1],x=rownames(table.shade),invert=T),]
				table.shade<-table.shade[grep(pattern=extra_fixed_effects[2],x=rownames(table.shade),invert=T),]
				table.shade<-table.shade[grep(pattern=extra_fixed_effects[3],x=rownames(table.shade),invert=T),]
			}else if(length(extra_fixed_effects)==4) {
				table.shade<-table.shade[grep(pattern=extra_fixed_effects[1],x=rownames(table.shade),invert=T),]
				table.shade<-table.shade[grep(pattern=extra_fixed_effects[2],x=rownames(table.shade),invert=T),]
				table.shade<-table.shade[grep(pattern=extra_fixed_effects[3],x=rownames(table.shade),invert=T),]
				table.shade<-table.shade[grep(pattern=extra_fixed_effects[4],x=rownames(table.shade),invert=T),]
			}

	table.shade$plant<-gsub("plant","",rownames(table.shade))
	table.shade$plant[1]<-refAccession
	print(table.shade)
	}
	# calculating absolute value
	if(refTreatment=="sun"){
		table.sun$Estimate[-1]<-table.sun$Estimate[-1] + table.sun$Estimate[1]	
		table.shade$Estimate[-1]<-table.shade$Estimate[-1] + table.shade$Estimate[1]
		table.shade$Estimate<-table.shade$Estimate + table.sun$Estimate

	} else if (refTreatment=="shade"){
		table.shade$Estimate[-1]<-table.shade$Estimate[-1] + table.shade$Estimate[1]	
		table.sun$Estimate[-1]<-table.sun$Estimate[-1] + table.sun$Estimate[1]
		table.sun$Estimate<-table.sun$Estimate + table.shade$Estimate
	}
	print("table.shade with absolute value is")
	print(table.shade)
	print("table.sun with absolute value is")
	print(table.sun)
	print("length(table.shade$Estimate) is")
	print(length(table.shade$Estimate))
	print("length(table.sun$Estimate) is")
	print(length(table.sun$Estimate))
	
	
	# select only Value, Std. Error, p-value columns (sun)
	table.sun<-table.sun[,c("plant","Estimate","Std. Error","p-value")]
	names(table.sun)<-c("plant",paste(phenotype,"_sun_mean","_r",refAccession,"_r",refTreatment,sep=""),paste(phenotype,"_sun_SE","_r",refAccession,"_r",refTreatment, sep=""),paste(phenotype,"_sun_pvalue","_r",refAccession,"_r",refTreatment,sep=""))
	# select only Value, Std. Error, p-value columns (shade)
	table.shade<-table.shade[,c("plant","Estimate","Std. Error","p-value")]
	names(table.shade)<-c("plant",paste(phenotype,"_shade_mean","_r",refAccession,"_r",refTreatment,sep=""),paste(phenotype,"_shade_SE","_r",refAccession,"_r",refTreatment, sep=""),paste(phenotype,"_shade_pvalue","_r",refAccession,"_r",refTreatment,sep=""))
	# merge
	table.summary<-merge(table.sun,table.shade,by="plant")
	## regression diagnosis # function of "modelcheck" is derived from Danielle's script (big_hyps_phenotype_analysis_simplified.R)
	## ggplot2
		rs <- residuals(temp.lme)
		fv <- fitted(temp.lme)
		a<-qplot(table.summary[,paste(phenotype,"_sun_mean","_r",refAccession,"_r",refTreatment,sep="")],table.summary[,paste(phenotype,"_shade_mean","_r",refAccession,"_r",refTreatment,sep="")] - table.summary[,paste(phenotype,"_sun_mean","_r",refAccession,"_r",refTreatment,sep="")],main="sun vs response") + labs(x="sun",y="response")
		b<-qplot(fv,rs,main="homogeneity")
		c<- qplot(sample=rs, stat="qq",main="normality (Q-Q plot)")
		d<-qplot(fv,sqrt(abs(rs)),main="?")		
		e<-qplot(plant,rs,data=data,geom="boxplot",main="independence") + facet_grid(treatment~.) # check independence
		library(grid)
		pdf(file=paste(phenotype,"modelcheckoutput.pdf",sep="_"),h=8,w=10.5)
			grid.newpage()
			grid.text(phenotype,vp=viewport(0.1,0.1,x=0.5,y=0.975))
			print(a, vp=viewport(0.3, 0.4, x=0.2, y=0.75))
			print(b, vp=viewport(0.3, 0.4, x=0.5, y=0.75))
			print(c, vp=viewport(0.3, 0.4, x=0.8, y=0.75))
			print(d, vp=viewport(0.3, 0.4, x=0.2, y=0.25))
		    print(e, vp=viewport(0.5, 0.4, x=0.65, y=0.25))
		dev.off()
	return(table.summary)
} # close model2c.lmer
	
object.name<-function(data) {
	# temp<-names(data)
	temp<-attributes(data)
	print(temp)
}
# update ggplot2 code (060816)
model2e.lmer<-function(phenotype,data,refAccession,refTreatment,formula,saveplot=F,unit="mm") { # specify level before running this function becasue releveling data and use it in lme did not work although data was releveled (in model1 function).
  # model2("hypocotyl",hyp_SAS_data2.healthy,"Col",lme(hypocotyl~plant*treatment+exp,random=~treatment|rep,data=hyp_SAS_data2.healthy))
  # able to work with refTreatment == "sun" or "shade" (052213)
  # requires library(lme4);library(lmerTest);library(ggplot2)
  
  # str(data)
  print(paste("phenotype is ",phenotype,".",sep=""))
  print(paste("reference (accession) is ",refAccession,".",sep=""))
  print(paste("reference (treatment) is ",refTreatment,".",sep=""))
  
  temp.lmer<-formula
  # print(temp.lme)
  # table<-as.data.frame(summary(temp.lme)$tTable)
  # table<-as.data.frame(summary(temp.lmer)@coefs)
  table<-as.data.frame(summary(temp.lmer)$coefficients)
  # calculate pvalue (require library(linguageR))
  #pvals<-pvals.fnc(templ.lmer,addPlot=TRUE,ndigits=3,nsim=1000) # for testing model2c.lmer nsim=1000, for real calc. set nsim=10000)
  #table$pvalue<-pvals # is this correct format?
  #### neeeds to work on "treatmentsun" for refAccession (052013). solved (052213)   
  if(refAccession=="Col") {
    rownames(table)<-sub("\\(Intercept\\)","plantCol",rownames(table)) 
    rownames(table)<-sub("^treatmentshade$","plantCol\\:treatmentshade",rownames(table))
    rownames(table)<-sub("^treatmentsun$","plantCol\\:treatmentsun",rownames(table))
  }
  else if(refAccession=="Ler") {
    rownames(table)<-sub("\\(Intercept\\)","plantLer",rownames(table)) 
    rownames(table)<-sub("^treatmentshade$","plantLer\\:treatmentshade",rownames(table))
    rownames(table)<-sub("^treatmentsun$","plantLer\\:treatmentsun",rownames(table))
  }
  else if(refAccession=="Ws") {
    rownames(table)<-sub("\\(Intercept\\)","plantWs",rownames(table)) 
    rownames(table)<-sub("^treatmentshade$","plantWs\\:treatmentshade",rownames(table))
    rownames(table)<-sub("^treatmentsun$","plantWs\\:treatmentsun",rownames(table))
  }
  print(table)
  # extract extra fixed effects (ohter than light, treatment)
  fixed_effects<-unique(sub("([[:print:]]+)\\:([[:print:]]+)","\\2",grep("\\:",rownames(table),value=T)))
  extra_fixed_effects<-grep(pattern="treatment",x=fixed_effects,invert=T,value=T)
  extra_fixed_effects<-c(extra_fixed_effects,unique(sub("([[:alpha:]]+)([[:digit:]]+)","\\1",rownames(table)[grep("plant",rownames(table),invert=T)[grep("plant",rownames(table),invert=T) %in% grep("treatment",rownames(table),invert=T)]])))
  
  print("fixed effects other than treatment and plant are:")
  print(extra_fixed_effects)
  # refTreatment is "sun" (needs to became more general expression to cope with various formula)
  if(refTreatment=="sun"){
    table.sun<-table[-grep(pattern="treatmentshade",x=rownames(table)),]
    if(length(extra_fixed_effects)==1) {
      table.sun<-table.sun[grep(pattern=extra_fixed_effects,x=rownames(table.sun),invert=T),]
    }else if(length(extra_fixed_effects)==2) {
      table.sun<-table.sun[grep(pattern=extra_fixed_effects[1],x=rownames(table.sun),invert=T),]
      table.sun<-table.sun[grep(pattern=extra_fixed_effects[2],x=rownames(table.sun),invert=T),]  			
    }else if(length(extra_fixed_effects)==3) {
      table.sun<-table.sun[grep(pattern=extra_fixed_effects[1],x=rownames(table.sun),invert=T),]
      table.sun<-table.sun[grep(pattern=extra_fixed_effects[2],x=rownames(table.sun),invert=T),]				
      table.sun<-table.sun[grep(pattern=extra_fixed_effects[3],x=rownames(table.sun),invert=T),]				
    }else if(length(extra_fixed_effects)==4) {
      table.sun<-table.sun[grep(pattern=extra_fixed_effects[1],x=rownames(table.sun),invert=T),]
      table.sun<-table.sun[grep(pattern=extra_fixed_effects[2],x=rownames(table.sun),invert=T),]				
      table.sun<-table.sun[grep(pattern=extra_fixed_effects[3],x=rownames(table.sun),invert=T),]				
      table.sun<-table.sun[grep(pattern=extra_fixed_effects[4],x=rownames(table.sun),invert=T),]
    }
    # table.sun<-table.sun[grep(pattern="plant",x=rownames(table.sun)),]
    print("table.sun is")
    print(table.sun)
    table.sun$plant<-gsub("plant","",rownames(table.sun))
    table.shade<-table[grep(pattern="treatmentshade",x=rownames(table)),]
    table.sun$plant<-gsub("plant","",rownames(table.sun))
    #
    table.shade$plant<-gsub("plant","",rownames(table.shade))
    table.shade$plant[1]<-refAccession 
    table.shade$plant<-gsub(":treatmentshade","",table.shade$plant)
    print("table.shade is")
    print(table.shade)	
  }
  else if(refTreatment =="shade") {
    
    table.sun<-table[grep(pattern="treatmentsun",x=rownames(table)),]
    # print(gsub("^treatmentsun","plantCol",rownames(table.sun)))
    
    # table.sun<-table.sun[grep(pattern="plant",x=rownames(table.sun)),]
    # table.sun<-table.sun[grep(pattern=extra_fixed_effects,x=rownames(table.sun),invert=T),]
    
    table.sun$plant<-gsub("plant","",rownames(table.sun))	
    table.sun$plant<-gsub(":treatmentsun","",table.sun$plant)
    print("table.sun is")
    print(table.sun)
    table.shade<-table[grep(pattern="treatmentsun",x=rownames(table),invert=T),]
    # print("table.shade is (before grep(pattern=extra_fixed_effects,x=rownames(table.shade),invert=T))")
    # print(table.shade)
    if(length(extra_fixed_effects)==1) {
      table.shade<-table.shade[grep(pattern=extra_fixed_effects,x=rownames(table.shade),invert=T),]
      
      # print("table.shade is (after grep(pattern=extra_fixed_effects,x=rownames(table.shade),invert=T))"	
      # print(table.shade)
      # print("grep result is")
      # print(grep(pattern=extra_fixed_effects,x=rownames(table.shade),invert=T))
    }else if(length(extra_fixed_effects)==2) {
      table.shade<-table.shade[grep(pattern=extra_fixed_effects[1],x=rownames(table.shade),invert=T),]
      table.shade<-table.shade[grep(pattern=extra_fixed_effects[2],x=rownames(table.shade),invert=T),]
    }else if(length(extra_fixed_effects)==3) {
      table.shade<-table.shade[grep(pattern=extra_fixed_effects[1],x=rownames(table.shade),invert=T),]
      table.shade<-table.shade[grep(pattern=extra_fixed_effects[2],x=rownames(table.shade),invert=T),]
      table.shade<-table.shade[grep(pattern=extra_fixed_effects[3],x=rownames(table.shade),invert=T),]
    }else if(length(extra_fixed_effects)==4) {
      table.shade<-table.shade[grep(pattern=extra_fixed_effects[1],x=rownames(table.shade),invert=T),]
      table.shade<-table.shade[grep(pattern=extra_fixed_effects[2],x=rownames(table.shade),invert=T),]
      table.shade<-table.shade[grep(pattern=extra_fixed_effects[3],x=rownames(table.shade),invert=T),]
      table.shade<-table.shade[grep(pattern=extra_fixed_effects[4],x=rownames(table.shade),invert=T),]
    }
    
    table.shade$plant<-gsub("plant","",rownames(table.shade))
    table.shade$plant[1]<-refAccession
    print(table.shade)
  }
  # calculating absolute value
  if(refTreatment=="sun"){
    table.sun$Estimate[-1]<-table.sun$Estimate[-1] + table.sun$Estimate[1]	
    table.shade$Estimate[-1]<-table.shade$Estimate[-1] + table.shade$Estimate[1]
    table.shade$Estimate<-table.shade$Estimate + table.sun$Estimate
    
  } else if (refTreatment=="shade"){
    table.shade$Estimate[-1]<-table.shade$Estimate[-1] + table.shade$Estimate[1]	
    table.sun$Estimate[-1]<-table.sun$Estimate[-1] + table.sun$Estimate[1]
    table.sun$Estimate<-table.sun$Estimate + table.shade$Estimate
  }
  print("table.shade with absolute value is")
  print(table.shade)
  print("table.sun with absolute value is")
  print(table.sun)
  print("length(table.shade$Estimate) is")
  print(length(table.shade$Estimate))
  print("length(table.sun$Estimate) is")
  print(length(table.sun$Estimate))
  
  
  # select only Value, Std. Error, p-value columns (sun)
  table.sun<-table.sun[,c("plant","Estimate","Std. Error","Pr(>|t|)")]
  names(table.sun)<-c("plant",paste(phenotype,"_sun_mean","_r",refAccession,"_r",refTreatment,sep=""),paste(phenotype,"_sun_SE","_r",refAccession,"_r",refTreatment, sep=""),paste(phenotype,"_sun_pvalue","_r",refAccession,"_r",refTreatment,sep=""))
  # select only Value, Std. Error, p-value columns (shade)
  table.shade<-table.shade[,c("plant","Estimate","Std. Error","Pr(>|t|)")]
  names(table.shade)<-c("plant",paste(phenotype,"_shade_mean","_r",refAccession,"_r",refTreatment,sep=""),paste(phenotype,"_shade_SE","_r",refAccession,"_r",refTreatment, sep=""),paste(phenotype,"_shade_pvalue","_r",refAccession,"_r",refTreatment,sep=""))
  # merge
  table.summary<-merge(table.sun,table.shade,by="plant")
  ## regression diagnosis # function of "modelcheck" is derived from Danielle's script (big_hyps_phenotype_analysis_simplified.R)
  ## ggplot2
  ####################
  modelcheck<-function(modelname,h,w, unit)	{
    ## ggplot2
    rs <- residuals(temp.lmer)
    fv <- fitted(temp.lmer)
    a<-qplot(table.summary[,paste(phenotype,"_sun_mean","_r",refAccession,"_r",refTreatment,sep="")],table.summary[,paste(phenotype,"_shade_mean","_r",refAccession,"_r",refTreatment,sep="")] - table.summary[,paste(phenotype,"_sun_mean","_r",refAccession,"_r",refTreatment,sep="")],main="sun vs response") + labs(x=paste("sun (",unit,")",sep=""),y=paste("response (",unit,")",sep=""))
        
    b<-qplot(fv,rs,main="homogeneity")
    ##c<- qplot(sample=rs, stat="qq",main="normality (Q-Q plot)")
    c <- ggplot(data.frame(x=rs,y=fv), aes(sample = y)) + stat_qq() + labs(title="normality (Q-Q plot)") # new version
    ##		d<-qplot(fv,sqrt(abs(rs)),main="?")		
    ##		e<-qplot(plant,rs,data=data,geom="boxplot",main="independence") + facet_grid(treatment~.) # check independence  
    library(grid)
    # drawing graph in pdf file			
    pdf(file=modelname,h=h,w=w)
    grid.newpage()
    grid.text(phenotype,vp=viewport(0.1,0.1,x=0.5,y=0.975))
    # ver2
    print(a, vp=viewport(0.45, 0.5, x=0.25, y=0.5))
#    print(b, vp=viewport(0.3, 0.4, x=0.5, y=0.5))
    print(c, vp=viewport(0.45, 0.5, x=0.7, y=0.5))
    dev.off()
  }
  #if(saveplot) 
    modelcheck(paste(phenotype,refAccession,"modelcheckoutput.pdf",sep="_"),h=8,w=10.5,unit=unit)
  #######
  return(table.summary)
} # close model2d.lmer

object.name<-function(data) {
  # temp<-names(data)
  temp<-attributes(data)
  print(temp)
}

### the end of model2d.lmer


## bar plot function (originally made by Julin)
###### plot.bar for ggplot2 ver 0.9.3
# for mixed effect model summary graph
plot.bar6 <- function(trait,data.input=data.plot,title="",ylabel,rownum=3,save.plot=T,show.plot=F,width=8,height=6,RFRratio=F) { 
  #make bar chart for a trait in Kazu's data
  if (title=="") title <- trait
  data.plot<-data.input[grep(paste("(^",trait,")\\w+(rCol)",sep=""),names(data.input),value=T)]
  print(data.plot)
  # names(data.plot)<-"diff" # Kazu added this line    
  data.plot$gene <- data.input$plant  
  # #Calculate difference between sun and shade to
  #sort plot
  #data.plot$dif <- (get(paste(trait,"shade",sep="_"),data.plot) -
  # # get(paste(trait,"sun",sep="_"),data.plot))
  data.plot$dif<-data.plot[,grep("sun_mean",names(data.plot))[grep("sun_mean",names(data.plot)) %in% grep("_rCol_rsun",names(data.plot))]]-data.plot[,grep("shade_mean",names(data.plot))[grep("shade_mean",names(data.plot)) %in% grep("_rCol_rsun",names(data.plot))]]    
  print("data.plot$dif is")
  print(data.plot$dif)     
  print("names(data.plot) are")
  print(names(data.plot))  
  # order according to response (data.plot$dif) # this is off in SA plot
  #data.plot <- data.plot[order(data.plot$dif),] 
  #print("ordered data is")
  print(data.plot)
  # #will need to preserve this ordering
  data.plot$rank <- 1:dim(data.plot)[1]  
  #subset and rearrange data
  data.plot.sun <- data.plot[c(grep("(_sun)",names(data.plot),value=T)[grep("(_sun)",names(data.plot)) %in% grep("(_rsun)",names(data.plot))],"gene","rank")]
  print("names(data.plot.sun) are")
  
  print(names(data.plot.sun))
  data.plot.shade <- data.plot[c(grep("(_shade)",names(data.plot),value=T)[grep("(_shade)",names(data.plot)) %in% grep("(_rshade)",names(data.plot))],"gene","rank")]
  print("names(data.plot.shade) are");print(names(data.plot.shade))
  data.plot.sun$trt <- "sun"
  data.plot.shade$trt <- "shade"
  names(data.plot.sun)[1:3] <- c("mean","SE","pvalue")
  names(data.plot.shade)[1:3] <- c("mean","SE","pvalue") # this pvalue is not "response" 
  # input response pvalue
  data.plot.shade$"pvalue"<-data.plot[,grep("(_shade_pvalue_rCol_rsun)",names(data.plot))]
  
  data.plot <- rbind(data.plot.sun,data.plot.shade)
  # if (length(data.plot)!=7) {
  # warning(paste("incorrect number of data columns for trait",trait,"not plotting"))
  # return()
  # }
  print("rearragend data.plot is");print(data.plot)
  data.plot$trt <- factor(data.plot$trt,levels=c("sun","shade"))
  data.plot$gene <- factor(data.plot$gene,levels=data.plot.sun$gene)
  data.plot$ymin=data.plot$mean - data.plot$SE
  data.plot$ymax=data.plot$mean + data.plot$SE
  
  # add significance label (p<0.05)
  data.plot$significance<-" "
  data.plot$significance[data.plot$pvalue<0.05]<-"*"
  # remove significance in Col in sun
  data.plot[data.plot$gene=="Col"&data.plot$trt=="sun","significance"]<-" "
  
  if(RFRratio) {
    data.plot$trt<-gsub("sun","high \nR/FR",as.character(data.plot$trt))
    data.plot$trt<-gsub("shade","low \nR/FR",as.character(data.plot$trt))
    data.plot$trt<-factor(data.plot$trt,levels=c("high \nR/FR","low \nR/FR"))
  } 
  
  print(data.plot)
  pl <- ggplot(data=data.plot)
  pl <- pl + geom_bar(mapping=aes(fill=trt,x=trt,y=mean),stat="identity")
  pl <- pl + facet_wrap(facets=~gene,nrow=rownum)
  pl <- pl + geom_errorbar(mapping=aes(x=trt,ymin=ymin,ymax=ymax))
  #  pl <- pl + opts(strip.text.x = theme_text(angle=90))
  pl <- pl + theme(strip.text.x = element_text(angle=90))
  #pl <- pl + theme(strip.text.x = element_text(angle=90,colour=data.plot$gene))
  pl <- pl + theme(strip.text.x = element_text(colour=data.plot$gene))
  #pl <- pl + opts(axis.title.x = theme_blank(), axis.title.y = theme_blank(),title=title)
  pl <-pl + theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90))
  pl <-pl + labs(title=title,y=ylabel)
  pl <- pl + geom_text(data=data.plot,aes(x=trt,y=ymax*1.05,label=factor(significance),color=trt))     # for significance
  pl
  if(save.plot) ggsave(filename=paste(title,"pdf",sep="."),
                       plot=pl,width=width,height=height)
  if(show.plot) pl
  return(pl) # added on 071715
}
### the end   



