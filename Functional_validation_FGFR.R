# Load packages
library(readxl)
library(reshape2)
library(ggplot2)
library(car)
library(dunn.test)
library(dplyr)

# Set directory
setwd('C:/Users/theod/Downloads')

## qPCR
# Load data
data<-read_xlsx('qPCR.xlsx')
data$logFC<-log2(data$FC/data$Control)
data$logFCF<-log2(data$`FC+F`/data$Control)
data<-data[,c(1,2,6,7)]

# Process data
data<-melt(data, id.vars=c('Sample', 'Analyte'))
data$variable<-factor(data$variable, levels=c("logFC", "logFCF"))

# Perform statistical analysis
results<-data.frame()

for(i in unique(data$Analyte)){
  data2<-data[data$Analyte==i,]
  
  groups<-unique(data2$variable)
  k<-length(groups)
  paired<-TRUE
  pair_id<-"Sample" 
  
  if(k < 2){
    print("Stop. Column of interest has fewer than 2 levels")
  }
  
  norm_p<-rep(NA, k)
  names(norm_p)<-groups
  
  if (k == 2) {
    
    if (paired) {
      wide<-reshape(data2,
                    idvar=pair_id,
                    timevar="variable",
                    direction="wide")
      
      if (nrow(wide) == 0) {
        stop("Stop. No complete pairs found")
      }
      
      differences<-wide[[3]] - wide[[5]]
      
      if(length(differences) >= 3) {
        norm_p<-shapiro.test(differences)$p.value
        normal<-norm_p > 0.05
      } else {
        normal <- FALSE
      }
      
      if (normal) {
        stat_p<-t.test(wide[[3]], wide[[5]], paired=TRUE, alternative='two.sided')
        tmp<-data.frame(parameter=i, test=stat_p$method, comparison=paste(colnames(wide)[2], 'vs', colnames(wide)[3], sep=' '), p.value=stat_p$p.value)
        results<-rbind(tmp, results)
        
      } else {
        stat_p<-wilcox.test(wide[[3]], wide[[5]], paired=TRUE, alternative='two.sided')
        tmp<-data.frame(parameter=i, test=stat_p$method, comparison=paste(colnames(wide)[2], 'vs', colnames(wide)[3], sep=' '), p.value=stat_p$p.value)
        results<-rbind(tmp, results)
      }
    }
  } 
}

results<-results[results$p.value <= 0.05,]

## LEGENDplex
# Load data
data<-read_xlsx('LEGENDplex.xlsx')

# Process data
data$logFC<-log2(data$FC/data$Control)
data$logFCF<-log2(data$`FC+F`/data$Control)
data<-data[,c(1,2,7,8)]

data<-melt(data, id.vars=c('Sample', 'Analyte'))
data$variable<-factor(data$variable, levels=c("logFC", "logFCF"))

# Perform statistical analysis
results<-data.frame()

for(i in unique(data$Analyte)){
  data3<-data[data$Analyte==i,]
  
  groups<-unique(data3$variable)
  k<-length(groups)
  paired<-TRUE
  pair_id<-"Sample" 
  
  if(k < 2){
    print("Stop. Column of interest has fewer than 2 levels")
  }
  
  norm_p<-rep(NA, k)
  names(norm_p)<-groups
  
  if (k == 2) {
    
    if (paired) {
      wide<-reshape(data3,
                    idvar=pair_id,
                    timevar="variable",
                    direction="wide")
      
      if (nrow(wide) == 0) {
        stop("Stop. No complete pairs found")
      }
      
      differences<-wide[[3]] - wide[[5]]
      
      if(length(differences) >= 3) {
        norm_p<-shapiro.test(differences)$p.value
        normal<-norm_p > 0.05
      } else {
        normal <- FALSE
      }
      
      if (normal) {
        stat_p<-t.test(wide[[3]], wide[[5]], paired=TRUE, alternative='two.sided')
        tmp<-data.frame(parameter=i, test=stat_p$method, comparison=paste(colnames(wide)[2], 'vs', colnames(wide)[4], sep=' '), p.value=stat_p$p.value)
        results<-rbind(tmp, results)
        
      } else {
        stat_p<-wilcox.test(wide[[3]], wide[[5]], paired=TRUE, alternative='two.sided')
        tmp<-data.frame(parameter=i, test=stat_p$method, comparison=paste(colnames(wide)[2], 'vs', colnames(wide)[4], sep=' '), p.value=stat_p$p.value)
        results<-rbind(tmp, results)
      }
    }
  } 
}

results<-results[results$p.value <= 0.05,]
