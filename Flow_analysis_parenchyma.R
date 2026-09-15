# Load packages
library(readxl)
library(ggplot2)
library(reshape2)
library(car)
library(dunn.test)

## SCS1

# Load data
data<-read_xlsx('6-1-21_Parenchyma statistics.xlsx')
colnames(data)<-c('LuT','Disease','Lin_percent','Granulocytes_percent','HLA.DR_myeloid cells_percent','Neutrophils_percent','HLA.DR_neutrophils_percent',
                  'Eosinophils_percent','CD14.monocytes_percent','Autofluorescence_percent','CD14.CD16_percent','AMs_percent','DCs_percent','Mast_percent','HLA.DR_mast_percent',
                  'Neutrophils','Eosinophils','CD14.monocytes','AMs','DCs','Mast')

data<-data[,c(1,2,16:21)]

#data$`HLA.DR_neutrophils_percent`<-as.numeric(data$`HLA.DR_neutrophils_percent`)
data[data$Disease=='Fibrosis EAA',]$Disease<-'Fibrosis'
data<-data[data$Disease %in% c("Tumor free",  "Emphysema", "Fibrosis"),]
data$Disease<-factor(data$Disease, levels=c("Tumor free",  "Emphysema", "Fibrosis"))

data<-melt(data, id.vars=c('Disease','LuT'))

# Perform statistical analysis
for (i in levels(data$variable)){
  df<-data[data$variable==i,]
  print(paste0('Test result for ',i))
  
  res.aov<-aov(value~Disease, df)
  var<-leveneTest(value~Disease, df)
  var<-var$`Pr(>F)`[1]
  resi<-shapiro.test(residuals(object=res.aov))
  resi<-resi$p.value
  
  if(var | resi <=0.05){
    dunn.test(x=df$value,g=df$Disease)
  }else{
    print(TukeyHSD(res.aov))
  }
}


## SCS2

# Load data
data<-read_xlsx('210115_flowDFG.xlsx', sheet=2)
colnames(data)<-c('LuT','Disease','Lin_percent','B.cells_percent','T.cells_percent','CD3.CD19_percent','CD4_percent','CD8_percent','NK_percent','ILC_percent',
                  'ILC2_percent','ILC3_percent','ILC1_percent','B.cells','T.cells','CD4','CD8','NK.cells','ILC')

data<-data[,c(1,2,14:18)]

data[data$Disease=='Fibrosis EAA',]$Disease<-'Fibrosis'
data<-data[data$Disease %in% c("Tumor free",  "Emphysema", "Fibrosis"),]
data$Disease<-factor(data$Disease, levels=c("Tumor free",  "Emphysema", "Fibrosis"))

data<-melt(data, id.vars=c('Disease','LuT'))

# Perform statistical analysis
for (i in levels(data$variable)){
  df<-data[data$variable==i,]
  print(paste0('Test result for ',i))
  
  res.aov<-aov(value~Disease, df)
  var<-leveneTest(value~Disease, df)
  var<-var$`Pr(>F)`[1]
  resi<-shapiro.test(residuals(object=res.aov))
  resi<-resi$p.value
  
  if(var | resi <=0.05){
    dunn.test(x=df$value,g=df$Disease)
  }else{
    print(TukeyHSD(res.aov))
  }
}
