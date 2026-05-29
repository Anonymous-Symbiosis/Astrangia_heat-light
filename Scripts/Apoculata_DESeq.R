#Astrangia poculata heat & light: reads to DESeq2
#Light environment shapes divergent coral thermal responses across symbiotic states
#This script was used for DESeq2 processing
#Note: In this script: apo = white corals, sym = brown corals



library(DESeq2)
library(tidyverse)


#----------Prepare countData----------
#import count data
countData <- read.table("../Data/holobiont_counts.txt")
#create column names
newColNames = c("ACC_10", "ACH_10", "ACC_11", "ACH_11", "ACC_12",
                "ACH_12", "ACC_13","ACH_13","ACC_14", "ACH_14",
                "ADC_15", "ADH_15", "ADC_16", "ADH_16", "ADC_17",
                "ADH_17", "ADC_18", "ADC_19", "ADH_19", "AHC_01",
                "AHH_01", "ADC_20", "ADH_20", "SHC_21", "SHH_21",
                "SHC_22", "SHH_22", "SHC_23", "SHH_23", "SHC_24",
                "SHH_24", "SHC_25", "SHH_25", "SHC_26", "SHH_26",
                "SCC_27", "SCH_27", "SCC_28", "SCH_28", "SCC_29",
                "SCH_29", "AHC_02", "AHH_02", "SCC_30", "SCC_31",
                "SCH_31", "SCC_32", "SCH_32", "SCC_33", "SCH_33",
                "SDC_34", "SDH_34", "SDC_35", "SDH_35", "SDC_36",
                "SDH_36", "SDC_37", "SDH_37", "SDC_38", "SDH_38",
                "SDC_39", "SDH_39", "AHC_03", "AHH_03", "ACC_40",
                "SCH_40", "AHC_04", "AHH_04", "AHH_05", "AHC_06",
                "AHH_06", "AHC_07", "AHH_07", "ACC_08", "ACH_08",
                "ACC_09", "ACH_09")
colnames(countData)=paste(newColNames)

new_order = sort(colnames(countData)) #reorder columns based on column names
countData <- countData[, new_order]
colnames(countData)

#separate host and symbiont counts 
head(grep("Sym.*", row.names(countData),ignore.case = T)) #identify where symbiont genes start
#832 833 834 835 836 837
new_order = sort(rownames(countData)) #reorder rows based on row names
countData <- countData[new_order, ]
countData <- countData[!grepl("^_", rownames(countData)), , drop = FALSE]
head(grep("Sym.*", row.names(countData),ignore.case = T)) #identify where symbiont genes start in new order
#45671 45672 45673 45674 45675 45676

#create countData for host and symbiont
countData_host <- countData[1:45670,]
countData_sym <- countData[45671:79126, ]

#overview of total raw counts
totalCounts=colSums(countData)
totalCounts_host=colSums(countData_host)
max(totalCounts_host) #699708
min(totalCounts_host) #79402
mean(totalCounts_host)
totalCounts_sym=colSums(countData_sym)
max(totalCounts_sym) #95119
min(totalCounts_sym) #964

#visualize total raw counts via barplots
barplot(totalCounts_host, col=c("deepskyblue","deepskyblue","deepskyblue",
                                "deepskyblue","deepskyblue","deepskyblue",
                                "deepskyblue","deepskyblue", "deepskyblue3",
                                "deepskyblue3", "deepskyblue3", "deepskyblue3",
                                "deepskyblue3", "deepskyblue3", "deepskyblue3",
                                "lightskyblue4", "lightskyblue4", "lightskyblue4",
                                "lightskyblue4", "lightskyblue4", "lightskyblue4",
                                "deepskyblue4", "deepskyblue4", "deepskyblue4",
                                "deepskyblue4", "deepskyblue4", "cadetblue1",
                                "cadetblue1", "cadetblue1", "cadetblue1", "cadetblue1",
                                "cadetblue1", "aquamarine",  "aquamarine", "aquamarine",
                                "aquamarine", "aquamarine", "aquamarine", "aquamarine",
                                "lightpink", "lightpink", "lightpink", "lightpink",
                                "lightpink", "lightpink", "lightpink", "lightsalmon",
                                "lightsalmon", "lightsalmon", "lightsalmon",
                                "lightsalmon", "lightsalmon", "lightsalmon", "brown1",
                                "brown1", "brown1", "brown1", "brown1", "brown1",
                                "coral3", "coral3",  "coral3", "coral3", "coral3",
                                "coral3", "hotpink1", "hotpink1", "hotpink1",
                                "hotpink1", "hotpink1", "hotpink1", "deeppink",
                                "deeppink", "deeppink", "deeppink", "deeppink", "deeppink"),
        ylab="Raw Counts", xlab="Apo and Sym Astrangia host in Response to Heat/Light Stress",
        main="Host")

barplot(totalCounts_sym, col=c("deepskyblue","deepskyblue","deepskyblue",
                               "deepskyblue","deepskyblue","deepskyblue",
                               "deepskyblue","deepskyblue", "deepskyblue3",
                               "deepskyblue3", "deepskyblue3", "deepskyblue3",
                               "deepskyblue3", "deepskyblue3", "deepskyblue3",
                               "lightskyblue4", "lightskyblue4", "lightskyblue4",
                               "lightskyblue4", "lightskyblue4", "lightskyblue4",
                               "deepskyblue4", "deepskyblue4", "deepskyblue4",
                               "deepskyblue4", "deepskyblue4", "cadetblue1",
                               "cadetblue1", "cadetblue1", "cadetblue1", "cadetblue1",
                               "cadetblue1", "aquamarine",  "aquamarine", "aquamarine",
                               "aquamarine", "aquamarine", "aquamarine", "aquamarine",
                               "lightpink", "lightpink", "lightpink", "lightpink",
                               "lightpink", "lightpink", "lightpink", "lightsalmon",
                               "lightsalmon", "lightsalmon", "lightsalmon",
                               "lightsalmon", "lightsalmon", "lightsalmon", "brown1",
                               "brown1", "brown1", "brown1", "brown1", "brown1",
                               "coral3", "coral3",  "coral3", "coral3", "coral3",
                               "coral3", "hotpink1", "hotpink1", "hotpink1",
                               "hotpink1", "hotpink1", "hotpink1", "deeppink",
                               "deeppink", "deeppink", "deeppink", "deeppink", "deeppink"),
        ylab="Raw Counts", xlab="In-hospite symbionts in Response to Heat/Light Stress",
        main="Symbiont")


#----------Prepare colData----------
coral=c(colnames(countData))
coralData <- data.frame(coral)
stat_treat=c("ACC", "ACC", "ACC", "ACC", "ACC", "ACC", "ACC", "ACC", "ACH", "ACH", "ACH", "ACH",
             "ACH", "ACH", "ACH", "ADC", "ADC", "ADC", "ADC", "ADC", "ADC", "ADH", "ADH", "ADH",
             "ADH", "ADH", "AHC", "AHC", "AHC", "AHC", "AHC", "AHC", "AHH", "AHH", "AHH", "AHH",
             "AHH", "AHH", "AHH", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCH", "SCH",
             "SCH", "SCH", "SCH", "SCH", "SCH", "SDC", "SDC", "SDC", "SDC", "SDC", "SDC", "SDH",
             "SDH", "SDH", "SDH","SDH", "SDH", "SHC", "SHC", "SHC", "SHC", "SHC", "SHC", "SHH",
             "SHH", "SHH", "SHH", "SHH", "SHH")
stat_treatData <- data.frame(stat_treat)

length(grep("A", stat_treat)) #39
length(grep("S", stat_treat)) #38
symstat <- rep(c("Apo","Sym"),times=c(39,38))
statData <- data.frame(symstat)
statData

length(grep("AC", stat_treat)) #apo control light: 15
length(grep("AD", stat_treat)) #apo no light: 11
length(grep("AH", stat_treat)) #apo control light: 13
length(grep("SC", stat_treat)) #sym control light 14 
length(grep("SD", stat_treat)) #sym control light 12
length(grep("SH", stat_treat)) #sym control light 12 
light <- rep(c("Control", "Dark", "High", "Control", "Dark", "High"),times=c(15,11,13,14,12,12))
lightData <- as.data.frame(light)
lightData

length(grep("ACC", stat_treat)) #8
length(grep("ACH", stat_treat)) #7
length(grep("ADC", stat_treat)) #6
length(grep("ADH", stat_treat)) #5
length(grep("AHC", stat_treat)) #6
length(grep("AHH", stat_treat)) #7
length(grep("SCC", stat_treat)) #7
length(grep("SCH", stat_treat)) #7
length(grep("SDC", stat_treat)) #6
length(grep("SDH", stat_treat)) #6
length(grep("SHC", stat_treat)) #6
length(grep("SHH", stat_treat)) #6
heat <- rep(c("Control", "Heated", "Control", "Heated", "Control", "Heated", "Control", "Heated","Control", "Heated","Control", "Heated"), times=c(8,7,6,5,6,7,7,7,6,6,6,6))
heatData <- as.data.frame(heat)
heatData

treat <- rep(c("CL_C", "CL_H", "NL_C", "NL_H", "HL_C", "HL_H","CL_C", "CL_H", "NL_C", "NL_H", "HL_C", "HL_H"), times=c(8,7,6,5,6,7,7,7,6,6,6,6))
treatData <- as.data.frame(treat)
treatData

genet <- sub(".*_", "", colnames(countData))
genetData <- as.data.frame(genet)
genetData

colData <- cbind(genetData,stat_treatData, statData, treatData, lightData, heatData)
row.names(colData) = coral

df <- read.csv("../Data/Metadata_Astrangia_Clean_May1.csv")
rownames(df) <- df$Sample_ID
colData$tank <- df[rownames(colData), "Tank"]

colData <- colData %>%
  mutate_at(vars(genet, tank, stat_treat, symstat, treat, light, heat), factor)


#----------DESeq----------
#Overview without specifying design for PCA and WGCNA
dds_host <- DESeqDataSetFromMatrix(countData = countData_host,
                                   colData = colData,
                                   design = ~1)
dds_host #dim: 45670 77
dds_host= DESeq(dds_host)
plotDispEsts(dds_host, main="Dispersion Uptake")

#Extract samples from control temperature and control light conditions for baseline comparison
countData_control <- countData_host[, grepl("^(ACC|SCC)_", colnames(countData_host))]
colData_control <- colData[grepl("^(ACC|SCC)_", rownames(colData)), ]
dds_host.symstat_control <- DESeqDataSetFromMatrix(countData = countData_control,
                                                   colData = colData_control,
                                                   design = ~symstat)
dds_host.symstat_control
dds_host.symstat_control= DESeq(dds_host.symstat_control)
plotDispEsts(dds_host.symstat_control, main="Dispersion Uptake")

#Treatment as the fixed effect for comparisons among treatments
dds_host.treat <- DESeqDataSetFromMatrix(countData = countData_host,
                                         colData = colData,
                                         design = ~treat)
dds_host.treat #dim: 45670 77
dds_host.treat= DESeq(dds_host.treat)
plotDispEsts(dds_host, main="Dispersion Uptake")

#Separate white and brown corals for subsequent analyses
countData_host.apo <- countData_host[,1:39]
countData_host.sym <- countData_host[,40:77]
colData.apo <- colData[1:39,]
colData.sym <- colData[40:77,]

dds_host.apo <- DESeqDataSetFromMatrix(countData = countData_host.apo,
                                       colData = colData.apo,
                                       design = ~treat)
dds_host.apo #dim: 45670 39
dds_host.apo= DESeq(dds_host.apo)

dds_host.sym <- DESeqDataSetFromMatrix(countData = countData_host.sym,
                                       colData = colData.sym,
                                       design = ~treat)
dds_host.sym #dim: 45670 38
dds_host.sym= DESeq(dds_host.sym)

save(countData_host, colData, countData_host.apo, countData_host.sym, colData.apo, colData.sym,
     dds_host, dds_host.symstat_control, dds_host.treat, dds_host.apo, dds_host.sym, file="dds.RData")

