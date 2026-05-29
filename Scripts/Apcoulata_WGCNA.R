#Astrangia poculata heat & light: WGCNA
#Light environment shapes divergent coral thermal responses across symbiotic states
#This script was used to generate Supplementary Fig. S2

library(DESeq2)
library(ggplot2)
library(WGCNA)
library(flashClust)


#----------Prepare filtered data for WGCNA----------
load("dds.RData")

res<- results(dds_host)
res3<-res[res$baseMean>3, ]
dim(res) #45670 6
dim(res3) #10736 6

vst <- vst(dds_host, blind=TRUE)
head(vst)
vst_wg=(assay(vst))
head(vst_wg)
nrow(vst_wg)
#45670
vstFiltered=(assay(vst))[(rownames((assay(vst))) %in% rownames(res3)),]
nrow(vstFiltered)
#10736

write.csv(vstFiltered,file="../Data/Astrangia_wgcna_allgenes.csv",quote=F,row.names=T)


#----------WGCNA----------
#Data input and cleaning
options(stringsAsFactors=FALSE)
allowWGCNAThreads()

dat=read.csv("../Data/Astrangia_wgcna_allgenes.csv")
head(dat) 
rownames(dat)<-dat$X
head(dat)
dat$X=NULL
head(dat)
names(dat)
nrow(dat)
#10736
datExpr0 = as.data.frame(t(dat))

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK #TRUE
#If false, run the script below
#if (!gsg$allOK)
#{if (sum(!gsg$goodGenes)>0)
#  printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse= ", ")));
#  if (sum(!gsg$goodSamples)>0)
#    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse=", ")))
#  datExpr0= datExpr0[gsg$goodSamples, gsg$goodGenes]
#}

dim(datExpr0) 
#10736  No change because there were never any outliers

#Outlier detection incorporated into trait measures 
traitData= read.csv("../Data/Phys_traits_WGCNA.csv", row.names=1)
dim(traitData) #79 15
head(traitData)
names(traitData)
traitData <- traitData[sort(rownames(traitData)), ]
rownames(traitData)
traitData<- traitData[(rownames(traitData) %in% rownames(datExpr0)),]
dim(traitData) #77 15

# Form a data frame analogous to expression data that will hold the clinical traits.
dim(datExpr0) #77 10736
rownames(datExpr0)
rownames(traitData)=rownames(datExpr0)
traitData$Sample= NULL 
# datTraits=allTraits
datTraits=traitData

table(rownames(datTraits)==rownames(datExpr0)) #should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraits)
head(datExpr0)

#sample dendrogram and trait heat map showing outliers
A=adjacency(t(datExpr0),type="signed")
# this calculates the whole network connectivity we choose signed because we care about direction of gene expression
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
thresholdZ.k=-2.5 # often -2.5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation where red indicates high values
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits))
datColors=data.frame(outlierC=outlierColor,traitColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,main="Sample dendrogram and trait heatmap")

#Remove outlying samples from expression and trait data
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
datExpr0=datExpr0[!remove.samples,]
datTraits=datTraits[!remove.samples,]
dim(datTraits) #75 15

save(datExpr0, datTraits, file="Astrangia_Samples_Traits_ALL.RData")


#----------Network construction and module detection-----------
options(stringsAsFactors = FALSE)
#enableWGCNAThreads use this in base R
allowWGCNAThreads() 
lnames = load(file="Astrangia_Samples_Traits_ALL.RData")

#Figure out proper SFT
# Choose a set of soft-thresholding powers
powers = c(seq(1,14,by=2), seq(15,30, by=0.5)); #may need to adjust these power values to hone in on proper sft value
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, networkType="signed", verbose = 2) #want smallest value, closest to 0.9 (but still under)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower=7 #smallest value to plateau at ~0.85
adjacency=adjacency(datExpr0, power=softPower,type="signed") #must change method type here too!!
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM= TOMsimilarity(adjacency,TOMType = "signed")
dissTOM= 1-TOM

geneTree= flashClust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
#each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes

minModuleSize=90 #we only want large modules
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
table(dynamicMods)
#1    2    3    4    5    6    7    8    9   10   11 
#3880 1275 1120  996  866  720  549  488  283  282  277 
dynamicColors= labels2colors(dynamicMods)
#plot dendrogram and colors underneath, pretty sweet
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")

#Merg modules whose expression profiles are very similar
#calculate eigengenes
MEList= moduleEigengenes(datExpr0, colors= dynamicColors,softPower = 5)
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")

save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_Astrangia_nomerge.RData")

lnames = load(file = "Network_Astrangia_nomerge.RData")
#plot
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

MEDissThres= 0.2
abline(h=MEDissThres, col="red")

merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors= merge$colors
mergedMEs= merge$newMEs

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

#save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_signed_0.2.RData")


#----------Relating modules to traits and finding important genes----------
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "Astrangia_Samples_Traits_ALL.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Network_signed_0.2.RData");
lnames

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
table(moduleColors)
#black        blue       brown       green greenyellow     magenta 
#549        1275        1120         866         277         283 
#pink      purple         red   turquoise      yellow 
#488         282         720        3880         996 

# Calculate Module Eigengenes and Correlations
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Prepare text matrix: display only where p-value < 0.05
textMatrix = ifelse(moduleTraitPvalue < 0.05,
                    paste(signif(moduleTraitCor, 2), "\n(",
                          signif(moduleTraitPvalue, 1), ")", sep = ""),
                    "")  # Empty string for p-values >= 0.05

dim(textMatrix) = dim(moduleTraitCor)

# Plot the heatmap
quartz()
sizeGrWindow(10, 6)
par(mar = c(10, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
                              xLabels = names(datTraits),
                              yLabels = names(MEs),  # Keep the module names
                              ySymbols = names(MEs),
                              colorLabels = FALSE,
                              colors = blueWhiteRed(50),
                              textMatrix = textMatrix,
                              setStdMargins = FALSE,
                              cex.text = 0.8,
                              zlim = c(-1, 1),
                              main = paste("Module-trait relationships (0.2)"))
dev.off()


#----------Calculate the number of genes in each module----------
moduleCounts = table(moduleColors)

#Define the order manually
custom_order = c("turquoise", "blue", "greenyellow", "magenta", "red", "pink", "yellow", "black", "purple","brown","green")

# Reorder moduleCounts according to custom_order
moduleCounts_ordered = moduleCounts[custom_order]
moduleCounts_ordered
#oduleColors
#turquoise        blue greenyellow     magenta         red        pink 
#3880        1275         277         283         720         488 
#yellow       black      purple       brown       green 
#996         549         282        1120         866 
barplot(moduleCounts_ordered, col = custom_order,
        border = "black", las = 2, ylim=c(0,4000))# Rotate x-axis labels vertically

midpoints <- barplot(moduleCounts_ordered, 
                     col = custom_order,
                     border = "black", 
                     las = 2, 
                     ylim = c(0, 4300))  # Extend y-axis

# Add text labels on top of each bar
text(x = midpoints, 
     y = moduleCounts_ordered, 
     labels = moduleCounts_ordered, 
     pos = 3,      # above the bars
     cex = 0.8)    # font size


#----------Prepare GO input files for modules of interest-----------
vsd <- read.csv("../Data/Astrangia_wgcna_allgenes.csv", row.names=1)
head(vsd)
options(stringsAsFactors=FALSE)
common_samples <- intersect(colnames(vsd), rownames(datExpr0))
vsd <- vsd[,common_samples]
data=t(vsd)
allkME =as.data.frame(signedKME(data, MEs))

whichModule="yellow" # name your color and execute to the end

length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file = file.path("..", "Data", "GO_Input", paste0(whichModule, "_fisher_0.2.csv")),quote = FALSE)
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file = file.path("..", "Data", "GO_Input", paste0(whichModule, "_kME_0.2.csv")),quote = FALSE)