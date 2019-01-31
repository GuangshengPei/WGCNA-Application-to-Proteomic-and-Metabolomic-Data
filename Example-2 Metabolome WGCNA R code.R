#=====================================================================================
#
#  Code chunk 1	Load Metabolome data for missing imputation
#
#=====================================================================================

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
# set workding directory where supplement file located 
workingDir = " ";
setwd(workingDir); 

# Load the mice package	(please install "mice" package before start)
# Load the RODBC package (please install "RODBC" package before start)
library("mice")
library("RODBC")
# link all data set in "Table S 2. Raw data of metabolomic example in WGCNA.xlsx "
channel=odbcConnectExcel2007("Suppl. Table S 2. Raw data of metabolomic example in WGCNA.xlsx")
sqlTables(channel)

# load metabolomic data set
MetData=sqlFetch(channel,"Metabolome")
rownames(MetData)=MetData[,1];MetData=MetData[,-1]
dim(MetData);names(MetData)

# determine whether data set are complete
complete.cases(MetData)			
md.pattern(MetData)	
# remove metabolites with too many missvalues (>90% samples)	
MetData[rowMeans(!is.na(MetData)) <0.9,] 

#replace the missing value by a minimum value of the metabolite quantified
Pre.MetData=MetData
Pre.MetData[is.na(Pre.MetData)]=0
Pre.MetData
Min.replace=apply(Pre.MetData,1,min)
Replace.matrix=matrix(rep(matrix(Min.replace),ncol(MetData)),nrow(MetData))
index=is.na(MetData)
MetData[index]=Replace.matrix[index]
MetData[index]

# save the data set after miss value imputation
write.csv(MetData,file="Example-2 Metabolome Table 1 Imputed express data.csv")

#=====================================================================================
#
#  Code chunk 2	Load WGCNA package 
#
#=====================================================================================

# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Read in the metabolic data set			
MetData = as.data.frame(t(MetData))
# Take a quick look at what is in the data set:
dim(MetData);names(MetData);

#=====================================================================================
#
#  Code chunk 3	Load traits data
#
#=====================================================================================

datTraits=sqlFetch(channel,"Trait")
rownames(datTraits)=datTraits[,1];datTraits=datTraits[,-1]
dim(datTraits);names(datTraits)

#=====================================================================================
#
#  Code chunk 4 SoftThreshold chosen
#
#=====================================================================================

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(MetData, powerVector = powers, verbose = 5)
# Plot the results:
pdf("Example-2 Metabolome Figure 1.pdf",10,5)
# sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.50,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#=====================================================================================
#
#  Code chunk 5	WGCNA construction
#
#=====================================================================================

softPower = 3;
net = blockwiseModules(MetData, power = softPower, networkType = "unsigned",
                       TOMType = "signed", minModuleSize = 5,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "Example-2 MetDataTOM", 
                       verbose = 3)
table(net$colors)

#=====================================================================================
#
#  Code chunk 6	Clustering dendrograms for metabolites
#
#=====================================================================================

# open a graphics window
pdf("Example-2 Metabolome Figure 2.pdf",12,9)
# sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = names(MetData), hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()				

#=====================================================================================
#
#  Code chunk 7 Modules analysis
#
#=====================================================================================

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
# save(MEs, moduleLabels, moduleColors, geneTree, file = "MetData-networkConstruction-auto.RData")

#=====================================================================================
#
#  Code chunk 8 Metabolite significance caculation
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(MetData);
nSamples = nrow(MetData);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(MetData, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#=====================================================================================
#
#  Code chunk 9	Module-trait relationships
#
#=====================================================================================

# open a graphics window
pdf("Example-2 Metabolome Figure 3.pdf",10,6)
#sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()		   
			   
#=====================================================================================
#
#  Code chunk 10 Correlated eigenmetabolite in modules associated with traits
#
#=====================================================================================

# Define variable Astaxanthin containing the Astaxanthin column of datTrait
Astaxanthin = as.data.frame(datTraits$Astaxanthin);
names(Astaxanthin) = "Astaxanthin"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(MetData, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(MetData, Astaxanthin, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Astaxanthin), sep="");
names(GSPvalue) = paste("p.GS.", names(Astaxanthin), sep="");

# Define variable Oxygenant containing the Oxygenant column of datTrait
Oxygenant = as.data.frame(datTraits$Oxygenant);
names(Oxygenant) = "Oxygenant"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(MetData, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance2 = as.data.frame(cor(MetData, Oxygenant, use = "p"));
GSPvalue2 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance2), nSamples));
names(geneTraitSignificance2) = paste("GS.", names(Oxygenant), sep="");
names(GSPvalue2) = paste("p.GS.", names(Oxygenant), sep="");

# Define variable Signal_transducer containing the Signal_transducer column of datTrait
Signal_transducer = as.data.frame(datTraits$Signal_transducer);
names(Signal_transducer) = "Signal_transducer"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(MetData, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance3 = as.data.frame(cor(MetData, Signal_transducer, use = "p"));
GSPvalue3 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance3), nSamples));
names(geneTraitSignificance3) = paste("GS.", names(Signal_transducer), sep="");
names(GSPvalue3) = paste("p.GS.", names(Signal_transducer), sep="");

#=====================================================================================
#
#  Code chunk 11	WGCNA results Output 
#
#=====================================================================================

# Create the starting data frame
geneInfo0 = data.frame(Metabolites = names(MetData),
                      moduleColor = moduleColors,
                      geneTraitSignificance,GSPvalue,
                      geneTraitSignificance2,GSPvalue2,
					  geneTraitSignificance3,GSPvalue3
					  )
# Order modules by their significance for Astaxanthin
modOrder = order(-abs(cor(MEs, Astaxanthin, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Astaxanthin));
geneInfo = geneInfo0[geneOrder, ]
#save WGCNA results
write.csv(geneInfo, file = "Example-2 Metabolome Table 2 WGCNA result.csv", row.names = FALSE)


#=====================================================================================
#
#  Code chunk 12	Network heatmap plot visulization
#
#=====================================================================================

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(MetData, power = softPower);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^3;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function

pdf("Example-2 Metabolome Figure 4.pdf",9,9)
#sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot")
dev.off()

#=====================================================================================
#
#  Code chunk 13	Metabolite enrichment analysis
#
#=====================================================================================

# Load the goseq package (please install "goseq" package before start)
library(goseq)
# load pathway annotation information 
Met2pathway=sqlFetch(channel,"Pathway")
# enrichment analysis for metabolite in blue module
pwf_blue=data.frame(geneInfo[,c(1,2)],blue=0,bias.data=100,pwf=0.2)
pwf_blue$blue[grep("blue",pwf_blue$moduleColor)]=1
pwf_blue=pwf_blue[,c(3:5)]
names(pwf_blue)=c("DEgenes","bias.data","pwf")
Enrich_blue=goseq(pwf=pwf_blue,gene2cat=Met2pathway);Enrich_blue$Module="blue"

# enrichment analysis for metabolite in green module
pwf_green=data.frame(geneInfo[,c(1,2)],green=0,bias.data=100,pwf=0.2)
pwf_green$green[grep("green",pwf_green$moduleColor)]=1
pwf_green=pwf_green[,c(3:5)]
names(pwf_green)=c("DEgenes","bias.data","pwf")
Enrich_green=goseq(pwf=pwf_green,gene2cat=Met2pathway);Enrich_green$Module="green"

# enrichment analysis for metabolite in turquoise module
pwf_turquoise=data.frame(geneInfo[,c(1,2)],turquoise=0,bias.data=100,pwf=0.2)
pwf_turquoise$turquoise[grep("turquoise",pwf_turquoise$moduleColor)]=1
pwf_turquoise=pwf_turquoise[,c(3:5)]
names(pwf_turquoise)=c("DEgenes","bias.data","pwf")
Enrich_turquoise=goseq(pwf=pwf_turquoise,gene2cat=Met2pathway);Enrich_turquoise$Module="turquoise"

# enrichment analysis for metabolite in black module
pwf_black=data.frame(geneInfo[,c(1,2)],black=0,bias.data=100,pwf=0.2)
pwf_black$black[grep("black",pwf_black$moduleColor)]=1
pwf_black=pwf_black[,c(3:5)]
names(pwf_black)=c("DEgenes","bias.data","pwf")
Enrich_black=goseq(pwf=pwf_black,gene2cat=Met2pathway);Enrich_black$Module="black"

# combine enrichment result in different modules 
Enrichment=rbind(Enrich_blue,Enrich_green,Enrich_turquoise,Enrich_black)
Enrich_module=Enrichment[Enrichment$over_represented_pvalue<0.05,]
Enrich_module=Enrich_module[,c(6,1:5)]
names(Enrich_module)=c("Module","Pathway","over_represented_pvalue","under_represented_pvalue","numMolInPath","numInPath")

head(Enrich_module)
# save the pathway enrichment results
write.csv(Enrich_module,"Example-2 Metabolome Table 3 Module enrichment results.csv",row.names=FALSE)

# Close the link to "Table S 1. Raw data of proteomic example in WGCNA.xlsx"
close(channel) 

# check interested metaboites
Met2pathway[Met2pathway$Pathway=="map01040 Biosynthesis of unsaturated fatty acids",1]
a=intersect(Met2pathway[Met2pathway$Pathway=="map01040 Biosynthesis of unsaturated fatty acids",1],row.names(geneInfo[geneInfo$moduleColor=="blue",]))

Met2pathway[Met2pathway$Pathway=="map00561 Glycerolipid metabolism",1]
intersect(Met2pathway[Met2pathway$Pathway=="map00561 Glycerolipid metabolism",1],row.names(geneInfo[geneInfo$moduleColor=="green",]))

Met2pathway[Met2pathway$Pathway=="map00564 Glycerophospholipid metabolism",1]
intersect(Met2pathway[Met2pathway$Pathway=="map00564 Glycerophospholipid metabolism",1],row.names(geneInfo[geneInfo$moduleColor=="green",]))

#=====================================================================================
#
#  Code chunk 14 Exporting Network to external visualization software
#
#=====================================================================================

# Cytoscape or VisANT plot and software
probes = names(MetData)
TOM=1-dissTOM
modules = c(moduleColors[!duplicated(moduleColors)])

for (i in seq(1, nrow(table(modules)))){
inModule=is.finite(match(moduleColors,modules[i]))
modProbes=probes[inModule]
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
	edgeFile=paste("Example-2 Metabolome network CytoEdge ",paste(modules[i],collapse="-"),".txt",sep=""),
	nodeFile=paste("Example-2 Metabolome network CytoNode ",paste(modules[i],collapse="-"),".txt",sep=""),
	weighted = TRUE, threshold = 0.02)	
vis = exportNetworkToVisANT(modTOM,
	file=paste("Example-2 Metabolome network VisANTInput ", modules[i],  sep=""),
	weighted=TRUE,threshold = 0.02)
}

#=====================================================================================
#
#  Code chunk 15 Exporting top hub Network to external visualization software
#
#=====================================================================================

# choose blue module
module = "blue";
inModule = (moduleColors==module);
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# choose the top 10 node with highest conncetion
nTop = 10;
IMConn = softConnectivity(MetData[, probes[is.finite(match(moduleColors,module))]]);
top = (rank(-IMConn) <= nTop)

# network file for cystoscape software 
cyt = exportNetworkToCytoscape(modTOM[top, top],
	edgeFile=paste("Example-2 Metabolome network CytoEdge ",paste(module,collapse="-"),"-top10.txt",sep=""),
	nodeFile=paste("Example-2 Metabolome network CytoNode ",paste(module,collapse="-"),"-top10.txt",sep=""),
	weighted = TRUE, threshold = 0.12)
# network file for visant software 	
vis = exportNetworkToVisANT(modTOM[top, top],
  file = paste("Example-2 Metabolome network VisANTInput ", module, "-top10.txt", sep=""),
  weighted = TRUE,threshold = 0.12)
  
# choose green module
module = "green";
inModule = (moduleColors==module);
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# choose the top 10 node with highest conncetion
nTop = 10;
IMConn = softConnectivity(MetData[, probes[is.finite(match(moduleColors,module))]]);
top = (rank(-IMConn) <= nTop)
cyt = exportNetworkToCytoscape(modTOM[top, top],
	edgeFile=paste("Example-2 Metabolome network CytoEdge ",paste(module,collapse="-"),"-top10.txt",sep=""),
	nodeFile=paste("Example-2 Metabolome network CytoNode ",paste(module,collapse="-"),"-top10.txt",sep=""),
	weighted = TRUE, threshold = 0.06)
	
vis = exportNetworkToVisANT(modTOM[top, top],
  file = paste("Example-2 Metabolome network VisANTInput ", module, "-top10.txt", sep=""),
  weighted = TRUE,threshold = 0.06)
  
# choose tuquoise module
module = "turquoise";
inModule = (moduleColors==module);
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# choose the top 10 node with highest conncetion
nTop = 10;
IMConn = softConnectivity(MetData[, probes[is.finite(match(moduleColors,module))]]);
top = (rank(-IMConn) <= nTop)

cyt = exportNetworkToCytoscape(modTOM[top, top],
	edgeFile=paste("Example-2 Metabolome network CytoEdge ",paste(module,collapse="-"),"-top10.txt",sep=""),
	nodeFile=paste("Example-2 Metabolome network CytoNode ",paste(module,collapse="-"),"-top10.txt",sep=""),
	weighted = TRUE, threshold = 0.2)
vis = exportNetworkToVisANT(modTOM[top, top],
  file = paste("Example-2 Metabolome network VisANTInput ", module, "-top10.txt", sep=""),
  weighted = TRUE,threshold = 0.2)
  
#=====================================================================================
#
#  Code chunk 16 Module membership vs. metabolite significance
#
#=====================================================================================

# choose blue module vs. Astaxanthin accumulation
module = "blue"								
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
# pdf(file = "Plots/MMvsGS-1.pdf", width=7, height=7);
par(mfrow = c(1,1));
verboseScatterplot(geneModuleMembership[moduleGenes, column],
geneTraitSignificance[moduleGenes, 1],
xlab = paste("Module Membership in", module, "module"),
ylab = "Metabolite significance for Astaxanthin accumulation",
main = paste("Module membership vs Astaxanthin accumulation\n"),
cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, col = module,xlim=c(-1,1),ylim=c(-1,1))
# If plotting into a file, close it
dev.off(); 

# choose green module vs. Oxygenant treatment
module = "green"								
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
# pdf(file = "Plots/MMvsGS-1.pdf", width=7, height=7);
par(mfrow = c(1,1));
verboseScatterplot(geneModuleMembership[moduleGenes, column],
geneTraitSignificance2[moduleGenes, 1],
xlab = paste("Module Membership in", module, "module"),
ylab = "Metabolite significance for Oxygenant response",
main = paste("Module membership vs Oxygenant response \n"),
cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, col = module,xlim=c(-1,1),ylim=c(-1,1))
# If plotting into a file, close it
dev.off(); 

# choose turquoise module vs. Signal transducer treatment
module = "turquoise"								
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
# pdf(file = "Plots/MMvsGS-1.pdf", width=7, height=7);
par(mfrow = c(1,1));
verboseScatterplot(geneModuleMembership[moduleGenes, column],
geneTraitSignificance3[moduleGenes, 1],
xlab = paste("Module Membership in", module, "module"),
ylab = "Metabolite significance for Signal transducer response",
main = paste("Module membership vs Signal transducer response\n"),
cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, col = module,xlim=c(-1,1),ylim=c(-1,1))
# If plotting into a file, close it
dev.off(); 
