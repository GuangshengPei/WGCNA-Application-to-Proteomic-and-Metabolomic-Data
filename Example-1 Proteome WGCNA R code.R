#=====================================================================================
#
#  Code chunk 1	Load Proteome data for missing imputation
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
# link all data set in "Table S 1. Raw data of proteomic example in WGCNA.xlsx"
channel=odbcConnectExcel2007("Suppl. Table S 1. Raw data of proteomic example in WGCNA.xlsx")
sqlTables(channel)

# load proteomic data set
ProData=sqlFetch(channel,"Proteome")
rownames(ProData)=ProData[,1];ProData=ProData[,-1]
dim(ProData);names(ProData)

# determine whether data set are complete
complete.cases(ProData)			
md.pattern(ProData)				
# remove peptides with too many missvalues (>90% samples)	
ProData[rowMeans(!is.na(ProData)) <0.9,]       
# Random simulation data 
imp=mice((ProData), m=10, seed=1)		
print(imp)
# linear regression
fit <- with(imp, lm(N48~ E24 + E48 + B24 + B48 + H24 + H48 + S24 + S48 + N24))
# regression result
pooled=pool(fit)
options(digits=3)
summary(pooled)

# predict the missvalue
Data.pre=ProData[is.na(ProData$N48),][,1:9]		
Data.pre=as.matrix(cbind(rep(1,93),Data.pre))
q=pooled$qbar                       
pre=Data.pre%*%q;pre                 
# replace the missing value by PMM imputation
index=is.na(ProData$N48)
ProData$N48[index]=pre               
ProData[index,]
# save the data set after miss value imputation
write.csv(ProData,file="Example-1 Proteome Table 1 Imputed express data.csv")
#=====================================================================================
#
#  Code chunk 2	Load WGCNA package 
#
#=====================================================================================

# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Read in the proteomic data set		
ProData = as.data.frame(t(ProData))
# Take a quick look at what is in the data set:
dim(ProData);names(ProData);

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
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(ProData, powerVector = powers, verbose = 5, networkType = "signed")
# Plot the results:
pdf("Example-1 Proteome Figure 1.pdf",10,5)
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
abline(h=0.85,col="red")
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

softPower = 29;
net = blockwiseModules(ProData, power = softPower, networkType = "signed",
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "Example-1 ProDataTOM", 
                       verbose = 3)
table(net$colors)

#=====================================================================================
#
#  Code chunk 6	Clustering dendrograms for peptides
#
#=====================================================================================

# open a graphics window
pdf("Example-1 Proteome Figure 2.pdf",12,9)
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()					

#=====================================================================================
#
#  Code chunk 7	Modules analysis
#
#=====================================================================================

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
# save(MEs, moduleLabels, moduleColors, geneTree, file = "ProData-networkConstruction-auto.RData")

#=====================================================================================
#
#  Code chunk 8	Peptide significance caculation
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(ProData);
nSamples = nrow(ProData);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(ProData, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#=====================================================================================
#
#  Code chunk 9	Module-trait relationships
#
#=====================================================================================

# open a graphics window
pdf("Example-1 Proteome Figure 3.pdf",5,6)
# sizeGrWindow(10,6)
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
#  Code chunk 10 Correlated eigenpeptide in modules associated with traits
#
#=====================================================================================

# Define variable Biofuel containing the Biofuel column of datTrait
Biofuel = as.data.frame(datTraits$Biofuel);
names(Biofuel) = "Biofuel"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(ProData, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(ProData, Biofuel, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Biofuel), sep="");
names(GSPvalue) = paste("p.GS.", names(Biofuel), sep="");

#=====================================================================================
#
#  Code chunk 11	WGCNA results Output 
#
#=====================================================================================

# Create the starting data frame
geneInfo0 = data.frame(Peptides = names(ProData),
                      moduleColor = moduleColors,
                      geneTraitSignificance,GSPvalue)
# Order modules by their significance for Biofuel
modOrder = order(-abs(cor(MEs, Biofuel, use = "p")));
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
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Biofuel));
geneInfo = geneInfo0[geneOrder, ]
#save WGCNA results
write.csv(geneInfo, file = "Example-1 Proteome Table 2 WGCNA result.csv", row.names = FALSE)

#=====================================================================================
#
#  Code chunk 12	Network heatmap plot visulization
#
#=====================================================================================

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(ProData, power = softPower, networkType = "signed" );
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function

pdf("Example-2 Proteome Figure 4.pdf",9,9)
# sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot")
dev.off()

#=====================================================================================
#
#  Code chunk 13	Protein enrichment analysis
#
#=====================================================================================

# Load the goseq package (please install "goseq" package before start)
library(goseq)
# load WGCNA results 
WGCNA=sqlFetch(channel,"Module")
rownames(WGCNA)=WGCNA[,1];WGCNA=WGCNA[,-1]

# load intereted module
Mod_red=WGCNA$red;names(Mod_red)=rownames(WGCNA);
Mod_pink=WGCNA$pink;names(Mod_pink)=rownames(WGCNA);
Mod_tan=WGCNA$tan;names(Mod_tan)=rownames(WGCNA);
Mod_magenta=WGCNA$magenta;names(Mod_magenta)=rownames(WGCNA);

# load length bias information for caculation then replace it by mean value  
Mybias.data=WGCNA$Length
pwf_red=nullp(DEgenes=Mod_red,bias.data=Mybias.data);pwf_red$pwf=mean(pwf_red$pwf)
pwf_pink=nullp(DEgenes=Mod_pink,bias.data=Mybias.data);pwf_pink$pwf=mean(pwf_pink$pwf)
pwf_tan=nullp(DEgenes=Mod_tan,bias.data=Mybias.data);pwf_tan$pwf=mean(pwf_tan$pwf)
pwf_magenta=nullp(DEgenes=Mod_magenta,bias.data=Mybias.data);pwf_magenta$pwf=mean(pwf_magenta$pwf)

# load pathway annotation information 
Met2pathway=sqlFetch(channel,"Pathway")

# pathway enrichment analysis for interested modules
Enrich_red=goseq(pwf=pwf_red,gene2cat=Met2pathway);Enrich_red$Module="red"
Enrich_pink=goseq(pwf=pwf_pink,gene2cat=Met2pathway);Enrich_pink$Module="pink"
Enrich_tan=goseq(pwf=pwf_tan,gene2cat=Met2pathway);Enrich_tan$Module="tan"
Enrich_magenta=goseq(pwf=pwf_magenta,gene2cat=Met2pathway);Enrich_magenta$Module="magenta"

# combine enrichment result in different modules 
Enrichment=rbind(Enrich_red, Enrich_pink,Enrich_tan,Enrich_magenta)
Enrich_module=Enrichment[Enrichment$over_represented_pvalue<0.05,]
Enrich_module=Enrich_module[,c(6,1:5)]
names(Enrich_module)=c("Module","Pathway","over_represented_pvalue","under_represented_pvalue","numMolInPath","numInPath")

head(Enrich_module)
# save the pathway enrichment results
write.csv(Enrich_module,"Example-1 Proteome Table 3 Module enrichment results.csv",row.names=FALSE)

#Close the link to "Table S 1. Raw data of proteomic example in WGCNA.xlsx"
close(channel)  

#=====================================================================================
#
#  Code chunk 14 Exporting Network to external visualization software
#
#=====================================================================================

# Cytoscape or VisANT plot and software
probes = names(ProData)
TOM=1-dissTOM
modules = c(moduleColors[!duplicated(moduleColors)])

for (i in seq(1, nrow(table(modules)))){
inModule=is.finite(match(moduleColors,modules[i]))
modProbes=probes[inModule]
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
	edgeFile=paste("Example-1 Proteome network CytoEdge ",paste(modules[i],collapse="-"),".txt",sep=""),
	nodeFile=paste("Example-1 Proteome network CytoNode ",paste(modules[i],collapse="-"),".txt",sep=""),
	weighted = TRUE, threshold = 0.02)
vis = exportNetworkToVisANT(modTOM,
	file=paste("Example-1 Proteome network VisANTInput ", modules[i],  sep=""),
	weighted=TRUE,threshold = 0.02)
}

#=====================================================================================
#
#  Code chunk 15 Exporting top hub Network to external visualization software
#
#=====================================================================================

# choose red module
module = c("red");
inModule = (moduleColors==module);
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# choose the top 10 node with highest conncetion
nTop = 10;
IMConn = softConnectivity(ProData[, probes[is.finite(match(moduleColors,module))]]);
top = (rank(-IMConn) <= nTop)

# network file for cystoscape software 
cyt = exportNetworkToCytoscape(modTOM[top, top],
	edgeFile=paste("Example-1 Proteome network CytoEdge ",paste(module,collapse="-"),"-top10.txt",sep=""),
	nodeFile=paste("Example-1 Proteome network CytoNode ",paste(module,collapse="-"),"-top10.txt",sep=""),
	weighted = TRUE, threshold = 0.21)
# network file for visant software 
vis = exportNetworkToVisANT(modTOM[top, top],
  file = paste("Example-1 Proteome network VisANTInput ", module, "-top10.txt", sep=""),
  weighted = TRUE,threshold = 0.21)
  
# choose pink module
module = c("pink");
inModule = (moduleColors==module);
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# choose the top 10 node with highest conncetion
nTop = 10;
IMConn = softConnectivity(ProData[, probes[is.finite(match(moduleColors,module))]]);
top = (rank(-IMConn) <= nTop)
cyt = exportNetworkToCytoscape(modTOM[top, top],
	edgeFile=paste("Example-1 Proteome network CytoEdge ",paste(module,collapse="-"),"-top10.txt",sep=""),
	nodeFile=paste("Example-1 Proteome network CytoNode ",paste(module,collapse="-"),"-top10.txt",sep=""),
	weighted = TRUE, threshold = 0.11)
	
vis = exportNetworkToVisANT(modTOM[top, top],
  file = paste("Example-1 Proteome network VisANTInput ", module, "-top10.txt", sep=""),
  weighted = TRUE,threshold = 0.11)
  
# choose tan module
module = c("tan");
inModule = (moduleColors==module);
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# choose the top 10 node with highest conncetion
nTop = 10;
IMConn = softConnectivity(ProData[, probes[is.finite(match(moduleColors,module))]]);
top = (rank(-IMConn) <= nTop)

cyt = exportNetworkToCytoscape(modTOM[top, top],
	edgeFile=paste("Example-1 Proteome network CytoEdge ",paste(module,collapse="-"),"-top10.txt",sep=""),
	nodeFile=paste("Example-1 Proteome network CytoNode ",paste(module,collapse="-"),"-top10.txt",sep=""),
	weighted = TRUE, threshold = 0.09)
vis = exportNetworkToVisANT(modTOM[top, top],
  file = paste("Example-1 Proteome network VisANTInput ", module, "-top10.txt", sep=""),
  weighted = TRUE,threshold = 0.09)

# choose megenta module 
module = c("magenta"); 
inModule = (moduleColors==module);
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# choose the top 10 node with highest conncetion
nTop = 10;
IMConn = softConnectivity(ProData[, probes[is.finite(match(moduleColors,module))]]);
top = (rank(-IMConn) <= nTop)
cyt = exportNetworkToCytoscape(modTOM[top, top],
	edgeFile=paste("Example-1 Proteome network CytoEdge ",paste(module,collapse="-"),"-top10.txt",sep=""),
	nodeFile=paste("Example-1 Proteome network CytoNode ",paste(module,collapse="-"),"-top10.txt",sep=""),
	weighted = TRUE, threshold = 0.12)
vis = exportNetworkToVisANT(modTOM[top, top],
  file = paste("Example-1 Proteome network VisANTInput ", module, "-top10.txt", sep=""),
  weighted = TRUE,threshold = 0.12)
  
#=====================================================================================
#
#  Code chunk 16 Module membership vs. peptide significance
#
#=====================================================================================

# choose red module vs. biofuels
module = "red"								
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
# pdf(file = "Plots/MMvsGS-1.pdf", width=7, height=7);
par(mfrow = c(1,1));
verboseScatterplot(geneModuleMembership[moduleGenes, column],
geneTraitSignificance[moduleGenes, 1],
xlab = paste("Module Membership in", module, "module"),
ylab = "Peptide significance for biofuel response",
main = paste("Module membership vs. biofuel response\n"),
cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, col = module,xlim=c(-1,1),ylim=c(-1,1))
# If plotting into a file, close it
dev.off(); 

# choose pink module vs. biofuels
module = "pink"								
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
# pdf(file = "Plots/MMvsGS-1.pdf", width=7, height=7);
par(mfrow = c(1,1));
verboseScatterplot(geneModuleMembership[moduleGenes, column],
geneTraitSignificance[moduleGenes, 1],
xlab = paste("Module Membership in", module, "module"),
ylab = "Peptide significance for biofuel response",
main = paste("Module membership vs. biofuel response\n"),
cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, col = module,xlim=c(-1,1),ylim=c(-1,1))
# If plotting into a file, close it
dev.off(); 
 
# choose tan module vs. biofuels 
module = "tan"								
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
# pdf(file = "Plots/MMvsGS-1.pdf", width=7, height=7);
par(mfrow = c(1,1));
verboseScatterplot(geneModuleMembership[moduleGenes, column],
geneTraitSignificance[moduleGenes, 1],
xlab = paste("Module Membership in", module, "module"),
ylab = "Peptide significance for biofuel response",
main = paste("Module membership vs. biofuel response\n"),
cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, col = module,xlim=c(-1,1),ylim=c(-1,1))
# If plotting into a file, close it

# choose magenta module vs. biofuels 
module = "magenta"								
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
# pdf(file = "Plots/MMvsGS-1.pdf", width=7, height=7);
par(mfrow = c(1,1));
verboseScatterplot(geneModuleMembership[moduleGenes, column],
geneTraitSignificance[moduleGenes, 1],
xlab = paste("Module Membership in", module, "module"),
ylab = "Peptide significance for biofuel response",
main = paste("Module membership vs. biofuel response\n"),
cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, col = module,xlim=c(-1,1),ylim=c(-1,1))
# If plotting into a file, close it