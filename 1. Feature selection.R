#### R Version 4.2.1
#### pacman Version 0.5.1
#### caret Version 6.3-932
#### dplyr Version 1.0.10
#### WGCNA Version 1.71
#### preprocessCor Version 1.58.0
#### extrafont Version 0.19
library(pacman)
library(dplyr)
library(WGCNA)
library(preprocessCore)
library(extrafont)
setwd("")
###################################################################################################################
train=read.table("data_for_FS.csv",sep=",")
colnames(train) <- train[1,]
rownames(train) <- train[,1]
train <- train[-1,-1]
datTraits <- as.data.frame(train[,2])
rownames(datTraits) <- rownames(train)
colnames(datTraits) <- c("cell viability")
datTraits=as.data.frame(lapply(datTraits,as.numeric))
#################################################
###################According to the signature of level3 for SD and MFC
signatures_level3 <- read.table("signatures_A549_24h_cp_level3_7495.csv",sep=",",header=0)
colnames(signatures_level3) <- signatures_level3[1,]
rownames(signatures_level3) <- signatures_level3[,1]
signatures_level3 <- signatures_level3[-1,-1]

#####################obtain of the level 3 data of compounds in training set
signatures_level3_train <- signatures_level3[which(signatures_level3$pert_id %in% train$pert),]
signatures_level3_train <- signatures_level3_train[,-1]
datTraits_level3 <- as.data.frame(signatures_level3_train[,c(0:3)])
rownames(datTraits_level3) <- rownames(signatures_level3_train)
colnames(datTraits_level3) <- c("cell viability","label0.5","label0.8")
signatures_level3_train_t <- as.data.frame(t(signatures_level3_train))
###########SD and MFC
mean <- as.data.frame(apply(signatures_level3_train_t,1,mean))
max <- as.data.frame(apply(signatures_level3_train_t,1,max))
var <-  as.data.frame(apply(signatures_level3_train_t,1,var))
SD <-  as.data.frame(apply(signatures_level3_train_t,1,sd))
MFC <- as.data.frame(max/mean)
colnames(var) <- c("var")
colnames(MFC) <- c("MFC")
colnames(SD) <- c("SD")
rownames(SD) <- rownames(signatures_level3_train_t)
rownames(MFC) <- rownames(signatures_level3_train_t)
#####Obtain of the data after removing the invariant genes
datExpr_t_2 <- cbind(SD,MFC,signatures_level3_train_t)
datExpr_train <- filter(datExpr_t_2,SD>1|MFC>2)
gene_index <- rownames(datExpr_train)
train_level5 <- train[,gene_index]
datExpr <-as.matrix(train_level5)
#####save the matrix after filtering 
train_matrix <- cbind(train[,c(1:5)],datExpr)
write.csv(train_matrix,"Train_variable_genes_matrix.csv")
######################################################################
#####################################################################
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
##sizeGrWindow(9, 5)
svg(filename = "Soft Threshold.svg",width = 15, height = 7)
font_import(pattern = "arial")
loadfonts()
par(family = "Arial", font = 2, cex.axis = 2,mfrow = c(1,2),mar = c(5,5,5,5));
cex1 = 2;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale free topology model fit",type="n", font = 2, cex.axis = 2,cex.lab=2,cex.main=2,
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red",lwd = 2)
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", font=2,cex.axis=2,cex.lab=2,cex.main=2,
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=0.90,col="red",lwd = 2)
axis(1, font = 2)
axis(2, font = 2)
dev.off()
write.csv(sft,"soft threshold data.csv")

cor <- WGCNA::cor
net = blockwiseModules(
  datExpr,
  power = sft$powerEstimate,
  maxBlockSize = 10000,
  TOMType = "unsigned", minModuleSize = 20,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "AS-green-FPKM-TOM",
  verbose = 3
)
cor<-stats::cor
table(net$colors) 
####################################
#############model##################
svg(filename = "WGCNA_model.svg",width = 12, height = 8)
mergedColors = labels2colors(net$colors)
table(mergedColors)
font_import(pattern = "arial")
loadfonts()
par(family = "Arial", font = 2, cex.axis = 1.8,mar = c(2,8,2,2),cex.lab=1.8,cex.main=1.8)
cex1 = 1.8;


# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",setLayout = TRUE, 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE,font=2, guideHang = 0.05,cex.colorLabels = 1.5,mar = c(5,8,5,5), cex.dendroLabels =1.5, 
                    cex.rowText =1.5
                    	)

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
dev.off()
####Module-Traits
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
svg(filename = "WGCNA_model_with_cell_viability.svg",width = 6, height = 8)
font_import(pattern = "arial")
loadfonts()
par(family = "Arial", font = 2, cex.axis = 1.5,mar = c(2,1,2,2),cex.lab=1.5,cex.main=1.5)
cex1 = 1.5;
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 8, 4, 8))
Labelname=c("blue","green","black","yellow","red","brown","turquoise","grey")
labeledHeatmap(Matrix = moduleTraitCor,xLabels = "Cell viability",yLabels = names(MEs),ySymbols = names(MEs),
               colors = blueWhiteRed(50),textMatrix = textMatrix,setStdMargins = FALSE,
               cex.text = 1.3,zlim = c(-1,1),cex.lab=1.2,
               main = paste("Module-trait relationships"))
dev.off()
###############################################################
#################################################################
# Define variable weight containing the weight column of datTrait
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
CB = as.data.frame(datTraits[,1]);
names(CB) = "Cell viability"
geneTraitSignificance = as.data.frame(cor(datExpr, CB, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(CB), sep="");
names(GSPvalue) = paste("p.GS.", names(CB), sep="")
module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;
svg(filename = "WGCNA_blue_Module membership vs. gene significance.svg",width = 6, height = 6)
mergedColors = labels2colors(net$colors)
table(mergedColors)
font_import(pattern = "arial")
loadfonts()
par(family = "Arial", font = 2, cex.axis = 1.8,cex.lab=1.8,cex.main=1.8)
cex1 = 1.8;
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for cell viability",
                   main = paste("Module membership vs. Gene significance\n"),
                   cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5,font = 2, col = module)
dev.off()
module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
svg(filename = "WGCNA_turquoise_Module membership vs. gene significance.svg",width = 6, height = 6)
mergedColors = labels2colors(net$colors)
table(mergedColors)
font_import(pattern = "arial")
loadfonts()
par(family = "Arial", font = 2, cex.axis = 1.8,cex.lab=1.8,cex.main=1.8)
cex1 = 1.8;
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for cell viability",
                   main = paste("Module membership vs. Gene significance\n"),
                   cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5,font = 2, col = module)
dev.off()
####################### blue turquoise
module = "blue"
column = match(module, modNames)
moduleGenes = moduleColors==module
module<-as.data.frame(dimnames(data.frame(datExpr))[[2]][moduleGenes])
names(module)="genename"
MM<-abs(geneModuleMembership[moduleGenes,column])
GS<-abs(geneTraitSignificance[moduleGenes, 1])
MMGS1<-as.data.frame(cbind(MM,GS))
rownames(MMGS1)=module$genename
write.csv(MMGS1, "MMGS_blue_ALL.csv")
module = "turquoise"
column = match(module, modNames)
moduleGenes = moduleColors==module
module<-as.data.frame(dimnames(data.frame(datExpr))[[2]][moduleGenes])
names(module)="genename"
MM<-abs(geneModuleMembership[moduleGenes,column])
GS<-abs(geneTraitSignificance[moduleGenes, 1])
MMGS2<-as.data.frame(cbind(MM,GS))
rownames(MMGS2)=module$genename
write.csv(MMGS2, "MMGS_turquoise_ALL.csv")
MMGS <- rbind(MMGS1,MMGS2)
write.csv(MMGS, "MMGS_blue_turquoise_ALL.csv")
###################################################################################
#############################BP####################################################
library(org.Hs.eg.db)
library(fgsea)
library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)
library(clusterProfiler)
library(stats)
matrix <- cbind(datTraits,datExpr)
###For 50% cell viability prediction
#matrix <- filter(matrix,matrix$`cell.viability`<0.5)
###For 80% cell viability prediction
matrix <- filter(matrix,matrix$`cell.viability`<0.8)
matrix <- matrix[,-1]
rows<-length(matrix[,1])
cols<-length(matrix[1,])
name<-colnames(matrix)
countMatrix<-data.frame(matrix(0,rows,cols))
names(countMatrix)<-c(name)
i=1
for(i in 1:rows){
  geneCount<-list()
  gene <-t(rbind(colnames(matrix),matrix[i,]))
  genelist <-as.double(gene[,2])
  names(genelist) <- as.character(gene[,1])
  genelist <- sort(genelist, decreasing = TRUE)
  ego2 <- gseGO(geneList     = genelist,
                ont          = "BP",  # "BP"、"MF"和"CC"或"ALL"
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = TRUE,
                seed = TRUE,
                by = "fgsea") 
  j=1
  if(is.null(ego2[1])){
    #print("no term")
    next()
  }
  while(!is.na(ego2[j]$ID)){
    
    a<-strsplit(as.character(ego2[j]$core_enrichment),'/')
    countMatrix[i,a[[1]]]=countMatrix[i,a[[1]]]+1
    j=j+1
    cat(i,'drug',j,'terms finished\n')
  }
}
countMatrix <- as.data.frame(t(countMatrix))
countMatrix$score <- apply(countMatrix,1,mean)
write.csv(countMatrix,"CountMatrix__cellviability0.8.csv")
BP_score <- as.data.frame(as.double(countMatrix$score))
BP1 <- as.data.frame(cbind(rownames(countMatrix),BP_score))
BP1 <- filter(BP1,BP1$`as.double(countMatrix$score)`>0)
BP <- as.double(BP1[,2])
names(BP) <- as.character(BP1[,1])
BP <- as.data.frame(sort(BP, decreasing = FALSE))
BP <- as.data.frame(cbind(rownames(BP),BP$`sort(BP, decreasing = FALSE)`))
colnames(BP) <- c("gene","score")
write.table(BP,"BP.txt",sep="\t",row.names = FALSE)

###################################################### MRMR cv50
library(mRMRe)
mrmr_feature<-train[,gene_index]
mrmr_feature$y<-train$cv0.5
target_indices = which(names(mrmr_feature)=='y')
for (m in which(sapply(mrmr_feature, class)!="numeric")){
  mrmr_feature[,m]=as.numeric(unlist(mrmr_feature[,m]))
} 
Data <- mRMR.data(data = data.frame(mrmr_feature))
#the counts of invariant genes:5135
mrmr=mRMR.ensemble(data = Data, target_indices = target_indices, 
                   feature_count = 5135, solution_count = 1)

score=as.data.frame(mrmr@scores[[as.character(mrmr@target_indices)]])

index2=mrmr@filters[[as.character(mrmr@target_indices)]]
train_feature <- mrmr_feature[,index2]
feature <- as.data.frame(colnames(train_feature))
feature_with_score <- cbind(feature,score)
write.csv(feature_with_score,"MRMR_cv0.5.csv")
###################################################### MRMR cv80
library(mRMRe)
mrmr_feature<-train[,gene_index]
mrmr_feature$y<-train$cv0.8
target_indices = which(names(mrmr_feature)=='y')
for (m in which(sapply(mrmr_feature, class)!="numeric")){
  mrmr_feature[,m]=as.numeric(unlist(mrmr_feature[,m]))
} 
Data <- mRMR.data(data = data.frame(mrmr_feature))
mrmr=mRMR.ensemble(data = Data, target_indices = target_indices, 
                   feature_count = 5135, solution_count = 1)
score=as.data.frame(mrmr@scores[[as.character(mrmr@target_indices)]])

index2=mrmr@filters[[as.character(mrmr@target_indices)]]
train_feature <- mrmr_feature[,index2]
feature <- as.data.frame(colnames(train_feature))
feature_with_score <- cbind(feature,score)

write.csv(feature_with_score,"MRMR_cv0.8.csv")


#####################################RRA###########################

setwd("C:\\Users\\wuyou\\Desktop\\20230607 machine leaning model of cell viability\\4_Feature selection\\RRA")
files=c("BP_CV50.txt", "MRMR_CV50.txt","WGCNA.txt")
upList=list()

allFCList=list()


for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.table(inputFile,header=T,sep = '\t',quote = '') # 注意文件读取
  header=unlist(strsplit(inputFile,"_"))
  upList[[header[1]]]=rev(as.vector(as.character(rt[,1])))
  fcCol=rt[,1:2]
  colnames(fcCol)=c("Gene",header[[1]])
  allFCList[[header[1]]]=fcCol
}
mergeLe=function(x,y){
  merge(x,y,by="Gene",all=T)}
newTab=Reduce(mergeLe,allFCList)
rownames(newTab)=newTab[,1]
newTab=newTab[,2:ncol(newTab)]
newTab[is.na(newTab)]=0
library(RobustRankAggreg)
upMatrix = rankMatrix(upList)
upAR1 = aggregateRanks(rmat=upMatrix,exact=TRUE	)
write.table(upAR1,file="CV50_up_exact_p_value.xls",sep="\t",quote=F,row.names=F)

############################################################
files=c("BP_CV80.txt", "MRMR_CV80.txt","WGCNA.txt")
upList=list()
allFCList=list()
for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.table(inputFile,header=T,sep = '\t',quote = '') # 注意文件读取
  header=unlist(strsplit(inputFile,"_"))
  
  upList[[header[1]]]=rev(as.vector(as.character(rt[,1])))
  fcCol=rt[,1:2]
  colnames(fcCol)=c("Gene",header[[1]])
  allFCList[[header[1]]]=fcCol
}
mergeLe=function(x,y){
  merge(x,y,by="Gene",all=T)}
newTab=Reduce(mergeLe,allFCList)
rownames(newTab)=newTab[,1]
newTab=newTab[,2:ncol(newTab)]
newTab[is.na(newTab)]=0
upMatrix = rankMatrix(upList)
upAR1 = aggregateRanks(rmat=upMatrix,exact=TRUE	)
write.table(upAR1,file="CV80_up_exact_p_value.xls",sep="\t",quote=F,row.names=F)


