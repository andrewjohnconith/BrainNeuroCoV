library(geomorph); library(tools); library(Morpho)

####NC Morphometrics####
setwd("/Users/Home/Dropbox/Amherst/Post-Doc/Muscles/Brain/Final/Writing/GitHub/Shape")

#Read in data and TPS file
AllData<-readland.tps("NC_20190809.tps", specID = "imageID", readcurves = F)
ReflectInds<-c("LFxTRC5563", "LFxTRC5526", "LFxTRC5547", "LFxTRC5552", "LFxTRC5509", "LFxTRC5089", "LFxTRC5613", "LFxTRC5606", "LFxTRC5086", "LFxTRC5569", "LFxTRC5614", "LFxTRC5093", "LFxTRC5109", "LFxTRC5108", "LFxTRC5601", "LFxTRC5589", "LFxTRC5098", "LFxTRC5527", "LFxTRC5095", "LFxTRC5118", "LFxTRC5582", "LFxTRC5099", "LFxTRC5105", "LFxTRC5591", "LFxTRC5115", "LFxTRC5620", "LFxTRC5087", "LFxTRC5603", "LFxTRC5578", "LFxTRC5110", "LFxTRC5009", "LFxTRC5088", "LFxTRC5106", "LFxTRC5595", "LFxTRC5583", "LFxTRC5117", "LFxTRC5112", "LFxTRC5510", "LFxTRC5084", "LFxTRC5113")
AllData[,1,ReflectInds]<-AllData[,1,ReflectInds]*-1

#Procrustes
Y.gpa<-gpagen(AllData)
gdf <- geomorph.data.frame(Y.gpa, species = dimnames(Y.gpa$coords)[[3]])
NCAllometry <- procD.lm(coords~Csize, data=gdf, iter=999)

summary(NCAllometry)

shape.resid <- arrayspecs(NCAllometry$residuals, p=dim(Y.gpa$coords)[1], k=dim(Y.gpa$coords)[2]) # size-adjusted residuals

#allometry-free shapes
NC.adj.shape <- shape.resid + array(Y.gpa$consensus, dim(shape.resid))

###Plot PCA
par(pty='s')
NeuroPC<-gm.prcomp(A = NC.adj.shape)
pc.plot<-plot(NeuroPC, axis1 = 1, axis2 = 2, pch=19)
#text(NeuroPC$x[,1:2], labels = rownames(NeuroPC$x))

#write.csv(NeuroPC$x, "NeuroPC_20220224.csv")


#Check NC LM number
lmk<- matrix(c(1,2,3,4,4,5,5,6,9,10,6,7,7,15,15,1,15,16,4,9,5,11,14,20,3,2,1,8,8,9,
               16,18,18,20,20,21,21,10,21,11,11,12,12,14,14,13,17,19,19,13,19,6,16,17,9,21),ncol=2, byrow = TRUE)

define.links(AllData[,,"LFxTRC5602"], ptsize = 1,  links = lmk)


##Modularity
  ##Rostrum bias
ModVec1<-c(rep("A",9),rep("B",5),"A",rep("B",6))
NC_ModTest_RB<-modularity.test(A = NC.adj.shape, partition.gp = ModVec1, iter = 9999)
summary(NC_ModTest_RB)
# Result implies modularity present

  ##Braincase bias
ModVec2<-c("A","A","A","A","B","B","B","A","B","B","B","B","B","B","A","B","B","B","B","B","B")
NC_ModTest_BcB<-modularity.test(A = NC.adj.shape, partition.gp = ModVec2, iter = 9999)
summary(NC_ModTest_BcB)

  ##Three modules
ModVec3<-c("A","A","A","A","C","C","C","A","C","B","B","B","B","B","A","B","B","B","B","B","B")
NC_ModTest_TM<-modularity.test(A = NC.adj.shape, partition.gp = ModVec3, iter = 9999)
summary(NC_ModTest_TM)


compare.CR(NC_ModTest_RB,NC_ModTest_TM,CR.null = F)

#  Effect sizes
####Examining modular support b/w hypotheses####
####compare.CR
#NOTE: more negative effects represent stronger modular signal! 
#  Effect sizes
#
#NC_ModTest_RB NC_ModTest_TM 
#-3.051926     -1.857668 

#Effect sizes for pairwise differences in CR effect size
#NC_ModTest_RB NC_ModTest_TM
#NC_ModTest_RB     0.0000000     0.5781602
#NC_ModTest_TM     0.5781602     0.0000000

#P-values
#NC_ModTest_RB NC_ModTest_TM
#NC_ModTest_RB      1.000000      0.563156
#NC_ModTest_TM      0.563156      1.000000

#Partition landmarks
# Two module partition
Rostrum.shape<-NC.adj.shape[c(1:9,15),,]
BrainCase.shape<-NC.adj.shape[c(5:7,9:14,16:21),,]

#Three module partition
Rostrum.shape.Only<-NC.adj.shape[c(1:4,8,15),,] #Rostrum
BrainCase.shape.Only<-NC.adj.shape[c(10:14,16:21),,] #Braincase
Central.shape.Only<-NC.adj.shape[c(5:7,9),,] #Center

#Rostrum
NeuroPCrm<-gm.prcomp(A = Rostrum.shape)
pc.plot<-plot(NeuroPCrm, axis1 = 1, axis2 = 2, pch=19)
#text(NeuroPCrm$x[,1:2], labels = rownames(NeuroPCrm$x))
#write.csv(NeuroPCrm$x, "NeuroPCrm_20220224.csv")

#Brain case
NeuroPCbc<-gm.prcomp(A = BrainCase.shape)
pc.plot<-plot(NeuroPCbc, axis1 = 1, axis2 = 2, pch=19)
#text(NeuroPCbc$x[,1:2], labels = rownames(NeuroPCbc$x))
#write.csv(NeuroPCbc$x, "NeuroPCbc_20220303.csv")


####Brain Morphometrics####
setwd("/Users/Home/Dropbox/Amherst/Post-Doc/Muscles/Brain/Final/Writing/GitHub/Shape")

#Read in data and TPS file
AllInds<-readland.tps(file = "Brain_20220224.tps", specID = "ID")

#Procrustes
SurfSlide<-seq(from = 58, to = 354, by = 1)
SlidersFile<-read.csv(file = "BrainCurveslide.csv")
Y.gpa<-gpagen(A = AllInds, surfaces = SurfSlide, curves = SlidersFile, ProcD = T)

#Allometry
gdf <- geomorph.data.frame(Y.gpa, species = dimnames(Y.gpa$coords)[[3]])
HybAllometry <- procD.lm(coords~Csize, data=gdf, iter=999)
summary(HybAllometry)
shape.resid <- arrayspecs(HybAllometry$residuals, p=dim(Y.gpa$coords)[1], k=dim(Y.gpa$coords)[2]) # size-adjusted residuals

#allometry-free shapes
adj.shape.brain <- shape.resid + array(Y.gpa$consensus, dim(shape.resid))
#save(adj.shape.brain, file="BrainAlloCorrect_20220125.rda")


###Plot PCA
par(pty='s')
BrainPC<-gm.prcomp(A = adj.shape.brain)
pc.plot<-plot(BrainPC, axis1 = 1, axis2 = 2, pch=19)
#text(BrainPC$x[,1:2], labels = rownames(BrainPC$x))
#write.csv(BrainPC$x[,1:2], "BrainPC_20220224.csv")



####2B-PLS####
##NC-Brain
CommonInds<-c(gsub( "_.*$", "", dimnames(adj.shape.brain)[[3]]), file_path_sans_ext(dimnames(NC.adj.shape)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(adj.shape.brain)[[3]]), file_path_sans_ext(dimnames(NC.adj.shape)[[3]])))]

#Brain setup
Brainnew<-adj.shape.brain[,,paste0(CommonInds,"_Brain")]
dimnames(Brainnew)[[3]]<-CommonInds

#NC setup
NCnew<-NC.adj.shape[,,CommonInds]
dimnames(NCnew)[[3]]<-CommonInds

BrainNCCorrs<-two.b.pls(A1 = Brainnew, A2 = NCnew, iter = 9999)
#plot(BrainNCCorrs, label = CommonInds)
plot(BrainNCCorrs)
summary(BrainNCCorrs)
#picknplot.shape(plot(BrainNCCorrs, label = CommonInds))

#For plotting and saving
#BrainNCCorrsPLOT<-cbind.data.frame(BrainNCCorrs$XScores[,1],BrainNCCorrs$YScores[,1])
#colnames(BrainNCCorrsPLOT)<-c("XBlock","YBlock")
#write.csv(BrainNCCorrsPLOT,"BrainNCCorrsPLOT_20220224.csv")

##Brain-Braincase
CommonInds<-c(gsub( "_.*$", "", dimnames(adj.shape.brain)[[3]]), file_path_sans_ext(dimnames(BrainCase.shape)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(adj.shape.brain)[[3]]), file_path_sans_ext(dimnames(BrainCase.shape)[[3]])))]

#Brain setup
Brainnew<-adj.shape.brain[,,paste0(CommonInds,"_Brain")]
dimnames(Brainnew)[[3]]<-CommonInds

#NC setup
NCnew<-BrainCase.shape[,,CommonInds]
dimnames(NCnew)[[3]]<-CommonInds

BrainBrainCaseCorrs<-two.b.pls(A1 = Brainnew, A2 = NCnew, iter = 9999)
plot(BrainBrainCaseCorrs)
#plot(BrainBrainCaseCorrs, label = CommonInds)
summary(BrainBrainCaseCorrs)
#picknplot.shape(plot(BrainBrainCaseCorrs, label = CommonInds))

#BrainBrainCaseCorrsPLOT<-cbind.data.frame(BrainBrainCaseCorrs$XScores[,1],BrainBrainCaseCorrs$YScores[,1])
#colnames(BrainBrainCaseCorrsPLOT)<-c("XBlock","YBlock")
#write.csv(BrainBrainCaseCorrsPLOT,"BrainBraincaseCorrsPLOT_20220224.csv")


##Rostrum-Brain
CommonInds<-c(gsub( "_.*$", "", dimnames(adj.shape.brain)[[3]]), file_path_sans_ext(dimnames(Rostrum.shape)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(adj.shape.brain)[[3]]), file_path_sans_ext(dimnames(Rostrum.shape)[[3]])))]

#Brain setup
Brainnew<-adj.shape.brain[,,paste0(CommonInds,"_Brain")]
dimnames(Brainnew)[[3]]<-CommonInds

#NC setup
NCnew<-Rostrum.shape[,,CommonInds]
dimnames(NCnew)[[3]]<-CommonInds

BrainRostrumCorrs<-two.b.pls(A1 = Brainnew, A2 = NCnew, iter = 9999)
plot(BrainRostrumCorrs)
#plot(BrainRostrumCorrs, label = CommonInds)
summary(BrainRostrumCorrs)
#picknplot.shape(plot(BrainRostrumCorrs, label = CommonInds))

#BrainRostrumCorrsPLOT<-cbind.data.frame(BrainRostrumCorrs$XScores[,1],BrainRostrumCorrs$YScores[,1])
#colnames(BrainRostrumCorrsPLOT)<-c("XBlock","YBlock")
#write.csv(BrainRostrumCorrsPLOT,"BrainRostrumCorrsPLOT_20220224.csv")



####Compare PLS####
##Assess strength of integration with the brain; rostrum vs braincase
compare.pls(BrainBrainCaseCorrs,BrainRostrumCorrs)

#Effect sizes
#BrainBrainCaseCorrs  BrainRostrumCorrs 
#2.093200             1.896705 

#Effect sizes for pairwise differences in PLS effect size
#                       BrainBrainCaseCorrs   BrainRostrumCorrs
#BrainBrainCaseCorrs             0.0000000           0.1451925
#BrainRostrumCorrs               0.1451925           0.0000000

#P-values
#                     BrainBrainCaseCorrs   BrainRostrumCorrs
#BrainBrainCaseCorrs           1.0000000         0.8845589
#BrainRostrumCorrs             0.8845589         1.0000000


####2B-PLS 2####
Rostrum.shape.Only<-NC.adj.shape[c(1:4,8,15),,] #Rostrum
BrainCase.shape.Only<-NC.adj.shape[c(10:14,16:21),,] #Braincase
Central.shape.Only<-NC.adj.shape[c(5:7,9),,] #Center

##Rostrum-Brain
CommonInds<-c(gsub( "_.*$", "", dimnames(adj.shape.brain)[[3]]), file_path_sans_ext(dimnames(Rostrum.shape.Only)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(adj.shape.brain)[[3]]), file_path_sans_ext(dimnames(Rostrum.shape.Only)[[3]])))]

#Brain setup
Brainnew<-adj.shape.brain[,,paste0(CommonInds,"_Brain")]
dimnames(Brainnew)[[3]]<-CommonInds

#NC setup
NCnew<-Rostrum.shape.Only[,,CommonInds]
dimnames(NCnew)[[3]]<-CommonInds

BrainRostrumCorrsONLY<-two.b.pls(A1 = Brainnew, A2 = NCnew, iter = 9999)
#plot(BrainRostrumCorrsONLY, label = CommonInds)
plot(BrainRostrumCorrsONLY)
summary(BrainRostrumCorrsONLY)
#picknplot.shape(plot(BrainRostrumCorrsONLY, label = CommonInds))

#For plotting and saving
#BrainRostrumCorrsONLYPLOT<-cbind.data.frame(BrainRostrumCorrsONLYPLOT$XScores[,1],BrainRostrumCorrsONLYPLOT$YScores[,1])
#colnames(BrainRostrumCorrsONLYPLOT)<-c("XBlock","YBlock")
#write.csv(BrainRostrumCorrsONLYPLOT,"BrainNCCorrsPLOT_20220224.csv")

##Brain-Braincase
CommonInds<-c(gsub( "_.*$", "", dimnames(adj.shape.brain)[[3]]), file_path_sans_ext(dimnames(BrainCase.shape.Only)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(adj.shape.brain)[[3]]), file_path_sans_ext(dimnames(BrainCase.shape.Only)[[3]])))]

#Brain setup
Brainnew<-adj.shape.brain[,,paste0(CommonInds,"_Brain")]
dimnames(Brainnew)[[3]]<-CommonInds

#NC setup
NCnew<-BrainCase.shape.Only[,,CommonInds]
dimnames(NCnew)[[3]]<-CommonInds

BrainBrainCaseCorrsONLY<-two.b.pls(A1 = Brainnew, A2 = NCnew, iter = 9999)
plot(BrainBrainCaseCorrsONLY)
#plot(BrainBrainCaseCorrsONLY, label = CommonInds)
summary(BrainBrainCaseCorrsONLY)
#picknplot.shape(plot(BrainBrainCaseCorrsONLY, label = CommonInds))

#BrainBrainCaseCorrsONLYPLOT<-cbind.data.frame(BrainBrainCaseCorrsONLYPLOT$XScores[,1],BrainBrainCaseCorrsONLYPLOT$YScores[,1])
#colnames(BrainBrainCaseCorrsONLYPLOT)<-c("XBlock","YBlock")
#write.csv(BrainBrainCaseCorrsONLYPLOT,"BrainBrainCaseCorrsONLYPLOT_20220224.csv")


##Central-Brain
CommonInds<-c(gsub( "_.*$", "", dimnames(adj.shape.brain)[[3]]), file_path_sans_ext(dimnames(Central.shape.Only)[[3]]))[duplicated(c(gsub( "_.*$", "", dimnames(adj.shape.brain)[[3]]), file_path_sans_ext(dimnames(Central.shape.Only)[[3]])))]

#Brain setup
Brainnew<-adj.shape.brain[,,paste0(CommonInds,"_Brain")]
dimnames(Brainnew)[[3]]<-CommonInds

#NC setup
NCnew<-Central.shape.Only[,,CommonInds]
dimnames(NCnew)[[3]]<-CommonInds

BrainCentralCorrsONLY<-two.b.pls(A1 = Brainnew, A2 = NCnew, iter = 9999)
plot(BrainCentralCorrsONLY)
#plot(BrainCentralCorrsONLY, label = CommonInds)
summary(BrainCentralCorrsONLY)
#picknplot.shape(plot(BrainCentralCorrsONLY, label = CommonInds))

#BrainCentralCorrsONLYPLOT<-cbind.data.frame(BrainCentralCorrsONLY$XScores[,1],BrainCentralCorrsONLY$YScores[,1])
#colnames(BrainCentralCorrsONLYPLOT)<-c("XBlock","YBlock")
#write.csv(BrainCentralCorrsONLYPLOT,"BrainCentralCorrsONLYPLOT_20220224.csv")


####Compare PLS 2####
##Assess strength of integration with the brain; rostrum vs braincase
compare.pls(BrainRostrumCorrsONLY,BrainBrainCaseCorrsONLY,BrainCentralCorrsONLY)

#Effect sizes
#BrainRostrumCorrsONLY BrainBrainCaseCorrsONLY   BrainCentralCorrsONLY 
#             2.195459                3.267109                1.399706 

#Effect sizes for pairwise differences in PLS effect size
#                         BrainRostrumCorrsONLY  BrainBrainCaseCorrsONLY BrainCentralCorrsONLY
#BrainRostrumCorrsONLY                0.000000                0.123997              1.710567
#BrainBrainCaseCorrsONLY              0.123997                0.000000              1.885305
#BrainCentralCorrsONLY                1.710567                1.885305              0.000000

#P-values
#                         BrainRostrumCorrsONLY BrainBrainCaseCorrsONLY BrainCentralCorrsONLY
#BrainRostrumCorrsONLY              1.00000000              0.90131766            0.08716103
#BrainBrainCaseCorrsONLY            0.90131766              1.00000000            0.05938864
#BrainCentralCorrsONLY              0.08716103              0.05938864            1.00000000
