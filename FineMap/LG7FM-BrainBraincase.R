library(qtl); library(dplyr); library(plotrix); library(ggplot2)

####Fine Mapping####
setwd("/Users/Home/Desktop/GitHub/FineMap/")

##Read in fine map
LG7 <- read.cross(format="csv",file="LG7_NcBrainCoV.csv",
                           na.strings="NA",genotypes=c("AA","AB","BB"),
                           alleles=c("A","B"),convertXdata=TRUE)

#Check sample sizes
apply(LG7$geno$`7`$data, 2, table)
colnames(LG7$pheno)

#Extract marker names w/ sample size greater than 10
ScafNames<-rownames(geno.table(LG7))[as.logical(apply(X = geno.table(LG7,7)[,3:5]>=10,MARGIN = 1,FUN = mean)==1)]
EffResult<-matrix(data = NA, nrow = length(ScafNames), ncol = 6)
  
#Calculate averge phenotypic effect for BrBC_XBlock
for (i in 1:length(ScafNames)){
  effect<-effectplot(LG7, pheno.col=2, mname1=ScafNames[i], var.flag=c("group"), draw = F)
  EffResult[i,]<-c(effect$Means, effect$SEs)
}

CombinedDataEffectE<-cbind.data.frame(ScafNames,as.numeric(sub("^[^_]*_", "", ScafNames)), EffResult,EffResult[,3]-EffResult[,1])
colnames(CombinedDataEffectE)<-c("Scaff", "Pos", "AA_mean", "AB_mean", "BB_mean", "AA_se", "AB_se", "BB_se", "AA.BB")

CombinedDataEffectEo<-CombinedDataEffectE[order(CombinedDataEffectE$Pos),]

#write.csv(CombinedDataEffectEo,"LG7_FineMapDataT_BrBCXBlock_Stringent_noaug.csv", quote = F, row.names = F)

#Calculate averge phenotypic effect for BCBr_YBlock
for (i in 1:length(ScafNames)){
  effect<-effectplot(LG7, pheno.col=3, mname1=ScafNames[i], var.flag=c("group"), draw = F)
  EffResult[i,]<-c(effect$Means, effect$SEs)
}

CombinedDataEffectE<-cbind.data.frame(ScafNames,as.numeric(sub("^[^_]*_", "", ScafNames)), EffResult,EffResult[,3]-EffResult[,1])
colnames(CombinedDataEffectE)<-c("Scaff", "Pos", "AA_mean", "AB_mean", "BB_mean", "AA_se", "AB_se", "BB_se", "AA.BB")

CombinedDataEffectEo<-CombinedDataEffectE[order(CombinedDataEffectE$Pos),]

#write.csv(CombinedDataEffectEo,"LG7_FineMapDataT_BCBrYBlock_Stringent_noaug.csv", quote = F, row.names = F)

##Combine the X-Y block datasets into a single dataframe defined my a 'Trait' Column
#See LG7_FineMapDataT_BrBcCombine.csv as an example

##Quick Plot Check##
#Read in files you just created, or skip this section and use the example files provided in the 'ggplot' section below
setwd("/Users/Home/Desktop/GitHub/FineMap/")
CombinedDataEffectEo<-read.csv("LG7_FineMapDataT_BCBrYBlock_Stringent_noaug.csv")

MapPP<-CombinedDataEffectEo
MapP<-MapPP[complete.cases(MapPP$AA.BB) & complete.cases(MapPP$AA_se) & complete.cases(MapPP$BB_se),]
plot(MapP$Pos, (MapP$AA.BB*-1), ylim=c(-0.01,0.016), type='l', lwd=2, lty=1, xlab = 'LG7', ylab = 'Avg. Pheno. Effect')

CombinedDataEffectEo<-read.csv("LG7_FineMapDataT_BrBCXBlock_Stringent_noaug.csv")
MapPP<-CombinedDataEffectEo
MapP<-MapPP[complete.cases(MapPP$AA.BB) & complete.cases(MapPP$AA_se) & complete.cases(MapPP$BB_se),]
lines(MapP$Pos, (MapP$AA.BB*-1), col="red", lwd=2)
abline(h=0, lty=3, lwd=1)

#BayesInt, Br-Bc XBlock
abline(v=50846240, col="purple") #scaffold_21_5242487 -> NC036786.1_50846240
abline(v=43819923, col="purple") #scaffold_0_5071706 NC036786.1_43819923

#BayesInt, Br-Bc XBlock
abline(v=47752984, col="blue") #scaffold_21_2195347 -> NC036786.1_47752984
abline(v=43819923, col="blue") #scaffold_0_5071706 NC036786.1_43819923



####Full map ggplot####
setwd("/Users/Home/Desktop/GitHub/FineMap/")
FullMap<-read.csv("LG7_FineMapDataT_BrBcCombine.csv")

####Line plot
MyLabels <- c(BrBC = "Brain [Braincase]", BCBr = "Braincase [Brain]")

FineMap_Plot <- ggplot(FullMap, aes(x=Pos, y=AA.BB, fill=Trait, col=Trait)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values=c("#0571b0","#ca0020")) +
  scale_fill_manual(values=c("#0571b0","#ca0020")) +
  facet_wrap(~ Trait, scales = "free", labeller=labeller(Trait = MyLabels), nrow = 1)

FineMap_PlotSE <- FineMap_Plot + theme_bw() + annotate("rect", xmin = 43819923, xmax = 50846240, ymin = -Inf, ymax = Inf, fill = "#404040", alpha = 0.2) +
  geom_ribbon(data = FullMap, alpha=0.25, aes(ymin=AA.BB-AA_se, ymax=AA.BB+AA_se), col=NA) +
  geom_ribbon(data = FullMap, alpha=0.25, aes(ymin=AA.BB-BB_se, ymax=AA.BB+BB_se), col=NA) +
  geom_line(size=0.9) +
  labs(x = "Marker position on LG7 (Mbs)", y = "Average Phenotypic Effect") +
  theme(legend.position = "none", aspect.ratio = 1, legend.key=element_blank(), axis.ticks = element_blank(), panel.border = element_rect(colour = "dark gray", fill=NA, size=1))

FineMap_PlotSE



####Reduced map plot####
setwd("/Users/Home/Desktop/GitHub/FineMap/")
MapPP<-read.csv("LG7_FineMapDataT_BCBrCombine.csv")
FullMap<-MapPP[complete.cases(MapPP$AA.BB) & complete.cases(MapPP$AA_se) & complete.cases(MapPP$BB_se),]

####Line plot
MyLabels <- c(BrBC = "Brain [Braincase]", BCBr = "Braincase [Brain]")

FineMap_Plot <- ggplot(FullMap, aes(x=UMD2a_start, y=AA.BB, fill=Trait, col=Trait)) +
  scale_color_manual(values=c("#0571b0","#ca0020")) +
  scale_fill_manual(values=c("#0571b0","#ca0020"))

FineMap_PlotRev <- FineMap_Plot + theme_bw() + annotate("rect", xmin = 43819923, xmax = 50846240, ymin = -Inf, ymax = Inf, fill = "#404040", alpha = 0.2) +
  geom_line(size=0.9) +
  labs(x = "Marker position on LG7 (Mbs)", y = "Average Phenotypic Effect") +
  scale_y_reverse() +
  theme(legend.position = "none", aspect.ratio = 1, legend.key=element_blank(), axis.ticks = element_blank(), panel.border = element_rect(colour = "dark gray", fill=NA, size=1))

FineMap_PlotRev


####Reduced map Fst Plot####
setwd("/Users/Home/Desktop/GitHub/FineMap/")
FullMap<-read.csv("LG7_FineMapDataT_BCBrCombine.csv")
MyLabels <- c(BrBC = "Brain (Braincase)", BCBr = "Braincase (Brain)")

#Insert a user defined alpha value
#To adjust the opacity range play with the
#'scale_alpha_continuous' and the 0 value below
FullMap$FstAlpha<-FullMap$Fst
FullMap$FstAlpha[which(FullMap$FstAlpha<0.6)]<-0

Zrange <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}
FullMap$FstAlpha[which(FullMap$FstAlpha>0.6)]<-Zrange(FullMap$FstAlpha[which(FullMap$FstAlpha>0.6)])

FineMap_Fst <- filter(FullMap, Fst > 0) %>%
  ggplot(aes(x=UMD2a_start, y=Fst))

FineMap_FstPoint <- FineMap_Fst + theme_bw() + annotate("rect", xmin = 43819923, xmax = 50846240, ymin = -Inf, ymax = Inf, fill = "#404040", alpha = 0.2) +
  geom_point(size=2, aes(alpha=FstAlpha)) +
  scale_alpha_continuous(range = c(0.2,1)) +
  labs(x = "Marker position on LG7 (Mbs)", y = "Fst") +
  theme(legend.position = "none", aspect.ratio = 1, legend.key=element_blank(), axis.ticks = element_blank(), panel.border = element_rect(colour = "dark gray", fill=NA, size=1))

FineMap_FstPointGene <- FineMap_FstPoint + 
  geom_hline(yintercept = 0.6, col="dark gray", linetype="dashed") #Standard Fst threshold

FineMap_FstPointGene
