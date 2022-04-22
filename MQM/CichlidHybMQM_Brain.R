library(qtl); library(dplyr); library(qtl2); library(qtl2pleio); library(ggplot2)

####NC Morphometrics####
setwd("/Users/Home/Desktop/GitHub/MQM/")

##Read in augmented data
augdata1 <- read.cross(format="csv",file="AugCrossF5_Covariation-20220407.csv",
                    na.strings="-",genotypes=c("AA","AB","BB"),
                    alleles=c("A","B"),convertXdata=TRUE)

####Br-Nc XBlock####
scan1<-mqmscan(augdata1,pheno.col=2,plot=T,model="dominance")
summary(scan1)

cofactorsindex<-NULL
cofactorsindex<-c(93,173,373)

#find.markerindex(augdata1,find.marker(augdata1,chr=3,pos=60))
#find.markerindex(augdata1,find.marker(augdata1,chr=5,pos=70))
#find.markerindex(augdata1,find.marker(augdata1,chr=11,pos=15))

homemadecofactors<-mqmsetcofactors(augdata1,cofactors=c(cofactorsindex))

homescan1<-mqmscan(cross=augdata1,cofactors=homemadecofactors,pheno.col=2,cofactor.significance=0.02,verbose=T,plot=T,model="dominance") #using a dominance model tests for both additive effects and dominance, i.e., AA+AB vs BB and BB+AB vs AA, default is additive only
summary(homescan1)

bayesint(homescan1,3,qtl.index=93,prob=0.95,lodcolumn=1, expandtomarkers=T)
bayesint(homescan1,5,qtl.index=173,prob=0.95,lodcolumn=1, expandtomarkers=T)

find.marker(augdata1,chr=3,pos=60)
A<-effectplot(F5All, pheno.col=2, mname1="scaffold_327_298291",var.flag=c("group"))
A
find.marker(augdata1,chr=5,pos=70)
B<-effectplot(F5All, pheno.col=2, mname1="scaffold_130_1262569",var.flag=c("group"))
B

NC1Results <- mqmpermutation(augdata1, pheno.col=2, scanfunction=mqmscan, cofactors=homemadecofactors, n.perm=100, multicore=T, plot=F)
sig.results <- mqmprocesspermutation(NC1Results)
summary(sig.results)
#mqmplot.permutations(NC1Results)


####Nc-Br YBlock####
scan1<-mqmscan(augdata1,pheno.col=3,plot=T,model="dominance")
summary(scan1)

cofactorsindex<-NULL
cofactorsindex<-c(237)

#find.markerindex(augdata1,find.marker(augdata1,chr=7,pos=20))

homemadecofactors<-mqmsetcofactors(augdata1,cofactors=c(cofactorsindex))

homescan1<-mqmscan(cross=augdata1,cofactors=homemadecofactors,pheno.col=3,cofactor.significance=0.02,verbose=T,plot=T,model="dominance") #using a dominance model tests for both additive effects and dominance, i.e., AA+AB vs BB and BB+AB vs AA, default is additive only
summary(homescan1)

bayesint(homescan1,7,qtl.index=237,prob=0.95,lodcolumn=1, expandtomarkers=T) #Use this one.

find.marker(augdata1,chr=7,pos=20)
A<-effectplot(F5All, pheno.col=3, mname1="scaffold_21_2195347",var.flag=c("group"))
A

NC2Results <- mqmpermutation(augdata1, pheno.col=3, scanfunction=mqmscan, cofactors=homemadecofactors, n.perm=100, multicore=T, plot=F)
sig.results <- mqmprocesspermutation(NC2Results)
summary(sig.results)
#mqmplot.permutations(NC2Results)


####Br-Ro XBlock####
scan1<-mqmscan(augdata1,pheno.col=4,plot=T,model="dominance")
summary(scan1)

cofactorsindex<-NULL
cofactorsindex<-c(93,178)

#find.markerindex(augdata1,find.marker(augdata1,chr=3,pos=60))
#find.markerindex(augdata1,find.marker(augdata1,chr=5,pos=75))

homemadecofactors<-mqmsetcofactors(augdata1,cofactors=c(cofactorsindex)) #Step 3: designate cofactors for multivariate scan.

homescan1<-mqmscan(cross=augdata1,cofactors=homemadecofactors,pheno.col=4,cofactor.significance=0.02,verbose=T,plot=T,model="dominance") #using a dominance model tests for both additive effects and dominance, i.e., AA+AB vs BB and BB+AB vs AA, default is additive only
summary(homescan1)

bayesint(homescan1,3,qtl.index=93,prob=0.95,lodcolumn=1, expandtomarkers=T)
bayesint(homescan1,5,qtl.index=178,prob=0.95,lodcolumn=1, expandtomarkers=T)

find.marker(augdata1,chr=3,pos=60)
A<-effectplot(F5All, pheno.col=4, mname1="scaffold_327_298291",var.flag=c("group"))
A
find.marker(augdata1,chr=5,pos=75)
B<-effectplot(F5All, pheno.col=4, mname1="scaffold_208_213442",var.flag=c("group"))
B

Rostrum1Results <- mqmpermutation(augdata1, pheno.col=4, scanfunction=mqmscan, cofactors=homemadecofactors, n.perm=100, multicore=T, plot=F)
sig.results <- mqmprocesspermutation(Rostrum1Results)
summary(sig.results)
#mqmplot.permutations(Rostrum1Results)


####Ro-Br YBlock####
scan1<-mqmscan(augdata1,pheno.col=5,plot=T,model="dominance")
summary(scan1)

cofactorsindex<-NULL
cofactorsindex<-c(237,421)

#find.markerindex(augdata1,find.marker(augdata1,chr=7,pos=20))
#find.markerindex(augdata1,find.marker(augdata1,chr=12,pos=19))

homemadecofactors<-mqmsetcofactors(augdata1,cofactors=c(cofactorsindex)) #Step 3: designate cofactors for multivariate scan.

homescan1<-mqmscan(cross=augdata1,cofactors=homemadecofactors,pheno.col=5,cofactor.significance=0.02,verbose=T,plot=T,model="dominance") #using a dominance model tests for both additive effects and dominance, i.e., AA+AB vs BB and BB+AB vs AA, default is additive only
summary(homescan1)

bayesint(homescan1,7,qtl.index=237,prob=0.95,lodcolumn=1, expandtomarkers=T)
bayesint(homescan1,12,qtl.index=421,prob=0.95,lodcolumn=1, expandtomarkers=T)

find.marker(augdata1,chr=7,pos=20)
A<-effectplot(F5All, pheno.col=5, mname1="scaffold_21_2195347",var.flag=c("group"))
A
find.marker(augdata1,chr=12,pos=19)
B<-effectplot(F5All, pheno.col=5, mname1="scaffold_9_6193445",var.flag=c("group"))
B

Rostrum2Results <- mqmpermutation(augdata1, pheno.col=5, scanfunction=mqmscan, cofactors=homemadecofactors, n.perm=100, multicore=T, plot=F)
sig.results <- mqmprocesspermutation(Rostrum2Results)
summary(sig.results)
#mqmplot.permutations(Rostrum2Results)


####Br-Bc XBlock####
scan1<-mqmscan(augdata1,pheno.col=6,plot=T,model="dominance")
summary(scan1)

cofactorsindex<-NULL
cofactorsindex<-c(237,560)

#find.markerindex(augdata1,find.marker(augdata1,chr=7,pos=20))
#find.markerindex(augdata1,find.marker(augdata1,chr=16,pos=50))

homemadecofactors<-mqmsetcofactors(augdata1,cofactors=c(cofactorsindex)) #Step 3: designate cofactors for multivariate scan.

homescan1<-mqmscan(cross=augdata1,cofactors=homemadecofactors,pheno.col=6,cofactor.significance=0.02,verbose=T,plot=T,model="dominance") #using a dominance model tests for both additive effects and dominance, i.e., AA+AB vs BB and BB+AB vs AA, default is additive only
summary(homescan1)

bayesint(homescan1,7,qtl.index=237,prob=0.95,lodcolumn=1, expandtomarkers=T)
bayesint(homescan1,16,qtl.index=560,prob=0.95,lodcolumn=1, expandtomarkers=T)

find.marker(augdata1,chr=7,pos=20)
A<-effectplot(F5All, pheno.col=6, mname1="scaffold_21_2195347",var.flag=c("group"))
A
find.marker(augdata1,chr=16,pos=50)
B<-effectplot(F5All, pheno.col=6, mname1="scaffold_86_1495785",var.flag=c("group"))
B

BraCase1Results <- mqmpermutation(augdata1, pheno.col=6, scanfunction=mqmscan, cofactors=homemadecofactors, n.perm=100, multicore=T, plot=F)
sig.results <- mqmprocesspermutation(BraCase1Results)
summary(sig.results)
#mqmplot.permutations(BraCase1Results)


####Bc-Br YBlock####
scan1<-mqmscan(augdata1,pheno.col=7,plot=T,model="dominance")
summary(scan1)

cofactorsindex<-NULL
cofactorsindex<-c(237)

#find.markerindex(augdata1,find.marker(augdata1,chr=7,pos=20))

homemadecofactors<-mqmsetcofactors(augdata1,cofactors=c(cofactorsindex)) #Step 3: designate cofactors for multivariate scan.

homescan1<-mqmscan(cross=augdata1,cofactors=homemadecofactors,pheno.col=7,cofactor.significance=0.02,verbose=T,plot=T,model="dominance") #using a dominance model tests for both additive effects and dominance, i.e., AA+AB vs BB and BB+AB vs AA, default is additive only
summary(homescan1)

bayesint(homescan1,7,qtl.index=237,prob=0.95,lodcolumn=1, expandtomarkers=T)

find.marker(augdata1,chr=7,pos=20)
A<-effectplot(F5All, pheno.col=7, mname1="scaffold_21_2195347",var.flag=c("group"))
A

BraCase2Results <- mqmpermutation(augdata1, pheno.col=7, scanfunction=mqmscan, cofactors=homemadecofactors, n.perm=100, multicore=T, plot=F)
sig.results <- mqmprocesspermutation(BraCase2Results)
summary(sig.results)
#mqmplot.permutations(BraCase2Results)



####Pleiotropy Analysis ptI####
F5All1 <- convert2cross2(F5All)

pmap <- insert_pseudomarkers(map=F5All1$gmap, step=1)
pr <- calc_genoprob(cross=F5All1, map=pmap, error_prob=0.002)
pp <- genoprob_to_alleleprob(pr)
kinship__loco <- calc_kinship(pp, "loco")


#Pleiotropy scan - LG7 - Br-Bc XBlock
BrBc_Xblock <- scan_pvl(probs = pp$`7`,
                           pheno = F5All1$pheno[,c(6,7)],
                           kinship = kinship__loco$`7`,
                           n_snp = 146)


out_lodsBrBcXblock <- BrBc_Xblock %>%
  calc_profile_lods() %>%
  add_pmap(pmap = F5All1$gmap$`7`)


BrBcXblock_Plot<-out_lodsBrBcXblock %>% 
  ggplot() + 
  geom_line(aes(x = marker_position, y = profile_lod, colour = trait), size=1.3, alpha = 1) +
  labs(x = "Marker position on LG7 (cM)", y = "LOD Score") +
  scale_colour_manual(values=c(pleiotropy="#404040",tr1="#ca0020",tr2="#0571b0"), name="Trace", labels=c("Pleiotropy", "Br-Bc", "Bc-Br")) +
  theme(legend.key=element_blank(), legend.position="bottom", axis.ticks = element_blank(),
        panel.border = element_rect(colour = "dark gray", fill=NA, size=1),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "dark gray"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "dark gray"),
        panel.background = element_rect(fill = "white", colour = "white", size = 1, linetype = "solid"))

MyRatio <- with(out_lodsBrBcXblock, diff(range(marker_position))/diff(range(profile_lod)))
BrBcXblock_PlotL <- BrBcXblock_Plot + coord_fixed(ratio=MyRatio)
BrBcXblock_PlotL

(mylrt_BrBc_Xblock <- calc_lrt_tib(BrBc_Xblock))
find_pleio_peak_tib(BrBc_Xblock, start_snp = 1)

#> (mylrt_BrBc_Xblock <- calc_lrt_tib(BrBc_Xblock))
#[1] 0.1392461
#> find_pleio_peak_tib(BrBc_Xblock, start_snp = 1)
#log10lik45 
#45 
#> pmap$`7`[45]
#scaffold_21_2195347 
#19.12355 


####Brain-braincase - Braincase-brain LG7 covariation bootstrap####
set.seed(20220406)
##
pleio_peak_index <- 45 #Br-Bc / Bc-Br
##

###############

# insert pseudomarkers
pmap <- insert_pseudomarkers(map=F5All1$gmap, step=1)
pr <- calc_genoprob(cross=F5All1, map=pmap, error_prob=0.002)
pp <- genoprob_to_alleleprob(pr)
kinship__loco <- calc_kinship(pp, "loco")

## ------------------------------------------------------------------------
phe <- F5All$pheno[,c(6,7)] #Br-Bc, Bc-Br
k <- kinship__loco$`7`
## ------------------------------------------------------------------------

###############
# simulate a phenotype
X1 <- pp$`7`[ , , pleio_peak_index]
## remove subjects with missing values of phenotype
is.na(phe[ , 1]) | is.na(phe[ , 2]) -> missing_indic
phe_nona <- phe[!missing_indic, ]
Xpre_nona <- X1[!missing_indic, ]
k_nona <- k[!missing_indic, !missing_indic]

##
gemma2::stagger_mats(Xpre_nona,Xpre_nona) -> X

calc_covs(pheno = phe_nona, kinship = k_nona) -> cc_out
(cc_out$Vg -> Vg)
(cc_out$Ve -> Ve)
# calculate Sigma
calc_Sigma(Vg = Vg, Ve = Ve,  kinship =  k_nona) -> Sigma
solve(Sigma) -> Sigma_inv
# calc Bhat 
B <- calc_Bhat(X = X, Sigma_inv = Sigma_inv, Y = phe_nona)
# Start loop
lrtBrBc_BcBr <- numeric()
for (i in 1:100){
  sim1(X = X, B = B, Sigma = Sigma) -> foo
  matrix(foo, ncol = 2, byrow = FALSE) -> Ysim
  rownames(Ysim) <- rownames(phe_nona)
  colnames(Ysim) <- c("Br_Bc", "Bc_Br")
  scan_pvl(probs = pp$`7`[!missing_indic, , ], pheno = Ysim, kinship = k_nona, start_snp = 1, n_snp = 146) -> loglik
  # in above call, s1 & nsnp come from command line args
  calc_lrt_tib(loglik) -> lrtBrBc_BcBr[i]
}

##LLRT
#mylrt_BrBc_Xblock
#0.1392461
#(pvalue <- mean(lrtBrBc_BcBr >= mylrt_BrBc_Xblock))
#p=0.54 #Indicates trait QTL are not separate, there is pleiotropy

run_num<-"Br-Bc_Bc-Br"
proc_num<-20220406

fn_out <- paste0("recla-boot_", run_num, "_", proc_num, ".csv")
write.csv(lrtBrBc_BcBr, fn_out)

save(lrtBrBc_BcBr, file="lrtBrBc_BcBr_20220406.rda")
