# install and load packages
require(limma);require(edgeR); require(EMMREML);require(gridExtra); require(grid); require(ggdendro); require(ggplot2);require(DESeq2)

##Load in the data
load("~/Dropbox/data_for_CA/dex_RNAseq.RData")

#what's in the package?
dim(counts) #a p x n matrix of counts for each gene, for 20 RNA-seq samples (16052 x 20)
dim(gene_lengths) #a p x 2 matrix giving the total exon length, in bp, for each gene in 'counts'

## Calculate RPKM for each gene and remove genes with a median RPKM<2
RPKM<-apply(counts,2,function(x){(as.numeric(x)*10^9)/(sum(as.numeric(x))*gene_lengths[rownames(counts),2])}); rownames(RPKM)=rownames(counts)
RPKM_filter=RPKM[which(apply(RPKM,1,function(a){quantile(a,0.5)})>=2),]
dim(RPKM_filter) ## print dimensions of RPKM_filter object

## Remove genes from the counts table that have a median RPKM < 2
counts_filter=counts[rownames(RPKM_filter),]
dim(counts_filter) ## print dimensions of counts_filter object

## Plot a few examples of the read counts and RPKMs for the first sample (column 1)
par(mfrow=c(2,1))
hist(log10(counts[,1]),breaks = 100,main="reads",xlab="log10(count)",col="steelblue4") ## there are a lot of genes with 0 (or very few) reads
hist(log10(counts_filter[,1]),breaks = 100,main="reads filtered at RPKM ≥ 2",xlab="log10(count)",col="steelblue4") ## our filter removed the lowly or not expressed genes

## Other normalization techniques
## voom normalization
v_norm=voom(calcNormFactors(DGEList(counts=counts)),plot=FALSE) #normalization is done with TMM using edgeR's calcNormFactors, then the data are transformed for linear modeling using voom 
v_norm=v_norm$E[rownames(RPKM_filter),] ## and filter to only include genes with a median RPKM ≥ 2

#distributions look good and comparable between samples, e.g.:
par(mfrow=c(1,1))
qqplot(v_norm[,1],v_norm[,2])
abline(0,1)

## PCA on the voom normalized matrix
pca=prcomp(t(v_norm)) 
qplot(pca$x[,1],pca$x[,2],col=metadata_dexetoh[rownames(pca$x),"Treatment"])+scale_colour_brewer("",palette = "Set1")+theme_bw(20)+geom_point(size=4)+ xlab("PC1") + ylab("PC2") ## Not much signal of treatment in these first 2 PCs, but it is strong in PC3:
qplot(pca$x[,2],pca$x[,3],col=metadata_dexetoh[rownames(pca$x),"Treatment"])+scale_colour_brewer("",palette = "Set1")+theme_bw(20)+geom_point(size=4)+ xlab("PC2") + ylab("PC3") ## Not much signal of treatment in these first 2 PCs, but it is strong in PC3:
qplot(pca$x[,1],pca$x[,2],col=metadata_dexetoh[rownames(pca$x),"ID"]) +scale_colour_brewer("",palette = "Paired")+theme_bw(20)+geom_point(size=4)+ xlab("PC1") + ylab("PC2") ## A strong individual effect on PC1, where the two samples from one individual are clear outliers (Ve12). Also, notice how the two samples from the same individual cluster together. 

## how much variance is explained by each PC?
summary(pca)


## modeling using limma
design=model.matrix(~metadata_dexetoh$Treatment+metadata_dexetoh$ID) ## control for ID, but you lose power because here it is included as a fixed effect with 10 levels. [NB: the design matrix can also be used to calculate precision weights in voom, which can be applied to limma modeling]

ptm <- proc.time()
limma=data.frame(eBayes(lmFit(v_norm,design)))
ptmlim <- proc.time()-ptm

## plot p-value distributions
par(mfrow=c(1,1))
hist(limma$p.value.metadata_dexetoh.Treatmentethanol,breaks = 100,main="effect of treatment",xlab="p-value",col="steelblue4") 

# 100 permutations, no enrichment over the background
perm_pval_limma=Reduce(c,lapply(1:100,function(perm){
   tr=Reduce(c,lapply(1:(nrow(metadata_dexetoh)/2),function(x) sample(c("dex","etoh")))) ## permute the reatment variable while keeping the paired nature of the data
   design=model.matrix(~tr+metadata_dexetoh$ID)
   limma_perm=as.data.frame(eBayes(lmFit(v_norm,design)))
   return(limma_perm$p.value.tretoh)
}))
hist(perm_pval_limma,breaks = 100,main="effect of permuted treatment",xlab="p-value",col="steelblue4") 

## another way to look at it is hierarchical clustering
dend1<-hclust(dist(t(v_norm)),met="ave")
hcdata<-dendro_data(dend1,type="rectangle")
hcdata$labels<-merge(x=hcdata$labels,y=metadata_dexetoh,by.x="label",by.y=0)
ggplot() +
  geom_segment(data=segment(hcdata), aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data = label(hcdata), aes(x=x, y=y, label=ID, colour = as.factor(Treatment), hjust=0), size=4) +
  geom_point(data = label(hcdata), aes(x=x, y=y), size=0.1, shape = 16) +
  coord_flip() +
  scale_y_reverse(expand=c(0.2, 0)) +
  scale_colour_manual("Treatment",values=c("firebrick3","dodgerblue3","darkorange2","forestgreen","darkorchid")) + 
  theme_dendro()

## You can clearly see that samples from the same individual always cluster together, rather than samples with the same treatment

## So we will need to control for relatedness and the paired nature of the samples
## Step 1, create a Z matrix for running EMMA
Z_matrix=matrix(nrow=nrow(metadata_dexetoh),ncol=length(unique(metadata_dexetoh$ID)))
colnames(Z_matrix)=unique(metadata_dexetoh$ID)
rownames(Z_matrix)=metadata_dexetoh$ID
for (r in 1:nrow(Z_matrix)) {
  for (c in 1:ncol(Z_matrix)) {
    if (rownames(Z_matrix)[r]==colnames(Z_matrix)[c]) {Z_matrix[r,c]=1}else {Z_matrix[r,c]=0}} }
rownames(Z_matrix)=rownames(metadata_dexetoh)

# check that your kin and Z matrices are in the same order
colnames(Z_matrix)==rownames(kin)

## model effects of Dex and in vivo vs. in vitro
ptm <- proc.time()
EMMA_model=t(apply(v_norm,1,function(y){
  design=model.matrix(~metadata_dexetoh$Treatment) ## we don't control for individual ID, that is controlled for in the Z_matrix
  emma=emmreml(y=y,X=design,Z=as.matrix(Z_matrix),K=as.matrix(kin),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
  p=emma$pvalbeta
  varb=emma$varbetahat
  b=emma$betahat
  return(c(b,varb,p[,"none"]))
}))
ptmmixed<-proc.time() - ptm

#compare runtimes
barplot(c(ptmlim[3],ptmmixed[3]),names=c('limma','emma'),horiz = T,xlab='runtime (seconds)',cex.axis=1.5,cex.names=1.2,cex.lab=1.5)

EMMA_model=as.data.frame(EMMA_model)
colnames(EMMA_model)=c('beta_intercept','beta_treatment','var_beta_intercept','var_beta_treatment','pval_intercept','pval_treatment')

## plot the pvalue distribution of the EMMA model 
par(mfrow=c(1,1))
hist(EMMA_model$pval_treatment,100,col="steelblue4",xlab="pval",main="effect of treatment")
##  and add in the limma p-value distribution
hist(limma$p.value.metadata_dexetoh.Treatmentethanol,100,col="lemonchiffon",add=T)

## 100 permutations for EMMA (don't do this in the class, it takes like 10min to run, but you can use the pre-loaded "perm_pval_emma" vector)
# perm_pval_emma=Reduce(c,lapply(1:100,function(perm){
#   tr=Reduce(c,lapply(1:(nrow(metadata_dexetoh)/2),function(x) sample(c("dex","etoh")))) ## permute the reatment variable while keeping the paired nature of the data
#   pvals=t(apply(v_norm,1,function(y){
#   design=model.matrix(~tr) ## we don't control for individual ID, that is controlled for in the Z_matrix
#   emma=emmreml(y=y,X=design,Z=as.matrix(Z_matrix),K=as.matrix(kin),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
#   p=emma$pvalbeta
#   return(c(p[2,"none"]))
# }))
#   return(pvals)}))


## There is an enrichment of low p-values here
hist(perm_pval_emma,breaks = 100,main="effect of permuted treatment",xlab="p-value",col="steelblue4",freq=F) 

## plot the pvalue distribution of the EMMA model to compare the permuted values against the observed data 
par(mfrow=c(1,2))
hist(EMMA_model$pval_treatment,100,col="steelblue4",xlab="pval",main="EMMA (mixed model)",freq=F,ylim=c(0,25))
##  and add in the permuted p-value distribution
hist(perm_pval_emma,100,col="lemonchiffon",add=T,freq=F)
legend("topright",c('observed','permuted'),fill=c('steelblue4','lemonchiffon'),bty='n')

hist(limma$p.value.metadata_dexetoh.Treatmentethanol,100,col="steelblue4",xlab="pval",main="limma (linear model)",freq=F,ylim=c(0,25))
##  and add in the permuted p-value distribution
hist(perm_pval_limma,100,col="lemonchiffon",add=T,freq=F)
legend("topright",c('observed','permuted'),fill=c('steelblue4','lemonchiffon'),bty='n')

## multiply Z and K matrices for running MACAU Poisson model:
ZK_mat=Z_matrix%*%as.matrix(kin)%*%t(Z_matrix)
## look at MACAU data
hist(macau_model$pvalue,100,col="steelblue4",xlab="pval",main="effect of treatment")
## and the permutations (depletion of low p-values for some reason...)
hist(perm_macau$pvalue,100,col="steelblue4",xlab="pval",main="permuted effect of treatment")

#limma, B-H correction
limcorrected<-unlist(lapply(limma$p.value.metadata_dexetoh.Treatmentethanol,function(x){((sum(perm_pval_limma<=x))/100)/(sum(limma$p.value.metadata_dexetoh.Treatmentethanol<=x))}))

#emma
emmcorrected<-unlist(lapply(EMMA_model$pval_treatment,function(x){((sum(perm_pval_emma<=x))/100)/(sum(EMMA_model$pval_treatment<=x))}))

#MACAU
maccorrected<-unlist(lapply(macau_model$pvalue,function(x){((sum(perm_macau$pvalue<=x))/100)/(sum(macau_model$pvalue<=x))}))

#Number of genes passing 10% FDR
fdrval<-0.1
barplot(c(sum(limcorrected<fdrval),sum(emmcorrected<fdrval),sum(maccorrected<fdrval)),names=c('limma','emma','macau'),ylab='no. genes FDR 10%',cex.axis=1.5,cex.lab=1.5,cex.names=1.5)

#Number of genes passing 1% FDR
fdrval<-0.01
barplot(c(sum(limcorrected<fdrval),sum(emmcorrected<fdrval),sum(maccorrected<fdrval)),names=c('limma','emma','macau'),ylab='no. genes FDR 10%',cex.axis=1.5,cex.lab=1.5,cex.names=1.5)

fdrval<-seq(0.01,0.2,by=0.01)
numgenes<-matrix(nrow=3,ncol=length(fdrval),NA)
for(i in 1:length(fdrval)){
  fdrthres<-fdrval[i]
  numgenes[1,i]<-sum(limcorrected<fdrthres)
  numgenes[2,i]<-sum(emmcorrected<fdrthres)
  numgenes[3,i]<-sum(maccorrected<fdrthres)
}

par(mfrow=c(1,1))
plot(numgenes[1,]~fdrval,type='l',xlab='fdr threshold',ylab='significant genes',lwd=5,frame.plot=F,cex.lab=1.5,cex.axis=1.5,col='steelblue4',ylim=c(0,2500),xlim=c(0,.2))
lines(numgenes[2,]~fdrval,type='l',lwd=5,col='purple4')
lines(numgenes[3,]~fdrval,type='l',col='goldenrod',lwd=5)
legend("bottomright",c('limma','emma','PMM'),fill=c('steelblue4','purple4','goldenrod'),bty='n',cex=1.5)