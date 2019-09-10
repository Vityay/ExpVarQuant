library(edgeR)
library(gamlss)
library(gamlss.tr)
library(parallel)
library(lmtest)
# library(aods3)

setwd("/media/yuri/Data/Dropbox/Figures_Mouse_BCV/Cor's Data/")

RNAseq <- read.table("20170306_RNASeq_Gertrud_rawcounts_0FPM.txt", header=T, sep="\t", row.names=1)
RNAseq <- RNAseq[grep("ENSMUSG", rownames(RNAseq)),]

samples <- read.csv("Phenotype.csv", header=T)
X <- paste(samples$mut, samples$Age, sep="_")
X <- X[match(samples$Library, colnames(RNAseq))]


RNAseq<- DGEList(RNAseq, group=X)
RNAseq <- calcNormFactors(RNAseq, method="upperquartile")

RNAseq$samples$lib.offset <- log(RNAseq$samples$lib.size*RNAseq$samples$norm.factors)
RNAseq$CPM <- sapply(seq_along(RNAseq$counts[1,]), function(i)  (RNAseq$counts[,i]*1e06)/(RNAseq$samples$lib.size[i]*RNAseq$samples$norm.factors[i])   )

ctr1 <- gamlss.control(n.cyc=100, gd.tol=1e05)
ctr2 <- glim.control(bf.cyc = 1e04)
warn1 <- function(w) NULL; warn2 <- function(e) NULL

RNAseq_gamlss <- mclapply(seq_along(RNAseq$counts[,1]), function(i) {
    cat(i, "\r")

    dat <- data.frame(x=RNAseq$samples$group,
                      mut = sapply(strsplit(as.character(RNAseq$samples$group), "_"), function(x) x[1]),
                      age = sapply(strsplit(as.character(RNAseq$samples$group), "_"), function(x) x[2]),
                      y=RNAseq$counts[i,],
                      lib.offset=RNAseq$samples$lib.offset)

    sink(file="/dev/null")

    m0.wt <- tryCatch( gamlss(y~0+age+offset(lib.offset), sigma.fo = ~ 1, data=subset(dat, mut=="wt"),  , method=RS(), control=ctr1, i.control= ctr2), warning= warn1, error= warn2)
    m1.wt <- tryCatch( gamlss(y~0+age+offset(lib.offset), sigma.fo = ~ 0+age, data=subset(dat, mut=="wt"), family=NBI(), method=RS(), control=ctr1, i.control= ctr2), warning= warn1, error= warn2)

    m0.kin <- tryCatch( gamlss(y~0+age+offset(lib.offset), sigma.fo = ~ 1, data=subset(dat, mut=="Kin"), family=NBI(), method=RS(), control=ctr1, i.control= ctr2), warning= warn1, error= warn2)
    m1.kin <- tryCatch( gamlss(y~0+age+offset(lib.offset), sigma.fo = ~ 0+age, data=subset(dat, mut=="Kin"), family=NBI(), method=RS(), control=ctr1, i.control= ctr2), warning= warn1, error= warn2)

    m0.young <- tryCatch( gamlss(y~0+mut+offset(lib.offset), sigma.fo = ~ 1, data=subset(dat, age=="young"), family=NBI(), method=RS(), control=ctr1, i.control= ctr2), warning= warn1, error= warn2)
    m1.young <- tryCatch( gamlss(y~0+mut+offset(lib.offset), sigma.fo = ~ 0+mut, data=subset(dat, age=="young"), family=NBI(), method=RS(), control=ctr1, i.control= ctr2), warning= warn1, error= warn2)

    m0.old <- tryCatch( gamlss(y~0+mut+offset(lib.offset), sigma.fo = ~ 1, data=subset(dat, age=="old"), family=NBI(), method=RS(), control=ctr1, i.control= ctr2), warning= warn1, error= warn2)
    m1.old <- tryCatch( gamlss(y~0+mut+offset(lib.offset), sigma.fo = ~ 0+mut, data=subset(dat, age=="old"), family=NBI(), method=RS(), control=ctr1, i.control= ctr2), warning= warn1, error= warn2)

    res <- data.frame(wt_young.mu = NA, wt_old.mu = NA, wt_young.ECV = NA, wt_old.ECV = NA, p.ECV.wt_old_to_young = NA,
                      Kin_young.mu = NA, Kin_old.mu = NA, Kin_young.ECV = NA, Kin_old.ECV = NA, p.ECV.Kin_old_to_young = NA,
                      p.ECV.wt_Kin_young = NA, p.ECV.wt_Kin_old = NA)

    if(!any(sapply(list(m0.wt, m1.wt), is.null))) {
      res$wt_young.mu  = exp(m1.wt$mu.coefficients + log(1e06))[2]
      res$wt_old.mu    = exp(m1.wt$mu.coefficients + log(1e06))[1]
      res$wt_young.ECV = sqrt(exp(m1.wt$sigma.coefficients[2]))
      res$wt_old.ECV   = sqrt(exp(m1.wt$sigma.coefficients[1]))
      res$p.ECV.wt_old_to_young = lrtest(m0.wt, m1.wt)$'Pr(>Chisq)'[2]
    }

    if(!any(sapply(list(m0.kin, m1.kin), is.null))) {
      res$Kin_young.mu  = exp(m1.kin$mu.coefficients + log(1e06))[2]
      res$Kin_old.mu    = exp(m1.kin$mu.coefficients + log(1e06))[1]
      res$Kin_young.ECV = sqrt(exp(m1.kin$sigma.coefficients[2]))
      res$Kin_old.ECV   = sqrt(exp(m1.kin$sigma.coefficients[1]))
      res$p.ECV.Kin_old_to_young = lrtest(m0.kin, m1.kin)$'Pr(>Chisq)'[2]
    }

    if(!any(sapply(list(m0.young, m1.young), is.null))) {
      res$p.ECV.wt_Kin_young = lrtest(m0.young, m1.young)$'Pr(>Chisq)'[2]
    }

    if(!any(sapply(list(m0.old, m1.old), is.null))) {
      res$p.ECV.wt_Kin_old = lrtest(m0.old, m1.old)$'Pr(>Chisq)'[2]
    }

    sink()

    return(res)

}, mc.cores=7)

RNAseq_gamlss <- do.call(rbind, RNAseq_gamlss)
rownames(RNAseq_gamlss) <- rownames(RNAseq$counts)

# let's polish a bit

for(i in c(3,4,8,9)) {
    RNAseq_gamlss[, i][RNAseq_gamlss[,i] > 3] <- NA
    RNAseq_gamlss[, i][RNAseq_gamlss[,i] < 1e-02] <- NA
}

RNAseq_gamlss_only_active <- RNAseq_gamlss[rowSums(RNAseq$CPM < 1)==0,]


boxplot(RNAseq_gamlss_only_active[, c(3,4,8,9)], outline=F)

cor(log(RNAseq_gamlss_only_active[,1]),  RNAseq_gamlss_only_active[,3], use="complete.obs")
cor(log(RNAseq_gamlss_only_active[,2]),  RNAseq_gamlss_only_active[,4], use="complete.obs")
cor(log(RNAseq_gamlss_only_active[,6]),  RNAseq_gamlss_only_active[,8], use="complete.obs")
cor(log(RNAseq_gamlss_only_active[,7]),  RNAseq_gamlss_only_active[,9], use="complete.obs")

rm(ctr1,ctr2,samples,warn1,warn2,X)


sum(p.adjust(RNAseq_gamlss_only_active$p.ECV.wt_old_to_young, "fdr") <= 0.1, na.rm=T)
sum(p.adjust(RNAseq_gamlss_only_active$p.ECV.Kin_old_to_young, "fdr") <= 0.1, na.rm=T)
sum(p.adjust(RNAseq_gamlss_only_active$p.ECV.wt_Kin_young, "fdr") <= 0.1, na.rm=T)
sum(p.adjust(RNAseq_gamlss_only_active$p.ECV.wt_Kin_old, "fdr") <= 0.1, na.rm=T)

library(gplots)


tmp1 <- as.character(na.omit(rownames(RNAseq_gamlss_only_active)[RNAseq_gamlss_only_active$p.ECV.wt_old_to_young <= 0.01]))
tmp2 <- as.character(na.omit(rownames(RNAseq_gamlss_only_active)[RNAseq_gamlss_only_active$p.ECV.Kin_old_to_young <= 0.01]))

venn(list(tmp1, tmp2))

tmp1 <- na.omit(RNAseq_gamlss_only_active[RNAseq_gamlss_only_active$p.ECV.wt_old_to_young  <= 0.01, c(8,9)] )
tmp2 <- na.omit(RNAseq_gamlss_only_active[RNAseq_gamlss_only_active$p.ECV.Kin_old_to_young  <= 0.01, c(8,9)] )

write.table(rownames(tmp1[tmp1[,2]-tmp1[,1] > 0, ]), "wt_old_young_higher_noise.txt", sep="\t", quote=F, row.names=F, col.names=F )
write.table(rownames(tmp1[tmp1[,2]-tmp1[,1] < 0, ]), "wt_old_young_lower_noise.txt", sep="\t", quote=F, row.names=F, col.names=F )

write.table(rownames(tmp2[tmp2[,2]-tmp2[,1] > 0, ]), "Kin_old_young_higher_noise.txt", sep="\t", quote=F, row.names=F, col.names=F )
write.table(rownames(tmp2[tmp2[,2]-tmp2[,1] < 0, ]), "Kin_old_young_lower_noise.txt", sep="\t", quote=F, row.names=F, col.names=F )


tmp1 <- as.character(na.omit(rownames(RNAseq_gamlss_only_active)[p.adjust(RNAseq_gamlss_only_active$p.ECV.wt_old_to_young, "fdr") <= 0.1]))
tmp2 <- as.character(na.omit(rownames(RNAseq_gamlss_only_active)[p.adjust(RNAseq_gamlss_only_active$p.ECV.Kin_old_to_young, "fdr") <= 0.1]))

venn(list(tmp1, tmp2))

tmp1 <- na.omit(RNAseq_gamlss_only_active[p.adjust(RNAseq_gamlss_only_active$p.ECV.wt_Kin_young, "fdr") <= 0.1,c(4,9)] )
tmp2 <- na.omit(RNAseq_gamlss_only_active[p.adjust(RNAseq_gamlss_only_active$p.ECV.wt_Kin_old, "fdr") <= 0.1,c(4,9)] )

boxplot(list(young=log2(tmp1[,1]/tmp1[,2]), old=log2(tmp2[,1]/tmp2[,2])))

tmp1 <- na.omit(RNAseq_gamlss_only_active[RNAseq_gamlss_only_active$p.ECV.wt_Kin_young <= 0.01, c(4,9)] )
tmp2 <- na.omit(RNAseq_gamlss_only_active[RNAseq_gamlss_only_active$p.ECV.wt_Kin_old <= 0.01, c(4,9)] )

write.table(rownames(tmp2[tmp2[,1]-tmp2[,2] > 0, ]), "Kin_wt_old_lower_noise.txt", sep="\t", quote=F, row.names=F, col.names=F )
write.table(rownames(tmp2[tmp2[,1]-tmp2[,2] < 0, ]), "Kin_wt_old_higher_noise.txt", sep="\t", quote=F, row.names=F, col.names=F )




rm(tmp1, tmp2)
save.image("Cor_gamlss.Rdat")



load("/media/yuri/Data/home1/MOUSE_NOISE_PROJECT/Figures_data_paper_07_Mar_2017/Figure1_NBI_HF.Rdat")
genes <- rownames(d.HF$counts)

BCV.NBI.Std      <- as.data.frame(do.call(cbind, lapply(gamlss.NBI_HF$NBI, function(x) {res <- x$BCV.Std; res[res < 1e-02 | res > 3] <- NA; res  }) ) )
rownames(BCV.NBI.Std ) <- genes

BCV.NBI.HF      <- as.data.frame(do.call(cbind, lapply(gamlss.NBI_HF$NBI, function(x) {res <- x$BCV.HF; res[res < 1e-02 | res > 3] <- NA; res  }) ) )
rownames(BCV.NBI.HF ) <- genes


load("/media/yuri/Data/Dropbox/Figures_Mouse_BCV/Cor's Data/Cor_gamlss.Rdat")

RNAseq_gamlss_only_active2 <- merge(BCV.NBI.Std, RNAseq_gamlss_only_active[,c(3,4,8,9)], by="row.names")
RNAseq_gamlss_only_active3 <- merge(BCV.NBI.HF, RNAseq_gamlss_only_active[,c(3,4,8,9)], by="row.names")

plot(hclust(as.dist(1-cor(RNAseq_gamlss_only_active2[,-1], use="pairwise.complete.obs") ), method="complete"), hang=-1)
plot(hclust(as.dist(1-cor(RNAseq_gamlss_only_active3[,-1], use="pairwise.complete.obs") ), method="complete"), hang=-1)
