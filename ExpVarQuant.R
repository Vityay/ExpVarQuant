library(gamlss)
library(edgeR)

setwd("path to gamlss_age.Rdat file")
load("gamlss_age.Rdat")

# Uncomment this code when working in Rstudio to prevent it from locking up
# sink(file = "diversionOfOutput.txt", append = TRUE) 

# Calculate library size and normalization factor (see edgeR manual).
Counts_edgeR <- calcNormFactors(Counts_edgeR, method="upperquartile")
# Calculate counts per million (CPM).
Counts_edgeR$CPM <- cpm.DGEList(Counts_edgeR)
# Calculate offset variable for each library as natural log of the product of library size and normalization factor.
ofs <- log(Counts_edgeR$samples$lib.size*Counts_edgeR$samples$norm.factors)
Counts_edgeR$samples$offset <- ofs

# Remove genes with 0 counts in any of the sample as analysis of zero inflated samples could introduce a significant bias.
idx <- !apply(Counts_edgeR$CPM,1, function(x) any(x == 0) )
Counts_edgeR$counts <- Counts_edgeR$counts[idx,]
Counts_edgeR$CPM <- Counts_edgeR$CPM[idx,]

# Remove lowly expressed genes as for these genes technical variations might have a substantial impact on gene noise estimation.
idx <- rowMeans(Counts_edgeR$CPM) > 1
Counts_edgeR$counts <- Counts_edgeR$counts[idx,]
Counts_edgeR$CPM <- Counts_edgeR$CPM[idx,]

# Remove idx and ofs variables from the workspace.
rm(idx, ofs)

# Estimate age effects on mean and overdispersion parameters with GAMLSS for each gene.
gene_i <- seq_along(Counts_edgeR$counts[,1])

# To try algorithm for just some genes, change gene_i variable. For example, set gene_i <- c(1:100) to estimate GAMLSS models for the first hundred genes.
gamlss_NB <- lapply(gene_i, function(i) {

# For each gene (i) a table containing: x - a factor of interest (age); y - RNA-seq. counts and offset (ofs) is created.
dat <- data.frame(
x = Counts_edgeR$samples$Age,
y = Counts_edgeR$counts[i,],
ofs = Counts_edgeR$samples$offset
)
# x is releveled to use young mice as a reference.
dat$x <- relevel(dat$x, ref = c("young"))

# Fit negative binomial (family = NBI()) GAMLSS model, which accounts for age effects on mean and overdispersion (non-Poisson noise).
# fo = y~0+x+offset(ofs) specifies model for mean and sigma.fo=~0+x for overdispersion, offset - offset(ofs) normalize counts to library size. sigma.start = 0.1 provides starting value for overdispersion estimation (default is 1). n.cyc – number of fitting algorithm cycles, see help(gamlss).
# In some cases, fitting of NB model may fail and tryCatch(..., warning= function(w) NULL, error= function(e) NULL) will return NULL as a result.
m0 <- tryCatch(
gamlss(fo = y ~ 0+x+offset(ofs), sigma.fo = ~ 0+x, data=dat,
family = NBI(), sigma.start = 0.1, n.cyc = 100),
warning= function(w) NULL, error= function(e) NULL
)

# Fit reduced model by omitting age factor from the estimation of overdispersion: sigma.fo = ~ 1. In essence, this model corresponds to the GLM model implemented in edgeR.
m1 <- tryCatch(
gamlss(fo = y ~ 0+x+offset(ofs), sigma.fo = ~ 1, data=dat,
family = NBI(), sigma.start = 0.1, n.cyc = 100),
warning= function(w) NULL, error= function(e) NULL
)

# Fit reduced model by omitting age factor from the estimation of mean: fo = y ~ offset(ofs).
m2 <- tryCatch(
gamlss(fo = y ~ offset(ofs), sigma.fo = ~ 0+x, data=dat,
family = NBI(), sigma.start = 0.1, n.cyc = 100),
warning= function(w) NULL, error= function(e) NULL
)

# Fit null model.
m3 <- tryCatch(
gamlss(fo = y ~ offset(ofs), sigma.fo = ~ 1, data=dat,
family = NBI(), sigma.start = 0.1, n.cyc = 100),
warning= function(w) NULL, error= function(e) NULL
)

# Create data frame res to store the results.
res <- data.frame(
cpm.young = NA,
cpm.old=NA,
LR.cpm = NA,
p_gamlss.cpm = NA,
p_glm.cpm = NA,
CV.young = NA,
CV.old = NA,
LR.cv = NA,
p.cv = NA
)

# Because fitting of the NB model may fail for some genes, check whether all models were fitted successfully. 
if(!any(sapply(list(m0,m1,m2,m3), is.null))) 
{
# Write GAMLSS estimations of gene’s mean (CPM) counts from the m0 model.
res$cpm.young = exp(m0$mu.coefficients+log(1e06))[[1]]
res$cpm.old = exp(m0$mu.coefficients+log(1e06))[[2]]

# Calculate log2 ratio for changes in CPMs between old and young mice.
res$LR.cpm = log2(exp(m0$mu.coefficients+log(1e06))[[2]]/exp(m0$mu.coefficients+log(1e06))[[1]])

	# GAMLSS log-likelihood ratio (LR) test for a significance of an age effect on gene’s mean (CPM) counts.
	# p_gamlss.cpm – p value of LR test statistic: D_μ=-2log⁡[L(μ_0,α_j  ┤|  X_ij)/L(μ_j,α_j  ┤|  X_ij), comparing m0 and m2 models.
res$p_gamlss.cpm = pchisq(2*(logLik(m0)-logLik(m2)), df=m0$df.fit-m2$df.fit, lower=F)[[1]]

# GLM log-likelihood ratio (LR) test for a significance of age effect on gene’s mean (CPM) counts.
# p_glm.cpm – p value of LR test statistic: D_(μ_GLM )=-2log⁡[L(μ_0,α_0  ┤|  X_ij)/L(μ_j,α_0  ┤|  X_ij), comparing m1 and m3 models.
res$p_glm.cpm = pchisq(2*(logLik(m1)-logLik(m3)), df=m1$df.fit-m3$df.fit, lower=F)[[1]]

	# Write GAMLSS estimations of gene’s non-Poisson noise from the m0 model: cv(μ)=√α.
res$CV.young  = sqrt(exp(m0$sigma.coefficients)[[1]])
res$CV.old = sqrt(exp(m0$sigma.coefficients)[[2]])

# Calculate log2 ratio for changes in cv(μ) between old and young mice.
res$LR.cv = log2(sqrt(exp(m0$sigma.coefficients)[[2]])/sqrt(exp(m0$sigma.coefficients)[[1]]))

# GAMLSS log-likelihood ratio (LR) test for a significance of an age effect on non-Poisson noise.
	# p.cv – p value of LR test statistic: D_α=-2log⁡[L(μ_j,α_0  ┤|  X_ij)/L(μ_j,α_j  ┤|  X_ij), comparing m0 and m1 models.
res$p.cv = pchisq(2*(logLik(m0)-logLik(m1)), df=m0$df.fit-m1$df.fit, lower=F)[[1]]
}
res
})

# Transform list gamlss_NB containing GAMLSS estimations for each gene to data frame
gamlss_NB <- do.call(rbind, gamlss_NB)
rownames(gamlss_NB) <- rownames(Counts_edgeR$counts)[gene_i]

# Because GAMLSS fitting might fail for some genes or yield inflated estimates of overdispersion, the results have to be cleaned.
# First, remove genes, for which GAMLSS model has failed.
gamlss_NB_clean <- na.omit(gamlss_NB)

# Second, remove genes, for which estimates of cv(μ)=√α were either inflated > 3 or close to Poisson < 10-3.
idx <- gamlss_NB_clean$CV.young > 3 | gamlss_NB_clean$CV.young < 1e-03 | gamlss_NB_clean$CV.old > 3 | gamlss_NB_clean$CV.old < 1e-03
gamlss_NB_clean <- gamlss_NB_clean[!idx,]

# Finally, calculate false discovery rates to account for multiple hypothesis testing.
gamlss_NB_clean$padj_gamlss.cpm <- p.adjust(gamlss_NB_clean$p_gamlss.cpm, "fdr")
gamlss_NB_clean$padj_glm.cpm <- p.adjust(gamlss_NB_clean$p_glm.cpm, "fdr")
gamlss_NB_clean$padj.cv <- p.adjust(gamlss_NB_clean$p.cv, "fdr")

# Save the results, which later can be loaded to R and used for further analysis
rm(idx, gene_i, gamlss_NB)
save.image("gamlss_age.Rdat")
