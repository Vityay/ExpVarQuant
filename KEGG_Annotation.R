library("org.Mm.eg.db")
library(KEGGREST)
library(h2o)

setwd("path to gamlss_age.Rdat file")
load("gamlss_age.Rdat")

# convert Ensembl genes IDs to Entrez gene IDs.
org.Mm <- as.list(org.Mm.egENSEMBL2EG)
org.Mm <- data.frame(EntrezID = unlist(org.Mm))
gamlss_NB_ENTREZ <- merge(org.Mm, gamlss_NB_clean, by="row.names")


# Retrieve KEGG pathways for mouse.
kegg_path <- keggLink("pathway", "mmu")
kegg_path <- split(kegg_path, kegg_path)
kegg_path <- lapply(kegg_path, names)
kegg_names <- sapply(seq_along(names(kegg_path)), function(i) {
cat(i, "\n")
res <- keggGet(names(kegg_path)[i])[[1]]$PATHWAY_MAP
})

# Map genes to KEGG pathways. If gene is present in a given pathway 1 is assigned, otherwise 0 is assigned.
KEGG_model_matrix <- sapply(seq_along(kegg_path), function(i) {
res <- rep(0, nrow(gamlss_NB_ENTREZ))
idx <- as.vector(na.omit(match(kegg_path[[i]], paste0("mmu:", gamlss_NB_ENTREZ$EntrezID))))
res[idx] <- 1
res
})
colnames(KEGG_model_matrix) <- kegg_names

# Remove pathways represented by less than 5 genes.
KEGG_model_matrix <- KEGG_model_matrix[,colSums(KEGG_model_matrix) >= 5]

# Remove genes which are not mapped to any of the KEGG pathways. Create a vector Y for dependent variable: Y=〖log〗_2 (〖cv(μ)〗_old/〖cv(μ)〗_young)  and a matrix X containing genes’ KEGG pathway annotation.
idx <- rowSums(KEGG_model_matrix) > 0
X <- KEGG_model_matrix[idx,]
Y <- gamlss_NB_ENTREZ$LR.cv[idx]

# Initialize h2o and fit ridge regression (alpha=0) KEGG model with 10-fold cross-validation (nfolds=10).
h2o.init()
df <- as.h2o(data.frame(Y=Y, X))
m1 <- h2o.glm(y=1, x=2:(ncol(X)+1), training_frame = df, family = "gaussian", nfolds = 10, alpha = 0,
lambda_search=TRUE, nlambdas=100, standardize = FALSE, remove_collinear_columns = FALSE, seed = 123)

# Predict 〖log〗_2 (〖cv(μ)〗_old/〖cv(μ)〗_young)  from the KEGG model (m1).
pred_m1 <- as.vector(predict(m1, newdata=df))

# Correlation between observed and predicted 〖log〗_2 (〖cv(μ)〗_old/〖cv(μ)〗_young).
cor_m1 <- cor(Y, pred_m1); cor_m1

# Extract KEGG model coefficients. The first coefficient corresponds to an intercept.
coef_m1 <- h2o.coef(m1)[-1]
coef_m1 <- coef_m1[coef_m1!=0]
coef_m1 <- sort(coef_m1)

# Plot the results.
par(mfrow=c(2,1))
# plot KEGG model vs observed values
plot(Y, pred_m1,  col="#00000088",
xlab="observed: log2(CV.old/CV.young)", ylab="predicted: KEGG model")
abline(lm(pred_m1~Y, data=data.frame(Y, pred_m1)), col=3, lwd=2)
legend("topleft", legend=paste("r =", round(cor_m1,3)))

# plot top 20 KEGG pathways associated with age-related increase in gene noise
par(mar=c(4,16,2,1))
barplot(coef_m1[(length(coef_m1)-19):length(coef_m1)], horiz=T, las=1,
cex.names=0.8, xlab="model coefficients")
mtext("top 20 KEGG pathways associated with age-related increase in gene noise",
side=3, adj=1, line=0, cex=1, font=2)