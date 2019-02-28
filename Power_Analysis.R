# Power analysis for non-Poisson noise
library(gamlss)
library(parallel)

# Function calculating likelihood ratio test for varying non-Poisson noise in Negative Binomial (NB) model
# N - sample size
# d - fold change for non-Poisson noise
# mu - average counts
# sigma - non-Poisson noise (coefficient of variation in counts - cv(X))
# offset - flag indicating whether to include (TRUE) or not (FALSE) offset into the models
# ofs - offset value

NB_func_sigma <- function(N, d, mu, sigma, offset = FALSE, ofs = log(1e07) ) {
  y1 <- rNBI(N, mu, sigma^2) # not that for gamlss paramitrization of NB, non-Poisson noise is estimated as cv^2(X).
  y2 <- rNBI(N, mu, (sigma*d)^2)
  x <- c(rep("a", N), rep("b", N)) # x - defines groups for which non-Poisson noise is compared.
  y <- c(y1, y2)
  dat <- data.frame(x=x, y=y, ofs = ofs)

  if(offset) {
    # m0 - reduced model, which does not account for the differences in non-Poisson noise. The model also includes offset.
    m0 <- gamlss(y ~ x + offset(ofs), data = dat, family=NBI(), mu.start = mu, sigma.start = sigma)
    # m1 - full model, which does  account for the differences in non-Poisson noise. The model also includes offset.
    m1 <- gamlss(y ~ x + offset(ofs), sigma.fo = ~x, data = dat, family=NBI(), mu.start = mu, sigma.start = sigma)
  } else {
    # m0 - reduced model, which does not account for the differences in non-Poisson noise. The model does not include offset.
    m0 <- gamlss(y ~ x, data = dat, family=NBI(), mu.start = mu, sigma.start = sigma)
    # m1 - full model, which does  account for the differences in non-Poisson noise. The model does not include offset.
    m1 <- gamlss(y ~ x, sigma.fo = ~x, data = dat, family=NBI(), mu.start = mu, sigma.start = sigma)
  }
  # likelihood ratio statistics
  stat <- 2*(logLik(m1)-logLik(m0))[[1]]
  # degrees of freedom
  df <- m1$df.fit-m0$df.fit
  # chi squared test
  p <- pchisq(stat, df, lower=F)

  # if p < 0.05, H0 hypothesis is rejected and H1 hypothesis is accepted: res <- 1
  # if p >= 0.05, H0 hypothesis is accepted and H1 hypothesis is rejected: res <- 0
  if(p < .05) {res <- 1} else {res <- 0}

  res
}

# Specify sample sizes
N_search <- c(5,10,15,20,25,30,40,50,100, 1000)
# Specify fold changes for non-Posson noise
d_search <- c(1/4, 1/3, 1/2, 2/3, 4/5, 5/4, 3/2, 2, 3, 4)

search <- expand.grid(params=list(N= N_search, d=d_search), KEEP.OUT.ATTRS = FALSE)
search <- search[order(search$d),]

# Specify the number of simulation runs
simNum <- 1000

# Specify offset
ofs <- c(0, log(1e06), log(1e08))
# Specify average counts
mu <- c(10, 100, 1000)
# Specify non-Posson noise
bcv <- c(0.1, 0.25, 0.5)

# power analysis for non-Poisson noise at different settings (varing offset, mu and sigma)
power_LR_sigma <- lapply(seq_along(ofs), function(O)
lapply(seq_along(mu), function(M)
  lapply(seq_along(bcv), function(CV) {
    cat("\n", paste(ofs[O], mu[M], bcv[CV], sep=";"), "\n")

    power_LR <- mclapply(c(1:nrow(search)), function(i) {
      cat(paste(search[i,], sep=";"), "\n")

      res <- lapply(c(1:simNum), function(j) {
          # cat(j, "\r")
          sink(file="/dev/null")
          res <- tryCatch(NB_func_sigma(N=search$N[i], d=search$d[i], mu=mu[M], sigma=bcv[CV], offset = TRUE, ofs = ofs[O]),
                   error = function(e) NULL, warning = function(w) NULL)
          sink()
          res
      })
      
      res <- unlist(res)
      res <- mean(res)
      res
    }, mc.cores=4)

    power_LR <- matrix(unlist(power_LR), ncol=10)
    rownames(power_LR) <- N_search
    colnames(power_LR) <- c("1/4", "1/3", "1/2", "2/3", "4/5", "5/4", "3/2", "2", "3", "4")

    res <- power_LR
  })
))

# Store information for further use
names(power_LR_sigma) <- paste("ofs:", c("NO", "1e06", "1e07"))
power_LR_sigma <- lapply(power_LR_sigma, function(x) { names(x) <- paste("mu:", mu); return(x) })
power_LR_sigma <- lapply(power_LR_sigma, function(x)
  lapply(x, function(x1) { names(x1) <- paste("bcv:", bcv); return(x1) })
)

#setwd("/media/Data/Downloads/PG_Review/")
save.image("power_LR_sigma.Rdat")


################################################################################
library(gamlss)
library(parallel)

# Function calculating likelihood ratio test for varying counts in Negative Binomial (NB) model
# N - sample size
# d - fold change for non-Poisson noise
# mu - average counts
# sigma - non-Poisson noise (coefficient of variation in counts - cv(X))
# offset - flag indicating whether to include (TRUE) or not (FALSE) offset into the models
# ofs - offset value

NB_func_mu <- function(N, d, mu, sigma, offset = FALSE, ofs = log(1e07) ) {
  y1 <- rNBI(N, mu, sigma^2)
  y2 <- rNBI(N, mu*d, sigma^2)
  x <- c(rep("a", N), rep("b", N))
  y <- c(y1, y2)
  dat <- data.frame(x=x, y=y, ofs = ofs)

  if(offset) {
    # m0 - reduced model, which does not account for the differences in average expression. The model also includes offset.
    m0 <- gamlss(y ~     offset(ofs), data = dat, family=NBI(), mu.start = mu, sigma.start = sigma)
    # m1 - full model, which does  account for the differences in averafe expression. The model also includes offset.
    m1 <- gamlss(y ~ x + offset(ofs), data = dat, family=NBI(), mu.start = mu, sigma.start = sigma)
  } else {
    # m0 - reduced model, which does  account for the differences in average expression. The model does not include offset.
    m0 <- gamlss(y ~ 1, data = dat, family=NBI(), mu.start = mu, sigma.start = sigma)
    # m1 - full model, which does  account for the differences in average expression. The model does not include offset.
    m1 <- gamlss(y ~ x, data = dat, family=NBI(), mu.start = mu, sigma.start = sigma)
  }
  # likelihood ratio statistics
  stat <- 2*(logLik(m1)-logLik(m0))[[1]]
  # degrees of freedom
  df <- m1$df.fit-m0$df.fit
  # chi squared test
  p <- pchisq(stat, df, lower=F)
  
  # if p < 0.05, H0 hypothesis is rejected and H1 hypothesis is accepted: res <- 1
  # if p >= 0.05, H0 hypothesis is accepted and H1 hypothesis is rejected: res <- 0
  
  if(p < .05) {res <- 1} else {res <- 0}

  res
}

# Specify sample sizes
N_search <- c(5,10,15,20,25,30,40,50,100, 1000)
# Specify fold changes for average expression
d_search <- c(1/4, 1/3, 1/2, 2/3, 4/5, 5/4, 3/2, 2, 3, 4)

search <- expand.grid(params=list(N= N_search, d=d_search), KEEP.OUT.ATTRS = FALSE)
search <- search[order(search$d),]

# Specify the number of simulation runs
simNum <- 1000

#specify offset
ofs <- c(0, log(1e06), log(1e08))
# specify average counts
mu <- c(10, 100, 1000)
# specify non-Poisson noise
bcv <- c(0.1, 0.25, 0.5)

# Power analysis for average expression with varying offset, mu and sigma
power_LR_mu <- lapply(seq_along(ofs), function(O)
lapply(seq_along(mu), function(M)
  lapply(seq_along(bcv), function(CV) {
    cat("\n", paste(ofs[O], mu[M], bcv[CV], sep=";"), "\n")

    power_LR <- mclapply(c(1:nrow(search)), function(i) {
      cat(paste(search[i,], sep=";"), "\n")

      res <- lapply(c(1:simNum), function(j) {
          # cat(j, "\r")
          sink(file="/dev/null")
          res <- tryCatch(NB_func_mu(N=search$N[i], d=search$d[i], mu=mu[M], sigma=bcv[CV], offset = TRUE, ofs = ofs[O]),
                   error = function(e) NULL, warning = function(w) NULL)
          sink()
          res
      })

      res <- unlist(res)
      res <- mean(res)
      res
    }, mc.cores=4)

    power_LR <- matrix(unlist(power_LR), ncol=10)
    rownames(power_LR) <- N_search
    colnames(power_LR) <- c("1/4", "1/3", "1/2", "2/3", "4/5", "5/4", "3/2", "2", "3", "4")

    res <- power_LR
  })
))

# Store information for further use
names(power_LR_mu) <- paste("ofs:", ofs)
power_LR_mu <- lapply(power_LR_mu, function(x) { names(x) <- paste("mu:", mu); return(x) })
power_LR_mu <- lapply(power_LR_mu, function(x)
  lapply(x, function(x1) { names(x1) <- paste("bcv:", bcv); return(x1) })
)

#setwd("/media/Data/Downloads/PG_Review/")
save.image("power_LR_mu.Rdat")

################################################################################
################################################################################
################################################################################
# Plotting power analysis results:

library(ggplot2)
library("egg")
library(reshape2)


#setwd("/media/Data/Downloads/PG_Review/")               # Set directory
load("power_LR_mu.Rdat")                                     # Load previously generated averages
load("power_LR_sigma.Rdat")                                  # Load previously generated sigmas


# Plot the results of the power-analysis
plot_power <- function(power_list, lowess_f = 0.1, rnd = 1) {

p <- lapply(power_list, function(x){
  lx <- apply(apply(x, 2, function(x) lowess(x, f=lowess_f)$y), 1, function(x) lowess(x, f=lowess_f)$y)
  lx <- t(lx)
  lx[lx > 1] <- 1
  lx[lx < 0] <- 0

  colnames(lx) <- colnames(x)
  rownames(lx) <- rownames(x)

  lx <- melt(lx)
  colnames(lx) <- c("N", "FC", "power")
  lx$N <- as.factor(lx$N)

  lx_round <- lx
  lx_round$power <- round(lx_round$power, rnd)

  res <- ggplot(data = lx, aes(x=N, y=FC, fill=power)) +
    geom_tile() + theme_bw() +
    scale_fill_gradient2(low = "#FF4444", high = "red", limit = c(0,1)) +
    geom_text(data = lx_round, aes(x=N, y=FC, label=power), color = "black", size = 3)

  res
})

res <- grid.arrange(grobs = p, nrow=3, ncol=3  )

}


## Write power analysis to PDFs
ggsave("power_LR_sigma_ofs0.pdf", plot_power(Reduce("c", power_LR_sigma[[1]]), lowess_f = 1), width=12, height= 9 )
ggsave("power_LR_sigma_ofs6.pdf", plot_power(Reduce("c", power_LR_sigma[[2]]), lowess_f = 1), width=12, height= 9 )
ggsave("power_LR_sigma_ofs8.pdf", plot_power(Reduce("c", power_LR_sigma[[3]]), lowess_f = 1), width=12, height= 9 )

ggsave("power_LR_mu_ofs0.pdf", plot_power(Reduce("c", power_LR_mu[[1]]), lowess_f = 1e-03), width=12, height= 9 )
ggsave("power_LR_mu_ofs6.pdf", plot_power(Reduce("c", power_LR_mu[[2]]), lowess_f = 1e-03), width=12, height= 9 )
ggsave("power_LR_mu_ofs8.pdf", plot_power(Reduce("c", power_LR_mu[[3]]), lowess_f = 1e-03), width=12, height= 9 )

################################################################################
################################################################################
################################################################################
