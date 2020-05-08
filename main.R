##############################################################################
# Project: Generalized Spatially Varying Coefficient Models
##############################################################################
rm(list = ls())

# library
# install.packages("devtools",repos="http://cran.r-project.org")
# library(devtools)
# install_github("funstatpackages/BPST")
# install_github("funstatpackages/gsvcm")

library('BPST') # basis matrix
library('mgcv') # horse shoe domain
library('MGLM') # kr function for design matrix
library('MASS') # generate negbinom response
library('gsvcm') # gsvcm

# choose family
family=poisson()

# Load matrices for triangles (Tr) and vertices (V)
# Select a pair of triangulations (Tr and V):
Tr = gsvcm::Tr0_horse
V = gsvcm::V0_horse

# Setup for simulation
n = 2000
nsim = 100
d = 2; r = 1

# Generate population data
ngrid = 0.05                                  # distance between grid points
all_pop = as.matrix(Datagenerator(family, ngrid)) # warnings() for NAs


# This contains population data only on the target domain.
pop.r=all_pop[!is.na(all_pop[,'m1']),]
N=nrow(pop.r)

# set up for smoothing parameters in the penalty term
lambda_start=0.0001
lambda_end=10
nlambda=10
lambda=exp(seq(log(lambda_start),log(lambda_end),length.out=nlambda))

# sample
set.seed(2020)
ind.s=sample(N,n,replace=FALSE)
data=as.matrix(pop.r[ind.s,])
y=data[,1]; beta0=data[,c(2:3)]; X=data[,c(4:5)]; S=data[,c(6:7)]

# fit the model:
mfit0 = fit.gsvcm(y, X, S, V, Tr, d, r, lambda, family)
y = mfit0$y; S = mfit0$S; X = mfit0$X
y_hat = predict(mfit0, X, S)

# k-fold cross-validation:
MSPE = cv.gsvcm(y, X, S, V, Tr, d = d, r = r, lambda, family)
mean(MSPE)

# GQLR test:
test_result =test.gsvcm(y, X, S, V, Tr, d, r, test_iter = 1, family = family, nB = 50)

# plot coefficients
plot(mfit0)

