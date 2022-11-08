# This script is used to remove the batch effect of the gene expression data
library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
# library(limma)
# setting the data: feature as columns, row as samples
pheno = pData(bladderEset)
edata = exprs(bladderEset)
# using model.matrix to generate object
mod = model.matrix(~as.factor(cancer), data=pheno) # cancer is the variable of interest, full model
mod0 = model.matrix(~1, data=pheno) # ~1 means no adjusted variables

# The sva function performs two different steps
## Firstly it identifies the number of latent factors that need to be estimated.
n.sv = num.sv(edata,mod,method="leek")

## Secondly we apply the sva function to estimate the surrogate variables:
svobj = sva(edata,mod,mod0,n.sv=n.sv) # here the full model and the adjusted model must be different

## Combat

# Applying the ComBat function to adjust for known batches
# The ComBat function adjusts for known batches. In order to use the function, you must have a known batch variable in your dataset.
# https://rdrr.io/bioc/sva/man/ComBat.html

pheno = pData(dat)
edata = exprs(dat)
batch = pheno$batch
mod = model.matrix(~as.factor(cancer), data=pheno)

# parametric adjustment
combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
