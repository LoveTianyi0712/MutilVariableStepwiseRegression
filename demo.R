rm(list = ls())
source("MultStepwiseReg.R")

X <- read.table("datax.txt")
Y <- read.table("datay.txt")

# randata <- matrix(rnorm(100000, 1, 1), 500, 200)
# datas$Y <- randata[,1:20]
# datas$X <- randata[,21:200]

datas <- list(X = X, Y = Y)
params <- list(
  ain = 0.10,
  aout = 0.10, 
  standardize = 1, 
  verbose = 2, 
  max_iter = 100)

result = MultStepwiseReg(datas = datas, params = params)
