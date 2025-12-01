source("../alternative.R")
library(argparse)

parser <- ArgumentParser(description='Simulation settings')
parser$add_argument("--D", type = "character")
parser$add_argument("--P", type = "character")
parser$add_argument("--dimr", type = "double")
parser$add_argument("--noiser", type = "double")
parser$add_argument("--lambda", type = "double")
parser$add_argument("--respath", type = "character")
parser$add_argument("--predict", action = "store_true")

args = parser$parse_args()

respath = args$respath
datapath = "/home/ubuntu/ebflow/data/"

###LOAD data
X_file_name = sprintf("X_%s_nbyp%.1f.rds",args$D,args$dimr)
dat = readRDS(paste0(datapath,X_file_name))
X = dat$X.train
y_file_name = sprintf("y_%s_%s_nbyp%.1f_invSNR%.1f.rds",args$P,args$D,args$dimr,args$noiser)
dat = readRDS(paste0(datapath,y_file_name))
Theta = dat$w.grid
w.true = dat$w.true
sigma = dat$sigma
theta = dat$theta.true
y = dat$y

max.iter = 1000
save.iter = 1

res_file_name = sprintf("VI_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f.rds",args$P,args$D,args$dimr,args$noiser,args$lambda)

if (file.exists(paste0(respath,res_file_name))) { quit() }

set.seed(args$seed)
t1 = Sys.time()
print(sprintf("RUN VI: %s", res_file_name))
out = VI(X, y, sigma, Theta, lambda=args$lambda, eta.w = 1.0,
                  w.true=w.true, max.iter=max.iter, save.iter=save.iter, verbose=TRUE)
print("FINISH")
if (args$predict) {
  theta.mean = out$hist[[max.iter+1]]$phi
  out[["theta.mean"]] = theta.mean
}
t2 = Sys.time()
print(t2-t1)
saveRDS(out, file=paste0(respath,res_file_name))

