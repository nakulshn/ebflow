source("../alternative.R")
library(argparse)

parser <- ArgumentParser(description='Simulation settings')
parser$add_argument("--D", type = "character")
parser$add_argument("--P", type = "character")
parser$add_argument("--dimr", type = "double")
parser$add_argument("--noiser", type = "double")
parser$add_argument("--T", type = "integer", help = "inner iters")
parser$add_argument("--lambda", type = "double")
parser$add_argument("--seed", type = "integer", help = "random seed")
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
y = dat$y

burn.iter = 200
max.iter = 10000
save.iter = 100

res_file_name = sprintf("Gibbs_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f_T%d_seed%d.rds",args$P,args$D,args$dimr,args$noiser,args$lambda,args$T,args$seed)

if (file.exists(paste0(respath,res_file_name))) { quit() }

set.seed(args$seed)
t1 = Sys.time()
print(sprintf("RUN Gibbs: %s", res_file_name))
out = Gibbs(X, y, sigma, Theta, lambda=args$lambda, eta.w = 1.0, gibbs.inner.iters = args$T,
                  burn.iters = burn.iter, w.true=w.true,
		  max.iter=max.iter, save.iter=save.iter, verbose=TRUE)
print("FINISH")
if (args$predict) {
  print(sprintf("Computing posterior mean"))
  w = out$w
  phi = out$hist[[max.iter+1]]$phi
  tmp.out = Gibbs(X, y, sigma, Theta, lambda=0, eta.w = 1.0, gibbs.inner.iters = 20000,
                    w.init = w, phi.init = phi, burn.iters = 0, w.true=w.true,
  		  max.iter=50000, save.iter=50, verbose=TRUE)
  theta.mean = rowMeans(sapply(seq(1,50001,50), function(it) tmp.out$hist[[it]]$phi))
  out[["theta.mean"]] = theta.mean
}
t2 = Sys.time()
print(t2-t1)
saveRDS(out, file=paste0(respath,res_file_name))

