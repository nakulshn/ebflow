## run_PolyaTree.R
##
## Skeleton runner for Polya-tree method, designed to mirror run_EBflow.R.
## Uses the same data files produced by gen_data.R and similar CLI arguments,
## but currently calls the simplified polya_tree_skeleton() routine.

source("../PolyaTree.R")
library(argparse)

parser <- ArgumentParser(description = "Simulation settings for Polya-tree method")

## Keep argument names consistent with run_EBflow.R so you can reuse job scripts.
parser$add_argument("--D",          type = "character")   # design name
parser$add_argument("--P",          type = "character")   # prior name
parser$add_argument("--dimr",       type = "double")      # n/p ratio
parser$add_argument("--noiser",     type = "double")      # inverse SNR
parser$add_argument("--seed",       type = "integer")
parser$add_argument("--respath",    type = "character")

args = parser$parse_args()

respath = args$respath
datapath = "/home/ubuntu/ebflow/data/"

###LOAD data
X_file_name = sprintf("X_%s_nbyp%.1f.rds",args$D,args$dimr)
dat = readRDS(paste0(datapath,X_file_name))
X = dat$X.train
y_file_name = sprintf("y_%s_%s_nbyp%.1f_invSNR%.1f.rds",args$P,args$D,args$dimr,args$noiser)
dat = readRDS(paste0(datapath,y_file_name))
grid = dat$w.grid
w.true = dat$w.true
sigma = dat$sigma
theta = dat$theta.true
y = dat$y

max.iter = 10000
save.iter = 1


## --- Set Polya-tree depth L, default from K = 2^L + 1 ---
K = 65
L = 6



res_file_name = sprintf("PolyaTree_%s_%s_nbyp%.1f_invSNR%.1f_seed%d.rds",args$P,args$D,args$dimr,args$noiser,args$seed)

if (file.exists(paste0(respath,res_file_name))) { quit() }

set.seed(args$seed)
t1 = Sys.time()
print(sprintf("RUN Polya Tree: %s", res_file_name))
pt_out = polya_tree_skeleton(X, y, sigma, grid, K=K, L=L,
		save.iter=save.iter, max.iter=max.iter, verbose=TRUE, w.true=w.true, theta.true=theta, print.iter=100)

print("FINISH Polya tree")
t2=Sys.time()
print(t2-t1)
saveRDS(pt_out, file=paste0(respath,res_file_name))