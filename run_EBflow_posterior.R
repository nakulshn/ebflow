source("../EBflow.R")
library(argparse)
library(pracma)

parser <- ArgumentParser(description='Simulation settings')
parser$add_argument("--D", type = "character")
parser$add_argument("--P", type = "character")
parser$add_argument("--dimr", type = "double")
parser$add_argument("--noiser", type = "double")
parser$add_argument("--eta.w", type = "character")
parser$add_argument("--eta.phi", type = "character")
parser$add_argument("--lambda", type = "double")
parser$add_argument("--precondition", action = "store_true")
parser$add_argument("--seed", type = "integer")
parser$add_argument("--respath", type = "character")

args = parser$parse_args()

respath = args$respath
datapath ="/home/ubuntu/ebflow/data/"

###LOAD data
X_file_name = sprintf("X_%s_nbyp%.1f.rds",args$D,args$dimr)
dat = readRDS(paste0(datapath,X_file_name))
X = dat$X.train; Xte = dat$X.test
y_file_name = sprintf("y_%s_%s_nbyp%.1f_invSNR%.1f.rds",args$P,args$D,args$dimr,args$noiser)
dat = readRDS(paste0(datapath,y_file_name))
res_file_name = sprintf("EBflow_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f_etaphi%s_etaw%s_seed%d.rds",args$P,args$D,args$dimr,args$noiser,args$lambda,args$eta.phi,args$eta.w,args$seed)
ebflow_out = readRDS(file=paste0(respath,res_file_name))


Theta = dat$w.grid
w.true = dat$w.true
sigma = dat$sigma
theta = dat$theta.true
y = dat$y

phi = ebflow_out$hist[[max.iter+1]]$phi

w = ebflow_out$w
phi = ebflow_out$hist[[max.iter+1]]$phi

set.seed(args$seed)
t1 = Sys.time()
tmp.out = gradient.flow.EB(X, y, sigma, Theta, lambda=0, eta.w=0, w.init=w, phi.init=phi,
		eta.phi=eta.phi[max.iter], tausq.scale=0.5, precondition=args$precondition,
		save.iter=1, max.iter=10000, verbose=TRUE, w.true=w.true)
phi = sapply(seq(1,10000), function(it) tmp.out$hist[[it]]$phi)
theta.mean = posterior.mean(w,phi,Theta,tmp.out$tau)


saveRDS(ebflow_out, file=paste0(respath,res_file_name))

