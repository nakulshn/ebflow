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

burn.iter = 200
max.iter = 10000
save.iter = 100

if (args$eta.phi == "decay") {
  eta.phi = logspace(0,-1,max.iter-burn.iter)
} else {
  eta.phi = rep(as.numeric(args$eta.phi),max.iter-burn.iter)
}
eta.phi = c(rep(1,burn.iter),eta.phi)

if (args$eta.w == "decay") {
  eta.w = logspace(-2,-3,max.iter-burn.iter)
} else {
  eta.w = rep(as.numeric(args$eta.w),max.iter-burn.iter)
}
eta.w = c(rep(0,burn.iter),eta.w)

if (args$precondition) {
  res_file_name = sprintf("EBflow_preconditioned_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f_etaphi%s_etaw%s_seed%d.rds",args$P,args$D,args$dimr,args$noiser,args$lambda,args$eta.phi,args$eta.w,args$seed)
} else {
  res_file_name = sprintf("EBflow_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f_etaphi%s_etaw%s_seed%d.rds",args$P,args$D,args$dimr,args$noiser,args$lambda,args$eta.phi,args$eta.w,args$seed)
}

if (file.exists(paste0(respath,res_file_name))) { quit() }

set.seed(args$seed)
t1 = Sys.time()
print(sprintf("RUN EBflow: %s", res_file_name))
ebflow_out = gradient.flow.EB(X, y, sigma, Theta, lambda=args$lambda, eta.w=eta.w,
		eta.phi=eta.phi, tausq.scale=0.5, precondition=args$precondition,
		save.iter=save.iter, max.iter=max.iter, verbose=TRUE, w.true=w.true)
print("FINISH EBflow")
if (args$predict) {
  print(sprintf("Computing posterior mean"))
  w = ebflow_out$w
  phi = ebflow_out$hist[[max.iter+1]]$phi
  tmp.out = gradient.flow.EB(X, y, sigma, Theta, lambda=0, eta.w=0, w.init=w, phi.init=phi,
  		eta.phi=eta.phi[max.iter], tausq.scale=0.5, precondition=args$precondition,
  		save.iter=50, max.iter=50000, verbose=TRUE, w.true=w.true)
  phi = sapply(seq(1,50001,50), function(it) tmp.out$hist[[it]]$phi)
  theta.mean = posterior.mean(w,phi,Theta,tmp.out$tau)
  ebflow_out[["theta.mean"]] = theta.mean
}
t2=Sys.time()
print(t2-t1)
saveRDS(ebflow_out, file=paste0(respath,res_file_name))

