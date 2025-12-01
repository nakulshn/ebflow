max.iter = 10000
plot.iter = 100

plot.TV = function(fname,prefix,w.true,iters,col,title) {
  pdf(file=fname,width=3,height=3)
  par(mar = c(4,4,4,4), cex=0.7)
  plot(iters,rep(0,length(iters)),type="n",xlim=c(0,max.iter),ylim=c(0,0.5),main=title,xlab="Iterations",ylab="TV")
  for (seed in 1:10) { 
    result.fname = sprintf("%s_seed%d.rds",prefix,seed)
    out = readRDS(result.fname)
    TV = sapply(iters,function(it) sum(abs(out$hist[[it]]$w-w.true)/2))
    lines(iters,TV,col=col,type="l")
  }
  dev.off()
}

prior = "gaussian"
lambda = 0.003
noiser = 0.5
dat = readRDS(sprintf("../data/y_%s_iid_nbyp0.5_invSNR%.1f.rds",prior,noiser))
w.true = dat$w.true
iters = seq(1,max.iter+1,by=plot.iter)
print("  EBflow")
for (etaphi in c("1.0","0.1")) {
  for (etaw in c("0.01","0.001")) {
    title = sprintf("EBflow, eta.phi=%s, eta.w=%s",etaphi,etaw)
    fname = sprintf("plots/EBflow_%s_etaphi%s_etaw%s.pdf",prior,etaphi,etaw)
    prefix = sprintf("results/EBflow_%s_iid_nbyp0.5_invSNR%.1f_lambda%.3f_etaphi%s_etaw%s",prior,noiser,lambda,etaphi,etaw)
    plot.TV(fname,prefix,w.true,iters,"red",title)
  }
}
title = sprintf("EBflow, eta.phi=decay, eta.w=decay")
fname = sprintf("plots/EBflow_%s_etaphidecay_etawdecay.pdf",prior)
prefix = sprintf("results/EBflow_%s_iid_nbyp0.5_invSNR%.1f_lambda%.3f_etaphidecay_etawdecay",prior,noiser,lambda)
plot.TV(fname,prefix,w.true,iters,"red",title)
print("  Langevin")
for (etaphi in c("1.0","0.1")) {
  for (T in c("100","1000")) {
    title = sprintf("Langevin MCEM, eta.phi=%s, T=%s",etaphi,T)
    fname = sprintf("plots/Langevin_%s_etaphi%s_T%s.pdf",prior,etaphi,T)
    prefix = sprintf("results/Langevin_%s_iid_nbyp0.5_invSNR%.1f_lambda%.3f_etaphi%s_T%s",prior,noiser,lambda,etaphi,T)
    plot.TV(fname,prefix,w.true,iters,"orange",title)
  }
}
print("  Gibbs")
for (T in c("10","100","1000")) {
  title = sprintf("Gibbs MCEM, T=%s",T)
  fname = sprintf("plots/Gibbs_%s_T%s.pdf",prior,T)
  prefix = sprintf("results/Gibbs_%s_iid_nbyp0.5_invSNR%.1f_lambda%.3f_T%s",prior,noiser,lambda,T)
  plot.TV(fname,prefix,w.true,iters,"blue",title)
}
