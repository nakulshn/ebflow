plot.w = function(fname,prefix,Theta,w.true,col,title,seeded=TRUE) {
  pdf(file=fname,width=3,height=3)
  par(mar = c(3,3,3,3), cex=0.7)
  plot(Theta,w.true,type="l",xlab="",ylab="",main=title)
  if (seeded) {
    for (seed in 1:10) { 
      out = readRDS(sprintf("%s_seed%d.rds",prefix,seed))
      lines(Theta,out$w,col=col,xlab="",ylab="")
    }
  } else {
    out = readRDS(sprintf("%s.rds",prefix))
    lines(Theta,out$w,col=col,xlab="",ylab="")
  }
  dev.off()
}

for (prior in c("gaussian","skew","cauchy","bimodal")) {
  print(prior)
  prior = "bimodal"
  if (prior == "gaussian") { lambda = 0.003; noiser = 0.5 }
  else { lambda = 0.001; noiser = 0.8 }
  dat = readRDS(sprintf("../data/y_%s_block02corr0.9_nbyp1.0_invSNR%.1f.rds",prior,noiser))
  Theta = dat$w.grid
  w.true = dat$w.true
  print("  EBflow")
  title = sprintf("EBflow")
  fname = sprintf("plots/EBflow_%s.pdf",prior)
  prefix = sprintf("results/EBflow_%s_block02corr0.9_nbyp1.0_invSNR%.1f_lambda%.3f_etaphidecay_etawdecay",prior,noiser,lambda)
  plot.w(fname,prefix,Theta,w.true,"red",title)
  print("  EBflow (preconditioned)")
  title = sprintf("EBflow (preconditioned)")
  fname = sprintf("plots/EBflow_preconditioned_%s.pdf",prior)
  prefix = sprintf("results/EBflow_preconditioned_%s_block02corr0.9_nbyp1.0_invSNR%.1f_lambda%.3f_etaphidecay_etawdecay",prior,noiser,lambda)
  plot.w(fname,prefix,Theta,w.true,"red",title)
  print("  Langevin")
  for (eta.phi in c(1.0,0.1)) {
    title = sprintf("Langevin (eta.phi=%.1f)",eta.phi)
    fname = sprintf("plots/Langevin_etaphi%.1f_%s.pdf",eta.phi,prior)
    prefix = sprintf("results/Langevin_%s_block02corr0.9_nbyp1.0_invSNR%.1f_lambda%.3f_etaphi%.1f_T100",prior,noiser,lambda,eta.phi)
    plot.w(fname,prefix,Theta,w.true,"orange",title)
  }
  print("  Gibbs")
  title = sprintf("Gibbs MCEM")
  fname = sprintf("plots/Gibbs_%s.pdf",prior)
  prefix = sprintf("results/Gibbs_%s_block02corr0.9_nbyp1.0_invSNR%.1f_lambda%.3f_T100",prior,noiser,lambda)
  plot.w(fname,prefix,Theta,w.true,"blue",title)
  print("  VI")
  title = sprintf("CAVI")
  fname = sprintf("plots/VI_%s.pdf",prior)
  prefix = sprintf("results/VI_%s_block02corr0.9_nbyp1.0_invSNR%.1f_lambda%.3f",prior,noiser,lambda)
  plot.w(fname,prefix,Theta,w.true,"green",title,seeded=FALSE)
  print("  PolyaTree")
  title = springf("Polya Tree")
  fname = sprintf("plots/PolyaTree_%s.pdf", prior)
  prefix = sprintf("results/PolyaTree_%s_block02corr0.9_nbyp1.0_invSNR%.1f", prior, noiser)
  plot.w(fname, prefix, Theta, w.true, "green", title, seeded=FALSE)

}

