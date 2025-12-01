compute.loglik = function(w,y,Theta,sigma) {
  n = length(y); y = c(y); K = length(Theta)
  mat = exp(-(Theta%o%rep(1,n)-rep(1,K)%o%y)^2/(2*sigma^2))
  return(mean(colSums(w*mat)))
}

print.table = function(tabs,design) {
  cat("\\hline\n")
  if (design == "identity") { tag = "Log-likelihood"; tag.2 = "Log-likelihood (unsmoothed)" }
  else { tag = "Prediction MSE"; tag.2 = "Prediction MSE (unsmoothed)" }
  cat(sprintf("\\multirow{6}{*}{%s} & %s & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n",
			  design,"EBflow",
			  tabs[["gaussian"]][["EBflow",tag]],tabs[["gaussian"]][["EBflow",tag.2]],
			  tabs[["skew"]][["EBflow",tag]],tabs[["skew"]][["EBflow",tag.2]],
			  tabs[["cauchy"]][["EBflow",tag]],tabs[["cauchy"]][["EBflow",tag.2]],
			  tabs[["bimodal"]][["EBflow",tag]],tabs[["bimodal"]][["EBflow",tag.2]]))
  for (method in c("EBflow.preconditioned", "Langevin.eta1.0", "Langevin.eta0.1", "Gibbs", "CAVI")) {
    cat(sprintf("& %s & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n",
  			  method, tabs[["gaussian"]][[method,tag]],tabs[["gaussian"]][[method,tag.2]],
  			  tabs[["skew"]][[method,tag]],tabs[["skew"]][[method,tag.2]],
  			  tabs[["cauchy"]][[method,tag]],tabs[["cauchy"]][[method,tag.2]],
  			  tabs[["bimodal"]][[method,tag]],tabs[["bimodal"]][[method,tag.2]]))
  
  }
}

methods = c("EBflow", "EBflow.preconditioned", "Langevin.eta1.0", "Langevin.eta0.1", "Gibbs", "CAVI")

outf = sprintf("tables/unsmoothed.tab")
sink(outf)
cat("Design & Method & \\multicolumn{2}{|c|}{gaussian} & \\multicolumn{2}{|c|}{skew} & \\multicolumn{2}{|c|}{cauchy} & \\multicolumn{2}{|c}{bimodal} \\\\\n")
cat("\\hline\n")
for (design in c("identity","iid","block02corr0.9","block10corr0.5")) {
  tabs = list()
  if (design == "identity") { cat("&& \\multicolumn{2}{|c|}{Log-lik} & \\multicolumn{2}{|c|}{Log-lik} & \\multicolumn{2}{|c|}{Log-lik} & \\multicolumn{2}{|c}{Log-lik} \\\\\n") }
  else if (design == "iid") { cat("\\hline\n"); cat("&& \\multicolumn{2}{|c|}{MSE} & \\multicolumn{2}{|c|}{MSE} & \\multicolumn{2}{|c|}{MSE} & \\multicolumn{2}{|c}{MSE} \\\\\n") }
  for (prior in c("gaussian","skew","cauchy","bimodal")) {
    if (prior == "gaussian") { noiser = 0.5; lambda = 0.003 }
    else { noiser = 0.8; lambda = 0.001 }
    fnames = c(sprintf("EBflow_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f_etaphi%s_etaw%s",prior,design,1.0,noiser,lambda,"decay","decay"),
    	           sprintf("EBflow_preconditioned_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f_etaphi%s_etaw%s",prior,design,1.0,noiser,lambda,"decay","decay"),
                 sprintf("Langevin_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f_etaphi%.1f_T%d",prior,design,1.0,noiser,lambda,1.0,100),
                 sprintf("Langevin_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f_etaphi%.1f_T%d",prior,design,1.0,noiser,lambda,0.1,100),
                 sprintf("Gibbs_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f_T%d",prior,design,1.0,noiser,lambda,100),
        	   sprintf("VI_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f",prior,design,1.0,noiser,lambda))
    fnames_unsmoothed = c(sprintf("EBflow_%s_%s_nbyp%.1f_invSNR%.1f_lambda0.000_etaphi%s_etaw%s",prior,design,1.0,noiser,"decay","decay"),
                 sprintf("EBflow_preconditioned_%s_%s_nbyp%.1f_invSNR%.1f_lambda0.000_etaphi%s_etaw%s",prior,design,1.0,noiser,"decay","decay"),
                 sprintf("Langevin_%s_%s_nbyp%.1f_invSNR%.1f_lambda0.000_etaphi%.1f_T%d",prior,design,1.0,noiser,1.0,100),
                 sprintf("Langevin_%s_%s_nbyp%.1f_invSNR%.1f_lambda0.000_etaphi%.1f_T%d",prior,design,1.0,noiser,0.1,100),
                 sprintf("Gibbs_%s_%s_nbyp%.1f_invSNR%.1f_lambda0.000_T%d",prior,design,1.0,noiser,100),
        	   sprintf("VI_%s_%s_nbyp%.1f_invSNR%.1f_lambda0.000",prior,design,1.0,noiser))
    dat = readRDS(sprintf("../data/y_%s_%s_nbyp%.1f_invSNR%.1f.rds",prior,design,1.0,noiser))
    Theta = dat$w.grid
    w.true = dat$w.true
    theta.true = dat$theta.true
    y = dat$y
    sigma = dat$sigma
    if (design != "identity") {
      dat = readRDS(sprintf("../data/X_%s_nbyp%.1f.rds",design,1.0))
      X.test = dat$X.test
    }
    tab = data.frame(matrix(ncol=2,nrow=0))
    if (design == "identity") { colnames(tab) = c("Log-likelihood","Log-likelihood (unsmoothed)") }
    else { colnames(tab) = c("Prediction MSE","Prediction MSE (unsmoothed)") }
    for (i in 1:length(methods)) {
      method = methods[i]
      stats = data.frame(matrix(ncol=2,nrow=1))
      colnames(stats) = colnames(tab)
      if (method != "CAVI") { out = readRDS(sprintf("results/%s_seed1.rds",fnames[i])) }
      else { out = readRDS(sprintf("results/%s.rds",fnames[i])) }
      if (design == "identity") { stats[["Log-likelihood"]] = compute.loglik(out$w,y,Theta,sigma) }
      else { stats[["Prediction MSE"]] = sum((X.test%*%(out$theta.mean-theta.true))^2)/sum((X.test%*%theta.true)^2) }
      if (method != "CAVI") { out = readRDS(sprintf("results/%s_seed1.rds",fnames_unsmoothed[i])) }
      else { out = readRDS(sprintf("results/%s.rds",fnames_unsmoothed[i])) }
      if (design == "identity") { stats[["Log-likelihood (unsmoothed)"]] = compute.loglik(out$w,y,Theta,sigma) }
      else { stats[["Prediction MSE (unsmoothed)"]] = sum((X.test%*%(out$theta.mean-theta.true))^2)/sum((X.test%*%theta.true)^2) }
      tab = rbind(tab,stats)
    }
    rownames(tab) = methods
    tabs[[prior]] = tab
  }
  print.table(tabs,design)
}
sink()
