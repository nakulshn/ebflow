compute.loglik = function(w,y,Theta,sigma) {
  n = length(y); y = c(y); K = length(Theta)
  mat = exp(-(Theta%o%rep(1,n)-rep(1,K)%o%y)^2/(2*sigma^2))
  return(mean(colSums(w*mat)))
}

methods = c("EBflow", "EBflow.preconditioned", "Langevin.eta1.0", "Langevin.eta0.1", "Gibbs", "CAVI", "PolyaTree")

print.table = function(tab,design) {
  cat("\\hline\n")
  if (design == "identity") { tag = "Log-likelihood" }
  else { tag = "Prediction MSE" }
  cat(sprintf("\\multirow{6}{*}{%s} & \\multirow{6}{*}{%.1f} & %s & %.3f & %.4f & %d & %.3f\\\\\n",
      		design,dimr,"EBflow",tab[["EBflow","TV"]],tab[["EBflow","TV.sd"]],
      		tab[["EBflow","Time (iters)"]],tab[["EBflow",tag]]))
  for (method in c("EBflow.preconditioned", "Langevin.eta1.0", "Langevin.eta0.1", "Gibbs", "CAVI", "PolyaTree")) {
    cat(sprintf("&& %s & %.3f & %.4f & %d & %.3f\\\\\n",
      		  method,tab[[method,"TV"]],tab[[method,"TV.sd"]],
      		  tab[[method,"Time (iters)"]],tab[[method,tag]]))
  }
}

# for (prior in c("gaussian","skew","cauchy","bimodal")) {
for (prior in c("bimodal")) {
  outf = sprintf("tables/%s.tab",prior)
  sink(outf)
  # for (design in c("identity","iid","block02corr0.9","block10corr0.5")) {
  for (design in c("block02corr0.9")) {
    # for (dimr in c(2.0,1.0,0.5)) {
    for (dimr in c(1.0)) {
      if (design == "identity" && dimr != 1.0) { next }
      if (design == "identity") { cat("Design & $n/p$ & Method & TV & TV.sd & Time (iters) & Log-likelihood \\\\\n") }
      else if (design == "iid" && dimr == 2.0) { cat("\\hline\nDesign & $n/p$ & Method & TV & TV.sd & Time (iters) & Prediction MSE \\\\\n") }
      if (prior == "gaussian") { noiser = 0.5; lambda = 0.003 }
      else { noiser = 0.8; lambda = 0.001 }
      fnames = c(sprintf("EBflow_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f_etaphi%s_etaw%s",prior,design,dimr,noiser,lambda,"decay","decay"),
      	           sprintf("EBflow_preconditioned_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f_etaphi%s_etaw%s",prior,design,dimr,noiser,lambda,"decay","decay"),
                   sprintf("Langevin_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f_etaphi%.1f_T%d",prior,design,dimr,noiser,lambda,1.0,100),
                   sprintf("Langevin_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f_etaphi%.1f_T%d",prior,design,dimr,noiser,lambda,0.1,100),
                   sprintf("Gibbs_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f_T%d",prior,design,dimr,noiser,lambda,100),
          	   sprintf("VI_%s_%s_nbyp%.1f_invSNR%.1f_lambda%.3f",prior,design,dimr,noiser,lambda),
          	   sprintf("PolyaTree_%s_%s_nbyp%.1f_invSNR%.1f",prior,design,dimr,noiser))
      dat = readRDS(sprintf("../data/y_%s_%s_nbyp%.1f_invSNR%.1f.rds",prior,design,dimr,noiser))
      Theta = dat$w.grid
      w.true = dat$w.true
      theta.true = dat$theta.true
      y = dat$y
      sigma = dat$sigma
      if (design != "identity") {
        dat = readRDS(sprintf("../data/X_%s_nbyp%.1f.rds",design,dimr))
        X.test = dat$X.test
      }
      tab = data.frame(matrix(ncol=4,nrow=0))
      if (design == "identity") { colnames(tab) = c("TV", "TV.sd", "Time (iters)", "Log-likelihood") }
      else { colnames(tab) = c("TV", "TV.sd", "Time (iters)", "Prediction MSE") }
      for (i in 1:length(methods)) {
        method = methods[i]
        if (method != "CAVI") {
          stats = data.frame(matrix(ncol=5,nrow=10))
          colnames(stats) = colnames(tab)
          for (seed in 1:10) {
            if (method != "PolyaTree") {
              out = readRDS(sprintf("results/%s_seed%d.rds",fnames[i],seed))
            } else {
              out = readRDS(sprintf("results_VIPT/%s_seed%d.rds",fnames[i],seed))
            }
            max.iter = 10000
            stats[[seed,"TV"]] = sum(abs(out$w-w.true)/2)


            #if method is Polya Tree, then we don't have out.hist
            #Instead we have out.w_samples
            #and w_guess at iteration it is computed as the mean of w_samples from burn_in to iter
            if (method == "PolyaTree") {
              w_samples_train = out$w_samples
              burn_in = 200
              TV.hist = sapply(seq(1,max.iter+1,100), function(it) sum(abs(colMeans(w_samples_train[(burn_in + 1):it, , drop = FALSE])-w.true)/2))
            } else {
              TV.hist = sapply(seq(1,max.iter+1,100), function(it) sum(abs(out$hist[[it]]$w-w.true)/2))
            }

            if (min(TV.hist) < 0.2) { iter = (min(which(TV.hist < 0.2))-1)*100 }
	    else { iter = max.iter }
            stats[[seed,"Time (iters)"]] = iter
            #time = out$hist[[iter+1]]$time; units(time) = "secs"
            #stats[[seed,"Time (seconds)"]] = as.numeric(time)
            if (design == "identity") { stats[[seed,"Log-likelihood"]] = compute.loglik(out$w,y,Theta,sigma) }
            else { stats[[seed,"Prediction MSE"]] = sum((X.test%*%(out$theta.mean-theta.true))^2)/sum((X.test%*%theta.true)^2) }
          }
        } else {
          stats = data.frame(matrix(ncol=5,nrow=1))
          colnames(stats) = colnames(tab)
          out = readRDS(sprintf("results/%s.rds",fnames[i]))
          max.iter = 1000
          stats[["TV"]] = sum(abs(out$w-w.true)/2)
          TV.hist = sapply(1:(max.iter+1), function(it) sum(abs(out$hist[[it]]$w-w.true)/2))
          if (min(TV.hist) < 0.2) { iter = min(which(TV.hist < 0.2))-1 }
          else { iter = max.iter }
          stats[["Time (iters)"]] = iter
          #time = out$hist[[iter+1]]$time; units(time) = "secs"
          #stats[["Time (seconds)"]] = time
          if (design == "identity") { stats[["Log-likelihood"]] = compute.loglik(out$w,y,Theta,sigma) }
          else { stats[["Prediction MSE"]] = sum((X.test%*%(out$theta.mean-theta.true))^2)/sum((X.test%*%theta.true)^2) }
        }
        stats.sum = data.frame(matrix(ncol=5,nrow=1))
        colnames(stats.sum) = colnames(tab)
        stats.sum[["TV"]] = mean(stats[["TV"]])
        if (method != "VI") { stats.sum[["TV.sd"]] = sd(stats[["TV"]]) }
        else { stats.sum[["TV.sd"]] = NA }
        stats.sum[["Time (iters)"]] = median(stats[["Time (iters)"]])
        if ((stats.sum[["Time (iters)"]] == 10000) || (stats.sum[["Time (iters)"]] == 1000 && method == "CAVI")) { stats.sum[["Time (iters)"]] = NA }
        #stats.sum[["Time (seconds)"]] = median(stats[["Time (seconds)"]])
        if (design == "identity") { stats.sum[["Log-likelihood"]] = mean(stats[["Log-likelihood"]]) }
        else { stats.sum[["Prediction MSE"]] = mean(stats[["Prediction MSE"]]) }
        tab = rbind(tab,stats.sum)
      }
      rownames(tab) = methods
      print.table(tab,design)
    }
  }
}

