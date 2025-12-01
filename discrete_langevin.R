library(pracma)
library(CVXR)

# compute the function [log N(0,tau^2)*g]'(phi)
dlog.conv <- function(centered, log.density, w, tau) {
  log.density.w = log.density + log(w)
  col.max = apply(log.density.w,2,max)
  density.w = t(exp(t(log.density.w)-col.max))
  num = colSums(-centered * density.w / tau^2)
  denom = colSums(density.w)
  return(num/denom)
}

# iteration on phi
langevin.step <- function(phi, centered, log.density, w, tau, XSinvX, XSinvy, eta.phi) {
  phi = phi - eta.phi * (XSinvX %*% phi - XSinvy)
  phi = phi + eta.phi * dlog.conv(centered, log.density, w, tau)
  phi = phi + sqrt(2*eta.phi) * rnorm(length(phi))
  return(c(phi))
}

# discrete difference approximation of second derivative
discrete.der = function(K,delta) {
  D = matrix(0,K-2,K)
  for (i in 1:(K-2)) {
    D[i,i] = 1/delta^2
    D[i,i+1] = -2/delta^2
    D[i,i+2] = 1/delta^2
  }
  return(D)
}

update_w = function(phi.save,Theta,tau,M) {
  K = length(Theta)
  p.tot = length(phi.save)
  centered = rep(1,K) %o% phi.save - Theta %o% rep(1,p.tot)
  log.density = dnorm(centered, sd=tau, log=TRUE)
  col.max = apply(log.density,2,max)
  density = t(exp(t(log.density)-col.max))
  w = Variable(K)
  obj = -mean(log(sum_entries(diag(w)%*%density,axis=2))) + sum((M%*%w)^2)
  const = list(sum(w) == 1, w>=0)
  prob = Problem(Minimize(obj),const)
  result = solve(prob)
  w.new = as.numeric(result$getValue(w))
  return(w.new)
}

langevin.EB.init <- function(X, y, sigma, Theta, eta.phi, max.iter, tausq.scale, w.init, phi.init) {
  n = nrow(X); p = ncol(X); K = length(Theta)
  # Set tau
  X.svd = svd(X,nv=p)
  max.sing = max(X.svd$d)
  if (is.null(tausq.scale)) {
    tausq.scale = 0.5
  } else if (tausq.scale >= 1.0) {
    stop("Sigma = sigma^2 I - tau^2 XX' is not positive definite")
  }
  tau = sqrt(tausq.scale) * sigma / max.sing
  # Set eta.phi
  if (is.null(eta.phi)) {
    eta.phi = logspace(0,-1,max.iter)
  }
  L = length(eta.phi)
  if (L < max.iter+1) {
    eta.phi = c(eta.phi,rep(eta.phi[L],max.iter+1-L))
  } else {
    eta.phi = eta.phi[1:(max.iter+1)]
  }
  # Precompute: X'Sigma^{-1}X
  #             X'Sigma^{-1}y
  D = X.svd$d^2/(rep(sigma^2,min(n,p))-tau^2*X.svd$d^2)
  if (p > n) { D = c(D,rep(0,p-n)) }
  XSinvX = X.svd$v%*%(D*t(X.svd$v))
  D2 = X.svd$d/(rep(sigma^2,min(n,p))-tau^2*X.svd$d^2)
  XSinvy = X.svd$v[,1:min(n,p)]%*%(D2*(t(X.svd$u)%*%y))
  eta.phi.scale = 1/(max(D)+1/tau^2)
  # Set w.init and phi.init
  if (is.null(w.init)) {
    w.init = rep(1/K,K)
  } else {
    w.init[which(w.init<1e-5)] = 1e-5
    w.init = w.init/sum(w.init)
  }
  if (is.null(phi.init)) {
    phi.init = rep(0,p)
  }
  return(list(tau=tau, eta.phi=eta.phi,
              XSinvX=XSinvX, XSinvy=XSinvy, eta.phi.scale=eta.phi.scale,
              w.init=w.init, phi.init=phi.init))
}

# X -- n x p design matrix 
# y -- response vector of length n
# sigma -- noise std dev
# Theta -- prior support points
# eta.w -- step size for prior update
# eta.phi -- step size scale for Langevin dynamics
#            Default: 0.1
# tausq.scale -- reparametrize by phi = theta + N(0,tau^2)
#        tau^2 = tausq.scale * sigma^2 / maximum eig of XX'
#        Default: tausq.scale = 0.5
# w.init -- initial prior weights
#           Default: uniform over prior support points Theta
# phi.init -- initialization for Langevin dynamics
#             Default: all 0
# save.iter -- save w and phi in increments of this many iterations
# max.iter -- total number of iterations to run
# verbose -- print progress every save.iter iterations
# w.true -- true w, used only to print progress
# rho -- smoothing parameter in L2, used only to print progress
langevin.EB <- function(X, y, sigma, Theta, lambda=1e-3,
        burn.iters=0, langevin.inner.iters=100,
	eta.phi=NULL, tausq.scale=NULL,
        w.init=NULL, phi.init=NULL, w.true=NULL,
        max.iter=10000, save.iter=100,
	print.iter=100, verbose=TRUE, plot=FALSE) {
  # Initialize default algorithm parameters
  K = length(Theta)
  p = ncol(X)
  inits = langevin.EB.init(X, y, sigma, Theta, eta.phi, max.iter,
          tausq.scale, w.init, phi.init)
  tau=inits$tau; eta.phi=inits$eta.phi
  XSinvX=inits$XSinvX; XSinvy=inits$XSinvy; eta.phi.scale=inits$eta.phi.scale
  w.init=inits$w.init; phi.init=inits$phi.init
  if (verbose) {
    print(sprintf("Langevin EB with sigma = %f, tau = %f", sigma, tau))
  }
  # Run gradient flow EB
  hist = list()
  w = w.init
  phi = phi.init
  # For smoothing spline
  delta = Theta[2]-Theta[1]
  D = discrete.der(K,delta)
  M = sqrt(lambda/(2*delta))*D
  phi.save = c()
  t.start = Sys.time()
  for (iter in 1:(max.iter+1)) {
    if ((iter %% save.iter == 1) | (save.iter == 1)) {
      eta.phi.curr = eta.phi[iter]
      t.curr = Sys.time()
      hist[[iter]] = list(w = w, phi = phi, time = t.curr - t.start)
      # Save history, print change in w, plot estimated w
    }
    if ((iter %% print.iter == 1) | (print.iter == 1)) {
      if (verbose) {
        if (is.null(w.true)) { err = NA }
        else { err = sum(abs(w-w.true))/2 }
        print(sprintf("Iteration: %d, eta.phi: %f, TV from truth: %f", iter, eta.phi.curr, err))
      }
      if (plot) {
        plot(Theta,w.true,type="l")
        lines(Theta,w,col="red")
      }
    }
    # Precompute Gaussian density at coordinates of phi - Theta
    if (max(abs(phi)) > 100 * max(abs(Theta))) {
      print("WARNING: Diverging iterates, consider reducing eta.phi")
    }
    centered = rep(1,K) %o% phi - Theta %o% rep(1,p)
    log.density = dnorm(centered, sd=tau, log=TRUE)
    # Take one Langevin step in phi
    phi = langevin.step(phi, centered, log.density, w, tau, XSinvX, XSinvy, eta.phi.curr*eta.phi.scale)
    phi.save = c(phi.save,phi)
    if (iter %% langevin.inner.iters == 0 && iter >= burn.iters) {
      if (length(phi.save) > 10000) {
        phi.save = sample(phi.save,size=10000,replace=FALSE)
      }
      try({w = update_w(phi.save,Theta,tau,M)}, silent=TRUE)
      w[which(w<1e-5)] = 1e-5
      w = w/sum(w)
      phi.save = c()
    }
  }
  return(list(w = w, hist = hist, sigma = sigma, tau = tau, langevin.inner.iters = langevin.inner.iters, burn.iters = burn.iters, eta.phi = eta.phi))
}

# Compute posterior mean estimate from w and samples of phi
posterior.mean = function(w,phi,Theta,tau) {
  p = nrow(phi); niters = ncol(phi); p.tot = niters*p; K = length(Theta)
  dim(phi) = p.tot
  centered = rep(1,K) %o% phi - Theta %o% rep(1,p.tot)
  log.density.w = dnorm(centered, sd=tau, log=TRUE) + log(w)
  col.max = apply(log.density.w,2,max)
  density.w = t(exp(t(log.density.w)-col.max))
  theta.mean = colMeans(Theta*density.w)/colMeans(density.w)
  dim(theta.mean) = c(p,niters)
  return(rowMeans(theta.mean))
}
