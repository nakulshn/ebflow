library(Matrix)
library(MASS)

# Generate weights for a discrete prior on [-3,3] with K grid points
# Set prior to:
#   gaussian -- normal(0,1)
#   bimodal -- 0.5 * normal(-1.5,0.5^2) + 0.5 * normal(1.5,0.5^2)
#   cauchy -- cauchy(0,0.6)
#   skew -- 0.33 * normal(-2,0.5^2) + 0.33 * normal(-1.5,1^2) + 0.33 * normal(0,2^2)
# all truncated to [-3,3] and discretized to the grid
#
# Return: w.grid -- support points for grid
#         w.true -- prior weight at each support point
gen_prior = function(K, prior) {
  w.grid = seq(from=-3, to=3, length=K)
  if (prior == 'gaussian') {
    w.true = dnorm(w.grid, mean=0, sd=1)
  } else if (prior == 'bimodal') {
    w.true = dnorm(w.grid, mean=-1.5, sd=0.5) + dnorm(w.grid, mean=1.5, sd=0.5)
  } else if (prior == 'cauchy') {
    w.true = dcauchy(w.grid,0,0.6)
  } else if (prior == 'skew') {
    w.true = dnorm(w.grid, mean=-2, sd=0.5) + dnorm(w.grid, mean=-1.5, sd=1) + dnorm(w.grid, mean=0, sd=2)
  }
  w.true = w.true/sum(w.true)
  return(list(w.grid=w.grid, w.true=w.true))
}

# Generate design matrix
#
# Return: X.train -- n x p design, columns having roughly mean 0, variance 1/p
gen_X = function(n, p, blocksize, corr, n.test=1000) {
  if (blocksize == 'identity') {
    X.train = diag(n)
    X.test = NA
  } else {
    block = matrix(corr,blocksize,blocksize)
    diag(block) = 1
    nblocks = ceiling(p/blocksize)
    X.stack = mvrnorm((n+n.test)*nblocks, mu=rep(0,blocksize), Sigma=block) / sqrt(n)
    X.train = matrix(0,n,p)
    X.test = matrix(0,n.test,p)
    for (l in 1:nblocks) {
      cols = ((l-1)*blocksize+1):(l*blocksize)
      X.train[,cols] = X.stack[((l-1)*(n+n.test)+1):((l-1)*(n+n.test)+n),]
      X.test[,cols] = X.stack[((l-1)*(n+n.test)+n+1):(l*(n+n.test)),]
    }
  }
  return(list(X.train=X.train, X.test=X.test))
}

# Generate all data and save to disk
set.seed(123)
priors = c("gaussian","skew","cauchy","bimodal")
designs = list(list(dname="identity",blocksize="identity",corr=NA),
        list(dname="iid",blocksize=1,corr=0),
        list(dname="block02corr0.9",blocksize=2,corr=0.9),
        list(dname="block10corr0.5",blocksize=10,corr=0.5))
dim.ratios = c(0.5,1,2)
noise.ratios = c(0.8,0.5)
K = 61
p = 1000

for (design in designs) {
  list2env(design,.GlobalEnv) # define dname, blocksize, and corr
  for (dim.ratio in dim.ratios) {
    if (dname == "identity" && dim.ratio != 1) { next }
    print(c(dname,dim.ratio))
    n = p * dim.ratio
    list2env(gen_X(n,p,blocksize,corr),.GlobalEnv) # define X.train and X.test
    fname = sprintf("data/X_%s_nbyp%.1f.rds", dname, dim.ratio)
    saveRDS(list(X.train=X.train,X.test=X.test), file=fname)
    for (prior in priors) {
      list2env(gen_prior(K,prior),.GlobalEnv) # define w.grid and w.true
      for (noise.ratio in noise.ratios) {
        y = list()
        theta.true = sample(w.grid, p, replace=TRUE, prob=w.true)
        sigma = sd(X.train %*% theta.true)*sqrt((1-noise.ratio)/noise.ratio)
        y = X.train %*% theta.true + rnorm(nrow(X.train), mean=0, sd=sigma)
        fname = sprintf("data/y_%s_%s_nbyp%.1f_invSNR%.1f.rds", prior, dname, dim.ratio, noise.ratio)
        saveRDS(list(w.grid=w.grid,w.true=w.true,theta.true=theta.true,y=y,sigma=sigma), file=fname)
      }
    }
  }
}

