library('microbenchmark')

x <- matrix(rnorm(10), 5)
w <- 1:5

cen <- colMeans(w*x)/sum(w)

## benchmark alternative rowsums
microbenchmark(
  colSums(w*x),
  (rep(1,nrow(x))%*%(w*x))[1,],
  times=1e4)

## benchmark matrix transpose
microbenchmark(t(x), aperm(x), times=1e5)

## benckmark centering columns of matrix
microbenchmark(
  scale(x, center=cen, scale=F),
  sweep(x, 2, cen),
  t(t(x) - cen),
  t.default(t.default(x) - cen),
  times=1e4)

x <- scale(x, center=cen, scale=F)

## t() is expensive; test alternative
microbenchmark(
  crossprod(t(t(x) - cen)*sqrt(w)),
  tcrossprod((t(x) - cen)*matrix(sqrt(w),dim(x)[2],dim(x)[1],byrow=T)),
  times=1e4)

## benckmark computing second moments
microbenchmark(
  cov.wt(x, w, center=F, method='ML')$cov*sum(w),
  crossprod(sqrt(w)*x),
  times=1e4)
