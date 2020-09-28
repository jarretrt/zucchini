## construct a symmetric tri-diagonal matrix
## dia - diagonal values
## off - off-diagonal values
tri_diag <- function(dia, off) {
  if(dia <= 0 || length(dia) != (length(off) + 1))
    stop('dia and off do not have conformal lengths')
  value <- diag(dia)
  rows <- row(value)
  cols <- col(value)
  value[cols == rows + 1] <- off
  value[cols == rows - 1] <- off
  value
}

## Compute Gauss-Hermite quadrature nodes (points) and
## weights using the Golub-Welsch algorithm; O(n^2)
## n - number of quadrature nodes (points)
## norm - whether to re-weight for computing integrals
##        with respect to standard normal density
gauss_hermite <- function(n, norm=TRUE) {
  # Golub-Welsch algorithm
  # 3-term recurrence coeffs
  beta <- sqrt(0.5 * (1:(n-1)))    
  # Jacobi matrix
  jaco <- tri_diag(rep(0,n), beta) 
  # Eigenvalue decomposition
  eige <- eigen(jaco)              
  # Hermite points
  indx <- order(eige$values)       
  x <- eige$values[indx]
  # Weights
  w <- sqrt(pi)*eige$vectors[1,indx]^2
  ## The weight function for Gauss-Hermite quadrature is 
  ## w(x) = exp(-x^2), see statmod::gauss.quad, for example
  ## For integrals with respect to standard normal density, need
  ## to adjust weights so that w(x) = (2*pi)^(1/2)*exp(-1/2*x^2)
  if(norm)
    w <- w * (2*pi)^(-1/2)*exp(x^2/2)
  return(list(points=x, weights=w))
}

## Compute multivariate Gauss-Hermite quadrature nodes and weights
## n     - number of points each dimension
## mu    - mean vector
## sigma - covariance matrix
mgauss_hermite <- function(n, mu, sigma) {
  if(!all(dim(sigma) == length(mu)))
    stop("mu and sigma have nonconformable dimensions")
  dm  <- length(mu)
  gh  <- gauss_hermite(n)
  #idx grows exponentially in n and dm
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh$points[idx],nrow(idx),dm)
  wts <- apply(matrix(gh$weights[idx],nrow(idx),dm), 1, prod)
  
  ## rotate, scale, translate points
  eig <- eigen(sigma) 
  rot <- eig$vectors %*% diag(sqrt(eig$values))
  pts <- t(rot %*% t(pts) + mu)
  return(list(points=pts, weights=wts))
}

