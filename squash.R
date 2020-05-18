## The coding in this R script is based primarily on the 
## theory and descriptions provided in the following article:
## William DuMouchel, Chris Volinsky, Theodore Johnson, 
## Corinna Cortes, and Daryl Pregibon. 1999. Squashing flat 
## files flatter. In Proceedings of the fifth ACM SIGKDD
## international conference on Knowledge discovery and data 
## mining (KDD '99). Association for Computing Machinery, 
## New York, NY, USA, 6-15. DOI:https://doi.org/10.1145/312129.312184

## questions:
##   - what happens if you treat indicators as numeric?
##   - logistic regresssion?

library('dplyr')
library('magrittr')

## use the data sphere method to partition quantitative data
## dat    - a data frame with only numeric variables
## layers - number of layers in partition
data_sphere_partition <- function(dat, layers=3) {
  
  if(!all(sapply(dat, is.numeric)))
    stop("'dat' may only contain numeric data")
  
  ## scale data 
  dsc <- dat %>% as.matrix %>% scale 
  
  ## and compute distances from center
  dst <- dsc %>% `^`(2) %>% rowSums %>% sqrt
  
  ## compute distance quantiles for layers
  lev <- quantile(dst, probs=seq(1/layers, 1-1/layers, 1/layers))
  
  ## partition into layers
  lay <- cut(dst, c(-Inf, lev, Inf)) %>% unclass
  
  ## partition into pyramids
  pmx <- cbind(-dsc %>% pmax(0), dsc %>% pmax(0))
  pyr1 <- pmx %>% apply(1, which.max)
  
  ## partition into sub-pyramids
  pmx[cbind(1:nrow(pmx),pyr1)] <- 0
  pyr2 <- pmx %>% apply(1, which.max)
  
  ## combine partition information
  prt <- 
    ifelse(lay == 1, paste(lay), 
      ifelse(lay == 2, paste(lay, pyr1, sep='.'),
        paste(lay, pyr1, pyr2, sep='.')))
  
  dat <- dat %>% mutate(`(partition)` = prt)
  attr(dat, 'scaled:center') <- attr(dsc, 'scaled:center')
  attr(dat, 'scaled:scale') <- attr(dsc, 'scaled:scale')
  attr(dat, 'layer:breaks') <- lev 
  return(dat)
}

## formula for computing reduced samples size
reduced_sample_size <- function(n, alpha)
  ceiling(max(1, alpha*log(n, 2)))

## compute up to fourth order moments for spedified data
## dat - a data frame with only numeric variables
## m   - reduced sample size (determines how many moments to return)
compute_moments_v1 <- function(dat, m = nrow(dat)) {

  ## if weight column missing, assume unit weights
  if(!('(weight)' %in% colnames(dat)))
    dat %<>% cbind('(weight)'= 1)
  
  ## extract weights
  wix <- match('(weight)', colnames(dat))
  wgt <- dat[,wix]
  
  ## compute weighted sums (for weighted mean)
  cen <- dat %>% apply(1, function(x) 
      x[wix] * x[-wix]) %>% rowSums
  
  ## remove weights  
  dsc <- dat[,-wix]
  
  ## center
  dsc <- dsc %>% scale(center=cen/sum(wgt), scale=F)
  
  ## rejoin weights
  dsc <- cbind(wgt, dsc)
  
  ## number of moments to return
  k <- m * ncol(dsc)
  
  ## compute weighted mixed moments up to fourth order
  ## x[1] are the weights, and x[-1] are the data
  ## this could be done in a 'for' loop to cut memory usage
  mmt <- c(
    sum(wgt),
    cen, 
    dsc %>% apply(1, function(x) 
      x[1] * (x[-1] %o% x[-1])) %>% rowSums,
    dsc %>% apply(1, function(x) 
      x[1] * (x[-1] %o% x[-1] %o% x[-1])) %>% rowSums,
    dsc %>% apply(1, function(x) 
      x[1] * (x[-1] %o% x[-1] %o% x[-1] %o% x[-1])) %>% rowSums)


  
  ## return first k moments and remove names
  return(unname(mmt[1:k]))
  
}

## compute up to second order moments for spedified data
## dat - a data frame with only numeric variables
compute_moments <- function(dat, ...) {

  ## if weight column missing, assume unit weights
  if(!('(weight)' %in% colnames(dat)))
    dat %<>% cbind('(weight)'= 1)
  
  ## extract weights
  wix <- match('(weight)', colnames(dat))
  wgt <- dat[,wix]
  dsc <- dat[,-wix]
  
  ## compute weighted sums (for weighted mean)
  cen <- colSums(wgt*dsc)

  ## center
  dsc <- dsc %>% scale(center=cen/sum(wgt), scale=F)
  
  ## compute weighted moments up to second order
  mmt <- c(sum(wgt), cen, crossprod(sqrt(wgt)*dsc))
  
  ## return moments and remove names
  return(unname(mmt))
}

## names for moments (not used)
collapse_names <- function(y) {
  expand.grid(dimnames(y)) %>% 
    apply(1, function(y0) {
      t0 <- table(y0)
      paste0(names(t0), 
        ifelse(t0 > 1, paste0('^',t0), ''),
        collapse='*')
    })  
}

## logistic transform for optimization variables
logi  <- function(x, ymin=0, ymax=1)
  ymin + (ymax-ymin) / (1 + exp(-x))

## inverse logistic transform for optimization variables
ilogi <- function(y, ymin=0, ymax=1)
   -log((ymax-ymin) / (y - ymin) - 1) 

## squash is a generic function
squash <- function(obj, ...)
  UseMethod('squash', obj)

## squash numeric data matrix (i.e., for an individual partition)
## dat   - numeric matrix
## alpha - degree of squashing (see 'reduced_sample_size') 
## ...   - arguments passed to 'optim'
squash.matrix <- function(dat, alpha=1, 
  method='CG', control=list(maxit=500), ...) {

  ## stop if not all data are numeric
  if(!is.numeric(dat))
    stop("'dat' must be a numeric matrix")
  
  ## full data sample size/dimension
  n <- nrow(dat)
  p <- ncol(dat)
  
  ## reduced sample size
  m <- reduced_sample_size(n, alpha)
 
  ## calculate moments for full data
  mmt <- compute_moments(dat, m)
  
  ## compute ranges (extended) for each variable
  extendrange <- function(x, f=0.05) {
    r <- range(x, na.rm=TRUE)
    r + c(-f, f) * (diff(r) + f)
  }
  rng <- dat %>% apply(2, extendrange) %>%
    cbind('(weight)'=c(0,n))

  ## vectors of min and max values
  ymn <- rep(rng[1,], rep(m,p+1)) %>% matrix(m)
  ymx <- rep(rng[2,], rep(m,p+1)) %>% matrix(m)

  ## created inital pseudo data for first partition
  dat_m <- dat[sample(n, m),] %>%
    cbind(`(weight)` = 1)
  
  ## store some info about the squashed data
  dmn_m <- dimnames(dat_m)
    
  opt <- optim(ilogi(dat_m, ymn, ymx), function(y) {
 
    ## reconstruct data
    dat_m <- matrix(logi(y, ymn, ymx),
        nrow=m, dimnames=dmn_m)
  
    ## compute moments for reduced data
    mmt_m <- compute_moments(dat_m)
  
    ## compute objective function; means and (co-)vars get priority
    u <- c(1000, rep(1000, p), diag(rep(500), p) + 500)
    if(length(u) < length(mmt_m))
      u <- c(u, rep(1, length(mmt_m)-length(u)))
    u <- u[1:length(mmt_m)]
    
    ## weighted sum of squares for moments
    sum(u * (mmt - mmt_m)^2)  

  }, method=method, control=control, ...)
  
  if(opt$convergence)
    warning("failed convergence ", opt$message)
  
  ## reconstruct reduced data
  dat_m <- opt$par %>% logi(ymn, ymx)
  
  return(dat_m)
    
}

## compute weighted mean and vcov
wtd_mean_vcov <- function(dat) {

  ## if weight column missing, assume unit weights
  if(!('(weight)' %in% colnames(dat)))
    dat[,'(weight)'] <- 1

  if(is.data.frame(dat))
    dat <- as.matrix(dat)

  if(!is.numeric(dat))
    stop("'dat' must contain only numeric data")
  
  wix <- which(colnames(dat) == '(weight)')

  wgt <- dat[,wix]

  cen <- dat %>% apply(1, function(x)
      x[wix]/sum(wgt) * x[-wix]) %>% rowSums

  dat[,-wix] <- scale(dat[,-wix], center=cen, scale=F)

  vcv <- dat %>% apply(1, function(x)
    x[wix]/sum(wgt) * x[-wix] %o% x[-wix]) %>% rowSums %>%
    matrix(length(cen), dimnames=list(names(cen),names(cen)))

  list(mean = cen, vcov = vcv)
}

missing_pattern <- function(dat) {
  
  if(is.data.frame(dat))
    dat <- as.matrix(dat)

  if(!is.numeric(dat))
    stop("'dat' must contain only numeric data")
  
  as.factor(apply(1-is.na(dat), 1, paste0, collapse=''))
  
}

## squash data frame with partitioning
## all combinations of categorical variables are treated as 
## distinct data partitions; missing values in categorical
## variables are treated as a distinct category; the missingness
## pattern (see 'missing_pattern') is treated as a categorical
## variable, which results in some partitions that are missing
## values for all records on one or more variables;
## dat    - data frame (with possibly categorical variables)
## alpha  - degree of squashing (see 'reduced_sample_size') 
## method, control - arguments passed to to 'optim'
squash.data.frame <- function(dat, alpha=1,
  method='CG', control=list(maxit=500)) {
  
  ## if any missing numeric data, create missing pattern
  num <- dat %>% select_if(is.numeric)
  if(any(is.na(num)))
    dat %<>% mutate('(missing)' = missing_pattern(num))
  
  # ## if all numeric data, shortcut
  # if(all(sapply(dat, is.numeric)))
  #   return(as.data.frame(squash(as.matrix(dat),
  #     alpha, method, control)))
  
  ## group by non-numeric variables
  dat %<>% group_by_if(~!is.numeric(.))
  
  ## squash numeric data by group
  dat_m <- group_modify(dat, function(x, y) {
    
    ## remove columns with all missing before squashing
    mis <- x %>% select_if(~all(is.na(.))) %>% distinct()
    obs <- x %>% select_if(~!all(is.na(.)))
    
    ## squash numeric data
    squ <- as.data.frame(squash(as.matrix(obs),
      alpha, method, control))
    
    ## reconstruct data
    merge(squ, mis) %>%
      select(names(x), '(weight)')
    })
  
  dat_m %<>%
    ungroup %>%
    select(names(dat), '(weight)')
  
  return(dat_m)
}

