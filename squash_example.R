library('dplyr')
source('squash.R')
data('iris')

## squash iris data
iris_squash <- squash(iris)

## plot original and squashed data
with(dat, plot(Sepal.Width, Sepal.Length, pch=20, col='gray'))
with(dat_m, points(Sepal.Width, Sepal.Length, pch=20, col='#800080'))

## data for prediction
pre <- expand.grid(
  Sepal.Width = seq(min(dat$Sepal.Width),
    max(dat$Sepal.Width), length.out=100),
  Species = unique(dat$Species))

## compare linear model fits using lm
fit_n <- lm(Sepal.Length ~ Sepal.Width * Species, data=dat)
pre_n <- predict(fit_n, newdata=pre, interval='confidence')
cbind(pre, pre_n) %>%
  group_by(Species) %>%
  group_walk(function(x, y) {
    lines(x$Sepal.Width, x$fit, col='gray', lwd=2)
    polygon(x=c(x$Sepal.Width, rev(x$Sepal.Width)),
            y=c(x$lwr, rev(x$upr)), border=NA, col='#00000022')
  })
fit_n %>% summary

## 'df.residual' needs to be corrected for this to work
fit_m <- lm(Sepal.Length ~ Sepal.Width * Species, data=dat_m,
  weights=dat_m$`(weight)`)
fit_m$df.residual <- with(fit_m, sum(weights)-length(coefficients))
pre_m <- predict(fit_m, newdata=pre, interval='confidence')
cbind(pre, pre_m) %>%
  group_by(Species) %>%
  group_walk(function(x, y) {
    lines(x$Sepal.Width, x$fit, col='#800080', lty=2, lwd=2)
    polygon(x=c(x$Sepal.Width, rev(x$Sepal.Width)),
            y=c(x$lwr, rev(x$upr)), border=NA, col='#80008022')
  })