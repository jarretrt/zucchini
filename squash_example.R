library('dplyr')
library('ggplot2')
library('magrittr')

source('squash.R')

data('iris')

## squash iris data
iris_squash <- squash(iris)

## plot original and squashed data
iris_all <- bind_rows(
  "Original"=iris,
  "Squashed"=iris_squash,
  .id='Dataset')

iris_ggp <- iris_all %>%
  ggplot(aes(x=Sepal.Width, y=Sepal.Length, col=Species)) +
  facet_wrap(~Dataset) +
  geom_point() + 
  theme(legend.position = 'top')
print(iris_ggp)

## data for prediction
rng <- range(iris$Sepal.Width)
pre <- expand.grid(
  Sepal.Width = seq(rng[1], rng[2], length.out=100),
  Species = unique(iris$Species))

## fit original data using lm
fit_n <- lm(Sepal.Length ~ Sepal.Width * Species, data=iris)
fit_n %>% summary
pre_n <- predict(fit_n, newdata=pre, interval='confidence', )
pre_n <- cbind(pre, pre_n) %>% mutate(Dataset='Original')

## fit squashed data using lm
## 'df.residual' needs to be corrected for this to work
fit_m <- lm(Sepal.Length ~ Sepal.Width * Species,
  data=iris_squash, weights=iris_squash$`(weight)`)
fit_m$df.residual <- with(fit_m, sum(weights)-length(coefficients))
pre_m <- predict(fit_m, newdata=pre, interval='confidence', )
pre_m <- cbind(pre, pre_m) %>% mutate(Dataset='Squashed')

## plot model fits for both data sets
pre <- bind_rows(pre_n, pre_m)
iris_ggp +
  geom_line(data=pre, aes(y=fit)) +
  geom_ribbon(data=pre, aes(y=fit, ymin=lwr, ymax=upr,
    fill=Species), alpha=0.5) +
  theme(legend.position = 'top')
