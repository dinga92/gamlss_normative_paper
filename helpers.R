
# install.packages('gamlss')
library(gamlss.dist)

params_to_scores <- function(y, params){
  # if the params has only 2 columns instead of expected 4
  # i.e. if the model was fitted using the gaussian family
  # instead of shash family, we add 0 to the other 2 columns
  # that corresponds to parameters of the gaussian distribution
  # as parametrized by shash distribution
  if(ncol(params) == 2){
    params <- cbind(params, rep(0, nrow(params)))
    params <- cbind(params, rep(log(1), nrow(params)))
  }

  # although the models were fitted using mgcv, we use gamlss to transform
  # the parameters of the distribution, because it's easier to use and it's
  # parametrized as other distributions in r
  log_densities <- gamlss.dist::dSHASHo2(y,
                                         mu=params[,1],
                                         # should be exp, because models were fitted
                                         # with the log link
                                         sigma=exp(params[,2]),
                                         nu=params[,3],
                                         # exp, because log link fink
                                         tau=exp(params[,4]),
                                         log=T)

  quantiles <- gamlss.dist::pSHASHo2(y, mu=params[,1],
                                     # exp, because log link fink
                                     sigma=exp(params[,2]),
                                     nu=params[,3],
                                     # exp, because log link fink
                                     tau=exp(params[,4]))

  # randomized residuals are quantiles transformed to corresponding
  # z-scores of standard normal
  z_randomized <- qnorm(quantiles, mean = 0, sd = 1)

  # censoring log densities that aren't between 1-99 percentiles
  log_densities_censored <- log_densities
  log_densities_censored[abs(z_randomized) > 2.326348] <- log(0.02)

  data.frame('log_densities' = log_densities,
             'log_densities_censored' = log_densities_censored,
             'quantiles' = quantiles,
             'z_randomized' = z_randomized)
}



# Skewness and kurtosis and their standard errors as implement by SPSS
#
# Reference: pp 451-452 of
# http://support.spss.com/ProductsExt/SPSS/Documentation/Manuals/16.0/SPSS 16.0 Algorithms.pdf
#
# See also: Suggestion for Using Powerful and Informative Tests of Normality,
# Ralph B. D'Agostino, Albert Belanger, Ralph B. D'Agostino, Jr.,
# The American Statistician, Vol. 44, No. 4 (Nov., 1990), pp. 316-321

calibration_descriptives <- function(x) {
  n <- length(x)
  m1 <- mean(x)
  m2 <- sum((x-m1)^2)
  m3 <- sum((x-m1)^3)
  m4 <- sum((x-m1)^4)
  s1 <-sd(x)
  skew <- n*m3/(n-1)/(n-2)/s1^3
  sdskew <- sqrt( 6*n*(n-1) / ((n-2)*(n+1)*(n+3)) )
  kurtosis <- (n*(n+1)*m4 - 3*m2^2*(n-1)) / ((n-1)*(n-2)*(n-3)*s1^4)
  sdkurtosis <- sqrt( 4*(n^2-1) * sdskew^2 / ((n-3)*(n+5)) )

  semean <- sqrt(var(x)/n)
  sesd <- s1/sqrt(2*(n-1))

  # TODO: find a better test
  # shapiro.test does not allow n > 5000 so we calculate it in a random
  # subsample instead, this shouldn't make a difference
  if (length(x) < 5000){
    W = shapiro.test(x)$statistic
  } else {
    W = shapiro.test(sample(x, 5000))$statistic
  }

  return(data.frame(t(c('mean'=m1,
                        'se.mean'=semean,
                        'sd'=s1,
                        'se.sd'=sesd,
                        'skew'=skew,
                        'se.skew'=sdskew,
                        'kurtosis'=kurtosis,
                        'se.kurtosis'=sdkurtosis,
                        'W' = W))))
}



params_to_quantiles <- function(params, quantiles){
  # if the predictions only have two columns, that means they are mu and sigma
  # coming from a gaussian model, therefore we add nu and tau == 0,
  # corresponding to nu and tau of a gaussian
  if(ncol(params) == 2){
    params <- cbind(params, matrix(0, nrow = nrow(params), ncol = 2))
  }
  names(params) <- c('predictions.mu', 'predictions.sigma',
                     'predictions.nu', 'predictions.tau')
  return(as.data.frame(sapply(quantiles,
                              function(q){gamlss.dist::qSHASHo2(q, params$predictions.mu,
                                                                exp(params$predictions.sigma),
                                                                params$predictions.nu,
                                                                exp(params$predictions.tau))})))
}
