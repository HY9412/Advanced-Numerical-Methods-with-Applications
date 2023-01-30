install.packages("rstan")
library(rstan)
library(bayesplot)

Ukr_data <- read.csv('TrendData.csv')
dat <- Ukr_data[c(-1,-2,-3)]
data <- list(y = dat, nrows = nrow(dat), ncolumns = ncol(dat))

Ukrmodel_fit <- stan(file = 'UkrL.stan',
                     data = data,
                     pars = c('psi','beta0','beta1','bet0','bet1'),
                     iter = 10000,
                     warmup = 5000,
                     chains = 3,
                     thin = 10,
                     algorithm = 'NUTS',
                     control = list(max_treedepth = 20),
                     cores = 6
)

posterior_extract <- as.matrix(Ukrmodel_fit)
# Define the vector of parameters
pars_bet0 <- c(paste0("bet0[", 1:27, "]"))
pars_bet1 <- c(paste0("bet1[", 1:27, "]"))
# Plot regional intercepts and regional slopes
bayesplot::mcmc_areas(Ukrmodel_fit, pars = pars_bet0, prob = 0.95)
bayesplot::mcmc_areas(Ukrmodel_fit, pars = pars_bet1, prob = 0.95)
                      

