  library(ggplot2)
  library(scales)
  library(RColorBrewer)
  library(ggrepel)
  library(data.table)
  library(bayesplot)
  library(loo)
  require(bridgesampling)
  require(brms)
  options(mc.cores = 4)

setwd('~/go_nuts/')
trials1 = fread('~/go_nuts/trials022319.csv', 
            sep = '\t',
            header = FALSE,
            col.names = c('start', 'surface', 'end'))
trials1[, .N, by = .(start, surface)]
trials2 = fread('~/go_nuts/trials022409.csv',
                sep = '\t',
                header = FALSE,
                col.names = c('start', 'surface', 'end'))
trials2[, .N, by = .(start, surface)]
trials3 = fread('~/go_nuts/trials022519.csv',
                sep = '\t',
                header = FALSE,
                col.names = c('start', 'surface', 'end'))
trials3[, .N, by = .(start, surface)]
trials4 = fread('~/go_nuts/trials02271408.csv',
                sep = '\t',
                header = FALSE,
                col.names = c('start', 'surface', 'end'))
trials4[, .N, by = .(start, surface)]

trials = rbind(
  trials1, 
  trials2, 
  trials3, 
  trials4
)
trials[, startf := factor(start, levels = c('T', 'H'))]
trials[, endf := factor(start, levels = c('T', 'H'))]
trials[, surfacef := factor(surface, levels = c('S', 'C'))]
trials[ , .N, by = .(startf, surfacef)]
View(trials)
trials[, start_heads := as.integer(ifelse(start == 'H', 
                        1, 
                        ifelse(start == 'T', 
                               0, 
                               NA)))]
dput(colnames(trials))
trials[, end_heads := as.integer(ifelse(end == 'H', 
                        1, 
                        ifelse(end == 'T', 
                               0, 
                               NA)))]
trials[, surface_cloth := as.integer(ifelse(surface  == 'C', 1,
                          ifelse(surface == 'S', 0,
                                 NA)))]
trial_count = trials[, .N, by = .(start_heads, surface_cloth, end_heads)]
trial_count2 = trial_count[, .(sumN = sum(N), end_heads_count = sum(ifelse(end_heads == 1, N, 0))), by = .(start_heads, surface_cloth)]
trial_count2

trials[, .(start_heads, surface_cloth, end_heads)]
bernoulli_prior <- brms::brm(formula = end_heads ~ start_heads + surface_cloth,
                           data = trials[, .(start_heads, surface_cloth, end_heads)],
                           family = 'bernoulli',
                           prior = set_prior('student_t(10, 0, 2)', class = 'b'),
                           chains = 4, 
                           cores = 4,
                           iter = 10000,
                           sample_prior = 'only')
summary(bernoulli_prior)
brms::marginal_effects(bernoulli_prior)
plot(bernoulli_prior)
pp = brms::pp_check(bernoulli_prior, nsamples = 100)
pp + theme_minimal()
pairs(bernoulli_prior)

(trials[, .(.N, 
            end_heads = sum(end_heads), 
            end_heads_rate = sum(end_heads)/.N), 
        keyby = .(surface_cloth, start_heads)]
)


bernoulli_fit <- brms::brm(
  formula = end_heads ~ 1 + intercept + start_heads + surface_cloth + surface_cloth : start_heads,
  data = trials[, .(start_heads, surface_cloth, end_heads)],
  family = 'bernoulli',
  prior = set_prior('student_t(10, 0, 2)', class = 'b'),
  chains = 4, 
  cores = 4,
  iter = 10000)

add_ic(bernoulli_fit) <- 'loo'
summary(bernoulli_fit)
brms::marginal_effects(bernoulli_fit) 
plot(bernoulli_fit)
pp = brms::pp_check(bernoulli_fit, nsamples = 100)
pp + theme_minimal()
pairs(bernoulli_fit)
# predictive_interval(bernoulli_fit)  # Not useful


bernoulli_s1 <- brms::brm(
  formula = end_heads ~ 0 + intercept + surface_cloth,
  data = trials[, .(start_heads, surface_cloth, end_heads)],
  family = 'bernoulli',
  prior = set_prior('student_t(10, 0, 2)', class = 'b'),
  chains = 4, 
  cores = 4,
  iter = 10000)

summary(bernoulli_s1)
brms::marginal_effects(bernoulli_s1) 
plot(bernoulli_s1)
pp = brms::pp_check(bernoulli_s1, nsamples = 100)
pp + theme_minimal()
pairs(bernoulli_s1)

bernoulli_sf <- brms::brm(
  formula = endf ~ 0 + intercept + surfacef,
  data = trials[, .(startf, surfacef, endf)],
  family = 'bernoulli',
  prior = set_prior('student_t(10, 0, 2)', class = 'b'),
  chains = 4, 
  cores = 4,
  iter = 10000)

summary(bernoulli_sf, prob = .89)
brms::marginal_effects(bernoulli_sf) 
plot(bernoulli_sf)
pp = brms::pp_check(bernoulli_sf, nsamples = 100)
pp + theme_minimal()
pairs(bernoulli_sf)
require(GGally)
require(ggmcmc)
S <- ggs(bernoulli_sf)
ggmcmc(S)
ggpl <- ggs_pairs(S, lower = list(continuous = 'density'))
ggpl


bernoulli_fitf <- brms::brm(
  formula = endf ~ 0 + intercept + startf + surfacef,
  data = trials[, .(startf, surfacef, endf)],
  family = 'bernoulli',
  prior = set_prior('student_t(10, 0, 1)', class = 'b'),
  chains = 4, 
  cores = 4,
  iter = 10000)

add_ic(bernoulli_fitf) <- 'loo'
summary(bernoulli_fitf)
brms::marginal_effects(bernoulli_fitf) 
plot(bernoulli_fitf)
pp = brms::pp_check(bernoulli_fitf, nsamples = 100)
pp + theme_minimal()
pairs(bernoulli_fitf)
# predictive_interval(bernoulli_fit)  # Not useful


fit1 <- brms::brm(formula = end_heads_count | trials(sumN) ~ start_heads + surface_cloth,
            data = trial_count2,
            family = binomial(),
            prior = set_prior('normal(0, 1)'),
            chains = 4, 
            cores = 4,
            iter = 4000)
summary(fit1)

brms::marginal_effects(fit1) 
plot(fit1)
pp = brms::pp_check(fit1, nsamples = 100)
pp + theme_minimal()
pairs(fit1)
