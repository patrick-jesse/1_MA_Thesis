# ******************************************************************** #
#                                                                      #
# Bayesian Linear Regression Model, Gibbs Sampler                      #
# Patrick Gyasi, MA Thesis, October 2023                               #  
# Objective: estimate the slope of the New Keynesian Phillips Curve    #
# for member states of the EA-11 country composition                   #
# pi_{i,t} = a*pi_{i,t+1} + b*x_{i,t-4} + e                                        #
#                                                                      #
# ******************************************************************** #

# idea: estimate and save the results of all models for each country within 
# a loop. models correspond to different measures of economic activity x_t

# to measure the duration of the code
start.time <- Sys.time()

library(tidyverse)
# library(stargazer)

load("C:/Users/patri/Desktop/1_MA_Thesis/data_NKPC.RData")

#------------------#
#####  set up  #####
#------------------#

postr_geo_list <- vector(mode = "list", length = length(country_names_ecb))
names(postr_geo_list) <- country_names_ecb

# to save the posterior mean, sd and credible intervals
postr.stats_slopeNKPC_list <- vector(mode = "list", length = length(country_names_ecb))
postr.stats_coeffiExpInfl_list <- vector(mode = "list", length = length(country_names_ecb))
postr.stats_sigma_sq_list <- vector(mode = "list", length = length(country_names_ecb))

# creating a list to save the results of each model with every iteration
# "model_i" will be used to access and subset elements of the list
model_i <- names(data_NKPC)[5:length(data_NKPC)]
postr_models_list <- vector(mode = "list", length = length(data_NKPC) - 4)
names(postr_models_list) <- model_i    # paste0("model_", 1:11)

# OLS estimation to cross-check with Bayesian results
ols_geo_list <- vector(mode = "list", length = length(country_names_ecb))
names(ols_geo_list) <- country_names_ecb

# a list to save the OLS results of each model
ols_models_list <- vector(mode = "list", length = length(data_NKPC) - 4)
names(ols_models_list) <- model_i      # paste0("model_", 1:11)

# list to save plots of summary statistics from the gtExtras package
summary_stats <- vector(mode ="list", length = length(country_names_ecb))
names(summary_stats) <- country_names_ecb

# --------------------------- #
# --- GIBBS SAMPLER SETUP --- #
# --------------------------- #

nsave <- 30000          # number of saved draws from the MCMC sampling
nburn <- 20000          # number of burn-ins (to eliminate initial value dependence)
ntot  <- nsave + nburn  # total number of interations
# P     <- ncol(X)        # number of regressors

# a dataframe to save the MCMC sampled estimates of the parameters of each 
# model: the slope of the NKPC, the coefficient of expected inflation and the 
# error variance sigma square.
# these estimates will be used to compute the posterior mean and plot the 
# posterior distribution.

postr_slopeNKPC <- matrix(data = NA, nrow = nsave, 
                          ncol = length(data_NKPC) - 4,
                          dimnames = list(NULL, model_i)
) |> as.data.frame()

postr_coeffiExpInfl <- matrix(data = NA, nrow = nsave, 
                              ncol = length(data_NKPC) - 4,
                              dimnames = list(NULL, model_i)
) |> as.data.frame()

postr_sigma_sq <- matrix(data = NA, nrow = nsave, 
                         ncol = length(data_NKPC) - 4,
                         dimnames = list(NULL, model_i)
) |> as.data.frame()

set.seed(123)

#--------------------------------------------------------------#
##### loop to estimate all models for each country at once #####
#--------------------------------------------------------------#

for(geo in country_names_ecb) { # unique(data_NKPC$country)
  
  y <- data_NKPC[data_NKPC$country == geo, "hicpX", drop = TRUE] # |> head()
  
  #----------------------------------------#
  #####  loop to estimate the models   #####
  #----------------------------------------#
  
  for(model in seq_along(postr_models_list)) {
    
    X <- cbind(data_NKPC[data_NKPC$country == geo, "exp_infl", drop = TRUE],
               data_NKPC[data_NKPC$country == geo, model+4, drop = TRUE]
    ) # |> head()
    
    # column names necessary so that the stargazer function can distinguish the
    # different models and create appropriate outputs
    colnames(X) <- c("exp_infl", model_i[model]) # paste0("model_", model))
    
    N <- nrow(X)  # number of observations
    
    # starting value. Gibbs sampling requires only one starting values, thus 
    # only defining a starting value for the error variance.
    sig2_draw <- 1
    
    # to store sampled values after burn-ins
    posterior_samples <- data.frame(coeffiExpInfl = numeric(nsave), # rep(0, nsave),
                                    slopeNKPC = NA,
                                    sigma_sq = NA
    )
    
    # defining constant quantities before the loop for Gibbs sampling.
    XX <- solve(t(X) %*% X)    # (X'X)^(-1)
    
    ols_models_list[[model]] <- lm(y ~ X - 1)
    # (OLS) beta_hat = (X'X)^(-1) * X'y:
    ols_coeff <- coefficients(ols_models_list[[model]])
    
    #---------------------------------#
    #####   MCMC: Gibbs sampling  #####
    #---------------------------------#
    
    for(i in 1:ntot) {
      
      # simulate beta (coefficients of the covariates) from its multivariate 
      # normal conditional:
      # to obtain draws from a multivariate normal distribution, we generate 
      # a kÃ—1 column vector (here: u) of independent N(0, 1) draws 
      # and then multiply this vector by the Cholesky factor (Lynch page 148)
      u <- t(rnorm(2, mean = 0, sd = 1))
      expl_vars_draw <- as.vector(ols_coeff + u %*% chol(sig2_draw * XX))
      
      # simulate sigma from its inverse gamma distribution
      resid <- y - X %*% expl_vars_draw
      SSR   <- t(resid) %*% resid
      # taking the inverse of the gamma distribution since R doesn't have an 
      # inverse gamma function.
      sig2_draw <- 1/rgamma(1, N/2, SSR/2)
      
      # store psoterior samples after burn-in period
      if(i > nburn) {
        posterior_samples[i - nburn, 1:2] <- expl_vars_draw
        posterior_samples[i - nburn, 3]   <- sig2_draw
      }
      
      # to print out the progress of the loop of the Gibbs sampler
      if(i %% 10000 == 0) {
        print(c(i, expl_vars_draw, sig2_draw))
      }
    }
    
    postr_models_list[[model]] <- posterior_samples
    
    postr_slopeNKPC[[model]]     <- posterior_samples$slopeNKPC
    postr_coeffiExpInfl[[model]] <- posterior_samples$coeffiExpInfl
    postr_sigma_sq[[model]] <- posterior_samples$sigma_sq
    
    postr_geo_list[[geo]] <- postr_models_list
    ols_geo_list[[geo]]   <- ols_models_list
    
  }
  
  # writexl::write_xlsx(postr_slopeNKPC,
  #                     paste0("4_results/quarterly_data/tables/1_slope_NKPC/",
  #                            geo, "_postr_slopeNKPC.xlsx")
  #                     )
  # writexl::write_xlsx(postr_coeffiExpInfl,
  #                     paste0("4_results/quarterly_data/tables/2_coeff_exp_infl/",
  #                            geo, "_postr_coeffiExpInfl.xlsx")
  #                     )
  # writexl::write_xlsx(postr_sigma_sq,
  #                     paste0("4_results/quarterly_data/tables/3_sigma_sq/",
  #                            geo, "_postr_sigma_sq.xlsx")
  #                     )
  
  # -------------------------- #
  # --- MCMC CONVERGENCE ----- #
  # -------------------------- #
  
  #---------------------------------------#
  #####  traceplot: slope of the NKPC #####
  #---------------------------------------#
  for(measure in model_i) {
    
    postr_slopeNKPC[measure] |>  # head()
      rename(plot_model = {{measure}}) |> # all_of(measure)) |>  # head()
      # selecting with an external variable is ambiguous because it is not clear 
      # whether you are refering to a data frame column with or an external object.
      # it displays an error message
      # using curly curly {{ }} is a remedy.
      ggplot(aes(x = 1:nsave, y = plot_model)) +
      geom_line() + 
      # geom_hline(aes(yintercept = coefficients(ols_models_list[[measure]])[2]
      # ), colour = "green") +
      theme_bw() + # theme_minimal
      labs(x         = "Iterations", y = "slope NKPC",
           # title     = paste0("trace plot, ", geo),
           subtitle  = paste0(geo, ", trace plot of the slope of the NKPC"),
           caption   = paste0("measure of economic activity: ", measure)
      )
    
    specification <- unlist(strsplit(measure, split = "_"))[1:2] |> 
      paste(collapse = "_")
    
    ggsave(paste0("4_results/quarterly_data/plots/1_traceplots_slopeNKPC/", 
                  geo, "_slopeNKPC_", specification, ".png"),
           width = 20, height = 12, units = "cm"
    )
  }
  
  
  #--------------------------------------------------------#
  #####  traceplot: coefficient of expected inflation  #####
  #--------------------------------------------------------#
  
  for(measure in model_i) {
    
    postr_coeffiExpInfl[measure] |>
      rename(plot_model = {{measure}}) |>
      ggplot(aes(x = 1:nsave, y = plot_model)) +
      geom_line() +
      # geom_hline(aes(yintercept = coefficients(ols_models_list[[measure]])[2]
      # ), colour = "green") +
      theme_bw() + # theme_minimal
      labs(x         = "Iterations", y = "coefficient expected inflation",
           # title     = paste0("trace plot, ", geo),
           subtitle  = paste0(geo, ", trace plot coefficient for expected inflation"),
           caption   = paste0("measure of economic activity: ", measure)
      )
    
    specification <- unlist(strsplit(measure, split = "_"))[1:2] |>
      paste(collapse = "_")
    
    ggsave(paste0("4_results/quarterly_data/plots/2_traceplots_exp_infl/",
                  geo, "_coeffiExpInfl_", specification, ".png"),
           width = 20, height = 12, units = "cm"
    )
  }
  
  #-------------------------------------#
  #####  traceplot: error variance  #####
  #-------------------------------------#
  
  for(measure in model_i) {
    
    postr_sigma_sq[measure] |>  # head()
      rename(plot_model = {{measure}}) |>  # head()
      ggplot(aes(x = 1:nsave, y = plot_model)) +
      geom_line() +
      theme_bw() + # theme_minimal
      labs(x         = "Iterations", y = "sigma square",
           # title     = paste0("trace plot, ", geo),
           subtitle  = paste0(geo, ", trace plot of the error variance"),
           caption   = paste0("measure of economic activity: ", measure)
      )
    
    specification <- unlist(strsplit(measure, split = "_"))[1:2] |>
      paste(collapse = "_")
    
    ggsave(paste0("4_results/quarterly_data/plots/3_traceplots_sigma_sq/",
                  geo, "_sigma_sq_", specification, ".png"),
           width = 20, height = 12, units = "cm"
    )
  }
  
  #----------------------------------------#
  ##### posterior analysis: slope NKPC #####
  #----------------------------------------#
  
  postr.stats_slopeNKPC <- data.frame(country = geo,
                                      specification = c("unemployment rate", 
                                                        "log(GDP)", "GDP growth"),
                                      posterior_mean = NA,
                                      posterior_sd = NA,
                                      quartile_low = NA,
                                      quartile_high = NA
  )
  
  # du berechnest spaltenweise die summary statistics in postr_slopeNKPC
  # und befÃ¼llst mit den Ergebnissen spaltenweiseweise postr.stats_slopeNKPC
  
  postr.stats_slopeNKPC$posterior_mean <- round(colMeans(postr_slopeNKPC), 3)
  postr.stats_slopeNKPC$posterior_sd <- round(apply(postr_slopeNKPC, 2, sd), 3)
  
  temp <- (length(postr.stats_slopeNKPC)-1):length(postr.stats_slopeNKPC)
  
  postr.stats_slopeNKPC[, temp] <- t(round(
    apply(postr_slopeNKPC, 2, quantile, prob = c(0.025,0.975)), 3))
  
  postr.stats_slopeNKPC <- postr.stats_slopeNKPC |> 
    mutate(
      "credible_interval_[95%]" = paste0("[", round(quartile_low, 3), ", ", 
                                         round(quartile_high, 3), "]")
    ) |> 
    mutate(
      quartile_low = NULL, 
      quartile_high = NULL
    )
  
  postr.stats_slopeNKPC_list[[geo]] <- postr.stats_slopeNKPC
  
  # writexl::write_xlsx(postr.stats_slopeNKPC, 
  #                     paste0("4_results/quarterly_data/tables/1_slope_NKPC/",
  #                            geo, "_postr_stats_slopeNKPC.xlsx")
  #                     )
  
  #---------------------------------------------------------------#
  ##### posterior analysis: coefficient of expected inflation #####
  #---------------------------------------------------------------#
  
  postr.stats_coeffiExpInfl <- data.frame(country = geo,
                                          specification = c("unemployment rate", 
                                                            "log(GDP)", "GDP growth"),
                                          posterior_mean = NA,
                                          posterior_sd = NA,
                                          quartile_low = NA,
                                          quartile_high = NA
  )
  
  postr.stats_coeffiExpInfl$posterior_mean <- round(
    colMeans(postr_coeffiExpInfl), 3)
  
  postr.stats_coeffiExpInfl$posterior_sd <- round(
    apply(postr_coeffiExpInfl, 2, sd), 3)
  
  postr.stats_coeffiExpInfl[, temp] <- t(round(
    apply(postr_coeffiExpInfl, 2, quantile, prob = c(0.025,0.975)), 3))
  
  postr.stats_coeffiExpInfl <- postr.stats_coeffiExpInfl |> 
    mutate(
      "credible_interval_[95%]" = paste0("[", round(quartile_low, 3), ", ", 
                                         round(quartile_high, 3), "]")
    ) |> 
    mutate(
      quartile_low = NULL, 
      quartile_high = NULL
    )
  
  
  postr.stats_coeffiExpInfl_list[[geo]] <- postr.stats_coeffiExpInfl
  
  # writexl::write_xlsx(postr.stats_coeffiExpInfl, 
  #                     paste0("4_results/quarterly_data/tables/2_coeff_exp_infl/",
  #                            geo, "_postr_stats_coeffiExpInfl.xlsx")
  #                     )
  
  #---------------------------------------------------------------#
  ##### posterior analysis: coefficient of the error variance #####
  #---------------------------------------------------------------#
  
  postr.stats_sigma_sq <- data.frame(country = geo,
                                     specification = c("unemployment rate", 
                                                       "log(GDP)", "GDP growth"),
                                     posterior_mean = NA,
                                     posterior_sd = NA,
                                     quartile_low = NA,
                                     quartile_high = NA
  )
  
  postr.stats_sigma_sq$posterior_mean <- round(
    colMeans(postr_sigma_sq), 3)
  
  postr.stats_sigma_sq$posterior_sd <- round(
    apply(postr_sigma_sq, 2, sd), 3)
  
  postr.stats_sigma_sq[, temp] <- t(round(
    apply(postr_sigma_sq, 2, quantile, prob = c(0.025,0.975)), 3))
  
  postr.stats_sigma_sq <- postr.stats_sigma_sq |> 
    mutate(
      "credible_interval_[95%]" = paste0("[", round(quartile_low, 3), ", ", 
                                         round(quartile_high, 3), "]")
    ) |> 
    mutate(
      quartile_low = NULL, 
      quartile_high = NULL
    )
  
  postr.stats_sigma_sq_list[[geo]] <- postr.stats_sigma_sq
  
  # writexl::write_xlsx(postr.stats_sigma_sq, 
  #                     paste0("4_results/quarterly_data/tables/3_sigma_sq/",
  #                            geo, "_postr_stats_sigma_sq.xlsx")
  #                     )
  # 
  #----------------------------------------------#
  #####  posterior distribution: slope NKPC  #####
  #----------------------------------------------#
  
  for(measure in model_i) {
    
    postr_slopeNKPC[measure] |>  # head()
      rename(plot_model = {{measure}}) |> # head()
      ggplot(aes(x = plot_model)) +
      geom_density() +
      geom_vline(aes(
        xintercept = coefficients(ols_models_list[[{{measure}}]])[2]
      ), linetype = "dotdash"
      ) +
      theme_bw() + # theme_minimal
      labs(x         = "slope of the NKPC", y = "density",
           title     = paste0(geo, ", posterior distribution of the slope of the NKPC"),
           subtitle  = paste0("measure of economic activity: ", measure) # ,
           #    caption   = "the dashed line corresponds to the OLS point estimate, 
           # which coincides with the posterior mean"
      )
    
    specification <- unlist(strsplit(measure, split = "_"))[1:2] |> 
      paste(collapse = "_")
    
    ggsave(paste0("4_results/quarterly_data/plots/1_posterior_density_slopeNKPC/", 
                  geo, "_slopeNKPC_", specification, ".png"),
           width = 20, height = 12, units = "cm")
  }
  
  
  postr_slopeNKPC |>
    pivot_longer(
      cols = everything(),
      names_to = "specifications",
      values_to = "slopeNKPC"
    ) |>
    ggplot(aes(x = slopeNKPC, y = specifications)) +
    geom_boxplot() + # (aes(colour = specifications)) +
    theme_bw() +
    labs(
      x = "Slope of the NKPC", 
      y = "measures of economic activity",
      title = paste0(geo, ", Slope of the NKPC"),
      subtitle = "Posterior distribution across all specifications"
    ) + 
    scale_y_discrete(label = c("GDP Growth YoY\nAdj. Quart.",
                               "log(GDP) Chain\nQuart.",
                               # "Unemployment Gap\nAdj. Quart.",
                               "Unemployment Rate\nAdj. Quart.")
    ) # + theme(axis.text.y = element_text(hjust = 1)) # [0,1]
  
  ggsave(paste0("4_results/quarterly_data/plots/1_boxplots_slopeNKPC/", 
                geo, "_slopeNKPC.png"),
         width = 20, height = 12, units = "cm"
  )
  
  
  #---------------------------------------------------------------------#
  #####  posterior distribution: coefficient of expected inflation  #####
  #---------------------------------------------------------------------#
  
  for(measure in model_i){
    
    postr_coeffiExpInfl[measure] |>  # head()
      rename(plot_model = {{measure}}) |> # head()
      ggplot(aes(x = plot_model)) +
      geom_density() +
      geom_vline(aes(
        xintercept = coefficients(ols_models_list[[{{measure}}]])[1]),
        linetype = "dotdash"
      ) +
      theme_bw() + # theme_minimal
      labs(x         = "coefficient expected inflation", y = "density",
           title     = paste0(geo, ", posterior distribution of the coefficient for expected inflation"),
           subtitle  = paste0("measure of economic activity: ", measure) # ,
           #    caption   = "the dashed line corresponds to the OLS point estimate,
           # which coincides with the posterior mean"
      )
    
    specification <- unlist(strsplit(measure, split = "_"))[1:2] |> 
      paste(collapse = "_")
    
    ggsave(paste0("4_results/quarterly_data/plots/2_posterior_density_exp_infl/",
                  geo, "_coeffiExpInfl_", specification, ".png"),
           width = 20, height = 12, units = "cm"
    )
  }
  
  
  postr_coeffiExpInfl |>
    pivot_longer(
      cols = everything(),
      names_to = "specifications",
      values_to = "expected_inflation"
    ) |>
    ggplot(aes(x = expected_inflation, y = specifications)) +
    geom_boxplot() + # (aes(colour = specifications)) +
    theme_bw() +
    labs(
      x = "coefficient of expected inflation", 
      y = "Measures of economic activity",
      title = paste0(geo, ", coefficient on expected inflation"),
      subtitle = "Posterior distribution across all specifications"
    ) + 
    scale_y_discrete(label = c("GDP Growth YoY\nAdj. Quart.",
                               "log(GDP) Chain\nQuart.",
                               # "Unemployment Gap\nAdj. Quart.",
                               "Unemployment Rate\nAdj. Quart.")
    ) # + theme(axis.text.y = element_text(hjust = 1)) # [0,1]
  
  ggsave(paste0("4_results/quarterly_data/plots/2_boxplots_exp_infl/", 
                geo, "_coeffiExpInfl.png"),
         width = 20, height = 12, units = "cm"
  )
  
  
  #----------------------------------------------------#
  #####  posterior distribution: error variance  #####
  #----------------------------------------------------#
  
  # head(y - fitted(ols_reg) == residuals(ols_reg))
  # manually computing the OLS estimate of the error variance sigma sq. so that
  # I can plot it on the posterior distributions to cross-check
  # write the formular from macroetrics here
  # sig2_OLS <- t(residuals(ols_reg)) %*% residuals(ols_reg) / (N - ncol(X))
  
  for(measure in model_i){
    
    ols_residuals <- residuals(ols_models_list[[measure]])
    sigma_sq.ols <- t(ols_residuals) %*% ols_residuals / (N - ncol(X))
    
    postr_sigma_sq[measure] |>  # head()
      rename(plot_model = {{measure}}) |> # head()
      ggplot(aes(x = plot_model)) +
      geom_density() +
      geom_vline(aes(xintercept = sigma_sq.ols), linetype = "dotdash") +
      theme_bw() + # theme_minimal
      labs(x         = "sigma square", y = "density",
           title     = paste0(geo, ", posterior distribution of the error variance"),
           subtitle  = paste0("measure of economic activity: ", measure) # ,
           #    caption   = "the dashed line corresponds to the OLS point estimate,
           # which coincides with the posterior mean"
      )
    
    measure <- unlist(strsplit(measure, split = "_"))[1:2] |>
      paste(collapse = "_")
    
    ggsave(paste0("4_results/quarterly_data/plots/3_posterior_density_sigma_sq/",
                  geo, "_sigma_sq_", measure, ".png"),
           width = 20, height = 12, units = "cm"
    )
    
  }
  
  postr_sigma_sq |>
    pivot_longer(
      cols = everything(),
      names_to = "specifications",
      values_to = "sigma_sq"
    ) |>
    ggplot(aes(x = sigma_sq, y = specifications)) +
    geom_boxplot() + # (aes(colour = specifications)) +
    theme_bw() +
    labs(
      x = "sigma_sq",
      y = "Measures of economic activity",
      title = paste0(geo, ", error variance"),
      subtitle = "Posterior distribution across all specifications"
    ) +
    scale_y_discrete(label = c("GDP Growth YoY\nAdj. Quart.",
                               "log(GDP) Chain\nQuart.",
                               # "Unemployment Gap\nAdj. Quart.",
                               "Unemployment Rate\nAdj. Quart.")
    ) # + theme(axis.text.y = element_text(hjust = 1)) # [0,1]
  
  ggsave(paste0("4_results/quarterly_data/plots/3_boxplots_sigma_sq/",
                geo, "_sigma_sq.png"),
         width = 20, height = 12, units = "cm"
  )
  
  #----------------------------#
  ##### summary statistics #####
  #----------------------------#
  
  # for(geo in country_names_ecb) {
  
  # data_NKPC[data_NKPC$country == geo, 3:8] |>
  #   gt_plt_summary()
  
  summary_stats[[geo]] <- data_NKPC |>
    filter(country == geo) |>
    select(!(country:period)) |>
    gtExtras::gt_plt_summary(title = paste0(geo, ", Summary Statistics"))
  
  # }
  
} # end of the loop estimating all models for each country at once

summary_stats$EA_11 <- data_NKPC |>
  select(!(country:period)) |>
  gtExtras::gt_plt_summary(title = "EA-11, Summary Statistics")

summary_stats

# #-------------------------------------------------#
# ##### outputs: posterior statistics #####
# #-------------------------------------------------#

postr.stats_slopeNKPC <- do.call(rbind, postr.stats_slopeNKPC_list) |> 
  arrange(country)

postr.stats_coeffiExpInfl <- do.call(rbind, postr.stats_coeffiExpInfl_list) |> 
  arrange(country)

postr.stats_sigma_sq <- do.call(rbind, postr.stats_sigma_sq_list) |> 
  arrange(country)

country_names <- c("Austria", "Belgium", "Germany", "Spain", "Finland",
                   "France", "Ireland", "Italy", "Luxembourg", "Netherlands", 
                   "Portugal"
)

for(i in seq_along(country_names_ecb)) {
  postr.stats_slopeNKPC$country <- gsub(country_names_ecb[i], 
                                        country_names[i], 
                                        postr.stats_slopeNKPC$country
  )
  
  postr.stats_coeffiExpInfl$country <- gsub(country_names_ecb[i], 
                                            country_names[i], 
                                            postr.stats_coeffiExpInfl$country
  )
  
  postr.stats_sigma_sq$country <- gsub(country_names_ecb[i], 
                                       country_names[i],
                                       postr.stats_sigma_sq$country
  )
}


postr.stats_slopeNKPC <- postr.stats_slopeNKPC |> 
  arrange(country) |> 
  rbind(EA.postr.stats_slopeNKPC)

postr.stats_slopeNKPC_gt <- postr.stats_slopeNKPC |> 
  gt::gt(groupname_col = "country", rowname_col = "specification") |> 
  tab_header(title = "Posterior statistics of the slope of the NKPC") |>
  # cols_label("95%-credible_interval" = "credible_interval [95%]") |> 
  gt::tab_options(
    data_row.padding = px(1),
    row_group.padding = px(2)
  ) |> 
  cols_align(align = "right", columns = "credible_interval_[95%]") |> 
  # gt::summary_rows(
  #   groups = TRUE,
  #   fns = list(
  #     "Minimum" = ~min(.),
  #     "Maximum" = ~max(.)
  #   ),
  #   formatter = fmt_number,
  #   decimals = 3
  # ) |> 
  gt::opt_stylize(style = 6, color = "gray")

gtsave(postr.stats_slopeNKPC_gt,
       "9_writing_process/postr.stats_slopeNKPC.tex")

# gtsave(postr.stats_slopeNKPC_gt, 
#        "4_results/quarterly_data/tables/1_slope_NKPC/postr.stats_slopeNKPC.png",
#        expand = 10)

writexl::write_xlsx(postr.stats_slopeNKPC,
                    "4_results/quarterly_data/tables/1_slope_NKPC/postr_stats_slopeNKPC.xlsx")



postr.stats_coeffiExpInfl <- postr.stats_coeffiExpInfl |> 
  arrange(country) |> 
  rbind(EA.postr.stats_coeffiExpInfl)

postr.stats_coeffiExpInfl_gt <- postr.stats_coeffiExpInfl |> 
  gt::gt(groupname_col = "country", rowname_col = "specification") |> 
  tab_header(title = "Posterior statistics of the coefficient on expected inflation") |> 
  gt::tab_options(
    data_row.padding = px(1),
    row_group.padding = px(2)
  ) |> 
  cols_align(align = "right", columns = "credible_interval_[95%]") |> 
  # gt::summary_rows(
  #   groups = TRUE,
  #   fns = list(
  #     "Minimum" = ~min(.),
  #     "Maximum" = ~max(.)
  #   ),
  #   formatter = fmt_number,
  #   decimals = 3
  # ) |> 
  gt::opt_stylize(style = 6, color = "gray")


gtsave(postr.stats_coeffiExpInfl_gt,
       "9_writing_process/postr.stats_coeffiExpInfl.tex")

# gtsave(postr.stats_coeffiExpInfl, 
#        "4_results/quarterly_data/tables/2_coeff_exp_infl/postr.stats_coeffiExpInfl.png",
#        expand = 10)

writexl::write_xlsx(postr.stats_coeffiExpInfl,
                    "4_results/quarterly_data/tables/2_coeff_exp_infl/postr.stats_coeffiExpInfl.xlsx")


postr.stats_sigma_sq <- postr.stats_sigma_sq |> 
  arrange(country) |> 
  rbind(EA.postr.stats_sigma_sq)

postr.stats_sigma_sq |> 
  gt::gt(groupname_col = "country", rowname_col = "specification") |> 
  tab_header(title = "Posterior statistics of the error variance") |> 
  gt::tab_options(
    data_row.padding = px(1),
    row_group.padding = px(2)
  ) |>
  cols_align(align = "right", columns = "credible_interval_[95%]") |> 
  # gt::summary_rows(
  #   groups = TRUE,
  #   fns = list(
  #     "Minimum" = ~min(.),
  #     "Maximum" = ~max(.)
  #   ),
  #   formatter = fmt_number,
  #   decimals = 3
  # ) |> 
  gt::opt_stylize(style = 6, color = "gray")


# gtsave(postr.stats_sigma_sq, 
#        "4_results/quarterly_data/tables/3_sigma_sq/postr.stats_sigma_sq.tex")

# gtsave(postr.stats_sigma_sq, 
#        "4_results/quarterly_data/tables/3_sigma_sq/postr.stats_sigma_sq.png",
#        expand = 10)

writexl::write_xlsx(postr.stats_sigma_sq,
                    "4_results/quarterly_data/tables/3_sigma_sq/postr.stats_sigma_sq.xlsx")



# to measure the duration of the code
print(Sys.time() - start.time)


# #-------------------------------------------------#
# ##### stargazer outputs: posterior statistics #####
# #-------------------------------------------------#
# 
# a <- stargazer::stargazer(postr.stats_slopeNKPC,
#           title = "Slope of the NKPC, AT",
#           summary = FALSE)
# 
# stargazer(postr.stats_coeffiExpInfl, 
#           title = "Expected Inflation, AT", 
#           summary = FALSE)
# 
# stargazer(postr.stats_sigma_sq, 
#           title = "Error Variance, AT", 
#           summary = FALSE)
# 
# summary statistics
# summary_stats <- as.data.frame(
#   data_NKPC[data_NKPC$country == geo, c(1, 3:8)]
#   )
# stargazer(summary_stats, title = "Summary Statistics AT")

data_description <- readxl::read_excel("9_writing_process/data.xlsx")
# , col_types = "character")

stargazer::stargazer(data_description, 
                     title = "Variable Description & Sources",
                     summary = FALSE)

##### OLS results #####
for(i in country_names_ecb) {
  print(i)
  for(model in model_i) {
    print(summary(ols_geo_list[[i]][[model]]))
  }
}

# #------------------------------------------#
# #####  stargazer outputs: OLS results  #####
# #------------------------------------------#

# summary(ols_models_list[["unemp_rate_adj"]])


stargazer(EA.ols_models_list, title = "EA-11 OLS Results", align = TRUE,
          dep.var.labels = "core inflation",
          covariate.labels = c("expected inflation", "unemployment rate",
                               "log(GDP)", "GDP growth"),
          no.space = TRUE
)

for(i in country_names_ecb) {
  
  stargazer(ols_geo_list[[i]], title = paste0(i, ", OLS Results"), 
            align = TRUE, dep.var.labels = "core inflation",
            covariate.labels = c("expected inflation", "unemployment rate",
                                 "log(GDP)", "GDP growth"),
            no.space = TRUE
  )
  
  # 
  #   print(i)
  #   for(model in model_i) {
  #     print(summary(ols_geo_list[[i]][[model]]))
  #   }
}

# c("expected inflation", "unemployment rate",
#   "unemployment gap", "log(GDP)", "log(GDP) gap",
#   "GDP growth")

##### checking if smaller economies have a larger slope coefficient than bigger #####
# economies.
# Large economies: Germany (24%), France(17%), Italy(12%), ~Spain(8%), (Netherlands, 6%)
# Small economies: Belgium (4%), Austria (3%), Ireland (3%), Finland (2%), Portugal (2%), Luxembourg
# Percentage refers to how much they contribute to the EU's gross domestic product.
# Source: https://www.weforum.org/agenda/2023/02/eu-countries-largest-economies-energy-gdp/ 

average_postr_mean <- mean(postr.stats_slopeNKPC$posterior_mean)

postr.stats_slopeNKPC[postr.stats_slopeNKPC$posterior_mean > average_postr_mean, 
                      c("country", "specification", "posterior_mean")]#, drop = FALSE]


postr.stats_slopeNKPC[postr.stats_slopeNKPC$posterior_mean <= average_postr_mean, 
                      c("country", "specification", "posterior_mean")]#, drop = FALSE]

range(postr.stats_slopeNKPC$posterior_mean[postr.stats_slopeNKPC$posterior_mean > 0.018])











