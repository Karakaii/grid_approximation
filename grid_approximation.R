# Temporary place to update my grid approximation function

#~############################################################################~#
# Grid Approximation function ----
#~############################################################################~#
# This is a swiss-knife function to do grid approximation with all sorts of
# distributions (adding them as I need them).

#~=======================================================~=
## Function parameters ----
#~=======================================================~=
grid_approximation <- function(
    # Grid information
  grid_length = 100000,
  grid_min,
  grid_max,
  # Prior and likelihood information
  prior,
  prior_type,
  likelihood,
  likelihood_type,
  # Graph information
  central_tendency = NULL,
  plot_x_lim = NULL,
  plot_x_lab = "unit",
  bayes_colors = c(
    "prior"="#ff5757", 
    "new data (likelihood)"="#4acdf1", 
    "posterior"="#6a0dad"
  ),
  custom_labels = NULL,
  seed = NULL
) {
  
  # Run a potential seed
  if(!is.null(seed)){set.seed(seed)}
  
  #~=======================================================~=
  ## Building the grid ----
  #~=======================================================~=
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ### Grid ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Building the grid
  # grid_n is the number of splits
  # need grid_n for the dnorm generation
  grid_data <- data.frame(
    grid_n = seq(from = grid_min, to = grid_max, length = grid_length)
  ) 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ### Disclaimer ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Can't use ifelse or case_when, both have issues
  # ifelse provides only one number
  # case_when calculates everything and then does the conditional, which means 
  # it breaks unnecessarily
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ### Prior ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # For the prior
  if(prior_type == "dnorm") {
    grid_data <- grid_data %>% mutate(
      prior = 
        dnorm(
          x    = grid_n,
          mean = prior[1],
          sd   = prior[2]
        )
    )
  }
  if(prior_type == "dbeta") {
    grid_data <- grid_data %>% mutate(
      prior = 
        dbeta(
          x    = grid_n,
          shape1 = prior[1],
          shape2 = prior[2]
        )
    )
  }
  if(prior_type == "data") {
    grid_data$prior <- density(
      prior,
      from = grid_min, to = grid_max,
      n = grid_length
    )$y
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ### Likelihood ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # For the likelihood
  if(likelihood_type == "dnorm") {
    grid_data <- grid_data %>% mutate(
      likelihood = 
        dnorm(
          x    = grid_n,
          mean = likelihood[1],
          sd   = likelihood[2]
        )
    )
  }
  if(likelihood_type == "dbeta") {
    grid_data <- grid_data %>% mutate(
      likelihood = 
        dbeta(
          x    = grid_n,
          shape1 = likelihood[1],
          shape2 = likelihood[2]
        )
    )
  }
  if(likelihood_type == "data") {
    grid_data$likelihood <- density(
      likelihood,
      from = grid_min, to = grid_max,
      n = grid_length
    )$y
  }
  
  #~=======================================================~=
  ## Calculating the posterior ----
  #~=======================================================~=
  # Then use bayes rule to calculate the posterior
  grid_data <- grid_data %>% mutate(
    unnormalised = likelihood * prior,
    posterior = unnormalised / sum(unnormalised)
  )
  # Confirm that the posterior approximation sums to 1
  sum(grid_data$posterior)
  
  # We then sample from the grid,
  # Here, the grid_n represents the posterior (the values sampled, based
  # on the weights given by the posterior)
  # So we extract it and rename it
  posterior_sample <- grid_data %>% 
    sample_n(size = grid_length, weight = posterior, replace = T) %>%
    select(grid_n) %>% 
    rename(posterior = grid_n)
  
  #~=======================================================~=
  ## Making the graph ----
  #~=======================================================~=
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ### Getting the data ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  # We prepare the data for the plot by extracting the data
  # for the prior, data, and posterior
  plot_data <- grid_data %>% 
    # Need to normalise the prior and likelihood for it all to look proportional
    mutate(
      prior = prior / sum(prior),
      likelihood = likelihood / sum(likelihood)
    ) %>% select(-unnormalised) %>% 
    # Pivot for the graph
    pivot_longer(
      cols = -grid_n,
      names_to = "distribution",
      values_to = "value"
    ) %>% 
    # Make likelihood easier to understand
    mutate(distribution = ifelse(
      distribution == "likelihood", 
      "new data (likelihood)", 
      distribution)
    )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ### Preparing the central tendencies ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Skip if central_tendency is NULL
  if(!is.null(central_tendency)) {
  
    # If there's only one central tendency specified, multiply it
    if(length(central_tendency) == 1) {
      central_tendency <- rep(central_tendency, 3)
    }
    
    #### Prior ----
    if(is.numeric(central_tendency[[1]])) {
      ct_prior <- central_tendency[[1]]
    } else {
      ct_prior <- case_when(
        central_tendency[[1]] == "median" ~ 
          case_when(
            prior_type == "dnorm" ~ prior[1],
            prior_type == "dbeta" ~ prior[1]/(prior[1]+prior[2]),
            prior_type == "data"  ~ median(prior)
          ),
        central_tendency[[1]] == "mean" ~       
          case_when(
            prior_type == "dnorm" ~ prior[1],
            prior_type == "dbeta" ~ qbeta(0.5, prior[1], prior[2]),
            prior_type == "data"  ~ mean(prior)
          )
      )
    }
    
    #### Likelihood ----
    if(is.numeric(central_tendency[[2]])) {
      ct_data <- central_tendency[[2]]
    } else {
      ct_data <- case_when(
        central_tendency[[2]] == "median" ~ 
          case_when(
            likelihood_type == "dnorm" ~ likelihood[1],
            likelihood_type == "dbeta" ~ likelihood[1]/
              (likelihood[1]+likelihood[2]),
            likelihood_type == "data"  ~ median(likelihood)
          ),
        central_tendency[[2]] == "mean" ~       
          case_when(
            likelihood_type == "dnorm" ~ likelihood[1],
            likelihood_type == "dbeta" ~ qbeta(0.5, likelihood[1], likelihood[2]),
            likelihood_type == "data"  ~ mean(likelihood)
          )
      )
    }
    
    #### Posterior ----
    if(is.numeric(central_tendency[[3]])) {
      ct_posterior <- central_tendency[[3]]
    } else {
      ct_posterior <- case_when(
        central_tendency[[3]] == "median" ~ median(posterior_sample$posterior),
        central_tendency[[3]] == "mean" ~ mean(posterior_sample$posterior)
      )
    }
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ### The graph ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Determine if there are custom labels for the distribution
  distribution_labels <- names(bayes_colors)
  if(!is.null(custom_labels)) {distribution_labels <- custom_labels}
  
  bayes_graph <- 
    ggplot(plot_data, aes(x = grid_n, y = value, color = distribution)) +
    
    # graph the distributions
    geom_line() +
    
    # General presentation
    xlab(plot_x_lab) + 
    ylab("density") +
    cowplot::theme_cowplot()  +
    theme(
      # Below specifies details of legends appearance
      legend.position      = c(1, 1),
      legend.justification = c("right", "top"),
      legend.margin        = margin(6, 6, 6, 6),
      legend.text          = element_text(size = 10),
      legend.text.align    = 0,
      legend.title         = element_text(size=11, face = "bold"),
      legend.key           = element_rect(linewidth = 0.4),
      legend.background    = element_blank(),
      legend.box.just      = "left",
      legend.box.background= element_rect(colour = "black"),
      # Modify the y axis
      axis.line.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y  = element_blank()
    ) + 
    scale_color_manual(
      name  = "Distributions", 
      values = bayes_colors,
      breaks = names(bayes_colors),
      labels = distribution_labels
      )
  
  if(!is.null(plot_x_lim)) {
    bayes_graph <- bayes_graph + coord_cartesian(xlim = plot_x_lim)
  }

  if(!is.null(central_tendency)) {
    bayes_graph <- bayes_graph +
      geom_vline(
        xintercept = ct_prior,
        linetype = "dashed",
        color = bayes_colors["prior"]
      ) +
      geom_vline(
        xintercept = ct_data,
        linetype = "dashed",
        color = bayes_colors["new data (likelihood)"]
      ) +
      geom_vline(
        xintercept = ct_posterior,
        linetype = "dashed",
        color = bayes_colors["posterior"]
      )
  }
  
  #~=======================================================~=
  ## Returning items ----
  #~=======================================================~=
  ## Return the graph and th posterior in a list
  return(list("posterior" = posterior_sample$posterior, "graph" = bayes_graph))
}

#~############################################################################~#
# Normal-normal ----
#~############################################################################~#

combine_normals <- function(mu1, sigma1, mu2, sigma2) {
  sigma1_sq <- sigma1^2
  sigma2_sq <- sigma2^2
  # All these different versions work
  mu_post <- (1/((1/sigma1_sq) + (1/sigma2_sq))) * ((mu1/sigma1_sq)+(mu2/sigma2_sq))
  mu_post <- mu2 * (sigma1_sq/(sigma1_sq+sigma2_sq)) + mu1 * (sigma2_sq/(sigma1_sq+sigma2_sq))
  mu_post <- (sigma1_sq * mu2 + sigma2_sq * mu1) / (sigma1_sq + sigma2_sq)
  sigma_post <- sqrt(1 / (1 / sigma1_sq + 1 / sigma2_sq))
  
  ci_lb <- mu_post - sigma_post*1.96
  ci_ub <- mu_post + sigma_post*1.96
  
  return(data.frame(mu_post, sigma_post, ci_lb, ci_ub))
}
