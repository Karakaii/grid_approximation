# I'm trying to make a flexible grid approximation function in order to 
# streamline the process of running grid approximation.

# The function takes in information about the grid you create.
# The grid_length represents how many parts the grid is made of.
# The more parts, the higher the resolution.
# This also represents the number of samples taken.
# grid_min and grid_max denote the start and end of the grid.

# The likelihood and the prior are entered in similar ways. An element with
# '_type' indicated how the function is to use the information entered.
# * if '_type' is "dnorm", then the prior or likelihood will need to provide an
#   array with two element, a mean and an SD.
# * if '_type' is "data", then the prior or likelihood will need to provide an
#   array with all the data, for which the density will be extracted using
#   density().

# central_tendency determines whether the mean or median is presented in
# the graph.

# plot_x_lim represent the limits (lower, upper) for the cartesian coordinates of the
# graph.

# The bayes_colors sets a default for the colours that users can change.
# Must have the names "prior", "data", and "posterior"

# The plot_x_lab allows users to set the x lab of the plot.

grid_approximation <- function(
  # Grid information
  grid_length = 10000,
  grid_min,
  grid_max,
  # Prior and likelihood information
  prior,
  prior_type,
  likelihood,
  likelihood_type,
  # Graph information
  central_tendency = "median",
  plot_x_lim,
  plot_x_lab = "unit",
  bayes_colors = c("prior"="#ff5757", "data"="#4acdf1", "posterior"="#d269ff")
) {
  
  # Building the grid
  # grid_n is the number of splits
  grid_data <- data.frame(
    grid_n = seq(from = grid_min, to = grid_max, length = grid_length)
  ) %>% mutate(
    prior = case_when(
      prior_type == "dnorm" ~ dnorm(
        x    = grid_n,
        mean = prior[1],
        sd   = prior[2]
      ),
      prior_type == "data" ~ density(
        prior,
        from = grid_min, to = grid_max,
        n = grid_length
      )$y
    ), 
    likelihood = case_when(
      likelihood_type == "dnorm" ~ dnorm(
        x    = grid_n,
        mean = likelihood[1],
        sd   = likelihood[2]
      ),
      likelihood_type == "data" ~ density(
        likelihood,
        from = grid_min, to = grid_max,
        n = grid_length
      )$y
    )
  ) %>% mutate(
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
  
  # We prepare the data for the plot by extracting the data
  # for the prior, data, and posterior
  plot_data <- rbind(
    grid_data %>% 
      select(grid_n, prior) %>% 
      rename(density = prior) %>%
      mutate(type = "prior", transparence=0.3),
    grid_data %>% 
      select(grid_n, likelihood) %>% 
      rename(density = likelihood) %>%
      mutate(type = "data", transparence=0.4),
    grid_data %>% 
      select(grid_n, unnormalised) %>% 
      rename(density = unnormalised) %>%
      mutate(type = "posterior", transparence=0.5)
  )
  
  # Determine the central tendencies
  ct_prior <- case_when(
    central_tendency == "median" ~ 
      case_when(
        prior_type == "dnorm" ~ prior[1],
        prior_type == "data"  ~ median(prior)
      ),
    central_tendency == "mean" ~       
      case_when(
        prior_type == "dnorm" ~ prior[1],
        prior_type == "data"  ~ mean(prior)
      )
  )
  ct_data <- case_when(
    central_tendency == "median" ~ 
      case_when(
        likelihood_type == "dnorm" ~ prior[1],
        likelihood_type == "data"  ~ median(prior)
      ),
    central_tendency == "mean" ~       
      case_when(
        likelihood_type == "dnorm" ~ prior[1],
        likelihood_type == "data"  ~ mean(prior)
      )
  )
  ct_posterior <- case_when(
    central_tendency == "median" ~ median(posterior_sample$posterior),
    central_tendency == "mean" ~ mean(posterior_sample$posterior)
  )
  
  ## Make the graph
  bayes_graph <- ggplot(plot_data, aes(x = grid_n, y = density, color = type)) +
    
    # graph the distributions
    geom_line() +
    
    # Graph central tendencies
    geom_point(
      y=0, x=ct_prior, size = 3,
      fill = NA, color = bayes_colors["prior"], shape = 21
    ) +
    geom_point(
      y=0, x=ct_data, size = 3,
      fill = NA, color = bayes_colors["data"], shape = 21
    ) +
    geom_point(
      y=0, x=ct_posterior, size = 3,
      fill = NA, color = bayes_colors["posterior"], shape = 21
    ) +
    
    # General presentation
    xlab(plot_x_lab) + ylab("density") +
    coord_cartesian(xlim = plot_x_lim) +
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
      legend.box.background= element_rect(colour = "black")
    ) + scale_color_manual(name  = "Distributions", values = bayes_colors)
  
  ## Return the graph and th posterior in a list
  return(list(posterior_sample$posterior, bayes_graph))
}
