This is my own R script trying to streamline the use of grid approximation in simple Bayesian analyses. 
This is based on the works of [Johnson et al., 2021](https://www.bayesrulesbook.com/chapter-6.html#grid-approximation) and [McElreath, 2020](https://xcelab.net/rm/statistical-rethinking/).

Requires cowplot and tidyverse.

To load in your R script use:
```{r}
devtools::source_url("https://raw.githubusercontent.com/Karakaii/grid_approximation/main/grid_approximation.R")
```