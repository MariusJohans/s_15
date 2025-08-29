
##=============================================================================
## 
## Artificial intelligence, health empowerment, and the 
## general practitioner system.
##
## Marius O. Johansen, 2025
##
## This code was written in  RStudio 2022.07.1+554 "Spotted Wakerobin"
## release. Functions and package compatibility may be subject to
## change.
##
##=============================================================================

library(foreign)
library(readr)
library(lavaan)
library(semPlot)
library(psych)
library(apaTables)
library(memisc)
library(car)
library(readxl)
library(knitr)
library(ggpubr)
library(e1071)
library(effsize)
library(ltm)
library(mice)
library(naniar)
library(finalfit)
library(writexl)
library(semTools)
library(MASS)
library(simsem)
library(simr)
library(ordinal)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)

Data <- read_excel("Data.xlsx", col_names = TRUE)
Data
Data$SR_T_num <- as.numeric(gsub("SR", "", Data$SR_T))

##=============================================================================
##
## Typical survival rate plot (Bergen)
##
##=============================================================================

simulate_relative_survival <- function(t, S0, t0, r = 0.07, noise_sd = 0.01) {
  delta_t <- t - t0
  base <- S0 * (1 - r)^delta_t
  noise <- rnorm(length(t), mean = 0, sd = noise_sd)
  survival <- base + noise
  pmin(pmax(survival, 0), 1)
}
t_vals <- seq(0, 30, by = 0.2)
bergen_curve <- data.frame(
  Municipality = "Bergen",
  t = t_vals,
  SR = simulate_relative_survival(t_vals, S0 = 0.477, t0 = 10.2, r = 0.07, 
                                  noise_sd = 0.001)
)
ggplot(bergen_curve, aes(x = t, y = SR)) +
  geom_line(color = "#0072B2", linewidth = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                     limits = c(0, 1)) +
  labs(
    ## title = "Standard survival rate plot",
    x = "Response time (minutes)",
    y = "Simulated survival rate"
  ) +
  theme_minimal(base_size = 14)

##=============================================================================
##
## Monte Carlo simulations
##
##=============================================================================

## Funksjon for aa simulere en kurve
simulate_one_curve <- function(S0, r = 0.07, noise_sd = 0.01, 
                               min_gain_max = 10) {
  minutes_gained <- seq(0, min_gain_max, by = 0.1)
  survival <- S0 * (1 / ((1 - r)^minutes_gained))
  noise <- rnorm(length(minutes_gained), mean = 0, sd = noise_sd)
  noisy_survival <- pmin(pmax(survival + noise, 0), 1)
  data.frame(
    MinutesGained = minutes_gained,
    SR = noisy_survival
    )
  }
## Monte Carlo: kjor mange simuleringer
simulate_mc <- function(S0, r = 0.07, noise_sd = 0.01, min_gain_max = 10, 
                        n_iter = 10000) {
  minutes_gained <- seq(0, min_gain_max, by = 0.1)
  ## kjor alle iterasjoner
  sims <- replicate(
    n_iter,
    simulate_one_curve(S0, r, noise_sd, min_gain_max)$SR,
    simplify = "matrix"
    )
  ## summer opp: gjennomsnitt + CI
  results <- data.frame(
    MinutesGained = minutes_gained,
    Mean = rowMeans(sims),
    Lower = apply(sims, 1, quantile, probs = 0.025),
    Upper = apply(sims, 1, quantile, probs = 0.975)
    )
  return(results)
}

set.seed(42)
bergen_mc <- simulate_mc(S0 = 0.477, r = 0.07,  noise_sd = 0.015, 
                         n_iter = 10000)
tokke_mc  <- simulate_mc(S0 = 0.2998, r = 0.07, noise_sd = 0.015,
                         n_iter = 10000)
luroy_mc  <- simulate_mc(S0 = 0.0925, r = 0.07, noise_sd = 0.015,
                         n_iter = 10000)
sorfo_mc  <- simulate_mc(S0 = 0.1939, r = 0.07, noise_sd = 0.015,
                         n_iter = 10000)
## Bergen
p1 <- ggplot(bergen_mc, aes(x = MinutesGained, y = Mean)) +
  geom_line(color = "#0072B2", linewidth = 1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "#0072B2") +
  scale_x_continuous(
    breaks = c(0, 1, 5, 10),
    labels = c("Current", "-1 min", "-5 min", "-10 min")
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = 
                       c(0, 1)) +
  labs(
    ## title = "Monte Carlo simulated survival curve (Bergen, n = 10,000)",
    title = "Bergen",
    x = "Earlier ambulance arrival",
    y = "Simulated survival rate"
  ) +
  theme_minimal(base_size = 14)

## Tokke
p2 <- ggplot(tokke_mc, aes(x = MinutesGained, y = Mean)) +
  geom_line(color = "grey44", linewidth = 1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "grey44") +
  scale_x_continuous(
    breaks = c(0, 1, 5, 10),
    labels = c("Current", "-1 min", "-5 min", "-10 min")
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = 
                       c(0, 1)) +
  labs(
    title = "Tokke",
    x = "Earlier ambulance arrival",
    y = "Simulated survival rate"
  ) +
  theme_minimal(base_size = 14)

## Luroy
p3 <- ggplot(luroy_mc, aes(x = MinutesGained, y = Mean)) +
  geom_line(color = "darkturquoise", linewidth = 1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "darkturquoise") +
  scale_x_continuous(
    breaks = c(0, 1, 5, 10),
    labels = c("Current", "-1 min", "-5 min", "-10 min")
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = 
                       c(0, 1)) +
  labs(
    title = "Luroy",
    x = "Earlier ambulance arrival",
    y = "Simulated survival rate"
  ) +
  theme_minimal(base_size = 14)

## Sorfold
p4 <- ggplot(sorfo_mc, aes(x = MinutesGained, y = Mean)) +
  geom_line(color = "darkseagreen4", linewidth = 1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "darkseagreen4") +
  scale_x_continuous(
    breaks = c(0, 1, 5, 10),
    labels = c("Current", "-1 min", "-5 min", "-10 min")
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = 
                       c(0, 1)) +
  labs(
    title = "Sorfold",
    x = "Earlier ambulance arrival",
    y = "Simulated survival rate"
  ) +
  theme_minimal(base_size = 14)

combined <- (p1 | p2) / (p3 | p4)
combined

library(spData)
install.packages('spDataLarge',
                 repos='https://nowosad.github.io/drat/', type='source')
library(spDataLarge)

library(spDataLarge)
class(world)
world$geometry[1]
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")  # "medium" gives reasonable detail



ggplot(data = world) +
  geom_sf() +
  xlab("") + ylab("") +
  coord_sf(xlim = c(2, 33), ylim = c(57, 72), expand = FALSE) +
  geom_point(aes(x = 5.320, y = 60.386), colour = "red", size = 3) +
  geom_point(aes(x = 5.363, y = 60.378), colour = "red", size = 3) +
  geom_point(aes(x = 6.087, y = 62.147), colour = "red", size = 3) +
  geom_point(aes(x = 5.326, y = 60.408), colour = "red", size = 3)




library(readxl)
library(ggplot2)
library(dplyr)
library(forcats)
library(nnet)
library(lmtest)
library(generalhoslem)
library(terra)
library(tmap)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(reshape) 
library(viridis)
library(maps)
library(spData)
library(tidyverse)
library(ggspatial)
library(gganimate)
library(gifski)
library(DescTools)
library(rcompanion)
library(vcd)
library(gtools)
library(reshape2)


##=============================================================================
##
## END
##
##=============================================================================