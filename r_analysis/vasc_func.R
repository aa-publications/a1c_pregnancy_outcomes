# script with functions for the vascular paper
# start by loading libraries
library(rms)
library(table1)
library(flextable)
library(ggpubr)
library(tidyverse)
library(extrafont)
loadfonts()

# now relevant functions

# make partial effects plot function
pp_outcome <- function(model, title, ybounds) {
  # have to split out t1d and t2d, since t2d estimate is only valid at 6
  df <- Predict(model, a1c_closest_to_supervision, dm_outcome, fun = plogis) %>%
    as_tibble()
  df_diab <- df %>% filter(dm_outcome == "Type 1 Diabetes")
  df_cont <- df %>% filter(dm_outcome == "No Type 1 Diabetes" & 
                             (round(a1c_closest_to_supervision, 2) <= 5.31 & round(a1c_closest_to_supervision, 2) >= 5.29))
  ggplot(df_diab) +
    geom_ribbon(aes(x = a1c_closest_to_supervision,
                    ymin = lower,
                    ymax = upper),
                alpha = 0.25,
                fill = "#D8AB4C") +
    geom_ribbon(aes(x = a1c_closest_to_supervision,
                    ymin = df_cont$lower,
                    ymax = df_cont$upper),
                alpha = 0.25,
                fill = "#000000") +
    geom_line(aes(x = a1c_closest_to_supervision, y = df_cont$yhat, color = "#000000")) +
    geom_line(aes(x = a1c_closest_to_supervision, y = yhat, color = "#D8AB4C")) +
    scale_y_continuous(name = "Adjusted Estimated Probability, %",
                       labels = scales::percent) +
    scale_x_continuous(name = expression("Hemoglobin"~A["1c"]~"near Conception, %")) +
    labs(caption = "", title = title) +
    scale_color_identity(guide = "legend",
                         name = "T1D Status",
                         labels = c("No", "Yes")) +
    coord_cartesian(ylim = ybounds, xlim = c(5, 10)) +
    theme_pubclean() +
    theme(legend.position = "right",
          legend.key = element_rect(fill = NA),
          text = element_text(family = "Arial"),
          plot.title = element_text(face = "bold"))
}

# get ors from models
# now make odds ratios from FP
ors_outcome <- function(model, control, comparators) {
  c <- contrast(model,
                list(a1c_closest_to_supervision = comparators, dm_outcome = "Type 1 Diabetes"),
                list(a1c_closest_to_supervision = control, dm_outcome = "No Type 1 Diabetes"))
  tibble(a1c = comparators,
         or = exp(c$Contrast),
         or_lb = exp(c$Lower),
         or_ub = exp(c$Upper)) %>%
    mutate(`OR (95% CI)` = paste0(round(or, 2), " (",
                                  round(or_lb, 2), ", ",
                                  round(or_ub, 2), ")"),
           A1c = factor(a1c))
}
