# set working directory
setwd("/Users/sackd/Library/CloudStorage/Box-Box/a1c_pregnancy_outcomes/vascular paper/")

# call functions and libraries
# first call some libraries
source("vasc_func.R")

# set seed for imputations
set.seed(1111)

# now pull in data and update dm_outcome to factor
data <- read_tsv("2021-10-29_predictors_outcomes.tsv", 
                 na = c("", "NA", "NaN")) %>%
  mutate(dm_outcome = factor(dm_outcome, levels = c("no_t1dm", "t1d"),
                             labels = c("No Type 1 Diabetes", "Type 1 Diabetes")),
         Race = factor(Race, levels = c("W", "B", "U", "A", "I"), 
                       labels = c("White", "Black", "Unknown", "Asian", "Native American")), # to improve confidence intervals
         Vascular = ifelse(Vascular == TRUE, 1, 0),
         Eclampsia = ifelse(`908.1` > 0 | `908.11` > 0 | `908.12`, 1, 0)) %>%
  rowwise() %>%
  mutate(a1c_closest_to_supervision = ifelse(dm_outcome == "No Type 1 Diabetes", rnorm(1, 5.3, 0.1), a1c_closest_to_supervision), # set A1c to a normal distribution in folks without T1D
         a1c_closest_to_delivery = ifelse(dm_outcome == "No Type 1 Diabetes", rnorm(1, 5.3, 0.1), a1c_closest_to_delivery)) # used this paper as reference https://care.diabetesjournals.org/content/32/5/828

# check A1c by T1D status
tapply(data$a1c_closest_to_supervision, data$dm_outcome, summary)
tapply(data$a1c_closest_to_delivery, data$dm_outcome, summary)

# check outcome status by dm_status to get sense of number of degrees of freedom to use
tapply(data$Vascular, data$dm_outcome, sum)
tapply(data$Eclampsia, data$dm_outcome, sum)

# add some columns for table 1
data <- data %>%
  mutate(fac.Vascular = factor(Vascular, levels = c(0, 1), labels = c("No", "Yes")),
         fac.Eclampsia = factor(Eclampsia, levels = c(0, 1), labels = c("No", "Yes")))

# add labels for table
label(data$a1c_closest_to_delivery) <- "Hemoglobin A1c Close to Delivery"
label(data$a1c_closest_to_supervision) <- "Hemoglobin A1c Close to Conception"
label(data$age_at_supervision_years) <- "Age"
units(data$age_at_supervision_years) <- "years"
label(data$Race) <- "Race"
label(data$DEP_INDEX) <- "Deprivation Index"
label(data$fac.Vascular) <- "Any Vascular Pathology"
label(data$fac.Eclampsia) <- "Pre-Eclampsia/Eclampsia/HELLP Syndrome"

# make a table of key variables
table1(~ a1c_closest_to_delivery + a1c_closest_to_supervision +
         age_at_supervision_years + Race + DEP_INDEX +
         fac.Vascular + fac.Eclampsia | dm_outcome,
       data = data, 
       render.continuous = c(.="Median [Q1, Q3]",
                             .= "Min, Max")) %>%
  t1flex() %>%
  save_as_docx("Table1", path = "Table1.docx")

# subset to only relevant covariates
data_vasc <- data %>%
  select(GRID, dm_outcome, a1c_closest_to_delivery, a1c_closest_to_supervision, age_at_supervision_years, 
         Race, DEP_INDEX, Vascular, Eclampsia)

# make datadist for plotting
dd <- datadist(data_vasc)
options(datadist = "dd")

# now do imputations
imp <- aregImpute(~ dm_outcome + a1c_closest_to_delivery + a1c_closest_to_supervision + age_at_supervision_years + Race +
                    DEP_INDEX + Vascular + Eclampsia,
                  data = data_vasc,
                  n.impute = 35)

# all vascular model
mod_vas <- fit.mult.impute(Vascular ~ dm_outcome*rcs(a1c_closest_to_supervision, 3) + Race + DEP_INDEX + age_at_supervision_years,
                           lrm, imp, data = data_vasc)
# anova table
anova(mod_vas)

# eclampsia model
mod_ec <- fit.mult.impute(Eclampsia ~ dm_outcome*rcs(a1c_closest_to_supervision, 3) + Race + DEP_INDEX + age_at_supervision_years,
                          lrm, imp, data = data_vasc)
anova(mod_ec)

# set control and comparators for odds ratios
control <- c(5.5)
comparators <- c(5.5, 6.5, 7.5, 8.5)

# create plots
# vascular plots
legend <- get_legend(pp_outcome(mod_vas, "A. Pregnant Patients with Vascular Complication", c(0, 0.3)) +
                       theme(legend.position = "bottom"))
fig4a <- pp_outcome(mod_vas, "A. Pregnant Patients with Vascular Complication", c(0, 0.3)) +
  theme(legend.position = "none")
fig4b <- ors_outcome(mod_vas, control, comparators) %>%
  ggplot() + 
  geom_vline(xintercept = 1) +
  geom_errorbarh(aes(y = fct_rev(A1c), xmin = or_lb, xmax = or_ub), height = 0.1, color = "#D8AB4C") +
  geom_point(aes(x = or, y = fct_rev(A1c)), color = "#D8AB4C") +
  geom_text(aes(x = 0.6, y = fct_rev(A1c), label = `OR (95% CI)`), hjust = 0) +
  annotate(geom = "text", x = 0.6, y = 4.5, label = "OR (95% CI)*",
           hjust = 0, fontface = "bold", family = "Arial") +
  geom_text(aes(x = 0.4, y = fct_rev(A1c), label = `A1c`), hjust = 0) +
  annotate(geom = "text", x = 0.35, y = 4.5, label = expression(bold("Hemoglobin"~A["1c"])),
           hjust = 0, family = "Arial") +
  labs(caption = expression("*Compared to"~A["1c"]~"of 5.5 in pregnant patients without T1D (OR = 1)"),
       title = "B. Odds Ratios of Vascular Complication") +
  ylab("") +
  scale_x_log10(name = "Odds Ratio (95% CI)") +
  theme_pubclean() +
  coord_cartesian(clip = "off") +
  theme(text = element_text(family = "Arial"),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 12),
        plot.title = element_text(face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = "dotted"),
        plot.caption = element_text(hjust = 0))

# eclampsia plots
fig4c <- pp_outcome(mod_vas, "C. Pregnant Patients with Pre-Eclampsia/Eclampsia/HELLP Syndrome", c(0, 0.3)) +
  theme(legend.position = "none")
fig4d <- ors_outcome(mod_ec, control, comparators) %>%
  ggplot() + 
  geom_vline(xintercept = 1) +
  geom_errorbarh(aes(y = fct_rev(A1c), xmin = or_lb, xmax = or_ub), height = 0.1, color = "#D8AB4C") +
  geom_point(aes(x = or, y = fct_rev(A1c)), color = "#D8AB4C") +
  geom_text(aes(x = 1.1, y = fct_rev(A1c), label = `OR (95% CI)`), hjust = 0) +
  annotate(geom = "text", x = 1.1, y = 4.5, label = "OR (95% CI)*",
           hjust = 0, fontface = "bold", family = "Arial") +
  geom_text(aes(x = 0.6, y = fct_rev(A1c), label = `A1c`), hjust = 0) +
  annotate(geom = "text", x = 0.35, y = 4.5, label = expression(bold("Hemoglobin"~A["1c"])),
           hjust = 0, family = "Arial") +
  labs(caption = expression("*Compared to"~A["1c"]~"of 5.5 in pregnant patients without T1D (OR = 1)"),
       title = "D. Odds Ratios of Pre-Eclampsia/Eclampsia/HELLP Syndrome") +
  ylab("") +
  scale_x_log10(name = "Odds Ratio (95% CI)") +
  theme_pubclean() +
  coord_cartesian(clip = "off") +
  theme(text = element_text(family = "Arial"),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 12),
        plot.title = element_text(face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = "dotted"),
        plot.caption = element_text(hjust = 0))

# now combine into one figure with a shared legend
fig4ab <- ggarrange(fig4a, fig4b, nrow = 1, widths = c(1.2, 1), align = "v")
fig4cd <- ggarrange(fig4c, fig4d, nrow = 1, widths = c(1.2, 1), align = "v")

# now combine them and add legend
fig4 <- ggarrange(fig4ab, fig4cd, as_ggplot(legend),
                  ncol = 1, heights = c(1, 1, 0.1))
# save as pdf
pdf("Figure4.pdf", width = 15, height = 10)
fig4
dev.off()
embedFonts("Figure4.pdf", fontpaths = "/Library/Fonts/Microsoft/Arial.ttf")

# save as tiff
tiff("Figure4.tiff", units = "in", width = 15, height = 10, res = 250)
fig4
dev.off()

# save as png
png("Figure4.png", units = "in", width = 15, height = 10, res = 250)
fig4
dev.off()