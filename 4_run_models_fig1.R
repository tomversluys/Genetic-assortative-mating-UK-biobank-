
##### THIS SRIPT CREATES SUBSETS AND RUNS MODELS ON ONE DATASET 
rm(list=ls())
library(gdata)
library(ggpubr)
library(dplyr)
library(ggthemes)
library(ggsci)
library(jtools)
library(interactions)
library(data.table)
library(stringr)
library(ggeffects)
library(tidyr)
library(ggplot2)
library(ggstance)
library(broom.mixed)
library(broom)
library(scales)
library(grid)
library(gratia)
library(mgcv)
library(magrittr)
#### LOAD PACKAGES
library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)
library(mgcv)
library(gratia)
library(purrr)
library(ggthemes)
library(ggpubr)
library(gridExtra)
library(jtools)

#### READ IN DATA
###############################################################################
main_couples_data <- fread("./NATURE_PAPER/DATA/couples_with_genetic_distances.csv")
# main_couples_data <- fread("<3_DATASET_HERE>") # load df output from previous script (3)

main_couples_data <- main_couples_data %>% mutate(ethnicity_female = as.factor(ethnicity_female), ethnicity_male = as.factor(ethnicity_male),
         female_previous_smoker = as.factor(female_previous_smoker)) %>%
  drop_na(distance_1:distance_40) %>%
  rowwise() %>%
  mutate(total_pregnancies = sum(number_of_live_births, number_of_miscarriages, number_of_stillbirths, na.rm = T),
         non_live_births = (total_pregnancies - number_of_live_births), 
         non_miscarriages = (total_pregnancies - number_of_miscarriages), 
         non_stillbirths = (total_pregnancies - number_of_stillbirths)) %>%
  drop_na(number_of_live_births, number_of_miscarriages, number_of_stillbirths) %>%
  data.table()

#### CREATE SUBSETS
###############################################################################
#### AGE CONTROLS
main_couples_age_controls <- main_couples_data %>% 
  drop_na(age_first_live_birth, age_last_live_birth) %>% 
  filter(age_first_live_birth > 10, age_last_live_birth > 10) %>% 
  mutate(age_female_recruitment = (age_first_live_birth + age_last_live_birth)/2) %>%
  data.table() 

#### WHITE UK
white_uk_data <- main_couples_data %>% 
  mutate(Ethnic_subset = ifelse(ethnicity_female %in% c("White") & 
                                  ethnicity_male %in% c("White"), 
                                "All white", "not")) %>%
  filter(born_uk == "yes") %>%
  filter(country == COUNTRY,
         Ethnic_subset == "All white")

#### WHITE UK AGE CONTROLS
white_uk_age_controls <- 
  white_uk_data %>% drop_na(age_first_live_birth, age_last_live_birth) %>% 
  filter(age_first_live_birth > 10, age_last_live_birth > 10) %>% 
  mutate(age_female_recruitment = (age_first_live_birth + age_last_live_birth)/2) %>%
  data.table() 

# Function to fit the GAM model for a given distance variable
###############################################################################
fit_gam <- function(distance_var, data) {
  model_formula <- update(base_formula, as.formula(paste("~ . + s(", distance_var, ', bs = "cr")')))
  look <- gam(model_formula,
              family = quasipoisson(),
              method = "REML",
              data = data)
  return(look)
}

# Function to create prediction data for a given data set, distance variable, and model
###############################################################################
create_pred_data <- function(data, distance_var, model) {
  pred_data <- data.frame(
    genetic_distance = seq(min(data[[distance_var]]), max(data[[distance_var]]), length.out = 100),
    old_bmi_female = mean(data$old_bmi_female, na.rm = TRUE),
    old_bmi_male = mean(data$old_bmi_male, na.rm = TRUE),
    age_female_recruitment = mean(data$age_female_recruitment, na.rm = TRUE),
    age_male_recruitment = mean(data$age_male_recruitment, na.rm = TRUE),
    old_qual_female = mean(data$old_qual_female, na.rm = TRUE),
    old_qual_male = mean(data$old_qual_male, na.rm = TRUE),
    minority_female = "no", 
    minority_male = "no",
    Number_of_minorities = 0,
    joint_income = mean(data$joint_income, na.rm = TRUE),
    health_female = mean(data$health_female, na.rm = TRUE),
    health_male = mean(data$health_male, na.rm = TRUE),
    female_previous_smoker = "0",
    age_first_live_birth = mean(data$age_first_live_birth, na.rm = TRUE),
    total_pregnancies = mean(data$total_pregnancies, na.rm = TRUE)
  )
  
  pred_data[[distance_var]] <- pred_data$genetic_distance
  preds <- predict(model, newdata = pred_data, type = "link", se.fit = TRUE)
  
  # Transform the predictions and confidence intervals back to the response scale
  ###############################################################################
  pred_data$fit <- exp(preds$fit)
  pred_data$lower <- exp(preds$fit - 1.96 * preds$se.fit)
  pred_data$upper <- exp(preds$fit + 1.96 * preds$se.fit)
  pred_data$peak <- which.max(pred_data$fit)
  
  fd <- gratia::derivatives(model, term = paste0("s(", distance_var, ")"), newdata = pred_data) %>%
    mutate(distance_40 = pred_data$genetic_distance)
  
  threshold <- 0.01
  pred_data <- pred_data %>% mutate(lower_d = fd$lower, upper_d = fd$upper,
                                    sig = ifelse((lower_d * upper_d) > 0, "sig", "not_sig"),
                                    line_type = ifelse(sig == "sig", "sig", "not_sig")) %>%
    select(genetic_distance, fit, upper, lower, line_type)
  
  pred_data <- pred_data %>%
    mutate(prev_fit = lag(fit),
           prev_genetic_distance = lag(genetic_distance),
           prev_line_type = lag(line_type)) %>%
    filter(!is.na(prev_fit))
  
  return(pred_data)
}


# PLOT DISTANCE VARIABLE FUNCTION FOR MAIN COUPLES
###############################################################################
plot_distance_variable <- function(distance_var, model, data, label) {
  
  pred_data <- create_pred_data(data, distance_var, model)
  
  # Calculate p-values for slopes
  original_p_value <- get_pvalue_distance(model, distance_var)
  
  pc_number <- gsub("distance_", "", distance_var)
  
  format_p_value <- function(p_value, min_exponent = 10) {
    if (p_value == 0.0) {
      return(sprintf("< 1.0e-%d", min_exponent))
    } else if (p_value < 0.001) {
      return(sprintf("%.1e", p_value)) # use scientific notation for smaller p-values
    } else {
      return(sprintf("%.3f", p_value))
    }
  }
  
  format_nature_p_value <- function(p_value) {
    formatted_p_value <- format_p_value(p_value)
    
    if (grepl("e", formatted_p_value)) {
      split_value <- strsplit(formatted_p_value, "e")[[1]]
      exponent <- as.integer(split_value[2])
      return(paste0(split_value[1], " \u00D7 10^", exponent))
    } else {
      return(formatted_p_value)
    }
  }
  
# Plotting predictions 
###############################################################################
  p <- ggplot() +
    geom_segment(data = pred_data, aes(color = "White UK", x = prev_genetic_distance, xend = genetic_distance, y = prev_fit, yend = fit, size = line_type), alpha = 0.6) +
    geom_ribbon(data = pred_data, aes(fill = "White UK", x = prev_genetic_distance, ymin = lower, ymax = upper), alpha = 0.2) +
    theme_few() +
    theme(plot.title = element_blank()) +
    scale_color_manual(name = "Sample", values = c("White UK" = "firebrick", "Age Control" = "steelblue")) +
    scale_fill_manual(name = "Sample", values = c("White UK" = "firebrick", "Age Control" = "steelblue")) +
    scale_size_manual(values = c("sig" = 1.5, "not_sig" = 0.5), guide = "none") +
    scale_linetype_manual(values = c("sig" = "solid", "not_sig" = "solid")) +
    guides(linetype = FALSE, colour = FALSE, fill = FALSE) +
    labs(x = paste("Log-genetic distance"))
    theme(legend.position = c(1, 1), legend.justification = c(1, 1), legend.box.just = "right",
          legend.key.size = unit(0.5, "cm"), legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
    theme(text = element_text(size = 18),
          axis.text = element_text(size = 16, color= "black"),
          axis.title = element_text(size = 18, color = "black"),
          legend.box.background = element_rect(colour = "black", size = 1),
          legend.position = c(0.95, 0.98),
          legend.justification = c(1, 1),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 13),
          legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.spacing.x = unit(0.1, "cm"),
          legend.spacing.y = unit(0.2, "cm"))
  
  return(p)
}

# base formulas
###############################################################################
# total_pregnancies
base_formula <-
  total_pregnancies ~ 
  old_qual_female +
  old_qual_male +
  joint_income +
  health_female +
  health_male +
  female_previous_smoker +
  s(old_bmi_female, bs = "cs") +
  old_bmi_male +
  Number_of_minorities +
  # minority_female + # factor alternative: comment in/out depending on preferred coding, as justified by theory
  # minority_male +
  age_male_recruitment +
  s(age_female_recruitment, bs = "cr")

distance_vars <- paste0("distance_", c(40)) # one
models <- map(distance_vars, fit_gam, data = main_couples_data)
labels <- paste("Log-genetic distance (",40, " PCs)", sep = "") # one
labels <- paste("data")

(plot <- mapply(plot_distance_variable, distance_vars, models, 
                        MoreArgs = list(data = main_couples_data, label = labels), SIMPLIFY = FALSE))

# number_of_live_births
base_formula <-
  number_of_live_births ~ # adjust response variable as desired
  old_qual_female +
  old_qual_male +
  joint_income +
  health_female +
  health_male +
  female_previous_smoker +
  s(old_bmi_female, bs = "cs") +
  old_bmi_male +
  Number_of_minorities +
  age_male_recruitment +
  s(age_female_recruitment, bs = "cr")

distance_vars <- paste0("distance_", c(40)) 
models <- map(distance_vars, fit_gam, data = main_couples_data) # change "main_couples_data" to preferred dataset
labels <- paste("Log-genetic distance (",40, " PCs)", sep = "") # one
labels <- paste("data")

(plot <- mapply(plot_distance_variable, distance_vars, models, 
                MoreArgs = list(data = main_couples_data, label = labels), SIMPLIFY = FALSE))


################################ RATIO MODELS
# library(gratia)
# Function to fit the GAM model for a given distance variable
fit_gam <- function(distance_var, data) {
  model_formula <- update(base_formula, as.formula(paste("~ . + s(", distance_var, ', bs = "cr")')))
  look <- gam(model_formula,
              family = binomial(),
              method = "REML",
              data = data)
  return(look)
}

# Function to create prediction data for a given data set, distance variable, and model
create_pred_data <- function(data, distance_var, model) {
  pred_data <- data.frame(
    genetic_distance = seq(min(data[[distance_var]]), max(data[[distance_var]]), length.out = 100),
    old_bmi_female = mean(data$old_bmi_female, na.rm = TRUE),
    old_bmi_male = mean(data$old_bmi_male, na.rm = TRUE),
    age_female_recruitment = mean(data$age_female_recruitment, na.rm = TRUE),
    age_male_recruitment = mean(data$age_male_recruitment, na.rm = TRUE),
    old_qual_female = mean(data$old_qual_female, na.rm = TRUE),
    old_qual_male = mean(data$old_qual_male, na.rm = TRUE),
    minority_female = "no", 
    minority_male = "no",
    Number_of_minorities = 0,
    joint_income = mean(data$joint_income, na.rm = TRUE),
    health_female = mean(data$health_female, na.rm = TRUE),
    health_male = mean(data$health_male, na.rm = TRUE),
    female_previous_smoker = "0",
    age_first_live_birth = mean(data$age_first_live_birth, na.rm = TRUE),
    total_pregnancies = mean(data$total_pregnancies, na.rm = TRUE)
  )
  
  table(main_couples_data$min)  
  pred_data[[distance_var]] <- pred_data$genetic_distance
  preds <- predict(model, newdata = pred_data, type = "link", se.fit = TRUE)
  
  
  pred_data$lower <- plogis(preds$fit - 1.96 * preds$se.fit) ### FOR BINOMIAL
  pred_data$upper <- plogis(preds$fit + 1.96 * preds$se.fit)
  pred_data$fit <- plogis(preds$fit)
  
  pred_data <- pred_data %>% mutate(fit = (1 - fit), 
                                    lower = (1 - lower),
                                    upper = (1 - upper))
  
  
  # Add the derivatives to the pred_data dataframe
  fd <- gratia::derivatives(model, term = paste0("s(", distance_var, ")"), newdata = pred_data) %>%
    mutate(distance_40 = pred_data$genetic_distance)
  
  threshold <- 0.01  # Adjust this value according to your requirements
  pred_data <- pred_data %>% mutate(lower_d = fd$lower, upper_d = fd$upper,
                                    sig = ifelse((lower_d * upper_d) > threshold, "sig", "not_sig"),
                                    line_type = ifelse(sig == "sig", "sig", "not_sig"))
  
  pred_data <- pred_data %>%
    mutate(prev_fit = lag(fit),
           prev_genetic_distance = lag(genetic_distance),
           prev_line_type = lag(line_type)) %>%
    filter(!is.na(prev_fit))
  
  return(pred_data)
}

get_pvalue_distance <- function(model, distance_var) {
  summary(model)$s.table[grepl(paste0("s\\(", distance_var), rownames(summary(model)$s.table)), "p-value"]
}

plot_distance_variable <- function(distance_var, model, age_control_model, data, age_control_data, label) {
  # Creating prediction data for the original data set
  pred_data <- create_pred_data(data, distance_var, model)
  
  # Creating prediction data for the age_control data set
  age_control_pred_data <- create_pred_data(age_control_data, distance_var, age_control_model)
  
  # Calculate p-values for both slopes
  original_p_value <- get_pvalue_distance(model, distance_var)
  age_control_p_value <- get_pvalue_distance(age_control_model, distance_var)
  
  pc_number <- gsub("distance_", "", distance_var)
  
  format_p_value <- function(p_value, min_exponent = 10) {
    if (p_value == 0.0) {
      return(sprintf("< 1.0e-%d", min_exponent))
    } else if (p_value < 0.001) {
      return(sprintf("%.1e", p_value)) # use scientific notation for smaller p-values
    } else {
      return(sprintf("%.3f", p_value))
    }
  }
  
  format_nature_p_value <- function(p_value) {
    formatted_p_value <- format_p_value(p_value)
    
    if (grepl("e", formatted_p_value)) {
      split_value <- strsplit(formatted_p_value, "e")[[1]]
      exponent <- as.integer(split_value[2])
      return(paste0(split_value[1], " \u00D7 10^", exponent))
    } else {
      return(formatted_p_value)
    }
  }
  
  
  # Plotting predictions for both data sets with a compact legend in the top right corner and p-values in the top center
  p <- ggplot() +
    
    geom_segment(data = pred_data, aes(color = "Original", x = prev_genetic_distance, xend = genetic_distance, y = prev_fit, yend = fit, size = line_type), alpha = 0.6) +
    geom_ribbon(data = pred_data, aes(fill = "Original", x = prev_genetic_distance, ymin = lower, ymax = upper), alpha = 0.2) +
    geom_segment(data = age_control_pred_data, aes(color = "Age Control", x = prev_genetic_distance, xend = genetic_distance, y = prev_fit, yend = fit, size = line_type), alpha = 0.6) +
    geom_ribbon(data = age_control_pred_data, aes(fill = "Age Control", x = prev_genetic_distance, ymin = lower, ymax = upper), alpha = 0.2) +
    
    scale_linetype_manual(values = c("sig" = "solid", "not_sig" = "dashed")) +
    
    theme_few() +
    labs(x = label) +
    theme(plot.title = element_blank()) +
    scale_color_manual(name = "", values = c("Original" = "firebrick", "Age Control" = "steelblue")) +
    scale_fill_manual(name = "", values = c("Original" = "firebrick", "Age Control" = "steelblue")) +
    scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed")) +
    guides(linetype = FALSE) +
    scale_size_manual(values = c("sig" = 1.5, "not_sig" = 0.5), guide = "none") +
    labs(x = paste("Log-genetic distance")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(breaks = c(3, 4, 5)) +
    theme(legend.position = c(1, 1), legend.justification = c(1, 1), legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6), legend.key.size = unit(0.5, "cm"), legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
    geom_rug(data = data, aes(x = distance_40), sides = "b", colour = "firebrick", alpha = 1) +
    geom_rug(data = age_control_data, aes(x = distance_40), sides = "b", colour = "steelblue", alpha = 0.1)
  return(p)
}

# # ANALYSIS FOR A SINGLE PC (40)
base_formula <-
  cbind(number_of_live_births, non_live_births) ~ 
  s(total_pregnancies, bs = "cs") +
  old_qual_female +
  old_qual_male +
  joint_income +
  health_female +
  health_male +
  female_previous_smoker +
  s(old_bmi_female, bs = "cs") +
  old_bmi_male +
  Number_of_minorities +
  # minority_female +
  # minority_male +
  age_male_recruitment +
  s(age_female_recruitment, bs = "cr")

distance_vars <- paste0("distance_", c(40)) # one
models <- map(distance_vars, fit_gam, data = main_couples_data)
age_control_models <- map(distance_vars, fit_gam, data = main_couples_age_controls)
labels <- paste("Log-genetic distance (",40, " PCs)", sep = "") # one
labels <- paste("data")
# #
(plots_births_rate <- mapply(plot_distance_variable, distance_vars, models, age_control_models, MoreArgs = list(data = main_couples_data, age_control_data = main_couples_age_controls, label = labels), SIMPLIFY = FALSE))

# # ANALYSIS FOR A SINGLE PC (40)
base_formula <-
  cbind(number_of_live_births, non_live_births) ~ 
  s(total_pregnancies, bs = "cs") +
  old_qual_female +
  old_qual_male +
  joint_income +
  health_female +
  health_male +
  female_previous_smoker +
  s(old_bmi_female, bs = "cs") +
  old_bmi_male +
  Number_of_minorities +
  # minority_female +
  # minority_male +
  age_male_recruitment +
  s(age_female_recruitment, bs = "cr")


distance_vars <- paste0("distance_", c(40)) # one
models <- map(distance_vars, fit_gam, data = white_uk_data)
age_control_models <- map(distance_vars, fit_gam, data = white_uk_age_controls)
labels <- paste("Log-genetic distance (",40, " PCs)", sep = "") # one
labels <- paste("data")
# #

(plots_births_rate_w <- mapply(plot_distance_variable, distance_vars, models, age_control_models, MoreArgs = list(data = white_uk_data, age_control_data = white_uk_age_controls, label = labels), SIMPLIFY = FALSE))




