library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(geosphere)
library(ggthemes)

###### NATURE FIGURES MAIN



########################################## NATURE FIGURES MAIN SUPPLEMENTARY 
########################################## 

### FIG 0 LOSS RATE AND AGE OF PREGNANCY
################################################################################
################################################################################

mod <- gam(cbind(non_live_births, number_of_live_births) ~ 
             s(distance_40, bs = "cs") +
             s(total_pregnancies, bs = "cs") +
             old_qual_female +
             old_qual_male +
             joint_income +
             health_female +
             female_previous_smoker +
             s(old_bmi_female, bs = "cs") +
             age_male +
             # minority_male +
             # minority_female +
             s(age_female_recruitment, bs = "cs"),
           family = "quasibinomial", 
           method = "REML", 
           # data = main_couples_age_controls)
          data = white_europe_age_controls)
summary(mod)

(mod_plot <- effect_plot(mod, pred = age_female_recruitment, interval = T, colors = "firebrick") + theme_few() +
  labs(x = "Average age of pregnancy", y = "Rate of pregnancy loss") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)))

ggsave("./newphd/graphics/plot.pdf",
       mod_plot, width = 18, height = 10,
       device='tiff', units = "cm", dpi=450)

### FIG 1 PREGNANCIES VS. LIVE BIRTHSM ETC
################################################################################
################################################################################
couples <- fread("./NATURE_PAPER/DATA/couples_with_genetic_distances.csv")

pregnancies_births <- gam(number_of_live_births ~ s(total_pregnancies) +
                            qual_female +
                            joint_income +
                            health_female +
                            female_previous_smoker +
                            s(bmi_female) +
                            age_male +
                            Number_of_minorities +
                            s(age_female_recruitment),
                          method = "REML", data = couples)
summary(pregnancies_births)
(preg <- effect_plot(pregnancies_births, pred = total_pregnancies, interval = T, colors = "firebrick") +
  scale_x_continuous(limits = c(0, 18), breaks = c(0:18)) +
  scale_y_continuous(limits = c(0, 7), breaks = c(0:7)) +
  labs(y = "Total children born", x = "Total pregnancies") +
    theme_few() +
    theme(text = element_text(size = 16),
          axis.text = element_text(size = 14, color= "black"),
          axis.title = element_text(size = 16, color = "black")))
  

pregnancies_miscarriages <- gam(total_complications ~ s(total_pregnancies) +
                            qual_female +
                            joint_income +
                            health_female +
                            female_previous_smoker +
                            s(bmi_female) +
                            age_male +
                            Number_of_minorities +
                            s(age_female_recruitment),
                          method = "REML", data = couples)
summary(pregnancies_miscarriages)
(miss <- effect_plot(pregnancies_miscarriages, pred = total_pregnancies, interval = T, colors = "firebrick") +
  scale_x_continuous(limits = c(0, 14), breaks = c(0:14)) +
  scale_y_continuous(limits = c(-0.1, 11), breaks = c(0:11)) +
  labs(y = "Total pregnancy losses", x = "Total pregnancies") +
    theme_few() +
    theme(text = element_text(size = 16),
          axis.text = element_text(size = 14, color= "black"),
          axis.title = element_text(size = 16, color = "black")))

(both <- ggarrange(preg, miss, labels = c("a", "b"), 
          ncol = 1))
ggsave("./newphd/graphics/plot.pdf",
       both, width = 24, height = 25,
       device='tiff', units = "cm", dpi=250)

### FIG 1 DISTRIBUTION OF GENETIC DISTANCES IN 2 MILLION RANDOM PAIRS
################################################################################
################################################################################

setwd("/Users/tomversluys/Documents")
data <- fread("./newphd/data/ukbb/ukbb_old_subset.csv")

random_pairs1 <- data %>% sample_n(2000000, replace = TRUE) %>% rename_with(tolower)
random_pairs2 <- data %>% sample_n(2000000, replace = TRUE) %>% rename_with(toupper)
random_pairs3 <- cbind(random_pairs1, random_pairs2) 
random_pairs3 <- random_pairs3 %>% mutate(newid = paste(1:2000000)) %>% arrange(newid) %>% data.table()

# calculate birth distance 
random_pairs3$birth_distance <- distHaversine(cbind(random_pairs3$placebirth_long, random_pairs3$placebirth_lat), 
                                     cbind(random_pairs3$PLACEBIRTH_LONG, random_pairs3$PLACEBIRTH_LAT))
random_pairs3 <- random_pairs3 %>% 
  mutate(birth_distance = as.numeric(birth_distance), birth_distance_km = birth_distance/1000)

# colnames(random_pairs3)
data <- random_pairs3 %>% arrange(newid) %>% select(pc1:pc40, PC1:PC40) %>% data.table() 

# # define function for measuring euclidean distance by column and row
product <- function(x, output){
  a = x[1:40]
  b = x[41:80] # first 5
  c = dist(rbind(a, b))
  return(c)
}

list <- data.frame(apply(data, 1, product)) %>% # first 5
  rename(distance_40 = apply.data..1..product.) %>%
  data.table()

dat1 <- cbind(list, random_pairs3) 
dat1 <- dat1 %>% mutate(log_distance = log(distance_40+0.5))
mean(dat1$log_distance, na.rm = T)
cor_coef <- cor.test(dat1$birth_distance_km, dat1$log_distance)

#### make plot
annotations <- data.frame(x = c(3.27, 3.4, 4.2, 4.5, 5.1, 5.3),
                          y = c(750, 1850, 750, 1850, 750, 1850), 
                          label = c("Born in England\n(G=3.27,\nN=1,225,412)",
                            "Born in different\n UK countries\n(G=3.4, N=463,005)",
                            "Born in England\n and Italy\n(G=4.2, N=5,083)",
                            "Born in England\n and Turkey\n(G=4.5, N=925)",
                            "Born in different\n African countries\n(G=5.1, N=1,012)",
                            "Born in Asia\n and Africa\n(G=5.3, N=2,145)"))

hist_dat <- fread("./NATURE_PAPER/DATA/couples_with_genetic_distances.csv")
(p <- ggplot(hist_dat, aes(x=distance_40)) +
    geom_histogram(binwidth = 0.05, colour = "#F39B7FFF", fill = "#F39B7FFF", alpha = 0.3) +
    geom_rect(aes(xmin = 3.237169, xmax = 3.241603, ymin = 0, ymax = Inf),  color='grey', linetype='solid', alpha=0.2) +
    geom_segment(aes(x = 3.239373, y = 0, xend = 3.239373, yend = Inf), linetype = "dashed", size = 1.2) +
    labs(y = "Number of couples", x = "Log-genetic distance (40 PCs)") +
    theme_few() + xlim(min(hist_dat$distance_40), max(hist_dat$distance_40)) +
    theme(text = element_text(size = 16), axis.text = element_text(size = 16)))

(p2 <- p + geom_segment(aes(x = 3.27, y = 0, xend = 3.27, yend = 750), linetype = "dashed") +
    geom_segment(aes(x = 3.4, y = 0, xend = 3.4, yend = 1850), linetype = "dashed") +
    geom_segment(aes(x = 4.2, y = 0, xend = 4.2, yend = 750), linetype = "dashed") +
    geom_segment(aes(x = 4.5, y = 0, xend = 4.5, yend = 1850), linetype = "dashed") +
    geom_segment(aes(x = 5.1, y = 0, xend = 5.1, yend = 750), linetype = "dashed") +
    geom_segment(aes(x = 5.3, y = 0, xend = 5.3, yend = 1850), linetype = "dashed") +
    geom_label(data = annotations, aes(x = x, y = y, label = paste(label)), size = 5, fill="white") + 
    theme(text = element_text(size = 16), axis.text = element_text(size = 14, color = "black"),
          axis.title = element_text(size = 16, color = "black")))

ggsave("./newphd/graphics/plot.pdf",
       p2, width = 24, height = 14,
       device='tiff', units = "cm", dpi=250)

### 1 SCREE PLOTS FOR PC VARIANCE EXPLAINED
################################################################################
################################################################################ 

data <- fread("./newphd/data/ukbb/ukbb_old_subset.csv")

data_pcs <- data_osb %>% dplyr::select(PC1:PC40) %>% drop_na()
variances <- apply(data_pcs, 2, var)
proportion_variance <- variances / sum(variances)
frame <- data_frame(PC = 1:40, y = proportion_variance)

# Plot the proportion of variance explained
(scree_plot <- ggplot(frame, aes(PC, y)) + 
    geom_point(color = "forestgreen", size = 3, alpha = 0.5) + 
    geom_line(colour = "forestgreen") + 
    labs(x = "PC number", y = "Proportion of Variance") + 
    theme_few() + 
    theme(text = element_text(size = 16),
           axis.text = element_text(size = 14, color = "black")))

ggsave("./newphd/graphics/plot1.pdf", width = 24, height = 14, device='tiff',
       units = "cm", dpi=250)

# # Find the elbow point
# diff_proportion_variance <- diff(proportion_variance)
# diff_between_diffs <- diff_proportion_variance[-1] - diff_proportion_variance[-length(diff_proportion_variance)]
# elbow_point <- which.max(diff_between_diffs) + 1
# points(elbow_point, proportion_variance[elbow_point], col = "red", pch = 19)

### 2 GENETIC DISTANCE 40 PCS AND BIRTH DISTANCE IN COUPLES
################################################################################
################################################################################ 
couples <- fread("./NATURE_PAPER/DATA/couples_with_genetic_distances.csv")
cor_coef <- cor.test(couples$birth_distance, couples$distance_40)
cor_coef <- cor_coef$estimate

(p <- ggplot(couples, aes(x=birth_distance, y=distance_40)) +
    geom_point(alpha=0.7, position=position_jitter(width=0.5, height=0.3), colour = "#D55E00") +
    geom_smooth(method="lm", se=FALSE, linetype="solid", color="black", size=1.5) +
    labs(
      title = "",
      x="Birth Distance (km)",
      y="Log-genetic distance (40 PCs)",
      caption=paste("Correlation Coefficient (r):", round(cor_coef, 2))) +
    theme_few() +
    theme(
      legend.position="none",
      plot.title = element_text(size = 20, face="bold", hjust = 0.5),
      text = element_text(size = 16),
      axis.text = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 16, color = "black"),
      plot.caption = element_text(size=16, face="bold")) +
    guides(size=FALSE))

ggsave("./newphd/graphics/plot1.pdf", p, width = 24, height = 14, device='tiff',
       units = "cm", dpi=250)

### 2 GENETIC DISTANCE KINSHIP ASSOCIATION
################################################################################
################################################################################ 

data <- fread("./newphd/data/ukbb/ukbb_old_subset.csv")
relatives <- data %>% 
  dplyr::select(ID, Ethnicity, PLaceBirth_lat, PLaceBirth_long, recruitment_age, PC1:PC40) %>% 
  drop_na() %>%
  data.table() # all 40
setwd("/Users/tomversluys/Documents")

rel <- fread("./newphd/data/ukbb/rel.csv")
rel1 <- rel %>% dplyr::select(1,3,4,5) %>% mutate(pairid = 1:length(ID1)) %>% rename(ID = ID1) %>%
  left_join(relatives, by = "ID") %>% 
  rename_with(tolower)
rel2 <- rel %>% dplyr::select(2) %>% mutate(pairid = 1:length(ID2)) %>% rename(ID = ID2) %>%
  left_join(relatives, by = "ID") %>%
  rename_with(toupper)

rel3 <- cbind(rel1, rel2) %>% arrange(pairid) %>% mutate(ethnic_match = ifelse(ethnicity == ETHNICITY, "yes", "no"),
                                                         age_diff = abs(recruitment_age - RECRUITMENT_AGE)) %>%
  data.table()

rel3 <- rel3 %>% dplyr::select(pc1:pc40, PC1:PC40, kinship, ibs0, ethnic_match, 
                               PLACEBIRTH_LAT, PLACEBIRTH_LONG, placebirth_lat, placebirth_long, 
                               age_diff) %>% filter(kinship >=-0)

# 3. Create a function to calculate Euclidean distance
euclidean_distance <- function(partner1, partner2) {
  sqrt(sum((partner1 - partner2)^2))
}

# 4. Calculate the Euclidean distances using an increasing number of principal components
distances <- data.frame(matrix(nrow = nrow(rel3), ncol = 0))
colnames(rel3)
for (i in 2:40) {
  pc_cols_partner1 <- paste0("pc", 1:i)
  pc_cols_partner2 <- paste0("PC", 1:i)
  distance_col_name <- paste0("distance_", i)
  distances[[distance_col_name]] <- apply(rel3, 1, function(row) {
    partner1_pcs <- as.numeric(row[pc_cols_partner1])
    partner2_pcs <- as.numeric(row[pc_cols_partner2])
    euclidean_distance(partner1_pcs, partner2_pcs)
  })
}

distances <- distances %>%
  mutate_at(vars(starts_with("distance_")), log)

# Bind the data frame with the Kinship column
distances <- cbind(distances, ibs0 = rel3$ibs0, kinship = rel3$kinship, 
                   placebirth_long = rel3$placebirth_long, 
                   placebirth_lat = rel3$placebirth_lat, 
                   PLACEBIRTH_LONG = rel3$PLACEBIRTH_LONG, 
                   PLACEBIRTH_LAT = rel3$PLACEBIRTH_LAT, 
                   age_diff = rel3$age_diff, 
                   ethnic_match = rel3$ethnic_match)
colnames(distances)
# calculate birth distance 
distances$birth_distance <- distHaversine(cbind(distances$placebirth_long, distances$placebirth_lat), 
                                          cbind(distances$PLACEBIRTH_LONG, distances$PLACEBIRTH_LAT))

# Convert the data frame to a data.table
distances_dt <- data.table(distances)

# 5. Calculate correlations between each distance variable and Kinship
correlations <- sapply(distances_dt[, 1:39], function(x) cor(x, distances_dt$ibs0, use = "complete.obs"))

# 6. Calculate 95% confidence intervals for the correlations
alpha <- 0.05
z_critical <- qnorm(1 - alpha/2)
correlations_ci <- sapply(correlations, function(r) {
  n <- nrow(distances_dt)
  z <- 0.5 * log((1 + r) / (1 - r))
  lower_bound <- tanh(z - z_critical / sqrt(n - 3))
  upper_bound <- tanh(z + z_critical / sqrt(n - 3))
  c(lower_bound, upper_bound)
})

# 7. Plot the correlations with error bars for the 95% confidence intervals
correlations_df <- data.frame(x = 2:40, y = correlations, lower = correlations_ci[1,], upper = correlations_ci[2,])

(cor_plot <- ggplot(correlations_df, aes(x, y)) +
    geom_point(color = "forestgreen", size = 3, alpha = 0.3) +
    # ylim(0, -1) +
    geom_line(color = "forestgreen") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, color = "forestgreen") +
    labs(x = "Number of PCs", y = "Kinship-genetic distance correlation (r)") +
    theme(text = element_text(size = 15),
          axis.text = element_text(size = 12, color= "black"),
          axis.title = element_text(size = 15, color = "black")) +
    theme_few())

ggsave("./newphd/graphics/plot.pdf",
       cor_plot, width = 18, height = 10,
       device='tiff', units = "cm", dpi=300)

######### GEOGRAPHICAL VALIDATION: BIRTH DISTANCE GENETIC DISTANCE CORRELATION
################################################################################
################################################################################ 
couples <- fread("./NATURE_PAPER/DATA/couples_with_genetic_distances.csv")

couples <- couples %>% drop_na(distance_2:distance_40, birth_distance) %>%
  data.table()

correlations <- sapply(couples[, 1:40], 
                       function(x) cor.test(x, couples$birth_distance, 
                                            use = "complete.obs"))

correlations <- lapply(2:40, function(i) {
  var_name <- paste("distance", i, sep = "_")
  test_cor <- cor.test(couples[[var_name]], couples$birth_distance)
  
  return(data.frame(
    distance_var = var_name,
    cor_coef = test_cor$estimate,
    p_value = test_cor$p.value,
    lower_ci = test_cor$conf.int[1],
    upper_ci = test_cor$conf.int[2]
  ))
})

# Combine results into a single data frame
results <- data.table(bind_rows(correlations)) %>%
  mutate(var = 2:40)

# Plot the results using ggplot
(geo_plots <- ggplot(results, aes(x = var, y = cor_coef)) +
    geom_point(colour = "firebrick", size = 3, alpha = 0.3) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, colour = "firebrick") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_aaas() +
    labs(x = "Number of PCs", y = "Birth distance-genetic distance correlation (r)") +
    # geom_hline(yintercept = 0, linetype = "dashed") +
    ylim(0.3, 0.4) +
    theme(text = element_text(size = 15),
  axis.text = element_text(size = 12, color= "black"),
  axis.title = element_text(size = 15, color = "black")) +
    theme_few())

ggsave("./newphd/graphics/plot.pdf",
       geo_plots, width = 18, height = 10,
       device='tiff', units = "cm", dpi=300)

### 3 DISTRIBUTION OF GENETIC DISTANCES
# birth distance analysis

### rank analysis
white <- lll %>% mutate(rank_10 = rank(log_distance), 
                        rank_40 = rank(log_distance_40))
cor.test(white$rank_10, white$rank_40)

low <- white %>% filter(rank_10 > 25200)
cor.test(low$rank_10, low$rank_40)

quantile(white$rank_10, 0.90)
quantile(white$log_distance_40, 0.90)

a <- white %>% filter(log_distance < 3)
b <- white %>% filter(log_distance_40 < 3.4)
c <- a %>% filter(id %in% b$id)

a <- c(2, 2)
b <- c(3, 3)
dist(rbind(a, b))
range(lll$e_distance_raw40)


# run relatives analysis 
#########
rel <- fread("./newphd/data/ukbb/rel.csv")
ll <- check2
rels_test <- ll %>% inner_join(rel, by = c("id" = "ID1", "ID" = "ID2"))
rels_test1 <- ll %>% inner_join(rel, by = c("id" = "ID2", "ID" = "ID1"))
rels_test2 <- rbind(rels_test, rels_test1) %>%
  rowwise() %>%
  mutate(total_pregnancies = sum(number_of_live_births, number_of_miscarriages, number_of_stillbirths, na.rm = T),
         non_live_births = (total_pregnancies - number_of_live_births), 
         non_miscarriages = (total_pregnancies - number_of_miscarriages), 
         non_stillbirths = (total_pregnancies - number_of_stillbirths)) 

test_mod <- gam(number_of_live_births ~ s(Kinship, bs = "cr") + 
                  s(age_female_recruitment, bs = "cr"),
                family = quasipoisson(),
                method = "REML",
                data = rels_test2)
summary(test_mod)
effect_plot(test_mod, pred = Kinship, plot.points = T, interval = T, colors = "firebrick") + 
  labs(y = "Total children born", x = "Kinship") +
  theme_few() +
  theme(text = element_text(size = 16),
        axis.text = element_text(size = 14, color= "black"),
        axis.title = element_text(size = 16, color = "black"))

test_mod1 <- gam(total_pregnancies ~ s(Kinship, bs = "cr") + 
                   s(age_female_recruitment, bs = "cr"),
                 family = quasipoisson(),
                 method = "REML",
                 data = rels_test2)
summary(test_mod1)
effect_plot(test_mod1, pred = Kinship, plot.points = T, interval = T, 
            colors = "firebrick") + 
  labs(y = "Total pregnancies", x = "Kinship") +
  theme_few() +
  theme(text = element_text(size = 16),
        axis.text = element_text(size = 14, color= "black"),
        axis.title = element_text(size = 16, color = "black"))

test_mod2 <- gam(total_complications ~ s(Kinship, bs = "cr") +
                   s(age_female_recruitment, bs = "cr"),
                 family = quasipoisson(),
                 method = "REML",
                 data = rels_test2)
summary(test_mod2)
effect_plot(test_mod2, pred = Kinship, plot.points = T, interval = T, 
            colors = "firebrick") + 
  labs(y = "Total losses", x = "Kinship") +
  theme_few() +
  theme(text = element_text(size = 16),
        axis.text = element_text(size = 14, color= "black"),
        axis.title = element_text(size = 16, color = "black"))

test_mod_p <- gam(cbind(number_of_live_births, non_live_births) ~ s(Kinship, bs = "cr") + 
                  s(total_pregnancies, bs = "cs") +
                  s(age_female_recruitment, bs = "cr"),
                family = quasibinomial(),
                method = "REML",
                data = rels_test2)

summary(test_mod_p)
effect_plot(test_mod_p, pred = Kinship, interval = T, colors = "firebrick") + 
  labs(y = "Total children born", x = "Kinship") +
  theme_few() +
  theme(text = element_text(size = 16),
        axis.text = element_text(size = 14, color= "black"),
        axis.title = element_text(size = 16, color = "black"))


