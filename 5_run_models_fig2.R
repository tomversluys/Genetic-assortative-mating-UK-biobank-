# Load Libraries
###############################################################################
library(mgcv)
library(dplyr)
library(jtools)
library(interactions)
library(geosphere)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gratia)
library(ggeffects)

# data <- fread("./newphd/data/ukbb/ukbb_old_subset.csv") 
data <- fread("<1_DATASET_HERE>") # load df output from previous script (1)

# Data selection and mutation
###############################################################################
data <- data %>%
  mutate(Ethnicity = as.factor(Ethnicity)) %>%
  mutate(minority = ifelse(! Ethnicity  %in% c(1, 1001, 1002, 1003), 1, 0)) %>%
  filter(OverallHealth >= 0, recruitment_age > 15, household_income > 0, 
         c(ChildFatherered >= 0 | is.na(ChildFatherered)), 
         c(number_of_live_births >= 0 | is.na(number_of_live_births)),
         c(number_of_stillbirths >= 0 | is.na(number_of_stillbirths)),
         c(number_of_miscarriages >= 0 | is.na(number_of_miscarriages))) %>%
  mutate(children_either = ifelse(Sex == 1, ChildFatherered, number_of_live_births),
         complications = number_of_miscarriages + number_of_stillbirths,
         NumberSexPart = ifelse(NumberSexPart < 0, NA, NumberSexPart),
         born_uk = ifelse(country %in% c(1, 2, 3, 4, 5), "yes", "no"),
         Ethnicity = as.factor(Ethnicity)) %>%
  drop_na(Heterozy) %>%
  rename_with(tolower) %>%
  drop_na(number_of_live_births, heterozy, minority, old_bmi, household_income, 
          old_qual, overallhealth, recruitment_age) %>%
  filter(is.na(missing_heterozygosity)) %>%
  data.table() 

# Fit model and print summary
###############################################################################
mod <- gam(children_either ~ s(heterozy, bs = "cs") + minority + s(old_bmi, bs = "cs") + household_income + old_qual + overallhealth + s(recruitment_age, bs = "cs"), METHOD = "REML",  family = "quasipoisson", data = data)
summary(mod)
effect_plot(mod, pred = heterozy, interval = T)

# Create predictions and plot
###############################################################################
newdata <- with(data, data.frame(heterozy = seq(min(heterozy), max(heterozy), length.out = 100),
                                 migration_distance = mean(migration_distance, na.rm = TRUE),
                                 minority = 0,
                                 old_bmi = mean(old_bmi),
                                 household_income = mean(household_income),
                                 old_qual = mean(old_qual),
                                 overallhealth = mean(overallhealth),
                                 recruitment_age = mean(recruitment_age)))

preds <- predict(mod, newdata = newdata, type = "link", se.fit = TRUE)
newdata$fit <- exp(preds$fit)
newdata$lower <- exp(preds$fit - 1.96 * preds$se.fit)
newdata$upper <- exp(preds$fit + 1.96 * preds$se.fit)
newdata$peak <- which.max(newdata$fit)
fd <- gratia::derivatives(mod, term = "heterozy", partial_match = TRUE, newdata = newdata)

# Modify newdata for plotting
###############################################################################
newdata <- newdata %>% mutate(lower_d = fd$lower, upper_d = fd$upper,
                              sig = ifelse((lower_d * upper_d) > 0, "sig", "not_sig"),
                              line_type = ifelse(sig == "sig", "sig", "not_sig")) %>%
  select(heterozy, fit, upper, lower, line_type)

newdata <- newdata %>%
  mutate(prev_fit = lag(fit),
         prev_heterozy = lag(heterozy),
         prev_line_type = lag(line_type)) %>%
  filter(!is.na(prev_fit))

# Plot the smooth with confidence intervals
###############################################################################
het <- data
(pp <- ggplot(newdata) +
    geom_rect(aes(xmin = 0.1888, xmax = 0.194, ymin = -Inf, ymax = Inf), 
              fill = "lightgray", alpha = 0.05) + 
    geom_rect(aes(xmin = 0.1888, xmax = 0.194, ymin = -Inf, ymax = Inf), 
              linetype = "dashed", color = "black", fill = NA, size = 0.5) +
    geom_segment(aes(x = prev_heterozy, xend = heterozy, y = prev_fit, yend = fit, size = line_type), color = "firebrick", alpha = 0.6) +
    geom_ribbon(aes(x = heterozy, ymin = lower, ymax = upper), fill = "firebrick", alpha = 0.2) +
    scale_size_manual(values = c("sig" = 1.5, "not_sig" = 0.5), guide = "none") +
    labs(x = "Heterozygosity", y = "Fertility") +
    theme_few() +
    scale_x_continuous(labels = scales::percent_format(accuracy = 0.5)) +
    theme(legend.position = c(1, 1), legend.box.just = "right", legend.key.size = unit(0.5, "cm"), 
          legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
    geom_rug(data = het, aes(x = heterozy), sides = "b", colour = "firebrick", alpha = 1))

# trio heterozy analysis
###############################################################################

# read in data
###############################################################################
data_path <- "<YOUR_PATH_HERE>"

# Load and preprocess data
###############################################################################
rel <- fread("<RELATIVES_DATASET_HERE>") # load df of relative pairs
# rel <- fread("./newphd/data/ukbb/rel.csv")

### how many feature in relatives?
###############################################################################
# rel_couple_test <- rel %>% mutate(idtest_1 = paste(ID1, ID2, sep = ""), 
#                                   idtest_2 = paste(ID2, ID1, sep = ""))
# 
# rel_couple_test1 <- main_couples_data %>% mutate(idtest_1 = paste(id_male, id_female, sep = ""), 
#                                                  idtest_2 = paste(id_female, id_male, sep = ""))
# 
# one <- rel_couple_test1 %>% filter(c(idtest_1 %in% rel_couple_test$idtest_1) | c(idtest_1 %in% rel_couple_test$idtest_2)) %>%
#   select(idtest_1, idtest_2, id_male, id_female, distance_40, ethnicity_female, ethnicity_male)
# 
# two <- one %>% inner_join(rel_couple_test, by = c("idtest_1" = "idtest_1")) %>%
#   select(-idtest_2.y) %>% rename(idtest_2 = idtest_2.x)
# three <- one %>% inner_join(rel_couple_test, by = c("idtest_1" = "idtest_2")) %>%
#   select(-idtest_1.y)
# 
# four <- rbind(two, three)

# give pair IDs
###############################################################################
check <- rel %>% mutate(pair_id = 1:length(ID1)) %>% group_by(pair_id)

check <- check %>% filter(between(Kinship, 0.177, 0.354), IBS0 < 0.0012 ) %>% arrange(pair_id) %>%
  data.table()

check1 <- check %>% mutate(id = ID1) %>% select(id, pair_id, Kinship) 

check2 <- check %>% mutate(id = ID2) %>% select(id, pair_id, Kinship) 

check3 <- rbind(check1, check2) %>% arrange(pair_id)

# data <- fread("./newphd/data/ukbb/ukbb_old_subset.csv") # this is original data
data <- fread("<1_DATASET_HERE>") # load df output from first script (1)

dem <- data %>% select(old_dob, Sex, ID) %>% rename(id = ID)
check3 <- check3 %>% inner_join(dem, by = "id")

check4 <- check3 %>% group_by(id) %>% mutate(cases = length(id)) %>% group_by(pair_id) %>%
  mutate(cases_pairs = length(pair_id)) %>% 
  filter(cases_pairs >=2) %>%
  data.table()

check5 <- check4 %>% mutate(type = ifelse(cases == 2, "child", "unknown")) %>% data.frame()

check5a <- check5 %>% group_by(pair_id) %>% mutate(child_not = length(type[type=="child"])) %>% 
  arrange(pair_id) %>%
  data.table()

check5b <- check5a %>% group_by(pair_id) %>% mutate(class_id = ifelse(type == "child", id, NA)) %>%
  arrange(pair_id, class_id) %>% mutate(class_id = first(class_id)) %>% 
  mutate(type = ifelse(type == "unknown", "parent", type)) %>% arrange(class_id, type) 

check5c <- check5b %>% group_by(class_id) %>% mutate(test = length(unique(id))) %>%
  select(-c(cases, cases_pairs, child_not, test)) %>% rename(ID = id, trio_id = class_id)


check5c <- check5c %>% group_by(trio_id) %>%
  arrange(trio_id, type) %>%
  filter( (old_dob[1]+old_dob[2]) > (old_dob[3]+old_dob[4])) %>%
  filter(!Sex[3] == Sex[4]) %>% data.table()

length(unique(check5c$trio_id)) # 1173 (vs 1066 in Bycroft et al)

#### merge with couples data
###############################################################################
data <- fread("<1_DATASET_HERE>") # load df output from first script (1)
# data <- fread("./newphd/data/ukbb/ukbb_old_subset.csv") 

couples <- data %>% select(ID, PC1:PC40)

parents <- check5c %>% filter(type == "parent") 
parents <- inner_join(parents, couples, by = "ID") %>% data.table()

parents1 <- parents[duplicated(trio_id)]
parents <- parents[order(trio_id)]
names(parents1) <- tolower(names(parents1))
parents2 <- parents[! ID %in% parents1$id]
parents2 <- parents2[order(trio_id)]
names(parents2) <- toupper(names(parents2))
combined <- cbind(parents1, parents2) %>% arrange(pair_id)

parents3 <- combined %>% arrange(pair_id) %>% select(pc1:pc40, PC1:PC40) %>% data.table() # all 40

# define function for measuring euclidean distance by column and row
###############################################################################
colnames(parents3)
product <- function(x, output){
  a = x[1:40]
  b = x[41:80] # all
  c = dist(rbind(a, b))
  return(c)
}
list_40 <- data.frame(apply(parents3, 1, product)) %>%  # all
  rename(e_distance_raw40 = apply.parents3..1..product.) %>%
  data.table()

parents3 <- combined %>% arrange(pair_id) %>% select(pc1:pc10, PC1:PC10) %>% data.table() # first 5

# define function for measuring euclidean distance by column and row
###############################################################################

parents4 <- cbind(list_40,
                  combined)
parents5 <- parents4 %>% drop_na(e_distance_raw40) %>% 
  mutate(distance_40 = log(e_distance_raw40 + 0.5))

parents5 <- parents5 %>% select(ID, distance_40, trio_id)

children <- check5c %>% filter(type == "child") %>% group_by(ID) %>% slice(1) %>% rename_with(tolower) %>%
  arrange(trio_id)


parents6 <- parents5 %>% filter(trio_id %in% children$trio_id) %>% rename_with(toupper) %>%
  arrange(TRIO_ID) %>% select(-c(ID))

children1 <- children %>% filter(trio_id %in% parents6$TRIO_ID) %>% rename(ID = id) %>%
  rename_with(toupper)


both <- inner_join(children1, parents6, by = "TRIO_ID") %>% rename_with(toupper) %>%
  select(-SEX)

# here, I merge parent-offspring trios with the main data
###############################################################################
data <- fread("<1_DATASET_HERE>") # load df output from first script (1)
# data <- fread("./newphd/data/ukbb/ukbb_old_subset.csv")
data1 <- data %>% 
  select(ID, Sex, OverallHealth, number_of_live_births, number_of_miscarriages, number_of_stillbirths,
                         recruitment_age, SmokingStatus, NewQualZScore, Ethnicity, household_income,
                         Home_long, Home_lat, BMIZScore, ChildFatherered, birth_weight_first_child,
                         Heterozy, HetPCcorrect,  old_qual, old_bmi, birth_continent,
                         age_first_live_birth,age_last_live_birth, migration_distance) %>%
  filter(OverallHealth >= 0, recruitment_age > 15, household_income > 0,
         c(ChildFatherered >= 0 | is.na(ChildFatherered)), c(number_of_live_births >= 0 | is.na(number_of_live_births)),
         c(number_of_stillbirths >= 0 | is.na(number_of_stillbirths)), c(number_of_miscarriages >= 0 | is.na(number_of_miscarriages))) %>%
  mutate(children_either = ifelse(Sex == 1, ChildFatherered, number_of_live_births),
         complications = number_of_stillbirths + number_of_miscarriages)

both2 <- inner_join(both, data1, by = "ID") %>% rename_with(tolower) %>% filter(distance_40 > 1) %>%
  mutate(log_het = log(heterozy)) %>% data.table()

# what is the effect of parental genetic distance on offspring het?
###############################################################################

hh <- gam(heterozy ~
            s(distance_40, bs = "cs") , METHOD = "REML", 
          data = both2)
der_het <- derivatives(hh)
ggpredict(hh, terms = "distance_40 [3.7, 5.2]")
summary(hh)

(p <- effect_plot(hh, pred = distance_40, interval = T, 
                  rug = T, rug.sides = "b", colors = "firebrick") + 
    theme_few() +
    labs(y = "Heterozygosity", x = "Log-genetic distance") +
    geom_rect(aes(ymin = 0.189, ymax = 0.194, xmin = -Inf, xmax = Inf), 
              fill = "lightgray", alpha = 0.05) + 
    geom_rect(aes(ymin = 0.189, ymax = 0.194, xmin = -Inf, xmax = Inf), 
              linetype = "dashed", color = "black", fill = NA, size = 0.5) + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.5), limits = c(0.18, 0.2)))





