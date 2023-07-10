# this script calculates the genetic distance and does further cleaning and returns a df for analysis
###############################################################################

rm(list=ls())
setwd("/Users/tomversluys/Documents")
# derived_combined <- fread("./NATURE_PAPER/DATA/couples.csv") 
derived_combined <- fread("<2_DATASET_HERE>") # load df output from previous script (2)

# Load required libraries
library(dplyr)
library(data.table)
library(tidyr)

# Preprocessing
# set each couple's pregnancy variables to female values
# set each male's chldren to male values
# set NAs to 0
###############################################################################
data_osb <- derived_combined %>%
  filter(sexcomposition == "OppositeSex") %>% 
  
  mutate(childfatherered = if_else(sex == 1, childfatherered, CHILDFATHERERED), 
         number_of_live_births = if_else(sex == 0, number_of_live_births, NUMBER_OF_LIVE_BIRTHS), 
         number_of_stillbirths = if_else(sex == 0, number_of_stillbirths, NUMBER_OF_STILLBIRTHS),
         number_of_miscarriages = if_else(sex == 0, number_of_miscarriages, NUMBER_OF_MISCARRIAGES),
         
         number_of_live_births = ifelse(is.na(number_of_live_births), 0, number_of_live_births), 
         number_of_stillbirths = ifelse(is.na(number_of_stillbirths), 0, number_of_stillbirths), 
         number_of_miscarriages = ifelse(is.na(number_of_miscarriages), 0, number_of_miscarriages), 
         
         total_pregnancies = sum(number_of_live_births, number_of_stillbirths, number_of_miscarriages, na.rm = T)) %>%
  filter(childfatherered == number_of_live_births) %>%
  filter(! number_of_live_births < 0) %>%
  filter(! number_of_stillbirths < 0) %>%
  filter(! number_of_miscarriages < 0) %>%
  filter(! total_pregnancies < 0,
         ! childfatherered < 0) %>%
  
  drop_na(number_of_live_births, number_of_stillbirths, number_of_miscarriages,
          total_pregnancies, childfatherered) %>%
  mutate(ethnicity = ifelse(ethnicity %in% c("1", "1001", "1002", "1003"), "White", ethnicity), 
         ETHNICITY = ifelse(ETHNICITY %in% c("1", "1001", "1002", "1003"), "White", ethnicity)) %>%
  data.table()

# Select the required columns for distance calculation
###############################################################################
data_osb1 <- data_osb %>% dplyr::select(pc1:pc40, PC1:PC40, id)

# Function to calculate Euclidean distance
###############################################################################
euclidean_distance <- function(partner1, partner2) {
  sqrt(sum((partner1 - partner2)^2))
}

# Calculate Euclidean distances for increasing number of principal components
###############################################################################
distances <- data.frame(matrix(nrow = nrow(data_osb1), ncol = 0))

for (i in 1:40) {
  pc_cols_partner1 <- paste0("pc", 1:i)
  pc_cols_partner2 <- paste0("PC", 1:i)
  distance_col_name <- paste0("distance_", i)
  distances[[distance_col_name]] <- apply(data_osb1, 1, function(row) {
    partner1_pcs <- as.numeric(row[pc_cols_partner1])
    partner2_pcs <- as.numeric(row[pc_cols_partner2])
    dist(rbind(partner1_pcs, partner2_pcs))
  })
}

# Find the minimum value across all distance variables
###############################################################################
min_value <- min(data_osb1[, sapply(distances, is.numeric)])

# Add a constant to ensure all values are non-negative
constant <- abs(min_value) + 1

# Add the constant to each distance variable and apply the log transformation
distances <- distances %>%
  mutate_at(vars(starts_with("distance_")), ~ log(. + constant))

check1 <- distances %>% cbind(data_osb) %>% data.table()


## filter and modify further
###############################################################################
check2 <- check1 %>%
  arrange(id) %>%
  filter(sexcomposition == "OppositeSex") %>%
  mutate(joint_income = c(household_income + HOUSEHOLD_INCOME)/2) %>%
  mutate(joint_income = ifelse(is.na(joint_income), mean(joint_income, na.rm = T), joint_income)) %>%
  mutate(smoking_status_female = ifelse(sex == 0, smokingstatus, SMOKINGSTATUS)) %>%
  mutate(female_previous_smoker = ifelse(smoking_status_female == 1, "yes", "no")) %>%
  mutate(health_female = ifelse(sex == 0, overallhealth, OVERALLHEALTH)) %>%
  mutate(health_male = ifelse(sex == 1, overallhealth, OVERALLHEALTH)) %>%
  mutate(old_qual_female = ifelse(sex == 0, old_qual, OLD_QUAL)) %>%
  mutate(old_qual_male = ifelse(sex == 1, old_qual, OLD_QUAL)) %>%
  mutate(old_bmi_female = ifelse(sex == 0, old_bmi, OLD_BMI)) %>%
  mutate(old_bmi_male = ifelse(sex == 1, old_bmi, OLD_BMI)) %>%
  mutate(centre = as.factor(centre)) %>%
  mutate(birth_lat_female = ifelse(sex == 0, placebirth_lat, PLACEBIRTH_LAT)) %>%
  mutate(birth_lat_male = ifelse(sex == 1, placebirth_lat, PLACEBIRTH_LAT)) %>%
  mutate(birth_long_female = ifelse(sex == 0, placebirth_long, PLACEBIRTH_LONG)) %>%
  mutate(birth_long_male = ifelse(sex == 1, placebirth_long, PLACEBIRTH_LONG)) %>%
  mutate(Number_of_minorities = case_when(ethnicity == "White" & ETHNICITY == "White" ~ 0,
                                          ethnicity == "White" & !ETHNICITY == "White"~ 1,
                                          !ethnicity == "White" & ETHNICITY == "White" ~ 1,
                                          !ethnicity == "White" & !ETHNICITY == "White" ~ 2)) %>%
  rowwise %>%
  mutate(total_complications = sum(number_of_miscarriages, 
                                   number_of_stillbirths, 
                                   na.rm = T)) %>%
  ungroup() %>%
  mutate(born_uk = if_else(between(country, 1, 5) & between(COUNTRY, 1, 5), "yes", "no")) %>%
  mutate(birth_location = ifelse(c(placebirth_long == PLACEBIRTH_LONG) & c(placebirth_lat==PLACEBIRTH_LAT), "same_birth", "different_birth")) %>%
  mutate(birth_continent_match = case_when(birth_continent == "Europe" & BIRTH_CONTINENT == "Europe" ~ "Europe",
                                           birth_continent == "Asia" & BIRTH_CONTINENT == "Asia" ~ "Asia",
                                           birth_continent == "North_America" & BIRTH_CONTINENT == "North_America" ~ "North_America",
                                           birth_continent == "Oceania" & BIRTH_CONTINENT == "Oceania" ~ "Oceania",
                                           birth_continent == "South_America" & BIRTH_CONTINENT == "South_America" ~ "South_America",
                                           birth_continent == "Africa" & BIRTH_CONTINENT == "Africa" ~ "Africa",
                                           TRUE ~ "Mismatch")) %>%
  mutate(age_female_recruitment = ifelse(sex == 0, recruitment_age, RECRUITMENT_AGE)) %>%
  mutate(age_male_recruitment = ifelse(sex == 1, recruitment_age, RECRUITMENT_AGE)) %>%
  mutate(age_male_abs = ifelse(sex == 1, old_dob, OLD_DOB)) %>%
  data.table()

check2 <- check2 %>% mutate(ethnicity_female = if_else(sex == 0, ethnicity, ETHNICITY)) %>%
  mutate(ethnicity_male = if_else(sex == 1, ethnicity, ETHNICITY)) %>%
  mutate(minority_female = if_else(ethnicity_female %in% c("White"), "no", "yes")) %>%
  mutate(minority_male = if_else(ethnicity_male %in% c("White"), "no", "yes")) %>%
  data.table()

# Create a long format data frame
###############################################################################
check2 <- check2 %>% 
  mutate(Ethnic_subset = as.factor(if_else(ethnicity == "White" & ETHNICITY == "White", "All white", "All ethnicities"))) %>%
  dplyr::select(id, ID, distance_1:distance_40, home_lat, home_long, birth_lat_female, birth_long_female,
                sexcomposition, number_of_miscarriages, number_of_stillbirths, 
                birth_lat_male, birth_long_male, Number_of_minorities, 
                health_female, female_previous_smoker,
                total_pregnancies, age_male_recruitment, health_male, 
                number_of_live_births, total_complications,ethnicity, ETHNICITY, 
                centre, age_female_recruitment, noinhouse, 
                NOINHOUSE, birth_continent, BIRTH_CONTINENT, joint_income, 
               country, COUNTRY, born_uk, household_income, ethnicity_male, ethnicity_female, minority_male, minority_female, 
                age_first_live_birth, age_last_live_birth, age_prim_women_birth_first_child, old_qual_male, old_qual_female, 
               old_bmi_male, old_bmi_female) %>% 
  mutate(total_pregnancies.fac = factor(total_pregnancies), age_female_recruitment = as.numeric(age_female_recruitment)) %>%
  arrange(id)

 # data
###############################################################################
check2 <- check2 %>% mutate(
                            ethnic_match = as.factor(ifelse(ethnicity == ETHNICITY, "same", "different")),
                            different_country = as.factor(ifelse(country == COUNTRY, "same", "different"))) %>%
  mutate(Ethnic_subset = as.factor(if_else(ethnicity == "White" & ETHNICITY == "White", "All white", "All ethnicities"))) %>%
  mutate(white_europe = as.factor(ifelse(birth_continent == "Europe" & ethnicity == "White" & ETHNICITY == "White",
                                         "yes", "no"))) %>%
  mutate(born_uk = as.factor(born_uk)) %>% filter(born_uk %in% c("no", "yes"))  %>%
  mutate(ethnic_country = as.factor(ifelse(different_country == "same" & ethnic_match == "same", "yes", "no"))) %>%
  mutate(country_match = ifelse(country == COUNTRY, "yes", "no")) %>%
  mutate(female_previous_smoker = ifelse(female_previous_smoker == "yes", 1, 0)) %>%
  data.table()

## save 
fwrite(check2, "<SAVE_TO_DIRECTORY>")
# fwrite(check2, "./NATURE_PAPER/DATA/couples_with_genetic_distances.csv")


