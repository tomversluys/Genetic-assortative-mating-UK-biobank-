# cleaning and preprocessing
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(ggplot2)
library(grid)
library(countrycode)
library(purrr)

# data
###############################################################################
derived_combined <- fread("<1_DATASET_HERE>")
# data_original <- fread("./newphd/data/ukbb/ukbb_old_subset.csv")  # this is original data
data <- data_original 

### begin couples identification
# filter data
data0 <- data %>%
  filter(
    AcType > 0 & OwnRent > 0) %>%
  drop_na(HomeE, HomeN)
colnames(data0)

# create grouping variable to match pairs
###############################################################################
data1<-data0%>%
  mutate(geogroup=paste(Centre, sep=""))%>%
  mutate(CoupleID=paste(geogroup, sep="_")) %>%
  group_by(CoupleID) %>%
  mutate(CoupleID = cur_group_id()) %>%
  mutate(N_total=length(CoupleID)) %>%
  mutate(N_males_group = length(Sex[Sex==1]), N_females_group = length(Sex[Sex==0])) %>%
  group_by(geogroup) %>%
  mutate(density_proxy = length(geogroup)) %>%
  dplyr::select(-c(geogroup, density_proxy)) %>% 
  data.table()
length(unique(data1$CoupleID))

data2 <- data1 %>%
  mutate(CoupleID = as.factor(CoupleID)) %>%
  group_by(CoupleID) %>%
  filter(n()>=50) %>%
  mutate(l = length(unique(ID))) %>%
  as.data.table()
colnames(data2)

# Define function for measuring euclidean distance by column and row
###############################################################################
product <- function(x, output){
  a = x[1:40]
  b = x[41:80] # all
  c = dist(rbind(a, b))
  c = log(c)
  return(c)
}

# # Create an empty list to store results
# mean_calculation_list <- vector("list", 5)

# Create an empty numeric vector to store mean distances
median_distances <- numeric()

colnames(sample_data1)
# Run loop 100 times
system.time(for (i in 1:1) {
  random1 <- map_dfr(seq_len(1910), ~ data2 %>%
                       group_by(CoupleID) %>%
                       sample_n(2, replace = F) %>%
                       mutate(newid = sample(1:2000000000, 1, replace = F)) %>%
                       mutate(sample_no = .x) %>%
                       ungroup)
  
  random1 <- random1 %>% mutate(newid = as.factor(newid)) %>%
    filter(newid %in% sample(levels(newid),42000))
  
  sample_data <- random1 %>% 
    rename_all(tolower) %>%
    arrange(newid) %>%
    group_by(newid) %>%
    filter(row_number()==1) %>%
    ungroup() %>%
    select(pc1:pc40) %>% 
    data.table()
  
  sample_data1 <- random1 %>% 
    rename_all(tolower) %>%
    arrange(newid) %>%
    group_by(newid) %>%
    filter(row_number()==2) %>% 
    ungroup() %>%
    select(pc1:pc40) %>% 
    data.table()
  
  sample_both <- cbind(sample_data, sample_data1)
  
  # Calculate distances
  distances <- apply(sample_both, 1, product)
  
  # Calculate mean distance for this iteration and store it in the vector
  median_distances[i] <- median(distances, na.rm = T)
})

median_distances_df <- median_distances_df
# median_distances_df <- fread("./NATURE_PAPER/DATA/median_distances_1000_random.csv")

# Create the ggplot

library(boot)
set.seed(123)  # Set a seed for reproducibility

median_stat <- function(data, indices) {
  # Select the data based on the indices
  resampled_data <- data[indices]
  
  # Calculate and return the median of the resampled data
  return(median(resampled_data))
}

check2_distance_40 <- main_couples_data$distance_40[!is.na(main_couples_data$distance_40)]  # Remove NA values
bs_result <- boot(check2_distance_40, median_stat, R = 1000)  # Run the bootstrap with 1000 replicates
ci <- boot.ci(bs_result, type = "perc", conf = 0.95)$percent
ci

set.seed(123)  # Set a seed for reproducibility
check2_distance_40 <- median_distances_df$median_distances[!is.na(median_distances_df$median_distances)]  # Remove NA values
bs_result <- boot(check2_distance_40, median_stat, R = 1000)  # Run the bootstrap with 1000 replicates
cir <- boot.ci(bs_result, type = "perc", conf = 0.95)$percent
cir

# Filter the data to keep only the values less than or equal to 3

# Calculate the density values
density_values <- density(median_distances_df$median_distances)
# Create a data frame with the density values
density_data <- data.frame(
  median_distances = density_values$x,
  density = density_values$y
)

density_data_filtered <- density_data %>%
  filter(median_distances <= 3.2314)

library(gridExtra)
(original1 <- ggplot(median_distances_df, aes(x = median_distances)) +
  geom_density(fill = "#0072B2", alpha = 0.8) +  # Change the fill color
  geom_ribbon(data = density_data_filtered, aes(y = density, ymin = 0, ymax = density),
              fill = "white", alpha = 1) +
  geom_vline(aes(xintercept = 3.226946), linetype = "dashed", size = 1.5, color = "darkgray") +  # Change line size and color
  geom_vline(aes(xintercept = 3.2314), linetype = "dashed", color = "#D55E00", size = 1.5) +  # Change line size and color
  geom_errorbarh(aes(xmin = ci[4], xmax = ci[5], y = 100)) +
  labs(x = "Median genetic distance", y = "Density") +
  theme_few() + theme(text = element_text(size = 16),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        panel.background = element_rect(fill = "white")))

(inset <- original1 + 
    xlim(3.205, 3.2505) +
    ylim(0, 350) +
    labs(x = "", y = "") +
    geom_label(aes(x = 3.2265, y = 150, label = "Random pairs"), size = 3.5, hjust = 1, vjust = 0, label.padding = unit(0.25, "lines"), label.size = 0.3, fill = "white") +
    geom_label(aes(x = 3.23123, y = 150, label = "True couples"), size = 3.5, hjust = -0.1, vjust = 0, label.padding = unit(0.25, "lines"), label.size = 0.3, fill = "white") +
    theme_few() + theme(text = element_text(size = 9),
                        axis.text = element_text(size = 9, color = "black"),
                        axis.title = element_text(size = 9, color = "black")))

(original1 <- original1 + 
    xlim(2.8, 6) +
    geom_segment(aes(x = 3.4, y = 350, xend = 4.5, yend = 500), linetype = "dashed", color = "gray30", alpha=0.2, size = 0.5) +
    geom_segment(aes(x = 3.4, y = 0, xend = 6, yend = 200), linetype = "dashed", color = "gray30", alpha=0.2 ,size = 0.5) +
    geom_rect(aes(xmin = 3.1, xmax = 3.4, ymin = 0, ymax = 350), 
                                   color='gray30', linetype='dashed', alpha=0, size = 0.5))

(original2 <- original1 + annotation_custom(grob = ggplotGrob(inset), 
                              xmin = 4.5, xmax = 6, ymin = 200, ymax = 500) +
  geom_rect(aes(xmin = 4.5, xmax = 6, ymin = 200, ymax = 500), 
            color='gray30', linetype='dashed', alpha=0, size = 0.5) + 
  ylim(0, 520))


ggsave("./newphd/graphics/plot1.pdf", original2,
       width = 20, height = 14,
       device='tiff', units = "cm", dpi=450)


p_value <- sum(median_distances_df$median_distances <= 3.231225) / length(median_distances_df$median_distances)
print(paste("P-value:", p_value))


p_inset
# Your original plot
p_main <- r_mating + xlim(2.58, 6.2) 
  
  # Inset: Create a zoomed version of your plot by adjusting the x and y scales
(p_inset <- r_mating + 
    xlim(3.2, 3.27))

# Adding the zoomed-in plot as an inset to the original plot
p_main <- p_main + 
  annotation_custom(grob = ggplotGrob(p_inset), xmin = 5, xmax = 6, ymin = 300, ymax = 500) +
  geom_rect(aes(xmin = 4.9, xmax = 6.1, ymin = 250, ymax = 610), color='black', linetype='dashed', alpha=0, size = 0.2)

print(p_main)
