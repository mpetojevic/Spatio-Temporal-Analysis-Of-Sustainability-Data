library(dplyr)
library(tidyr)
library(stringr)
library(rgeoboundaries)
library(class)
library(plotly)
library(dplyr)
library(tools)  
library(leaflet)
library(rnaturalearth)
library(sf)
library(reshape2)
library(dtw)
library(dbscan)
library(RColorBrewer)
library(mclust)
library(cluster)
library(factoextra)
library(spdep)
library(patchwork)
library(tmap)
library(nnet)
library(randomForest)

################################################### LOAD DATA IN WIDE AND LONG FORMAT ##########################################################################

getwd()  # Check the current working directory
setwd("~/Uni/2024S/Bachelor/Spatio-Temporal-Analysis-Of-Sustainability-Data/Shiny-App") 
getwd()
data <- read.csv("sustainability_data_central_europe.csv")

time_series_data <- data %>%
  select(shapeName, starts_with("water_percentage"), starts_with("tree_percentage"),
         starts_with("flooded_vegetation_percentage"), starts_with("crops_percentage"),
         starts_with("built_area_percentage"), starts_with("bare_ground_percentage"),
         starts_with("snow_ice_percentage"), starts_with("clouds_percentage"),
         starts_with("rangeland_percentage"), starts_with("burned_area"),
         starts_with("CO2_total"), starts_with("PM25_total"),
         starts_with("TPC_total"), starts_with("NMHC_total"),
         starts_with("OC_total"), starts_with("CH4_total"),
         starts_with("SO2_total"), starts_with("BC_total")) %>%
  pivot_longer(cols = -shapeName,
               names_to = "parameter",
               values_to = "value") %>%
  mutate(year = as.integer(str_extract(parameter, "\\d{4}")),
         parameter = str_remove(parameter, "_\\d{4}")) %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  group_by(shapeName, year) %>%
  summarize(across(everything(), mean, na.rm = TRUE), .groups = 'drop')


time_series_data <- data %>%
  select(shapeName, shapeGroup, total_area_km2) %>%
  left_join(time_series_data, by = "shapeName")

time_series_data <- time_series_data %>% mutate(country = case_when(
  shapeGroup == "AUT" ~ "Austria",
  shapeGroup == "CZE" ~ "Czechia",
  shapeGroup == "HUN" ~ "Hungary",
  shapeGroup == "DEU" ~ "Germany",
  shapeGroup == "POL" ~ "Poland",
  shapeGroup == "SVK" ~ "Slovakia",
  shapeGroup == "CHE" ~ "Switzerland",
  shapeGroup == "SVN" ~ "Slovenia",
  TRUE ~ NA_character_
)) 

central_europe <- rbind(gb_adm1("Austria"), gb_adm1("Germany"), 
                        gb_adm1("Czech Republic"), gb_adm1("Poland"), 
                        gb_adm1("Slovakia"), gb_adm1("Hungary"), 
                        gb_adm1("Switzerland"), gb_adm1("Slovenia"))

data_with_state_boarders <- central_europe %>% 
  select(shapeName, geometry) %>%
  left_join(data, by = "shapeName")

data_with_state_boarders <- data_with_state_boarders %>% mutate(country = case_when(
  shapeGroup == "AUT" ~ "Austria",
  shapeGroup == "CZE" ~ "Czechia",
  shapeGroup == "HUN" ~ "Hungary",
  shapeGroup == "DEU" ~ "Germany",
  shapeGroup == "POL" ~ "Poland",
  shapeGroup == "SVK" ~ "Slovakia",
  shapeGroup == "CHE" ~ "Switzerland",
  shapeGroup == "SVN" ~ "Slovenia",
  TRUE ~ NA_character_
)) 

combined_data <- read.csv("combined_treecover_loss_data.csv") %>%
  mutate(country = ifelse(country == "Czech Republic", "Czechia", country))  %>%
  group_by(country) %>%
  summarise(tree_loss = sum(umd_tree_cover_loss__ha, na.rm = TRUE))

countries <- ne_countries(scale = "medium", returnclass = "sf")
central_europe_countries <- c("Austria", "Czechia", "Germany", "Hungary", 
                              "Poland", "Slovakia", "Slovenia", "Switzerland")
geo_data <- countries %>%
  left_join(combined_data, by = c("name" = "country")) %>%
  filter(name %in% central_europe_countries) 


geo_data_centroids <- st_centroid(geo_data)

tree_loss_palette <- colorRampPalette(brewer.pal(9, "YlOrRd"))

# Plot with binning color scales for tree loss
ggplot() +
  geom_sf(data = geo_data, aes(fill = tree_loss), color = "black", size = 0.5) +  
  geom_text(data = geo_data_centroids, aes(label = paste(name, "\n -", round(tree_loss, 2), "ha"),
                                           geometry = geometry), 
            stat = "sf_coordinates", size = 3, color = "black", check_overlap = TRUE) + 
  scale_fill_gradientn(colors = tree_loss_palette(10),  
                       breaks = pretty(geo_data$tree_loss, n = 5)) +  
  labs(
    title = "Total Tree Cover Loss Across Central Europe Between 2001 and 2023",
    fill = "Tree Loss (ha)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

############################################################### SPATIAL MORAN I ##################################################################################


nb <- poly2nb(data_with_state_boarders$geometry, queen = TRUE)  # Use 'queen' contiguity for neighbors
lw <- nb2listw(nb, style = "W")  # Spatial weights matrix, row-standardized

# List of the years columns in wide dataframe
year_columns <- c("tree_percentage_2018", "tree_percentage_2019", "tree_percentage_2020", 
                  "tree_percentage_2021", "tree_percentage_2022", "tree_percentage_2023")


#Moran scatter plots

par(mfrow = c(3, 2))
for (year_col in year_columns) {
   tree_percentage <- data_with_state_boarders[[year_col]]
   moran.plot(tree_percentage, lw, main = paste("Moran's I Scatterplot for", year_col))
}


# Function to calculate Moran's I for a given year
calculate_morans_i_for_year <- function(data, value_column, lw) {
  tree_percentage <- data[[value_column]]
  morans_i <- moran.test(tree_percentage, lw)
  return(morans_i)
}


# Store results in a list
morans_i_results <- list()

# Loop through the columns and calculate Moran's I for each year
for (year_col in year_columns) {
  morans_i_results[[year_col]] <- calculate_morans_i_for_year(
    data_with_state_boarders, 
    value_column = year_col, 
    lw = lw
  )
}

morans_i_results

# Extract Moran's I statistic and p-values for each year
morans_i_data <- data.frame(
  year = year_columns,  # List of years
  morans_i_stat = sapply(morans_i_results, function(res) res$estimate["Moran I statistic"]),
  p_value = sapply(morans_i_results, function(res) res$p.value)
)

print(morans_i_data)

library(scales)
ggplot(morans_i_data, aes(x = as.factor(year), y = morans_i_stat, fill = year)) +
  geom_bar(stat = "identity", color = "black") + 
    geom_text(aes(label = paste0("Moran's I: ", round(morans_i_stat, 3))), 
            vjust = -0.5, size = 3.5) +
    geom_text(aes(label = paste0("p-Value: ", scientific(p_value))), 
            vjust = 1.5, color = "white", size = 3) +  
  labs(title = "Moran's I for Tree Percentage by Year (2018-2023)",
       x = "Year",
       y = "Moran's I Statistic") +
  theme_minimal() +
  theme(legend.position = "none")

str(morans_i_data)

# No bars
ggplot(morans_i_data, aes(y = as.factor(year), x = morans_i_stat)) +
  geom_segment(aes(y = as.factor(year), yend = as.factor(year), 
                   x = 0, xend = morans_i_stat), color = "black", size = 1) +
  geom_point(aes(x = morans_i_stat, y = as.factor(year)), color = "black", size = 3) + 
  geom_text(aes(label = paste0("Moran's I: ", round(morans_i_stat, 3))), 
            hjust = -0.3, size = 3.5) +  # Label next to the lines
  geom_text(aes(label = paste0("p-Value: ", scientific(p_value))), 
            vjust = -1, hjust = 1.2, color = "black", size = 3) +  
  labs(title = "Moran's I for Tree Percentage by Year (2018-2023)",
       y = "Year",
       x = "Moran's I Statistic") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, max(morans_i_data$morans_i_stat) + 0.05))  




#################################################### MONTE CARLO GLOBAL MORAN APPROACH #############################################################################

monte_carlo_results <- list()

for (year_col in year_columns) {
  tree_percentage <- data_with_state_boarders[[year_col]]
  monte_carlo_results[[year_col]] <- moran.mc(tree_percentage, lw, nsim = 999)
}

print(monte_carlo_results)
monte_carlo_results[["tree_percentage_2020"]]


monte_carlo_data <- data.frame()

# Loop through each year and collect the simulated Moran's I results
for (year_col in year_columns) {
  monte_carlo_year_data <- data.frame(
    year = year_col,
    simulated_moran = monte_carlo_results[[year_col]]$res,
    observed_moran = monte_carlo_results[[year_col]]$statistic,
    p_value = monte_carlo_results[[year_col]]$p.value
  )
  monte_carlo_data <- rbind(monte_carlo_data, monte_carlo_year_data)
}

monte_carlo_data$year <- gsub("tree_percentage_", "", monte_carlo_data$year)

ggplot(monte_carlo_data, aes(x = simulated_moran)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  geom_vline(aes(xintercept = observed_moran), color = "red", linetype = "dashed", size = 1) +
  facet_wrap(~ year, ncol = 1) +  # Create one row per year, or adjust ncol to change layout
  geom_text(aes(x = Inf, y = Inf, label = paste0("p-value: ", round(p_value, 4))),
            vjust = 2, hjust = 2, size = 3, color = "black", inherit.aes = FALSE) +  # Position and adjust p-value label
  labs(title = "Monte Carlo Simulations for Moran's I (2018-2023)",
       x = "Simulated Moran's I",
       y = "Frequency") +
  theme_minimal()



############################################################## LOCAL MORAN ########################################################################################

calculate_local_morans_i_for_year <- function(data, value_column, lw) {
  tree_percentage <- data[[value_column]]
  lmoran <- localmoran(tree_percentage, lw, alternative = "greater")
  
  return(data.frame(
    Ii = lmoran[, "Ii"],                  # Local Moran's I statistic
    EIi = lmoran[, "E.Ii"],               # Expected value of Ii
    VarIi = lmoran[, "Var.Ii"],           # Variance of Ii
    ZIi = lmoran[, "Z.Ii"],               # Z-scores
    p_value = lmoran[, "Pr(z > E(Ii))"]   # P-values for greater test
  ))
}

local_morans_i_results <- list()

for (year_col in year_columns) {
  local_morans_i_results[[year_col]] <- calculate_local_morans_i_for_year(
    data_with_state_boarders, 
    value_column = year_col, 
    lw = lw
  )
}

for (i in seq_along(year_columns)) {
  year_col <- year_columns[i]
  year_result <- local_morans_i_results[[year_col]]
  data_with_state_boarders[[paste0("Ii_", year_col)]] <- year_result$Ii
  data_with_state_boarders[[paste0("ZIi_", year_col)]] <- year_result$ZIi
  data_with_state_boarders[[paste0("p_value_", year_col)]] <- year_result$p_value
}

year_plots_list <- list()

for (year_col in year_columns) {
  
  year_label <- gsub("tree_percentage_", "", year_col)
  
  # 1. Tree Percentage for the selected year
  p1 <- tm_shape(data_with_state_boarders) +
    tm_polygons(col = year_col, title = paste("Tree Percentage (", year_label, ")", sep = ""), style = "quantile") +
    tm_layout(legend.outside = TRUE)
  
  # 2. Local Moran's I for the selected year
  p2 <- tm_shape(data_with_state_boarders) +
    tm_polygons(col = paste0("Ii_", year_col), title = "Local Moran's I", style = "quantile") +
    tm_layout(legend.outside = TRUE)
  
  # Commented this out bc Z scores give us the same information as p values!
  # 3. Z-scores for the selected year (Categorized into Negative, No, Positive SAC)
  #p3 <- tm_shape(data_with_state_boarders) +
   # tm_polygons(col = paste0("ZIi_", year_col), title = "Z-score",
    #            breaks = c(-Inf, 1.65, Inf),   
     #           labels = c("Z < 1.65", "Z >= 1.65"),  
      #          palette = c("blue", "red")) +  # Blue for Z < 1.65, Red for Z >= 1.65
    #tm_layout(legend.outside = TRUE)
  
  # 4. p-values for the selected year
  p4 <- tm_shape(data_with_state_boarders) +
    tm_polygons(col = paste0("p_value_", year_col), title = "p-value",
                breaks = c(-Inf, 0.05, Inf), 
                palette = c("red", "white")) +  
    tm_layout(legend.outside = TRUE)
  
  combined_plot <- tmap_arrange(p1, p2, p4, ncol = 3)
  
  year_plots_list[[year_label]] <- combined_plot
}

tmap_mode("view")
year_plots_list[["2018"]]
year_plots_list[["2019"]]
year_plots_list[["2020"]]
year_plots_list[["2021"]]
year_plots_list[["2022"]]
year_plots_list[["2023"]]


############################################################ TWO SIDED LOCAL MORAN #############################################################################

calculate_local_morans_i_for_year_two_sided <- function(data, value_column, lw) {
  tree_percentage <- data[[value_column]]
  
  lmoran <- localmoran(tree_percentage, lw, alternative = "two.sided")
  
  return(data.frame(
    Ii = lmoran[, "Ii"],                   # Local Moran's I statistic
    EIi = lmoran[, "E.Ii"],                # Expected value of Ii
    VarIi = lmoran[, "Var.Ii"],            # Variance of Ii
    ZIi = lmoran[, "Z.Ii"],                # Z-scores
    p_value = lmoran[, "Pr(z != E(Ii))"]   # P-values for greater test
  ))
}

local_morans_i_results_two_sided <- list()

for (year_col in year_columns) {
  local_morans_i_results_two_sided[[year_col]] <- calculate_local_morans_i_for_year_two_sided(
    data_with_state_boarders, 
    value_column = year_col, 
    lw = lw
  )
}

for (i in seq_along(year_columns)) {
  year_col <- year_columns[i]
  year_result <- local_morans_i_results_two_sided[[year_col]]
  data_with_state_boarders[[paste0("Ii_ts_", year_col)]] <- year_result$Ii
  data_with_state_boarders[[paste0("ZIi_ts_", year_col)]] <- year_result$ZIi
  data_with_state_boarders[[paste0("p_value_ts_", year_col)]] <- year_result$p_value
}

sac_plots_list <- list()

for (year_col in year_columns) {
  
  year_label <- gsub("tree_percentage_", "", year_col)

  sac_plot <- ggplot(data_with_state_boarders) +
    geom_sf(aes(fill = cut(.data[[paste0("ZIi_ts_", year_col)]],
                           breaks = c(-Inf, -1.96, 1.96, Inf),
                           labels = c("Negative SAC", "No SAC", "Positive SAC"))),
            color = "black") +
    scale_fill_manual(
      values = c("Negative SAC" = "blue", 
                 "No SAC" = "white", 
                 "Positive SAC" = "red")
    ) +
    labs(title = paste("Z-score (", year_label, ")", sep = ""), fill = "Spatial Auto-correlation") +
    theme_minimal()
  
    sac_plots_list[[year_label]] <- sac_plot
}

# Arrange with patchwork
combined_sac_plots <- wrap_plots(sac_plots_list, ncol = 2)  # Arrange in a 2-column grid

print(combined_sac_plots)


############################################################ CLUSTERS BASED ON TWO SIDED LOCAL MORAN ##############################################################

#2018

mp <- moran.plot(as.vector(scale(data_with_state_boarders$tree_percentage_2018)), lw, ylab = "Spatially Lagged Scaled tree Percentage for Year 2018", xlab = "Scaled Tree Percentage in 2018")
head(mp)

outliers <- which(mp$is_inf) 
text(mp$x[outliers], mp$wx[outliers], 
     labels = data_with_state_boarders$shapeName[outliers], 
     pos = 4, col = "blue", cex = 0.8) 


look_up <- data_with_state_boarders %>%
  mutate(x = mp$x, wp = mp$wx) %>%
  select(shapeName, tree_percentage_2018, x, wp, p_value_ts_tree_percentage_2018) %>%
  filter(p_value_ts_tree_percentage_2018 < 0.05)


data_with_state_boarders$quadrant_2018 <- NA
# high-high
data_with_state_boarders[(mp$x >= 0 & mp$wx >= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2018 <= 0.05), "quadrant_2018"]<- 1
# low-low
data_with_state_boarders[(mp$x <= 0 & mp$wx <= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2018 <= 0.05), "quadrant_2018"]<- 2
# high-low
data_with_state_boarders[(mp$x >= 0 & mp$wx <= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2018 <= 0.05), "quadrant_2018"]<- 3
# low-high
data_with_state_boarders[(mp$x <= 0 & mp$wx >= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2018 <= 0.05), "quadrant_2018"]<- 4
# non-significant
data_with_state_boarders[(data_with_state_boarders$p_value_ts_tree_percentage_2018 > 0.05), "quadrant_2018"] <- 5

tmap_mode("view")
tm_shape(data_with_state_boarders) + 
  tm_fill(col = "quadrant_2018", title = "",
          breaks = c(1, 2, 3, 4, 5, 6),
          palette = c("red", "blue", "lightpink", "skyblue2", "white"),
          labels = c("High-High", "Low-Low", "High-Low", "Low-High", "Non-significant"),
          alpha = 0.7) +
  tm_borders(alpha = 0.5) +
  tm_layout(frame = FALSE, title = "Clusters", legend.outside = TRUE)

#2019

mp <- moran.plot(as.vector(scale(data_with_state_boarders$tree_percentage_2019)), lw)
head(mp)

look_up <- data_with_state_boarders %>%
  mutate(x = mp$x, wp = mp$wx) %>%
  select(shapeName, tree_percentage_2019, x, wp, p_value_ts_tree_percentage_2019) %>%
  filter(p_value_ts_tree_percentage_2019 < 0.05)


data_with_state_boarders$quadrant_2019 <- NA
# high-high
data_with_state_boarders[(mp$x >= 0 & mp$wx >= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2019 <= 0.05), "quadrant_2019"]<- 1
# low-low
data_with_state_boarders[(mp$x <= 0 & mp$wx <= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2019 <= 0.05), "quadrant_2019"]<- 2
# high-low
data_with_state_boarders[(mp$x >= 0 & mp$wx <= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2019 <= 0.05), "quadrant_2019"]<- 3
# low-high
data_with_state_boarders[(mp$x <= 0 & mp$wx >= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2019 <= 0.05), "quadrant_2019"]<- 4
# non-significant
data_with_state_boarders[(data_with_state_boarders$p_value_ts_tree_percentage_2019 > 0.05), "quadrant_2019"] <- 5

tmap_mode("view")
tm_shape(data_with_state_boarders) + 
  tm_fill(col = "quadrant_2019", title = "",
          breaks = c(1, 2, 3, 4, 5, 6),
          palette = c("red", "blue", "lightpink", "skyblue2", "white"),
          labels = c("High-High", "Low-Low", "High-Low", "Low-High", "Non-significant"),
          alpha = 0.7) +
  tm_borders(alpha = 0.5) +
  tm_layout(frame = FALSE, title = "Clusters", legend.outside = TRUE)


#2020

mp <- moran.plot(as.vector(scale(data_with_state_boarders$tree_percentage_2020)), lw)
head(mp)

data_with_state_boarders$quadrant_2020 <- NA
# high-high
data_with_state_boarders[(mp$x >= 0 & mp$wx >= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2020 <= 0.05), "quadrant_2020"]<- 1
# low-low
data_with_state_boarders[(mp$x <= 0 & mp$wx <= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2020 <= 0.05), "quadrant_2020"]<- 2
# high-low
data_with_state_boarders[(mp$x >= 0 & mp$wx <= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2020 <= 0.05), "quadrant_2020"]<- 3
# low-high
data_with_state_boarders[(mp$x <= 0 & mp$wx >= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2020 <= 0.05), "quadrant_2020"]<- 4
# non-significant
data_with_state_boarders[(data_with_state_boarders$p_value_ts_tree_percentage_2020 > 0.05), "quadrant_2020"] <- 5

tmap_mode("view")
tm_shape(data_with_state_boarders) + 
  tm_fill(col = "quadrant_2020", title = "",
          breaks = c(1, 2, 3, 4, 5, 6),
          palette = c("red", "blue", "lightpink", "skyblue2", "white"),
          labels = c("High-High", "Low-Low", "High-Low", "Low-High", "Non-significant"),
          alpha = 0.7) +
  tm_borders(alpha = 0.5) +
  tm_layout(frame = FALSE, title = "Clusters", legend.outside = TRUE)


#2021

mp <- moran.plot(as.vector(scale(data_with_state_boarders$tree_percentage_2021)), lw)
head(mp)

data_with_state_boarders$quadrant_2021 <- NA
# high-high
data_with_state_boarders[(mp$x >= 0 & mp$wx >= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2021 <= 0.05), "quadrant_2021"]<- 1
# low-low
data_with_state_boarders[(mp$x <= 0 & mp$wx <= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2021 <= 0.05), "quadrant_2021"]<- 2
# high-low
data_with_state_boarders[(mp$x >= 0 & mp$wx <= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2021 <= 0.05), "quadrant_2021"]<- 3
# low-high
data_with_state_boarders[(mp$x <= 0 & mp$wx >= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2021 <= 0.05), "quadrant_2021"]<- 4
# non-significant
data_with_state_boarders[(data_with_state_boarders$p_value_ts_tree_percentage_2021 > 0.05), "quadrant_2021"] <- 5

tmap_mode("view")
tm_shape(data_with_state_boarders) + 
  tm_fill(col = "quadrant_2021", title = "",
          breaks = c(1, 2, 3, 4, 5, 6),
          palette = c("red", "blue", "lightpink", "skyblue2", "white"),
          labels = c("High-High", "Low-Low", "High-Low", "Low-High", "Non-significant"),
          alpha = 0.7) +
  tm_borders(alpha = 0.5) +
  tm_layout(frame = FALSE, title = "Clusters", legend.outside = TRUE)

#2022

mp <- moran.plot(as.vector(scale(data_with_state_boarders$tree_percentage_2022)), lw)
head(mp)

data_with_state_boarders$quadrant_2022 <- NA
# high-high
data_with_state_boarders[(mp$x >= 0 & mp$wx >= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2022 <= 0.05), "quadrant_2022"]<- 1
# low-low
data_with_state_boarders[(mp$x <= 0 & mp$wx <= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2022 <= 0.05), "quadrant_2022"]<- 2
# high-low
data_with_state_boarders[(mp$x >= 0 & mp$wx <= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2022 <= 0.05), "quadrant_2022"]<- 3
# low-high
data_with_state_boarders[(mp$x <= 0 & mp$wx >= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2022 <= 0.05), "quadrant_2022"]<- 4
# non-significant
data_with_state_boarders[(data_with_state_boarders$p_value_ts_tree_percentage_2022 > 0.05), "quadrant_2022"] <- 5

tmap_mode("view")
tm_shape(data_with_state_boarders) + 
  tm_fill(col = "quadrant_2022", title = "",
          breaks = c(1, 2, 3, 4, 5, 6),
          palette = c("red", "blue", "lightpink", "skyblue2", "white"),
          labels = c("High-High", "Low-Low", "High-Low", "Low-High", "Non-significant"),
          alpha = 0.7) +
  tm_borders(alpha = 0.5) +
  tm_layout(frame = FALSE, title = "Clusters", legend.outside = TRUE)

#2023

mp <- moran.plot(as.vector(scale(data_with_state_boarders$tree_percentage_2023)), lw)
head(mp)

data_with_state_boarders$quadrant_2023 <- NA
# high-high
data_with_state_boarders[(mp$x >= 0 & mp$wx >= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2023 <= 0.05), "quadrant_2023"]<- 1
# low-low
data_with_state_boarders[(mp$x <= 0 & mp$wx <= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2023 <= 0.05), "quadrant_2023"]<- 2
# high-low
data_with_state_boarders[(mp$x >= 0 & mp$wx <= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2023 <= 0.05), "quadrant_2023"]<- 3
# low-high
data_with_state_boarders[(mp$x <= 0 & mp$wx >= 0) & (data_with_state_boarders$p_value_ts_tree_percentage_2023 <= 0.05), "quadrant_2023"]<- 4
# non-significant
data_with_state_boarders[(data_with_state_boarders$p_value_ts_tree_percentage_2023 > 0.05), "quadrant_2023"] <- 5

tmap_mode("view")
tm_shape(data_with_state_boarders) + 
  tm_fill(col = "quadrant_2023", title = "",
          breaks = c(1, 2, 3, 4, 5, 6),
          palette = c("red", "blue", "lightpink", "skyblue2", "white"),
          labels = c("High-High", "Low-Low", "High-Low", "Low-High", "Non-significant"),
          alpha = 0.7) +
  tm_borders(alpha = 0.5) +
  tm_layout(frame = FALSE, title = "Clusters", legend.outside = TRUE)



############################################################# TEMPORAL AUTOCORRELATION ############################################################################


library(ggplot2)
ggplot(time_series_data, aes(x = year, y = tree_percentage, group = shapeName, color = shapeName)) +
  geom_line() +
  labs(title = "Tree Percentage Over Time", x = "Year", y = "Tree Percentage") +
  theme_minimal()

ggplot(time_series_data, aes(x = year, y = tree_percentage, group = shapeName, color = shapeName)) +
  geom_line() +
  labs(title = "Tree Percentage Over Time by Country", x = "Year", y = "Tree Percentage") +
  facet_wrap(~ country) +  
  theme_minimal() +
  theme(legend.position = "none")  


data_summary <- time_series_data %>%
  group_by(year) %>%
  summarize(mean_tree_percentage = mean(tree_percentage))

ggplot(data_summary, aes(x = year, y = mean_tree_percentage)) +
  geom_line(color = "blue") +
  labs(title = "Average Tree Percentage Over Time", x = "Year", y = "Mean Tree Percentage") +
  theme_minimal()

# Lagged correlation
time_series_data <- time_series_data %>%
  arrange(shapeName, year) 

# Lagged variable for tree percentage
time_series_data <- time_series_data %>%
  group_by(shapeName) %>%
  mutate(lag_tree_percentage = lag(tree_percentage))

# Correlation between tree_percentage and lag_tree_percentage
correlation_results <- time_series_data %>%
  group_by(shapeName) %>%  
  summarise(correlation = cor(tree_percentage, lag_tree_percentage, use = "complete.obs"))  

print(correlation_results)

# Scatter plot to visualize relationship between tree_percentage and lagged values
ggplot(time_series_data, aes(x = lag_tree_percentage, y = tree_percentage, color = shapeName)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  
  labs(title = "Tree Percentage vs Lagged Tree Percentage",
       x = "Lagged Tree Percentage",
       y = "Tree Percentage") +
  theme_minimal() +
  theme(legend.position = "none")  

# Bar plot showing correlation for each region
ggplot(correlation_results, aes(x = reorder(shapeName, correlation), y = correlation, fill = correlation)) +
  geom_col() +
  coord_flip() +  
  labs(title = "Correlation between Tree Percentage and Lagged Tree Percentage by Region",
       x = "Region",
       y = "Correlation") +
  scale_fill_gradient2(low = "darkred", mid = "yellow", high = "blue", midpoint = 0) +  
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),  
    axis.title = element_text(size = 14), 
    plot.title = element_text(size = 16),  
    panel.grid.major = element_line(color = "gray", size = 0.3)  
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")  

######################################## POSITIVE TEMP AUTOCORRELATION - LOSS OR GROWTH ? ##########################################################################

# Filter regions with negative correlation
positive_corr <- correlation_results %>%
  filter(correlation > 0)

# Slope of tree percentage over time for each region
trend_results <- time_series_data %>%
  filter(shapeName %in% positive_corr$shapeName) %>%
  group_by(shapeName) %>%
  summarise(slope = coef(lm(tree_percentage ~ year))[2])  # Extract the slope from linear regression

# Merge the slope with the correlation results
positive_corr_trend <- positive_corr %>%
  left_join(trend_results, by = "shapeName")

ggplot(positive_corr_trend, aes(x = reorder(shapeName, slope), y = slope, fill = slope)) +
  geom_col() +
  coord_flip() +  # Flip axes for better readability
  labs(title = "Tree Growth or Loss in Regions with Positive Correlation",
       x = "Region",
       y = "Slope of Tree Percentage Over Time") +
  scale_fill_gradient2(low = "darkred", mid = "yellow", high = "darkgreen", midpoint = 0) +  # Color gradient for slopes
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),  
    axis.title = element_text(size = 14), 
    plot.title = element_text(size = 16),  
    panel.grid.major = element_line(color = "gray", size = 0.3)
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")  

ggplot(time_series_data %>% filter(shapeName == "Győr-Moson-Sopron"), 
       aes(x = year, y = tree_percentage, group = shapeName, color = shapeName)) +
  geom_line() +
  labs(title = "Tree Percentage Over Time for Győr-Moson-Sopron", x = "Year", y = "Tree Percentage") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(15, 18)) 
ggsave("gyor_plot.png", width = 6, height = 3) 


ggplot(time_series_data %>% filter(shapeName == "Thüringen"), 
       aes(x = year, y = tree_percentage, group = shapeName, color = shapeName)) +
  geom_line(size = 1) +
  labs(title = "Tree Percentage Over Time for Thüringen", x = "Year", y = "Tree Percentage") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(28, 33))

ggsave("thuringen_plot.png", width = 6, height = 3) 

time_series_data_with_geometry <- data_with_state_boarders %>% 
  select(shapeName, geometry) %>%
  left_join(time_series_data, by = "shapeName")

# correlation with crops percentage, rangeland percentage, built area percentage
correlation_results <- time_series_data_with_geometry %>%
  group_by(shapeName) %>%
  summarise(
    cor_crops = cor(tree_percentage, crops_percentage, use = "complete.obs"),
    cor_rangeland = cor(tree_percentage, rangeland_percentage, use = "complete.obs"),
    cor_built_area = cor(tree_percentage, built_area_percentage, use = "complete.obs")
  )

library(tmap)
tmap_mode("view")
tm_shape(correlation_results) +
  tm_polygons(col = "cor_built_area", 
              palette = "-RdBu",  # Red to blue palette
              title = "Correlation: Tree vs Built Area",
              style = "quantile", 
              n = 5) + 
  tm_layout(legend.outside = TRUE)

tm_shape(correlation_results) +
  tm_polygons(col = "cor_crops", 
              palette = "-RdBu",  # Red to blue palette
              title = "Correlation: Tree vs Crops Percentage",
              style = "quantile", 
              n = 5) + 
  tm_layout(legend.outside = TRUE)

tm_shape(correlation_results) +
  tm_polygons(col = "cor_rangeland", 
              palette = "-RdBu",  # Red to blue palette
              title = "Correlation: Tree vs Rangeland Percentage",
              style = "quantile", 
              n = 5) + 
  tm_layout(legend.outside = TRUE)

# Plot for crops_percentage
ggplot(correlation_results, aes(x = reorder(shapeName, cor_crops), y = cor_crops, fill = cor_crops)) +
  geom_col() +
  coord_flip() +  # Flip the axes for better readability
  labs(title = "Correlation between Tree Percentage and Crops Percentage by Region",
       x = "Region",
       y = "Correlation (Tree vs Crops)") +
  scale_fill_gradient2(low = "darkred", mid = "white", high = "blue", midpoint = 0) +
  theme_minimal()

# Plot for rangeland_percentage
ggplot(correlation_results, aes(x = reorder(shapeName, cor_rangeland), y = cor_rangeland, fill = cor_rangeland)) +
  geom_col() +
  coord_flip() +  # Flip the axes for better readability
  labs(title = "Correlation between Tree Percentage and Rangeland Percentage by Region",
       x = "Region",
       y = "Correlation (Tree vs Rangeland)") +
  scale_fill_gradient2(low = "darkred", mid = "white", high = "blue", midpoint = 0) +
  theme_minimal()

# Plot for built_area_percentage
ggplot(correlation_results, aes(x = reorder(shapeName, cor_built_area), y = cor_built_area, fill = cor_built_area)) +
  geom_col() +
  coord_flip() +  # Flip the axes for better readability
  labs(title = "Correlation between Tree Percentage and Built Area Percentage by Region",
       x = "Region",
       y = "Correlation (Tree vs Built Area)") +
  scale_fill_gradient2(low = "darkred", mid = "white", high = "blue", midpoint = 0) +
  theme_minimal()



######################################################## LINEAR REGRESSION ###########################################################################

model <- lm(tree_percentage ~ built_area_percentage + crops_percentage + rangeland_percentage + water_percentage
            + flooded_vegetation_percentage + bare_ground_percentage + snow_ice_percentage + burned_area
            + CO2_total + PM25_total + TPC_total + NMHC_total + OC_total + CH4_total + SO2_total
            + BC_total, data = time_series_data)

plot(model, which = 1)  # not the best fit
time_series_data$predicted_lm <- predict(model)

ggplot(time_series_data, aes(x = predicted_lm, y = tree_percentage)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Actual Tree Percentage Linear Model",
       x = "Predicted Tree Percentage",
       y = "Actual Tree Percentage") +
  theme_minimal() # bad fitting


############################################################# LINEAR TRAIN + TEST ################################################################################

train_data <- time_series_data %>% filter(year < 2023)
test_data <- time_series_data %>% filter(year == 2023)

# Prepare data matrices
lm_model <- lm(tree_percentage ~ built_area_percentage + crops_percentage + rangeland_percentage + water_percentage
               + flooded_vegetation_percentage + bare_ground_percentage + snow_ice_percentage + burned_area
               + CO2_total + PM25_total + TPC_total + NMHC_total + OC_total + CH4_total + SO2_total + BC_total,
               data = train_data)

# Make predictions on the test data (2023)
test_data$predicted_lm <- predict(lm_model, newdata = test_data)

rmse_lm <- sqrt(mean((test_data$tree_percentage - test_data$predicted_lm)^2))
print(paste("RMSE for Linear Model on test data (2023):", rmse_lm))

ggplot(test_data, aes(x = predicted_lm, y = tree_percentage)) +
  geom_point(color = "gray", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Actual Tree Percentage (Linear Model, Test Data 2023)",
       x = "Predicted Tree Percentage",
       y = "Actual Tree Percentage") +
  theme_minimal()


####################################################### RANDOM FOREST REFGRESSION ################################################################################

model_rf <- randomForest(tree_percentage ~ built_area_percentage + crops_percentage + rangeland_percentage + water_percentage
                         + flooded_vegetation_percentage + bare_ground_percentage + snow_ice_percentage + burned_area
                         + CO2_total + PM25_total + TPC_total + NMHC_total + OC_total + CH4_total + SO2_total
                         + BC_total, data = time_series_data)
summary(model_rf)
varImpPlot(model_rf, main = "Variable Importance in Random Forest")


time_series_data$predicted_rf <- predict(model_rf)

ggplot(time_series_data, aes(x = predicted_rf, y = tree_percentage)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Actual Tree Percentage Random Forest",
       x = "Predicted Tree Percentage",
       y = "Actual Tree Percentage") +
  theme_minimal() # bad fitting


######################################################### RF TRAIN + TEST ################################################################################

train_data <- time_series_data %>% filter(year < 2023)
test_data <- time_series_data %>% filter(year == 2023)

# Fit Random Forest model on training data
model_rf <- randomForest(tree_percentage ~ built_area_percentage + crops_percentage + rangeland_percentage + water_percentage
                         + flooded_vegetation_percentage + bare_ground_percentage + snow_ice_percentage + burned_area
                         + CO2_total + PM25_total + TPC_total + NMHC_total + OC_total + CH4_total + SO2_total
                         + BC_total, data = train_data, ntree = 500)  # Increase ntree for better results

test_data$predicted_rf <- predict(model_rf, newdata = test_data)

rmse_rf <- sqrt(mean((test_data$tree_percentage - test_data$predicted_rf)^2))
print(paste("RMSE for Random Forest on test data (2023):", rmse_rf))

ggplot(test_data, aes(x = predicted_rf, y = tree_percentage)) +
  geom_point(alpha = 0.5, color = "gray") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Actual Tree Percentage (Random Forest, Test Data 2023)",
       x = "Predicted Tree Percentage",
       y = "Actual Tree Percentage") +
  theme_minimal()

# Plot variable importance
varImpPlot(model_rf, main = "Variable Importance in Random Forest")



######################################################### NNET REGRESSION #################################################################################

nn_model <- nnet(tree_percentage ~ built_area_percentage + crops_percentage + rangeland_percentage + water_percentage +
                   flooded_vegetation_percentage + bare_ground_percentage + snow_ice_percentage + burned_area +
                   CO2_total + PM25_total + TPC_total + NMHC_total + OC_total + CH4_total + SO2_total + BC_total,
                 data = time_series_data, size = 10, linout = TRUE)

time_series_data$predicted_nnet <- predict(nn_model, time_series_data)

ggplot(time_series_data, aes(x = predicted_nnet, y = tree_percentage)) +
  geom_point(color = "gray", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +  # Reference line (y = x)
  labs(title = "Predicted vs Actual Tree Percentage (Neural Network)",
       x = "Predicted Tree Percentage",
       y = "Actual Tree Percentage") +
  theme_minimal()


######################################################## NNET TEST TRAIN ##################################################################################

train_data <- time_series_data %>% filter(year < 2023)
test_data <- time_series_data %>% filter(year == 2023)

nn_model <- nnet(tree_percentage ~ built_area_percentage + crops_percentage + rangeland_percentage + water_percentage +
                   flooded_vegetation_percentage + bare_ground_percentage + snow_ice_percentage + burned_area +
                   CO2_total + PM25_total + TPC_total + NMHC_total + OC_total + CH4_total + SO2_total + BC_total,
                 data = train_data, size = 10, linout = TRUE, maxit = 500)

# Make predictions on the test data (2023)
test_data$predicted_nnet <- predict(nn_model, test_data)

rmse_nnet <- sqrt(mean((test_data$tree_percentage - test_data$predicted_nnet)^2))
print(paste("RMSE for Neural Network on test data (2023):", rmse_nnet))

ggplot(test_data, aes(x = predicted_nnet, y = tree_percentage)) +
  geom_point(color = "gray", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +  # Reference line (y = x)
  labs(title = "Predicted vs Actual Tree Percentage (Neural Network, Test Data 2023)",
       x = "Predicted Tree Percentage",
       y = "Actual Tree Percentage") +
  theme_minimal()



########################################################## XGBOOST REGRESSION ######################################################################

library(xgboost)

# Convert to matrix format for xgboost
data_matrix <- model.matrix(tree_percentage ~ built_area_percentage + crops_percentage + rangeland_percentage + water_percentage +
                              flooded_vegetation_percentage + bare_ground_percentage + snow_ice_percentage + burned_area +
                              CO2_total + PM25_total + TPC_total + NMHC_total + OC_total + CH4_total + SO2_total + BC_total,
                            data = time_series_data)

xgb_model <- xgboost(data = data_matrix, label = time_series_data$tree_percentage, nrounds = 100, objective = "reg:squarederror")

time_series_data$predicted_xgboost <- predict(xgb_model, data_matrix)

#time_series_data$residuals_xgboost <- time_series_data$tree_percentage - time_series_data$predicted_tree_percentage

ggplot(time_series_data, aes(x = predicted_xgboost, y = tree_percentage)) +
  geom_point(color = "gray", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +  # Reference line (y = x)
  labs(title = "Predicted vs Actual Tree Percentage (XGBoost)",
       x = "Predicted Tree Percentage",
       y = "Actual Tree Percentage") +
  theme_minimal()

importance_matrix <- xgb.importance(feature_names = colnames(data_matrix), model = xgb_model)

print(importance_matrix)

#xgb.plot.importance(importance_matrix, top_n = 10, measure = "Gain")

# If you want a custom ggplot for better control over styling
importance_df <- as.data.frame(importance_matrix)

ggplot(importance_df, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Feature Importance for Tree Percentage Prediction",
       x = "Features",
       y = "Importance (Gain)") +
  theme_minimal()


######################################### TEST XGBOOST ##########################################################################################################

train_data <- time_series_data %>% filter(year < 2023)
test_data <- time_series_data %>% filter(year == 2023)

# Prepare data matrices
train_matrix <- model.matrix(tree_percentage ~ built_area_percentage + crops_percentage + rangeland_percentage + water_percentage +
                               flooded_vegetation_percentage + bare_ground_percentage + snow_ice_percentage + burned_area +
                               CO2_total + PM25_total + TPC_total + NMHC_total + OC_total + CH4_total + SO2_total + BC_total,
                             data = train_data)

test_matrix <- model.matrix(tree_percentage ~ built_area_percentage + crops_percentage + rangeland_percentage + water_percentage +
                              flooded_vegetation_percentage + bare_ground_percentage + snow_ice_percentage + burned_area +
                              CO2_total + PM25_total + TPC_total + NMHC_total + OC_total + CH4_total + SO2_total + BC_total,
                            data = test_data)

# Train the model on the training data (2018-2022)
xgb_model <- xgboost(data = train_matrix, label = train_data$tree_percentage, nrounds = 100, objective = "reg:squarederror")

# Predict on the test data (2023)
test_data$predicted_xgboost <- predict(xgb_model, test_matrix)

# Calculate RMSE for the test set
rmse <- sqrt(mean((test_data$tree_percentage - test_data$predicted_xgboost)^2))
print(paste("RMSE on test data (2023):", rmse))

ggplot(test_data, aes(x = predicted_xgboost, y = tree_percentage)) +
  geom_point(color = "gray", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Actual Tree Percentage (Test Data, 2023)",
       x = "Predicted Tree Percentage",
       y = "Actual Tree Percentage") +
  theme_minimal()


######################################################## CORRPLOT ##############################################################################################

library(corrplot)

corr_data <- time_series_data %>%
  select(tree_percentage, built_area_percentage, crops_percentage, rangeland_percentage, water_percentage,
         flooded_vegetation_percentage, bare_ground_percentage, snow_ice_percentage, burned_area, CO2_total,
         PM25_total, TPC_total, NMHC_total, OC_total, CH4_total, SO2_total, BC_total)

corr_matrix <- cor(corr_data, use = "complete.obs")

corrplot(corr_matrix, method = "circle", type = "lower", tl.col = "black")



####################################################### MODEL COMPARISON #######################################################################################

rsquared <- function(actual, predicted) {
  cor(actual, predicted) ^ 2
}

# R² for linear model on test data
rsquared_lm <- rsquared(test_data$tree_percentage, test_data$predicted_lm)

# R² for nnet model on test data
rsquared_nnet <- rsquared(test_data$tree_percentage, test_data$predicted_nnet[,1])

# R² for Random Forest model on test data
rsquared_rf <- rsquared(test_data$tree_percentage, test_data$predicted_rf)

# R² for XGBoost model on test data
rsquared_xgb <- rsquared(test_data$tree_percentage, test_data$predicted_xgboost)

# Print R² values for comparison
cat("R² - Linear Model:", rsquared_lm, "\n")
cat("R² - Random Forest:", rsquared_rf, "\n")
cat("R² - NNet Model:", rsquared_nnet, "\n")
cat("R² - XGBoost:", rsquared_xgb, "\n")


###################################################### RMSE  ######################################################################################

rmse_lm <- sqrt(mean((test_data$tree_percentage - test_data$predicted_lm)^2))
rmse_rf <- sqrt(mean((test_data$tree_percentage - test_data$predicted_rf)^2))
rmse_nnet <- sqrt(mean((test_data$tree_percentage - test_data$predicted_nnet)^2))
rmse <- sqrt(mean((test_data$tree_percentage - test_data$predicted_xgboost)^2))

cat("RMSE for Linear Model on test data (2023):", rmse_lm, "\n")
cat("RMSE for Random Forest on test data (2023):", rmse_rf, "\n")
cat("RMSE for Neural Network on test data (2023):", rmse_nnet, "\n")
cat("RMSE on test data (2023):", rmse, "\n")


###################################################### PREDICT TREE CHANGE ############################################################################

land_use_vars <- c("tree_percentage", "built_area_percentage", "crops_percentage", 
                   "rangeland_percentage", "water_percentage", "flooded_vegetation_percentage",
                   "bare_ground_percentage", "snow_ice_percentage")

calculate_yearly_change <- function(data, variable, start_year, end_year) {
  for (year in start_year:end_year) {
    col_current <- paste0(variable, "_", year)
    col_base <- paste0(variable, "_", year - 1)
    col_change <- paste0(variable, "_change_", year - 1, "_", year)
    
    data[[col_change]] <- data[[col_current]] - data[[col_base]]
  }
  return(data)
}

start_year <- 2019  
end_year <- 2023 

for (variable in land_use_vars) {
  data <- calculate_yearly_change(data, variable, start_year, end_year)
}


data_long <- data %>%
  select(shapeName, contains("_change_")) %>%  
  pivot_longer(cols = contains("_change_"), names_to = "Variable", values_to = "Change") %>%
  # Extract the second year (2019 from tree_percentage_change_2018_2019)
  mutate(Year = as.numeric(sub(".*_change_\\d{4}_(\\d{4})$", "\\1", Variable)))


data_long <- data_long %>%
  mutate(ChangeType = sub("_change_\\d{4}_\\d{4}$", "", Variable)) %>%  
  select(-Variable) %>%  
  mutate(ChangeType = as.factor(ChangeType)) 

data_long <- data_long %>%
  pivot_wider(names_from = ChangeType, values_from = Change)


################## XG BOOST

train_data <- data_long %>% filter(Year < 2023)
test_data <- data_long %>% filter(Year == 2023)

train_matrix <- model.matrix(tree_percentage ~ built_area_percentage + crops_percentage + rangeland_percentage + water_percentage +
                               flooded_vegetation_percentage + bare_ground_percentage + snow_ice_percentage, 
                             data = train_data)

test_matrix <- model.matrix(tree_percentage ~ built_area_percentage + crops_percentage + rangeland_percentage + water_percentage +
                              flooded_vegetation_percentage + bare_ground_percentage + snow_ice_percentage, 
                            data = test_data)

xgb_model <- xgboost(data = train_matrix, label = train_data$tree_percentage, nrounds = 100, objective = "reg:squarederror", verbose = 0)

test_data$predicted_xgboost <- predict(xgb_model, test_matrix)

ggplot(test_data, aes(x = predicted_xgboost, y = tree_percentage)) +
  geom_point(color = "gray", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Actual Tree Percentage Change (2023)",
       x = "Predicted Tree Percentage Change",
       y = "Actual Tree Percentage Change") +
  theme_minimal()


test_data$residuals <- test_data$tree_percentage - test_data$predicted_xgboost

summary(test_data$residuals)

threshold <- 1 

outliers <- test_data %>% filter(abs(residuals) > threshold)

print(outliers)

library(ggrepel)

ggplot(test_data, aes(x = predicted_xgboost, y = tree_percentage)) +
  geom_point(color = "gray", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_point(data = outliers, aes(x = predicted_xgboost, y = tree_percentage), color = "red", size = 2) +  
  geom_text_repel(data = outliers, aes(label = shapeName), size = 3, color = "blue", max.overlaps = 10) +  
  coord_cartesian(xlim = c(min(test_data$predicted_xgboost) - 0.5, max(test_data$predicted_xgboost) + 0.5), 
                  ylim = c(min(test_data$tree_percentage) - 0.5, max(test_data$tree_percentage) + 0.5)) +  
  labs(title = "Predicted vs Actual Tree Percentage Change (2023) with Outliers Highlighted",
       x = "Predicted Tree Percentage Change",
       y = "Actual Tree Percentage Change") +
  theme_minimal()


########################### PREDICT 2024

test_data_2024 <- test_data  

test_matrix_2024 <- model.matrix(tree_percentage ~ built_area_percentage + crops_percentage + rangeland_percentage + water_percentage +
                                   flooded_vegetation_percentage + bare_ground_percentage + snow_ice_percentage, 
                                 data = test_data_2024)

test_data_2024$predicted_xgboost_2024 <- predict(xgb_model, test_matrix_2024)

print(test_data_2024$predicted_xgboost_2024)

ggplot(test_data_2024, aes(x = shapeName, y = predicted_xgboost_2024)) +
  geom_point(color = "blue", size = 3) +
  labs(title = "Predicted Tree Percentage Change for 2024",
       x = "State/Region",
       y = "Predicted Tree Percentage Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


l <- train_data %>% select (shapeName, tree_percentage, Year)


m <- test_data_2024 %>% select(shapeName, tree_percentage, Year)

combined_data <- rbind(l, m)

test_data_2024$Year <- 2024
test_data_2024$tree_percentage <- test_data_2024$predicted_xgboost_2024

k <- test_data_2024 %>% select(shapeName, tree_percentage, Year)

ccc <- rbind(combined_data, k)


ggplot(ccc, aes(x = Year, y = tree_percentage, group = shapeName, color = shapeName)) +
  geom_line() + 
  geom_point(size = 2) + 
  labs(title = "Tree Percentage Change: Actual (2019-2023) and Predicted (2024)",
       x = "Year", y = "Tree Percentage Change") +
  theme_minimal() +
  theme(legend.position = "none")


###################################### TRY LINEAR

train_data <- data_long %>% filter(Year < 2023)
test_data <- data_long %>% filter(Year == 2023)

lm_model <- lm(tree_percentage ~ built_area_percentage + crops_percentage + rangeland_percentage + water_percentage +
                 flooded_vegetation_percentage + bare_ground_percentage + snow_ice_percentage, 
               data = train_data)
test_data$predicted_lm <- predict(lm_model, newdata = test_data)

ggplot(test_data, aes(x = predicted_lm, y = tree_percentage)) +
  geom_point(color = "black", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Actual Tree Percentage Change (2023) using Linear Model",
       x = "Predicted Tree Percentage Change",
       y = "Actual Tree Percentage Change") +
  theme_minimal()


################################################# TRY RANDOM FOREST -> worse

train_data <- data_long %>% filter(Year < 2023)
test_data <- data_long %>% filter(Year == 2023)

rf_model <- randomForest(tree_percentage ~ built_area_percentage + crops_percentage + rangeland_percentage + water_percentage +
                           flooded_vegetation_percentage + bare_ground_percentage + snow_ice_percentage, 
                         data = train_data, 
                         ntree = 500, 
                         importance = TRUE)

test_data$predicted_rf <- predict(rf_model, test_data)

ggplot(test_data, aes(x = predicted_rf, y = tree_percentage)) +
  geom_point(color = "black", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Actual Tree Percentage Change (2023) using Random Forest",
       x = "Predicted Tree Percentage Change",
       y = "Actual Tree Percentage Change") +
  theme_minimal()

####################################################### SIMULATIONS WITH XGBOOST ####################################################################

simulate_tree_percentage <- function(test_data, model, built_area_change, crops_change, pollution_reduction) {
  simulated_data <- test_data
  
  simulated_data$built_area_percentage <- simulated_data$built_area_percentage * (1 - built_area_change)
  simulated_data$crops_percentage <- simulated_data$crops_percentage * (1 - crops_change)
  simulated_data$CO2_total <- simulated_data$CO2_total * (1 - pollution_reduction)
  
  simulated_data_matrix <- model.matrix(tree_percentage ~ built_area_percentage + crops_percentage + rangeland_percentage + water_percentage +
                                          flooded_vegetation_percentage + bare_ground_percentage + snow_ice_percentage + burned_area +
                                          CO2_total + PM25_total + TPC_total + NMHC_total + OC_total + CH4_total + SO2_total + BC_total,
                                        data = simulated_data)
  
  simulated_data$predicted_tree_percentage <- predict(model, simulated_data_matrix)
  
  return(simulated_data)
}

# Reduce built area by 10%, crops by 10%, and pollution by 20%
simulated_scenario <- simulate_tree_percentage(test_data, xgb_model, built_area_change = 0.10, crops_change = 0.10, pollution_reduction = 0.20)

comparison <- test_data %>%
  select(shapeName, tree_percentage, predicted_xgboost) %>%
  mutate(simulated_tree_percentage = simulated_scenario$predicted_tree_percentage)

# Calculate the difference between simulated and predicted values
comparison <- comparison %>%
  mutate(difference = simulated_tree_percentage - predicted_xgboost)

ggplot(comparison, aes(x = predicted_xgboost, y = simulated_tree_percentage, color = difference)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_gradient2(low = "red", mid = "pink", high = "blue", midpoint = 0) +
  labs(title = "Simulated vs Predicted Tree Percentage (Test Data 2023)",
       x = "Current Predicted Tree Percentage",
       y = "Simulated Optimized Tree Percentage",
       color = "Change") +
  theme_minimal()


help <- regional_comparison %>% filter(difference < 0)

help <- help %>% inner_join(correlation_results, "shapeName") 
help <- help %>% select (-cor_rangeland)

# Summarize the total improvements and declines
summary_by_region <- comparison %>%
  summarise(
    total_improved = sum(difference > 0),
    total_declined = sum(difference < 0),
    avg_improvement = mean(difference[difference > 0], na.rm = TRUE),
    avg_decline = mean(difference[difference < 0], na.rm = TRUE)
  )

print(summary_by_region)

# Regional improvement
ggplot(regional_comparison, aes(x = avg_predicted_tree_percentage, y = avg_simulated_tree_percentage, color = difference)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_gradient2(low = "red", mid = "yellow", high = "blue", midpoint = 0) +
  labs(title = "Simulated vs Predicted Tree Percentage (by Region)",
       x = "Average Predicted Tree Percentage",
       y = "Average Simulated Tree Percentage",
       color = "Change") +
  theme_minimal()


help_table <- time_series_data %>%
  filter(shapeName == "Somogy") %>%
  select (shapeName, country, year, tree_percentage, built_area_percentage, crops_percentage, rangeland_percentage)
help_table



###################################################################################### SPATIO-TEMPORAL CLUSTERING #########################################################################################################


polygon_distances <- st_distance(data_with_state_boarders$geometry)
spatial_distances <- as.matrix(polygon_distances)

temporal_data <- data %>% 
  select(starts_with("tree_percentage")) %>%  
  as.matrix()

# Compute the Euclidean distance between time series (temporal data)
temporal_distances <- as.matrix(dist(temporal_data)) 

normalize <- function(mat) {
  mat / sqrt(sum(mat^2))  # Frobenius norm
}

spatial_distances_norm <- normalize(spatial_distances)
temporal_distances_norm <- normalize(temporal_distances)

spatial_distances <- as.numeric(spatial_distances)
spatial_distances <- matrix(spatial_distances, nrow = nrow(spatial_distances_norm), ncol = ncol(spatial_distances_norm))

temporal_distances <- as.numeric(temporal_distances)
temporal_distances <- matrix(temporal_distances, nrow = nrow(temporal_distances_norm), ncol = ncol(temporal_distances_norm))


spatial_distances_norm <- normalize(spatial_distances)
temporal_distances_norm <- normalize(temporal_distances)

alpha <- 0.5  # 0 only spatial clustering 1 only feature clustering
combined_distances <- alpha * temporal_distances_norm + (1 - alpha) * spatial_distances_norm


gap_stat <- clusGap(combined_distances, FUN = pam, K.max = 10, B = 50)  # B = number of bootstraps
fviz_gap_stat(gap_stat)
optimal_clusters <- maxSE(f = gap_stat$Tab[, "gap"], SE.f = gap_stat$Tab[, "SE.sim"])

set.seed(123)
pam_result <- pam(combined_distances, k = optimal_clusters)  

summary(pam_result)

data_with_state_boarders$cluster <- pam_result$clustering

ggplot(data_with_state_boarders) +
  geom_sf(aes(fill = as.factor(cluster)), color = "black", size = 0.5) +
  labs(
    title = "Spatio-Temporal Clustering Map, Tree Percentage Alpha = 0.5",
    fill = "Clusters"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


# FOR LOOP FOR ALPHA FROM 0.1 TO 0.9
plots_list <- list()
for (alpha in seq(0.1, 0.9, by = 0.1)) {
  
  combined_distances <- alpha * temporal_distances_norm + (1 - alpha) * spatial_distances_norm
  
  gap_stat <- clusGap(combined_distances, FUN = pam, K.max = 10, B = 50)  # B = number of bootstraps
  
  optimal_clusters <- maxSE(f = gap_stat$Tab[, "gap"], SE.f = gap_stat$Tab[, "SE.sim"])
  print(paste("alpha: ", alpha, " opt b: ", optimal_clusters))
  
  set.seed(123)
  pam_result <- pam(combined_distances, k = optimal_clusters)
  
  data_with_state_boarders$cluster <- pam_result$clustering
  
  p <- ggplot(data_with_state_boarders) +
    geom_sf(aes(fill = as.factor(cluster)), color = "black", size = 0.5) +
    labs(
      title = paste("alpha =", alpha, "nclust =", optimal_clusters),
      fill = "Clusters"
    ) +
    theme_minimal()
  
  plots_list[[paste("alpha_", alpha, sep = "")]] <- p
}


plots <- wrap_plots(plots_list, ncol = 3, nrows = 3) 
plots

data$cluster <- pam_result$clustering

time_series_data <- time_series_data %>% select(-c(cluster))

time_series_data <- data %>%
  select(shapeName, cluster) %>%
  left_join(time_series_data,  by = "shapeName")

ggplot(time_series_data, aes(x = year, y = tree_percentage, group = shapeName, color = as.factor(cluster))) +
  geom_line(alpha = 0.5) +  # Use lines to show the trends over time for each region
  labs(title = "Tree Percentage Trends by Region and Cluster for alpha = 0.5", 
       x = "Year", y = "Tree Percentage", color = "Cluster") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for clarity
    plot.title = element_text(hjust = 0.5)
  )



###################################################### ALL VARIABLES FOR CLUSTERING ############################################################################

polygon_distances <- st_distance(data_with_state_boarders$geometry)
spatial_distances <- as.matrix(polygon_distances)

temporal_data <- data %>% 
  select(starts_with(c("tree_percentage", "water_percentage", "crops_percentage", "built_area_percentage", "rangeland_percentage", "bare_ground", "snow_ice", "flooded"))) %>%  
  mutate(across(everything(), as.numeric)) %>%  
  as.matrix()

# Compute the Euclidean distance between time series (temporal data)
temporal_distances <- as.matrix(dist(temporal_data)) 

normalize <- function(mat) {
  mat / sqrt(sum(mat^2))  # Frobenius norm
}

spatial_distances_norm <- normalize(spatial_distances)
temporal_distances_norm <- normalize(temporal_distances)

spatial_distances <- as.numeric(spatial_distances)
spatial_distances <- matrix(spatial_distances, nrow = nrow(spatial_distances_norm), ncol = ncol(spatial_distances_norm))

temporal_distances <- as.numeric(temporal_distances)
temporal_distances <- matrix(temporal_distances, nrow = nrow(temporal_distances_norm), ncol = ncol(temporal_distances_norm))


spatial_distances_norm <- normalize(spatial_distances)
temporal_distances_norm <- normalize(temporal_distances)

alpha <- 0.5  # 0 only spatial clustering 1 only feature clustering
combined_distances <- alpha * temporal_distances_norm + (1 - alpha) * spatial_distances_norm


gap_stat <- clusGap(combined_distances, FUN = pam, K.max = 10, B = 50)  # B = number of bootstraps
fviz_gap_stat(gap_stat)
optimal_clusters <- maxSE(f = gap_stat$Tab[, "gap"], SE.f = gap_stat$Tab[, "SE.sim"])

set.seed(123)
pam_result <- pam(combined_distances, k = optimal_clusters)  


data_with_state_boarders$cluster <- pam_result$clustering

ggplot(data_with_state_boarders) +
  geom_sf(aes(fill = as.factor(cluster)), color = "black", size = 0.5) +
  labs(
    title = "Spatio-Temporal Clustering Map for Alpha = 0.5 And Important Land Use Features",
    fill = "Clusters"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


######################################################### CLUSTERING VISUALIZATION MULTIPLE VARS ##############################################################

library(GGally)

data$cluster <- pam_result$clustering
#time_series_data <- time_series_data %>% select(-c(cluster))

time_series_data <- data %>%
  select(shapeName, cluster) %>%
  left_join(time_series_data,  by = "shapeName")


time_series_data <- time_series_data %>%
  mutate(across(c(tree_percentage, water_percentage, crops_percentage, built_area_percentage, rangeland_percentage), 
                as.numeric))  # Ensure numeric type

time_series_data <- time_series_data %>%
  filter(!is.na(tree_percentage) & !is.na(water_percentage) & !is.na(crops_percentage) & 
           !is.na(built_area_percentage) & !is.na(rangeland_percentage))  # Remove rows with NAs

time_series_data <- time_series_data %>%
  mutate(cluster = as.factor(cluster))

generate_parallel_plot <- function(year_data, year) {
  ggparcoord(
    data = year_data, 
    match(c('tree_percentage', 'water_percentage', 'crops_percentage', 'built_area_percentage', 'rangeland_percentage'),
          names(year_data)
    ),
    groupColumn = "cluster",  # Group by cluster
    scale = "globalminmax",  
    alphaLines = 0.6  
  ) +
    labs(title = paste("Parallel Coordinate Plot for Clusters in Year", year),
         color = "Cluster") +
    theme_minimal() +
    theme(legend.position = "bottom")  
}

years <- 2018:2023
plots_list <- list()

for (year in years) {
  year_data <- time_series_data %>% filter(year == !!year)  
  
  if(nrow(year_data) > 0) {  
    plot <- generate_parallel_plot(year_data, year)  
    plots_list[[paste0("plot_", year)]] <- plot  
  }
}


plots_list[["plot_2018"]]
plots_list[["plot_2019"]]
plots_list[["plot_2020"]]
plots_list[["plot_2021"]]
plots_list[["plot_2022"]]
plots_list[["plot_2023"]]


# Feature importance for clusters

data$cluster <- pam_result$clustering

time_series_data <- time_series_data %>% select(-c(cluster))

time_series_data <- data %>%
  select(shapeName, cluster) %>%
  left_join(time_series_data,  by = "shapeName")

library(xgboost)

importance_list <- list()
for (clust in unique(data_with_state_boarders$cluster)) {
  
  subset_data <- time_series_data %>% filter(cluster == clust)
  
  model <- xgboost(data = model.matrix(~ crops_percentage + built_area_percentage + water_percentage + rangeland_percentage + flooded_vegetation_percentage + 
                                         bare_ground_percentage +  snow_ice_percentage - 1, subset_data),
                   label = subset_data$tree_percentage, nrounds = 50, objective = "reg:squarederror", verbose = 0)
  
  importance <- xgb.importance(feature_names = colnames(model.matrix(~ crops_percentage + built_area_percentage + water_percentage  + rangeland_percentage 
                                                                     + flooded_vegetation_percentage + 
                                                                       bare_ground_percentage +  snow_ice_percentage - 1, subset_data)),
                               model = model)
  
  importance$Cluster <- as.factor(clust)
  
  importance_list[[paste("Cluster", clust)]] <- importance
}

importance_df <- do.call(rbind, importance_list)

ggplot(importance_df, aes(x = reorder(Feature, Gain), y = Gain, fill = Cluster)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ Cluster, scales = "free") +
  labs(title = "Feature Importance for Each Cluster", x = "Features", y = "Importance (Gain)") +
  theme_minimal() +
  theme(legend.position = "none")




specific_cluster_data <- time_series_data %>%
  filter(cluster == 2)

# Plot the tree percentage trends for just this cluster
ggplot(specific_cluster_data, aes(x = year, y = tree_percentage, group = shapeName, color = shapeName)) +
  geom_line() +  # Plot the individual region trends within the cluster
  labs(title = "Tree Percentage Trends in Cluster 2", 
       x = "Year", y = "Tree Percentage", color = "Region") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

help <- data_with_state_boarders %>%
  select(shapeName, shapeGroup, cluster, starts_with("tree_percentage"))


########################################################## STDBSAN ATTEMPT - no sence in results ################################################################################

stdbscan = function(spatial_distances_norm, 
                    temporal_distances_norm, 
                    eps, 
                    eps2, 
                    minpts, 
                    cldensity = TRUE) { 
  
  countmode = 1:nrow(spatial_distances_norm)  # Number of points
  seeds = TRUE
  
  # Use the precomputed spatial and temporal distances
  data_spatial <- spatial_distances_norm
  data_temporal <- temporal_distances_norm
  n <- nrow(data_spatial)
  
  # Initialize clustering labels and other structures
  classn <- cv <- integer(n)  # Cluster labels for each point
  isseed <- logical(n)        # To mark seed points
  cn <- integer(1)            # Cluster number
  
  for (i in 1:n) {
    if (i %in% countmode)
      #cat("Processing point ", i, " of ", n, ".\n")
      unclass <- (1:n)[cv < 1]  # Points that have not been classified yet
    
    if (cv[i] == 0) {  # If point i has not been visited yet
      # Find neighbors based on both spatial and temporal distances
      reachables <- intersect(unclass[data_spatial[i, unclass] <= eps],  
                              unclass[data_temporal[i, unclass] <= eps2])
      
      # If there are not enough points to form a cluster, mark as noise
      if (length(reachables) + classn[i] < minpts) {
        cv[i] <- (-1)  # Mark as noise (cv = -1)
      } else {
        # Create a new cluster
        cn <- cn + 1                   
        cv[i] <- cn
        isseed[i] <- TRUE
        reachables <- setdiff(reachables, i)  # Remove i from reachables
        unclass <- setdiff(unclass, i)        # Remove i from unclassified
        
        # Update class count for the neighbors
        classn[reachables] <- classn[reachables] + 1
        
        # Continue clustering with the reachable neighbors
        while (length(reachables)) {
          cv[reachables] <- cn  # Assign cluster label
          ap <- reachables      # Store current reachables
          reachables <- integer()
          
          for (i2 in seq(along = ap)) {
            j <- ap[i2]
            
            # Find neighbors of j based on spatial and temporal distances
            jreachables <- intersect(unclass[data_spatial[j, unclass] <= eps], 
                                     unclass[data_temporal[j, unclass] <= eps2])
            
            if (length(jreachables) + classn[j] >= minpts) {
              isseed[j] <- TRUE
              cv[jreachables[cv[jreachables] < 0]] <- cn  # Assign cluster to previously noise points
              reachables <- union(reachables, jreachables[cv[jreachables] == 0])
            }
            classn[jreachables] <- classn[jreachables] + 1  # Update class count
            unclass <- setdiff(unclass, j)  # Remove j from unclassified
          }
        }
      }
    }
    if (!length(unclass))
      break
  }
  
  # Assign noise points to cluster 0
  if (any(cv == (-1))) {
    cv[cv == (-1)] <- 0
  }
  
  # Return the clustering result
  out <- list(cluster = cv, eps = eps, minpts = minpts, density = classn)
  
  # Add seeds (core points) to the result if applicable
  if (seeds && cn > 0) {
    out$isseed <- isseed
  }
  
  class(out) <- "stdbscan"
  return(out)
}


hist(as.numeric(spatial_distances_norm), breaks = 50, main = "Histogram of Spatial Distances", xlab = "Spatial Distance")

# Visualize temporal distances
hist(as.numeric(temporal_distances_norm), breaks = 50, main = "Histogram of Temporal Distances", xlab = "Temporal Distance")

result <- stdbscan(
  spatial_distances_norm = spatial_distances_norm,      # Normalized spatial distances
  temporal_distances_norm = temporal_distances_norm,    # Normalized temporal distances
  eps = 0.001,                                           # Spatial distance threshold
  eps2 = 0.002,                                          # Temporal distance threshold
  minpts = 5,                                          # Minimum points for forming a dense cluster
  cldensity = TRUE                                      # Optional: to track density of clusters
)


data_with_state_boarders$st_dbscan_cluster <- result$cluster

ggplot(data_with_state_boarders) +
  geom_sf(aes(fill = as.factor(st_dbscan_cluster)), color = "black", size = 0.5) +
  labs(
    title = "Spatio-Temporal Clustering Map STDBSCAN",
    fill = "ST-DBSCAN Clusters"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )




####################################### FOLLOWING CODE IS OT TO BE USED (these are just some approaches I tried, but didn't use) #####################################################

centroids <- st_centroid(data_with_state_boarders$geometry)
spatial_coords <- st_coordinates(centroids)

data_with_state_boarders$longitude <- spatial_coords[, 1]
data_with_state_boarders$latitude <- spatial_coords[, 2]

gg <- ggplot(data_with_state_boarders, aes(x = longitude, y = latitude, color = factor(cluster))) +
  geom_point(size = 3) +
  labs(title = "Spatio-Temporal Clustering", x = "Longitude", y = "Latitude", color = "Cluster") +
  theme_minimal()

gg

# Calculate the mean tree percentage for each state across the years 2018–2023
data_with_state_boarders <- data_with_state_boarders %>%
  rowwise() %>%  # Apply the operation row-wise
  mutate(mean_tree_percentage = mean(c_across(starts_with("tree_percentage_")), na.rm = TRUE))


tree_percentage_long <- data_with_state_boarders %>%
  select(shapeName, cluster, starts_with("tree_percentage")) %>%  # Select relevant columns
  pivot_longer(cols = starts_with("tree_percentage"),
               names_to = "year",
               names_prefix = "tree_percentage_",  # Remove prefix for clarity
               values_to = "tree_percentage") %>% 
  mutate(year = as.numeric(year))  # Convert year column to numeric

# Plot tree percentage trends by cluster
ggplot(tree_percentage_long, aes(x = year, y = tree_percentage, color = factor(cluster), group = shapeName)) +
  geom_line(size = 1) +
  facet_wrap(~ cluster) +  # Separate plots for each cluster
  labs(title = "Tree Percentage Trends by Cluster (2018-2023)",
       x = "Year",
       y = "Tree Percentage",
       color = "Cluster") +
  theme_minimal()


cluster_summary <- data_with_state_boarders %>%
  group_by(cluster) %>%
  summarize(mean_tree_percentage = mean(mean_tree_percentage, na.rm = TRUE))

print(cluster_summary)

ggplot(tree_percentage_long, aes(x = year, y = tree_percentage, color = factor(cluster), group = shapeName)) +
  geom_line(size = 0.5, alpha = 0.5) +  # Individual state trends
  stat_summary(aes(group = factor(cluster)), fun = "mean", geom = "line", size = 1.5, color = "black") +  # Mean trend
  facet_wrap(~ cluster) +
  labs(title = "Tree Percentage Trends by Cluster (2018-2023)",
       x = "Year",
       y = "Tree Percentage",
       color = "Cluster") +
  theme_minimal()



stdbscan_with_features <- function(x,
                                   y,
                                   time,
                                   features,
                                   eps,
                                   eps2,
                                   eps_features, 
                                   minpts) {
  
  countmode <- 1:length(x)
  seeds <- TRUE
  
  data_spatial <- as.matrix(dist(cbind(y, x)))  
  print(data_spatial)
  data_temporal <- as.matrix(dist(time))        
  data_features <- as.matrix(dist(features))   
  

  print("LOOK")
  #print(data_features)
  
  n <- nrow(data_spatial)
  classn <- cv <- integer(n)
  isseed <- logical(n)
  cn <- integer(1)
  
  for (i in 1:n) {
    unclass <- (1:n)[cv < 1]
    
    if (cv[i] == 0) {
      # Combine weighted spatial, temporal, and feature distances to find reachable points
      reachables_spatial <- unclass[data_spatial[i, unclass] <= eps]
      reachables_temporal <- unclass[data_temporal[i, unclass] <= eps2]
      reachables_features <- unclass[data_features[i, unclass] <= eps_features]
      
      # Print reachable points based on feature distance
      print(paste("Observation", i, "feature-reachable points:", paste(reachables_features, collapse = ", ")))
      
      # Combine all reachables
      reachables <- intersect(intersect(reachables_spatial, reachables_temporal), reachables_features)
      
      if (length(reachables) + classn[i] < minpts){
        cv[i] <- (-1)  # Mark as noise
      }
      else {
        cn <- cn + 1
        cv[i] <- cn
        isseed[i] <- TRUE
        reachables <- setdiff(reachables, i)
        unclass <- setdiff(unclass, i)
        classn[reachables] <- classn[reachables] + 1
        
        while (length(reachables)) {
          cv[reachables] <- cn  # Assign cluster number
          ap <- reachables
          reachables <- integer()
          
          for (i2 in seq(along = ap)) {
            j <- ap[i2]
            
            jreachables <- intersect(
              intersect(unclass[data_spatial[j, unclass] <= eps], 
                        unclass[data_temporal[j, unclass] <= eps2]),
              unclass[data_features[j, unclass] <= eps_features]
            )
            
            if (length(jreachables) + classn[j] >= minpts) {
              isseed[j] <- TRUE
              cv[jreachables[cv[jreachables] < 0]] <- cn
              reachables <- union(reachables, jreachables[cv[jreachables] == 0])
            }
            classn[jreachables] <- classn[jreachables] + 1
            unclass <- setdiff(unclass, j)
          }
        }
      }
    }
    if (!length(unclass)) break
  }
  
  if (any(cv == (-1))) {
    cv[cv == (-1)] <- 0
  }
  
  out <- list(cluster = cv, eps = eps, eps2 = eps2, eps_features = eps_features, minpts = minpts, density = classn)
  if (seeds && cn > 0) {
    out$isseed <- isseed
  }
  class(out) <- "stdbscan"
  return(out)
}



  time_series_data_with_geometry <- data_with_state_boarders %>% 
    select(shapeName, geometry) %>%
    left_join(time_series_data, by = "shapeName")
  
  time_series_data_with_geometry <- time_series_data_with_geometry %>%
    arrange(shapeGroup, shapeName, year)
  
  centroids <- st_centroid(time_series_data_with_geometry$geometry)
  spatial_coords <- st_coordinates(centroids)
  
  selected_variables <- c("tree_percentage")
  print("HELLO")
  print(selected_variables)

  temporal_data <- time_series_data %>%
    select(all_of(selected_variables))
  
  print("Selected vars:")
  print(head(temporal_data))
  
  scaled_temporal_data <- scale(temporal_data)
  #print(head(scaled_temporal_data))
  
  scaled_spatial_coords <- scale(spatial_coords)
  time_data <- time_series_data$year
  scaled_time_data <- scale(time_data)
  

  result <- stdbscan_with_features(
    x = scaled_spatial_coords[, 1],           # long
    y = scaled_spatial_coords[, 2],           # lat
    time = scaled_time_data,                  # time (year)
    features = scaled_temporal_data,   # only selected variables (features) for feature distance
    eps = 2, #input$eps,                   # spatial distance threshold
    eps2 = 2, #input$eps2,                 # temp distance threshold
    eps_features = 0.03,      # feat distance threshold (adjust as needed)
    minpts = 20              # min points for a cluster
  )
  
  time_series_data_with_geometry$st_dbscan_cluster <- result$cluster
  time_series_data_with_geometry <- time_series_data_with_geometry %>%
    mutate(st_dbscan_cluster = as.factor(st_dbscan_cluster))
  
  time_series_data_with_geometry$longitude <- spatial_coords[, 1]
  time_series_data_with_geometry$latitude <- spatial_coords[, 2]
  
  # Before scaling
  hist(time_series_data_with_geometry$tree_percentage, breaks = 20, main = "Tree Percentage Distribution (Original)")
  
  # After scaling
  scaled_tree_percentage <- scale(time_series_data_with_geometry$tree_percentage)
  hist(scaled_tree_percentage, breaks = 20, main = "Tree Percentage Distribution (Scaled)")
  
  cluster_color_mapping <- get_color_palette_clusters(time_series_data_with_geometry$st_dbscan_cluster)
  
  gg <- ggplot(time_series_data_with_geometry, aes(x = longitude, y = latitude, color = st_dbscan_cluster)) +
    geom_point(aes(
      text = paste(
        "Cluster:", st_dbscan_cluster,
        "<br>State:", shapeName,
        "<br>Tree %:", round(tree_percentage, 2),
        "<br>Longitude:", round(longitude, 2),
        "<br>Latitude:", round(latitude, 2)
      )
    ), size = 3) +
    scale_color_manual(values = cluster_color_mapping) + 
    labs(title = "ST-DBSCAN Clustering With Features", x = "Longitude", y = "Latitude", color = "Cluster") +
    theme_minimal()
  
  ggplotly(gg, tooltip = "text")
  
  ggplot(time_series_data_with_geometry, aes(x = st_dbscan_cluster, y = tree_percentage, fill = st_dbscan_cluster)) +
    geom_boxplot() +
    labs(title = "Distribution of Tree Percentage Across Clusters",
         x = "Cluster", y = "Tree Percentage") +
    theme_minimal()
  
  data_spatial <- as.matrix(dist(cbind(time_series_data_with_geometry$longitude, time_series_data_with_geometry$latitude)))  
  data_temporal <- as.matrix(dist(time_series_data_with_geometry$year))        
  data_features <- as.matrix(dist(scaled_temporal_data))   
  
  plot_heatmap <- function(distance_matrix, title) {
    dist_melt <- melt(as.matrix(distance_matrix))
    ggplot(dist_melt, aes(Var1, Var2, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = "blue", high = "red") +
      labs(title = title, x = "Observations", y = "Observations", fill = "Distance") +
      theme_minimal()
  }
  
  # Visualize spatial distance matrix
  spatial_heatmap <- plot_heatmap(data_spatial, "Spatial Distance Matrix")
  spatial_heatmap
  
  # Visualize temporal distance matrix
  temporal_heatmap <- plot_heatmap(data_temporal, "Temporal Distance Matrix")
  temporal_heatmap
  
  # Visualize feature distance matrix
  feature_heatmap <- plot_heatmap(data_features, "Feature Distance Matrix")
  feature_heatmap
  
  
  
  pal <- colorFactor(palette = unname(cluster_color_mapping), domain = time_series_data_with_geometry$st_dbscan_cluster)
  
  leaflet(data = time_series_data_with_geometry) %>%
    addTiles() %>%
    addPolygons(
      fillColor = ~pal(st_dbscan_cluster),  
      fillOpacity = 0.7,
      color = "black",                    
      weight = 1,                          
      highlightOptions = highlightOptions(weight = 2, color = "#666", fillOpacity = 0.7, bringToFront = TRUE),
      popup = ~paste0(
        "<div style='font-size: 14px;'><b>State:</b> ", shapeName, "<br>",
        "<b>Cluster:</b> ", st_dbscan_cluster, "<br>",
        "<b>Tree percentage:</b> ", round(tree_percentage, 2), "<br>"
      ), 
      layerId = ~shapeName
    ) %>%
    addLegend(
      pal = pal,
      values = ~st_dbscan_cluster,    
      title = "ST-DBSCAN Clusters",
      position = "bottomright"
    ) %>%
    setView(lng = 15, lat = 48, zoom = 4) %>%  
    addControl("<strong>Spatio-Temporal Clustering Map</strong>", 
               position = "bottomleft", 
               className = "map-title")
  
  

  
  time_series_data_with_geometry$longitude <- spatial_coords[, 1]
  time_series_data_with_geometry$latitude <- spatial_coords[, 2]
  
  

  gg <- ggplot(time_series_data_with_geometry, aes(x = longitude, y = latitude, color = st_dbscan_cluster)) +
    geom_point(aes(
      text = paste(
        "Cluster:", st_dbscan_cluster,
        "<br>State:", shapeName,
        "<br>Tree %:", round(tree_percentage, 2),
        "<br>Longitude:", round(longitude, 2),
        "<br>Latitude:", round(latitude, 2)
      )
    ), size = 3) +
    labs(title = "ST-DBSCAN Clustering With Features", x = "Longitude", y = "Latitude", color = "Cluster") +
    theme_minimal()
  
  ggplotly(gg, tooltip = "text")
  
  pal <- colorFactor(palette = unname(cluster_color_mapping), domain = time_series_data_with_geometry$st_dbscan_cluster)
  
  leaflet(data = time_series_data_with_geometry) %>%
    addTiles() %>%
    addPolygons(
      fillColor = ~pal(st_dbscan_cluster),  
      fillOpacity = 0.7,
      color = "black",                    
      weight = 1,                          
      highlightOptions = highlightOptions(weight = 2, color = "#666", fillOpacity = 0.7, bringToFront = TRUE),
      popup = ~paste0(
        "<div style='font-size: 14px;'><b>State:</b> ", shapeName, "<br>",
        "<b>Cluster:</b> ", st_dbscan_cluster, "<br>",
        "<b>Tree percentage:</b> ", round(tree_percentage, 2), "<br>"
      ), 
      layerId = ~shapeName
    ) %>%
    addLegend(
      pal = pal,
      values = ~st_dbscan_cluster,    
      title = "ST-DBSCAN Clusters",
      position = "bottomright"
    ) %>%
    setView(lng = 15, lat = 48, zoom = 4) %>%  
    addControl("<strong>Spatio-Temporal Clustering Map</strong>", 
               position = "bottomleft", 
               className = "map-title")
  
  
  stdbscan <- function (spatial_distances, 
                        time,
                        eps, 
                        eps2, 
                        minpts, 
                        cldensity = TRUE) { 
    
    countmode <- 1:nrow(spatial_distances)  # Correct this to reflect the size of the data
    seeds <- TRUE
    
    data_spatial <- spatial_distances  # Precomputed spatial distance matrix
    data_temporal <- as.matrix(dist(time))  # Compute temporal distance matrix from 'time' data
    n <- nrow(data_spatial)
    
    classn <- cv <- integer(n)
    isseed <- logical(n)
    cn <- integer(1)
    
    for (i in 1:n) {
      if (i %in% countmode) {
        unclass <- (1:n)[cv < 1]  # Unclassified points
        
        if (cv[i] == 0) {
          reachables <- intersect(unclass[data_spatial[i, unclass] <= eps],  unclass[data_temporal[i, unclass] <= eps2])
          
          if (length(reachables) + classn[i] < minpts)
            cv[i] <- -1  # Mark as noise                    
          else {
            cn <- cn + 1                   
            cv[i] <- cn
            isseed[i] <- TRUE
            reachables <- setdiff(reachables, i)
            unclass <- setdiff(unclass, i)       
            classn[reachables] <- classn[reachables] + 1
            
            while (length(reachables)) {
              cv[reachables] <- cn           
              ap <- reachables                           
              reachables <- integer(0)
              
              for (i2 in seq_along(ap)) {
                j <- ap[i2]
                
                jreachables <- intersect(unclass[data_spatial[j, unclass] <= eps], unclass[data_temporal[j, unclass] <= eps2])
                
                if (length(jreachables) + classn[j] >= minpts) {
                  isseed[j] <- TRUE
                  cv[jreachables[cv[jreachables] < 0]] <- cn
                  reachables <- union(reachables, jreachables[cv[jreachables] == 0])
                }
                classn[jreachables] <- classn[jreachables] + 1
                unclass <- setdiff(unclass, j)
              }
            }
          }
        }
        if (!length(unclass)) break
      }
    }
    
    # If any points were marked as noise, keep them marked as noise (cv = 0)
    if (any(cv == -1)) {
      cv[cv == -1] <- 0
    }
    
    out <- list(cluster = cv, eps = eps, minpts = minpts, density = classn)
    if (seeds && cn > 0) {
      out$isseed <- isseed
    }
    class(out) <- "stdbscan"
    return(out)
  }
  
  
  print(dist(temporal_data[1:10, ]))
 
  #scaled_spatial_coords <- scale(spatial_coords)  # Assuming longitude and latitude are in this matrix
  scaled_feature_data <- scale(temporal_data)     # Assuming feature data is in this matrix
  
  polygon_distances <- st_distance(data_with_state_boarders$geometry)
  spatial_distances <- as.matrix(polygon_distances)
  
  spatial_distances <- as.numeric(spatial_distances)
  spatial_distances <- matrix(spatial_distances, nrow = nrow(spatial_distances_norm), ncol = ncol(spatial_distances_norm))
  
  result <- stdbscan(
    spatial_distances = spatial_distances,
    time = scaled_feature_data,     # Feature data (tree percentage, etc.)
    eps = 1.5,                             # Spatial distance threshold
    eps2 = 5,                             # Feature distance threshold
    minpts = 15
    )
  
  print(result$cluster)
  

  time_series_data_with_geometry$st_dbscan_cluster <- result$cluster
  time_series_data_with_geometry <- time_series_data_with_geometry %>%
    mutate(st_dbscan_cluster = as.factor(st_dbscan_cluster))
  
  time_series_data_with_geometry$longitude <- spatial_coords[, 1]
  time_series_data_with_geometry$latitude <- spatial_coords[, 2]
  
  #dta <- time_series_data_with_geometry %>%
   # filter(year == 2018)
  
  
  gg <- ggplot(time_series_data_with_geometry, aes(x = longitude, y = latitude, color = st_dbscan_cluster)) +
    geom_point(aes(
      text = paste(
        "Cluster:", st_dbscan_cluster,
        "<br>State:", shapeName,
        "<br>Tree %:", round(tree_percentage, 2),
        "<br>Longitude:", round(longitude, 2),
        "<br>Latitude:", round(latitude, 2)
      )
    ), size = 3) +
    labs(title = "ST-DBSCAN Clustering With Features", x = "Longitude", y = "Latitude", color = "Cluster") +
    theme_minimal()
  
  ggplotly(gg, tooltip = "text")
  
 
  
  leaflet(data = time_series_data_with_geometry) %>%
    addTiles() %>%
    addPolygons(
      fillColor = ~pal(st_dbscan_cluster),  
      fillOpacity = 0.7,
      color = "black",                    
      weight = 1,                          
      highlightOptions = highlightOptions(weight = 2, color = "#666", fillOpacity = 0.7, bringToFront = TRUE),
      popup = ~paste0(
        "<div style='font-size: 14px;'><b>State:</b> ", shapeName, "<br>",
        "<b>Cluster:</b> ", st_dbscan_cluster, "<br>",
        "<b>Tree percentage:</b> ", round(tree_percentage, 2), "<br>"
      ), 
      layerId = ~shapeName
    ) %>%
    addLegend(
      pal = pal,
      values = ~st_dbscan_cluster,    
      title = "ST-DBSCAN Clusters",
      position = "bottomright"
    ) %>%
    setView(lng = 15, lat = 48, zoom = 4) %>%  
    addControl("<strong>Spatio-Temporal Clustering Map</strong>", 
               position = "bottomleft", 
               className = "map-title")
  
  tree_percentage_and_cluster <- time_series_data_with_geometry %>% 
    select(shapeName, tree_percentage, year, longitude, latitude, st_dbscan_cluster)
  
  
  # Choose minpts (you can change this based on your data)
  minpts <- 15
  
  # Compute the distance matrix for spatial data
  distance_matrix_spatial <- as.matrix(dist(cbind(spatial_coords[, 1], spatial_coords[, 2])))
  
  # Sort the distances for each point and pick the minpts-th nearest neighbor distance
  k_distance <- apply(distance_matrix_spatial, 1, function(row) sort(row)[minpts])
  
  # Plot the k-distance values to find the "elbow"
  plot(sort(k_distance), type = "l", main = "K-distance plot", xlab = "Points sorted by distance", ylab = paste(minpts, "-th nearest neighbor distance"))
  

  # Visualize the feature distance distribution
  distance_matrix_feature <- as.matrix(dist(temporal_data))  # Feature data like tree_percentage
  
  hist(distance_matrix_feature, breaks = 50, main = "Feature Distance Distribution", xlab = "Feature Distance")
  
  
  #MORGAN I SPACE
  library(spdep)
  # Create neighbors based on the geometry column (using "queen" contiguity)
  nb <- poly2nb(data_with_state_boarders$geometry, queen = TRUE)
  
  # Create spatial weights matrix from neighbors list
  lw <- nb2listw(nb, style = "W")  # Style "W" for row-standardized weights

  # Extract the variable of interest (e.g., tree percentage)
  tree_percentage <- data_with_state_boarders$tree_percentage_2021
  
  # Calculate Moran's I
  morans_i <- moran.test(tree_percentage, lw)
  print(morans_i)
  
  # Calculate spatial lag
  tree_percentage_lag <- lag.listw(lw, tree_percentage)
  
  # Create a data frame for plotting
  moran_plot_data <- data.frame(
    tree_percentage = tree_percentage,
    tree_percentage_lag = tree_percentage_lag
  )
  
  # Plot Moran's I scatterplot
  ggplot(moran_plot_data, aes(x = tree_percentage, y = tree_percentage_lag)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +  # Linear fit line
    labs(title = "Moran's I Plot for Tree Percentage",
         x = "Tree Percentage",
         y = "Spatial Lag of Tree Percentage") +
    theme_minimal()
  
  
  #MORAN I TIME
  
  create_temporal_weights <- function(data, year_column) {
    # Extract unique years
    years <- unique(data[[year_column]])
    n <- length(years)
    
    # Create an empty temporal weights matrix
    temporal_weights <- matrix(0, nrow = n, ncol = n, dimnames = list(years, years))
    
    # Define neighbors as consecutive years (1-step lag)
    for (i in 2:n) {
      temporal_weights[i, i-1] <- 1  # Previous year is a neighbor
      temporal_weights[i-1, i] <- 1  # Next year is a neighbor
    }
    
    # Convert to a listw object for Moran's I calculation
    lw_temporal <- mat2listw(temporal_weights, style = "W")  # Row-standardized weights
    return(lw_temporal)
  }
  
  # Assuming the year column in your long data is "year"
  lw_temporal <- create_temporal_weights(time_series_data, "year")
  
  morans_i_temporal <- moran.test(time_series_data$tree_percentage, lw_temporal)
  
  # View the results
  print(morans_i_temporal)
  
