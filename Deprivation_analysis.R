library(tidyverse)
library(sf)
library(randomForest)
library(tmap)
library(corrplot)
library(caret)
library(vip)
library(pdp)
library(ggplot2)
library(gridExtra)
library(spdep)
library(GWmodel)
library(RColorBrewer)
library(scales)
library(lmtest)
library(sp)
library(dplyr)
library(viridis)

## ----- Setup and Data Loading ---------
# Load Edinburgh SIMD data zone shapefile - with rural classification
simd <- st_read("/Users/mariamuwais/VisualAnalytics/rural_classification/rural_classification_shp.shp")

# Change NAs in RuralUrban column  'Urban'
simd$RuralUrban[simd$RuralUrban == "" | is.na(simd$RuralUrban)] <- "Urban"
head(simd, 5)

# Load 2020 SIMD data and join with shapefile
simd2020 <- read_csv("/Users/mariamuwais/VisualAnalytics/SIMD_Edinburgh_2020csv.csv")
simd <- simd %>%
  left_join(simd2020, by = c("DataZone" = "Data_Zone"))

# Check data types  and NA values
summary(simd) #All columns needed ofr analysis (2020 SIMD rank, 2020 domain ranks) are numeric

cat("NA count per column:\n")
print(colSums(is.na(simd))) # 0 NA values in 2020 SIMD rank and 2020 domain ranks need replacing 

# Key domain rank columns (lower rank = more deprived in SIMD)
domains <- c("SIMD2020v2_Income_Domain_Rank", "SIMD2020_Employment_Domain_Rank", "SIMD2020_Health_Domain_Rank",
             "SIMD2020_Education_Domain_Rank", "SIMD2020_Access_Domain_Rank", "SIMD2020_Crime_Domain_Rank", "SIMD2020_Housing_Domain_Rank")
domain_labels <- c("Income", "Employment", "Health",
                   "Education", "Access", "Crime", "Housing")

# Seperate Rural and Urban Subsets
rural <- simd %>% filter(RuralUrban == "Rural")
urban <- simd %>% filter(RuralUrban == "Urban")

## ------------------------------------------------- GEOGRAPHICALLY WEIGHTED REGRESSION --------------------------------------------------
# 1. Check projection to BNG CRS (for metre-based distances)
simd_bng <- st_transform(simd, crs = 27700)
names(simd_bng)

# 2. Convert to SpatialPointsDataFrame via centroids (GWmodel requires sp-class objects with point geometry)
centroids_sf <- st_centroid(simd_bng)
centroids_sp <- as(centroids_sf, "Spatial")

# 3. Define GWR formula
formula_gwr <- SIMD2020v2_Rank ~ SIMD2020v2_Income_Domain_Rank + SIMD2020_Employment_Domain_Rank + 
  SIMD2020_Health_Domain_Rank + SIMD2020_Education_Domain_Rank + SIMD2020_Access_Domain_Rank +
  SIMD2020_Crime_Domain_Rank + SIMD2020_Housing_Domain_Rank 

# 4. Global OLS baseline
ols_model <- lm(formula_gwr, data = centroids_sp@data)
summary(ols_model)

# 5. Select optimal adaptive bandwidth (number of nearest neighbours used per local regression)
bw_adapt <- bw.gwr(
  formula = formula_gwr,
  data = centroids_sp,
  approach = "AICc",
  kernel = "bisquare",
  adaptive = TRUE
)

# 6. Run GWR
gwr_model <- gwr.basic(
  formula = formula_gwr,
  data = centroids_sp,
  bw = bw_adapt,
  kernel = "bisquare",
  adaptive = TRUE,
  F123.test = TRUE 
)

print(gwr_model)

# 7. Attach GWR results back to the polygon sf object
simd_bng$local_R2 <- gwr_model$SDF$Local_R2
simd_bng$coef_income <- gwr_model$SDF$SIMD2020v2_Income_Domain_Rank
simd_bng$coef_employ <- gwr_model$SDF$SIMD2020v2_Employment_Domain_Rank
simd_bng$coef_health <- gwr_model$SDF$SIMD2020v2_Health_Domain_Rank
simd_bng$coef_educ <- gwr_model$SDF$SIMD2020v2_Education_Domain_Rank
simd_bng$coef_access <- gwr_model$SDF$SIMD2020v2_Access_Domain_Rank
simd_bng$coef_crime <- gwr_model$SDF$SIMD2020v2_Crime_Domain_Rank
simd_bng$coef_housing <- gwr_model$SDF$SIMD2020v2_Housing_Domain_Rank
simd_bng$gwr_residuals <- gwr_model$SDF$residual

# 8. Map Local R-squared
ggplot(simd_bng) +
  geom_sf(aes(fill = local_R2), colour = NA) +
  scale_fill_viridis_c(option = "magma", name = "Local R Squared") +
  labs(title = "GWR Local Squared",
       subtitle = "SIMD 2020 - Overall Rank, Edinburgh Data Zone") + 
  theme_minimal()

ggsave("gwr_localr2.png", width = 8, height = 7, dpi = 300)

# 9. Map Local coefficients for each domain
plot_coef <- function(col, label) {
  ggplot(simd_bng) +
    geom_sf(aes(fill = .data[[col]]), colour = NA) +
    scale_fill_viridis_c(option = "plasma", name = "Coefficient") + 
    labs(title = paste("GWR Local Coefficient -", label),
         subtitle = "SIMD 2020, Edinburgh Data Zones") +
    theme_minimal()
}

coef_cols <- c("coef_income", "coef_employ", "coef_health", "coef_educ", "coef_access", "coef_crime", "coef_housing")
coef_labels <- c("Income", "Employment", "Health", "Education", "Access", "Crime", "Housing")

for (i in seq_along(coef_cols)) {
  p <- plot_coef(coef_cols[i], coef_labels[i])
  ggsave(paste0("gwr_coef_", tolower(coef_labels[i]), ".png"),
         plot = p, width = 8, height = 7, dpi = 300)
}

# 10. Map GWR Residuals
# Extract results to sf for mapping 
results_df <- as.data.frame(gwr_model$SDF)


#Rename coeficient columns to short domain labels
coef_cols_raw <- paste0(domains)
for (i in seq_along(domains)) {
  names(results_df)[names(results_df) == domains[i]] <- domain_labels[i]
}

shp_results <- simd_bng %>%
  bind_cols(results_df%>%
              select(all_of(domain_labels),
                     Local_R2 = "Local_R2",
                     Residual = "residual"))


# 11. Map a single coefficient 
map_coefficient <- function(data, col, title, subtitle = NULL) {
  ggplot(data) + 
    geom_sf(aes(fill = .data[[col]]), colour ="white", linewidth = 0.08) +
    scale_fill_distiller(
      palette = "RdYlBu",
      direction = -1,
      name = "Coefficient",
      labels = number_format(accuracy = 0.01)
    ) + 
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      plot.subtitle = element_text(hjust = 0.5, colour = 'grey40', size = 8),
      legend.key.height = unit(0.8, "cm")
    )
}

# 12. Map local coefficnets for each domain 
coef_maps <- lapply(seq_along(domain_labels), function(i) {
  map_coefficient(
    shp_results,
    col = domain_labels[i],
    title = domain_labels[i],
    subtitle = 'Local GWR coefficient'
  )
})

#Arrange all 7 domain maps in a grid
n_cols <- 4
n_rows <- ceiling(length(domain_labels)/ n_cols)
p_coefs <- arrangeGrob(
  grobs = coef_maps,
  ncol = n_cols,
  top = "GWR Local Coefficients of SIMD Domains\nRed = stronger positive effect on deprivation rank"
)

ggsave("gwr_domain+coefficients.png", p_coefs,
       width = 14, height = n_rows * 4.5, dpi = 150)

# 13. Map GWR Residuals 
p_resid <- ggplot(shp_results) +
  geom_sf(aes(fill = Residual), colour = "white", linewidth = 0.08) +
  scale_fill_distiller(
    palette   = "RdBu",
    direction = -1,
    name      = "Residual",
    labels    = comma
  ) +
  labs(
    title    = "GWR Residuals of Domain Ranks ",
    subtitle = "Red = model under-predicts deprivation | Blue = model over-predicts"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey40", size = 9)
  )

ggsave("gwr_residuals.png", p_resid, width = 8, height = 8, dpi = 150)

# 14. Map — OLS vs GWR residuals comparison 
simd_bng$OLS_Residuals <- residuals(ols_model) #Save global residuals for mapping

p_ols_resid <- ggplot(simd_bng) +
  geom_sf(aes(fill = OLS_Residuals), colour = "white", linewidth = 0.08) +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       name = "Residual", labels = comma) +
  labs(title = "OLS Residuals",
       subtitle = "Global model — same coefficients everywhere") +
  theme_void(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, colour = "grey40", size = 8))

p_gwr_resid2 <- ggplot(shp_results) +
  geom_sf(aes(fill = Residual), colour = "white", linewidth = 0.08) +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       name = "Residual", labels = comma) +
  labs(title = "GWR Residuals",
       subtitle = "Local model — coefficients vary by location") +
  theme_void(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, colour = "grey40", size = 8))

p_resid_compare <- arrangeGrob(
  p_ols_resid, p_gwr_resid2,
  ncol = 2,
  top  = "OLS vs GWR Residuals — Domain Ranks"
)

ggsave("gwr_ols_vs_gwr_residuals.png", p_resid_compare,
       width = 14, height = 7, dpi = 150)


# 15. Coefficient range plot 
# Shows how much each domain's effect varies across Edinburgh
coef_ranges <- map_dfr(domain_labels, function(lbl) {
  vals <- shp_results[[lbl]]
  data.frame(
    Domain  = lbl,
    Min     = min(vals, na.rm = TRUE),
    Q25     = quantile(vals, 0.25, na.rm = TRUE),
    Median  = median(vals, na.rm = TRUE),
    Q75     = quantile(vals, 0.75, na.rm = TRUE),
    Max     = max(vals, na.rm = TRUE),
    IQR     = IQR(vals, na.rm = TRUE),
    Global  = coef(ols_model)[domains[match(lbl, domain_labels)]]
  )
})

p_coef_range <- ggplot(coef_ranges,
                       aes(x = reorder(Domain, IQR))) +
  geom_linerange(aes(ymin = Min, ymax = Max),
                 colour = "grey70", linewidth = 0.8) +
  geom_linerange(aes(ymin = Q25, ymax = Q75),
                 colour = "#2C7BB6", linewidth = 2) +
  geom_point(aes(y = Median), colour = "#2C7BB6", size = 3) +
  geom_point(aes(y = Global), colour = "#D7191C",
             size = 3, shape = 18) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  coord_flip() +
  labs(
    title    = "GWR Coefficient Ranges by Domain",
    subtitle = "Blue bar = IQR across data zones | Red diamond = global OLS coefficient\nWider bar = more spatial variation in that domain's effect",
    x        = NULL,
    y        = "Coefficient value"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave("gwr_coefficient_ranges.png", p_coef_range,
       width = 9, height = 6, dpi = 150)

# 16. Export results 
output_csv <- shp_results %>%
  st_drop_geometry() %>%
  select(DataZone, Intermediate_Zone, SIMD2020v2_Rank,
         all_of(domain_labels), Local_R2, Residual) %>%
  arrange(DataZone)

write.csv(output_csv, "gwr_results.csv", row.names = FALSE)

## ----------------------------------------------------RANDOM FOREST ANALYSIS------------------------------------------------------------
# 1. split train and test data 
target_col <- "SIMD2020v2_Rank"

model_df <- simd %>%
  select(all_of(c(target_col, domains)))  %>%
  drop_na()

set.seed(42)
train_idx  <- createDataPartition(model_df[[target_col]], p = 0.75, list = FALSE)
train_data <- model_df[train_idx, ]
test_data  <- model_df[-train_idx, ]

cat("Training rows:", nrow(train_data), "| Test rows:", nrow(test_data), "\n\n")

# 2. Fit Random Forest (regression) 
set.seed(42)
rf_model <- randomForest(
  x        = train_data[, domains],
  y        = train_data[[target_col]],
  ntree    = 500,           # number of trees
  mtry     = floor(sqrt(length(domains))),  # features per split
  importance = TRUE,
  keep.forest = TRUE
)

print(rf_model)

# 3. Evaluate on test set 
preds      <- predict(rf_model, newdata = test_data[, domains])
actuals    <- test_data[[target_col]]

rmse       <- sqrt(mean((preds - actuals)^2))
mae        <- mean(abs(preds - actuals))
r_squared  <- cor(preds, actuals)^2

cat("\n── Test Set Performance ──────────────────────────────\n")
cat(sprintf("  RMSE      : %.1f  (ranks — out of ~6,976 Scotland-wide)\n", rmse))
cat(sprintf("  MAE       : %.1f\n", mae))
cat(sprintf("  R-squared : %.3f\n", r_squared))

# 4. Feature importance plot 
imp <- importance(rf_model, type = 1)   # %IncMSE
imp_df <- data.frame(
  Feature    = rownames(imp),
  IncMSE     = imp[, 1]
) %>% arrange(desc(IncMSE))

p_imp <- ggplot(imp_df, aes(x = reorder(Feature, IncMSE), y = IncMSE)) +
  geom_col(fill = "#2C7BB6", alpha = 0.85) +
  coord_flip() +
  labs(
    title    = "Feature Importance: Edinburgh SIMD 2020 Random Forest",
    subtitle = "% Increase in MSE when feature is permuted (higher = more important)",
    x        = NULL,
    y        = "% Increase in MSE"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave("feature_importance.png", p_imp, width = 10, height = 8, dpi = 150)

# 5. Domain-level importance (grouped) 
domain_map <- c(
  SIMD2020v2_Income_Domain_Rank      = "Income",
  SIMD2020_Employment_Domain_Rank    = "Employment",
  SIMD2020_Health_Domain_Rank        = "Health",
  SIMD2020_Education_Domain_Rank     = "Education",
  SIMD2020_Access_Domain_Rank        = "Access",
  SIMD2020_Crime_Domain_Rank         = "Crime",
  SIMD2020_Housing_Domain_Rank       = "Housing",
  income_rate = "Income", income_count = "Income",
  employment_rate = "Employment", employment_count = "Employment",
  CIF = "Health", ALCOHOL = "Health", DRUG = "Health",
  SMR = "Health", DEPRESS = "Health", LBWT = "Health", EMERG = "Health",
  Attendance = "Education", Attainment = "Education",
  no_qualifications = "Education", not_participating = "Education",
  University = "Education",
  crime_count = "Crime", crime_rate = "Crime",
  overcrowded_count = "Housing", nocentralheating_count = "Housing",
  overcrowded_rate = "Housing", nocentralheating_rate = "Housing",
  drive_petrol = "Access", drive_GP = "Access", drive_post = "Access",
  drive_primary = "Access", drive_retail = "Access", drive_secondary = "Access",
  PT_GP = "Access", PT_post = "Access", PT_retail = "Access",
  broadband = "Access"
)
domain_imp <- imp_df %>%
  mutate(Domain = domain_map[Feature]) %>%
  group_by(Domain) %>%
  summarise(Total_IncMSE = sum(IncMSE), .groups = "drop") %>%
  arrange(desc(Total_IncMSE))

p_domain <- ggplot(domain_imp, aes(x = reorder(Domain, Total_IncMSE),
                                   y = Total_IncMSE, fill = Domain)) +
  geom_col(show.legend = FALSE, alpha = 0.85) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Domain-Level Importance (Edinburgh SIMD 2020)",
    subtitle = "Sum of %IncMSE across all indicators within each domain",
    x        = NULL,
    y        = "Cumulative % Increase in MSE"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("domain_importance.png", p_domain, width = 7, height = 5, dpi = 150)
print(domain_imp)

# 6. Correlation matrix of domain ranks 
cor_matrix <- cor(model_df[, domains], use = "complete.obs")
colnames(cor_matrix) <- rownames(cor_matrix) <-
  c("Income", "Employment", "Health", "Education", "Access", "Crime", "Housing", "Overall")

png("domain_correlation.png", width = 700, height = 650, res = 120)
corrplot(cor_matrix, method = "color", type = "lower",
         addCoef.col = "black", number.cex = 0.75,
         tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("#D7191C", "white", "#2C7BB6"))(200),
         title = "Domain Rank Correlations — Edinburgh SIMD 2020",
         mar = c(0, 0, 2, 0))
dev.off()
cat("Saved: domain_correlation.png\n")


# 7. Summary 
cat("\n══ SUMMARY ══════════════════════════════════════════════\n")
cat("Top 5 most important features:\n")
print(head(imp_df, 5))
cat("\nMost important domain overall:", domain_imp$Domain[1], "\n")
cat(sprintf("Model R² on unseen data: %.3f\n", r_squared))

## ------------------------------------------- REWEIGHTING SIMD ------------------------------------------
# 1. Function to extract domain importances from random forest 
get_rf_importance <- function(data, domains, target = "SIMD2020v2_Rank") {
  
  df <- data %>%
    st_drop_geometry() %>%
    select(all_of(c(target, domains))) %>%
    drop_na()
  
  rf_model <- randomForest(
    x = df[, domains],
    y = df[[target]],
    ntree = 500,
    importance = TRUE #model calculates how much each domain contributes to predictive accuracy 
  )
  
  importance(rf_model, type = 1) %>%   # %IncMSE (percentage increase in mean squared error) importance
    as.data.frame() %>%
    rownames_to_column("domain") %>%
    rename(importance = `%IncMSE`) %>%
    mutate(importance = pmax(importance, 0))  # floors any negative importance values at 0
}

rural_importance <- get_rf_importance(rural, domains)
urban_importance <- get_rf_importance(urban, domains)

# 2. Normalise importances to sum to 1 to become the weights
# Convert raw importance scores into proportional weights that sum to 1, making them directly usable as index weights
rural_weights <- rural_importance %>%  # divide each domains importance by the total importance across all domains
  mutate(weight = importance / sum(importance))

urban_weights <- urban_importance %>%
  mutate(weight = importance / sum(importance))
print(rural_weights)
print(urban_weights)

# 3. Compare weights 
bind_rows(
  rural_weights %>% mutate(context = "rural"),
  urban_weights %>% mutate(context = "urban")
) %>%
  ggplot(aes(x = reorder(domain, weight), y = weight, fill = context)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_x_discrete(labels = c(
    "SIMD2020v2_Income_Domain_Rank" = "Income",
    "SIMD2020v2_Employment_Domain_Rank" = "Employment",
    "SIMD2020v2_Health_Domain_Rank" = "Health",
    "SIMD2020v2_Education_Domain_Rank" = "Education",
    "SIMD2020v2_Access_Domain_Rank" = "Access",
    "SIMD2020v2_Crime_Domain_Rank" = "Crime",
    "SIMD2020v2_Housing_Domain_Rank" = "Housing",
  )) +
  labs(title = "Domain weights by rural/urban context",
       x = "Domain", y = "Normalised importance weight") +
  theme_minimal()

# 4. Construct the Reweighted Index
# Replace SIMD weights with RF-derived weights per context (rural or urban)
compute_weighted_rank <- function(data, weights_df) {
  w <- setNames(weights_df$weight, weights_df$domain)   # Pull weights as named vector
  data %>%
    st_drop_geometry() %>%
    select(data_zone, all_of(domains)) %>%
    mutate(
      # Weighted sum of domain ranks
      # Note: ranks are already normalised within SIMD (1 = most deprived)
      # So lower weighted sum = more deprived
      reweighted_score = rowSums(                       #rowSums() adds the weighted domain values together to produce a single composite score per data zone
        across(all_of(domains), ~ .x * w[cur_column()]) #multiplies each domain column by its corresponding weight 
      )
    ) %>%
    mutate(
      # Convert score to rank (1 = most deprived)
      reweighted_rank = rank(reweighted_score)
    ) %>%
    select(data_zone, reweighted_score, reweighted_rank)
}

rural_reweighted <- compute_weighted_rank(rural, rural_weights)
urban_reweighted <- compute_weighted_rank(urban, urban_weights)
print(rural_reweighted)

# 5. Combine and join back to spatial data
reweighted_all <- bind_rows(rural_reweighted, urban_reweighted)

simd_bng <- simd_bng %>%
  left_join(reweighted_all, by = "data_zone")

# 6. Compare Standard vs Reweighted Rankings 
simd_bng <- simd_bng %>%
  mutate(
    rank_change = SIMD2020v2_Rank - reweighted_rank,
    # Positive = zone is MORE deprived under new index than SIMD suggested
    direction = case_when(
      rank_change > 50  ~ "More deprived (underestimated by SIMD)",
      rank_change < -50 ~ "Less deprived (overestimated by SIMD)",
      TRUE              ~ "Broadly stable"
    )
  )

names(simd_bng)
# 7. Summarise rank changes by rural/urban
simd_bng %>%
  st_drop_geometry() %>%
  group_by(RuralUrban) %>%
  summarise(
    mean_rank_change = mean(rank_change, na.rm = TRUE),
    sd_rank_change   = sd(rank_change, na.rm = TRUE),
    n_underestimated = sum(rank_change > 50, na.rm = TRUE),
    n_overestimated  = sum(rank_change < -50, na.rm = TRUE)
  )

# 8. Map the Results
tmap_mode("plot")
# Map 1: Standard SIMD rank
m1 <- tm_shape(simd_bng) +
  tm_fill("SIMD2020v2_Rank",
          palette = "RdYlBu",
          n = 5,
          style = "quantile",
          title = "Standard SIMD Rank") +
  tm_borders(alpha = 0.3) +
  tm_layout(title = "Standard SIMD")

# Map 2: Reweighted rank
m2 <- tm_shape(simd_bng) +
  tm_fill("reweighted_rank",
          palette = "RdYlBu",
          n = 5,
          style = "quantile",
          title = "Reweighted Rank") +
  tm_borders(alpha = 0.3) +
  tm_layout(title = "Rurally-Sensitive SIMD")

# Map 3: Rank change 
m3 <- tm_shape(simd_bng) +
  tm_fill("rank_change",
          palette = "-RdBu",
          midpoint = 0,
          style = "quantile",
          title = "Rank Change") +
  tm_borders(alpha = 0.3) +
  tm_layout(title = "Change: Positive = more deprived under new index")

tmap_arrange(m1, m2, m3, ncol = 3)

