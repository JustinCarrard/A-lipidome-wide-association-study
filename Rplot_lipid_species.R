#import_library
library(readxl)
library(car)
library(broom)
library(dplyr)
library(tidyr)
library(ggdendro)
library(egg)
library(ggplot2)
library(reshape2)
library(emmeans)

#install.packages("writexl")
library("writexl")

#install.packages("data.table")
library(data.table)

# Paths
getwd()
setwd("/Volumes/nuxyde32/Uni-Basel/DSBG/Forschung/COmPLETE Health/Metabolomics/Pilotstudie/Book chapter/R") # Main path

data_path <- "./data" # Path for data
graphics_path <- "./output/graphics" # Path for graphics
text_path <- "./output/text" # Path for text-output

#import_data
dat <- as.data.frame(read_excel(paste0(data_path, "/", "data_book_chapter.xlsx")))

head(dat)
str(dat)

#remove wrong date from sampling time column
dat$"Sampling time" = gsub("1899-12-31","",dat$"Sampling time")

#Group sampling time into 5 categories
dat$Sampling_time_cat <- cut(
  chron::times(dat$"Sampling time")
  , breaks = chron::times(c(
    "08:00:00"
    , "10:00:00"
    , "12:00:00"
    , "14:00:00"
    , "16:00:00"
    , "23:59:59"
  ))
  , labels = c("8-9.59","10-11.59","12-13.59","14-15.59","16-18.59")
  , include.lowest = TRUE
  , right = TRUE
)

dat[, c("Sampling time", "Sampling_time_cat")]

#log2_transformation_lipids
(variables_to_transform <- names(dat)[24:237]) # Name of variables to transform with log2

for (i in variables_to_transform) {
  dat[, paste0(i, "_log2")] <- log2(dat[, i])
}

names(dat)
head(dat)

#z_standardization_dependent_variables
(variables_to_standardize <- names(dat)[c(6:23, 239:452)])

for (i in variables_to_standardize) {
  dat[, paste0(i, "_std")] <- scale(dat[, i], center = TRUE, scale = TRUE)
}

names(dat)
head(dat)

# Convert Decade, Sex and Statins to factors
dat$Decade <- factor(dat$Decade)
dat$Sex <- factor(dat$Sex)
dat$Statins <- factor(dat$Statins)
dat$Sampling_time_cat <- factor(dat$Sampling_time_cat)

#218_univariate_multiple_linear_regressions
my_lms <- lapply(dat[,c(471:684)], function(x) lm(x ~ dat$Decade * dat$Sex + dat$Statins + dat$PBF_std + dat$Total_PA_std + dat$HbA1c_std + dat$Sampling_time_cat + dat$Fasting_std, data = dat))

#extract Decade results
res_frame_Decade <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {
  
  sum_tmp <- summary(my_lms[[i]])
  
  res_frame_Decade$estimate[which(res_frame_Decade$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Estimate"]["dat$Decade1"]
  res_frame_Decade$std.error[which(res_frame_Decade$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Std. Error"]["dat$Decade1"]
  res_frame_Decade$statistic[which(res_frame_Decade$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "t value"]["dat$Decade1"]
  res_frame_Decade$p.value[which(res_frame_Decade$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Pr(>|t|)"]["dat$Decade1"]
  res_frame_Decade$response[which(res_frame_Decade$term %in% names(my_lms)[i])] <- "Aged"
  
  rm(sum_tmp)
  
}

#extract Sex results
res_frame_Sex <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {
  
  sum_tmp <- summary(my_lms[[i]])
  
  res_frame_Sex$estimate[which(res_frame_Sex$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Estimate"]["dat$Sex1"]
  res_frame_Sex$std.error[which(res_frame_Sex$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Std. Error"]["dat$Sex1"]
  res_frame_Sex$statistic[which(res_frame_Sex$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "t value"]["dat$Sex1"]
  res_frame_Sex$p.value[which(res_frame_Sex$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Pr(>|t|)"]["dat$Sex1"]
  res_frame_Sex$response[which(res_frame_Sex$term %in% names(my_lms)[i])] <- "Male"
  
  rm(sum_tmp)
  
}

#extract Statins results
res_frame_Statins <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {
  
  sum_tmp <- summary(my_lms[[i]])
  
  res_frame_Statins$estimate[which(res_frame_Statins$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Estimate"]["dat$Statins1"]
  res_frame_Statins$std.error[which(res_frame_Statins$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Std. Error"]["dat$Statins1"]
  res_frame_Statins$statistic[which(res_frame_Statins$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "t value"]["dat$Statins1"]
  res_frame_Statins$p.value[which(res_frame_Statins$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Pr(>|t|)"]["dat$Statins1"]
  res_frame_Statins$response[which(res_frame_Statins$term %in% names(my_lms)[i])] <- "Statins intake"
  
  rm(sum_tmp)
  
}

#extract PBF results
res_frame_PBF <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {
  
  sum_tmp <- summary(my_lms[[i]])
  
  res_frame_PBF$estimate[which(res_frame_PBF$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Estimate"]["dat$PBF_std"]
  res_frame_PBF$std.error[which(res_frame_PBF$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Std. Error"]["dat$PBF_std"]
  res_frame_PBF$statistic[which(res_frame_PBF$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "t value"]["dat$PBF_std"]
  res_frame_PBF$p.value[which(res_frame_PBF$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Pr(>|t|)"]["dat$PBF_std"]
  res_frame_PBF$response[which(res_frame_PBF$term %in% names(my_lms)[i])] <- "Body fat (%)"
  
  rm(sum_tmp)
  
}

#extract HbA1c results
res_frame_HbA1c <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {
  
  sum_tmp <- summary(my_lms[[i]])
  
  res_frame_HbA1c$estimate[which(res_frame_HbA1c$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Estimate"]["dat$HbA1c_std"]
  res_frame_HbA1c$std.error[which(res_frame_HbA1c$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Std. Error"]["dat$HbA1c_std"]
  res_frame_HbA1c$statistic[which(res_frame_HbA1c$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "t value"]["dat$HbA1c_std"]
  res_frame_HbA1c$p.value[which(res_frame_HbA1c$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Pr(>|t|)"]["dat$HbA1c_std"]
  res_frame_HbA1c$response[which(res_frame_HbA1c$term %in% names(my_lms)[i])] <- "HbA1c"
  
  rm(sum_tmp)
  
}

#extract Total PA results
res_frame_Total_PA <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {
  
  sum_tmp <- summary(my_lms[[i]])
  
  res_frame_Total_PA$estimate[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Estimate"]["dat$Total_PA_std"]
  res_frame_Total_PA$std.error[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Std. Error"]["dat$Total_PA_std"]
  res_frame_Total_PA$statistic[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "t value"]["dat$Total_PA_std"]
  res_frame_Total_PA$p.value[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Pr(>|t|)"]["dat$Total_PA_std"]
  res_frame_Total_PA$response[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- "Daily physical activity"
  
  rm(sum_tmp)
  
}

#Combine_data_frames
Overall <- rbind(
  res_frame_Decade
  , res_frame_Sex
  , res_frame_Statins
  , res_frame_PBF
  , res_frame_HbA1c
  , res_frame_Total_PA
)

#Adjust p-values
Overall <- Overall[order(Overall$p.value),]
Overall$BH <- p.adjust(Overall$p.value, method = "BH")

#Categorise BH p-values
Overall$BH_cat <- NA

Overall$BH_cat[Overall$BH > 0.05] <- "> 0.05"
Overall$BH_cat[Overall$BH <= 0.05 & Overall$BH > 0.01] <- "≤ 0.05"
Overall$BH_cat[Overall$BH <= 0.01 & Overall$BH > 0.001] <- "≤ 0.01"
Overall$BH_cat[Overall$BH <= 0.001 & Overall$BH > 0.0001] <- "≤ 0.001"
Overall$BH_cat[Overall$BH <= 0.0001] <- "≤ 0.0001"

Overall$BH_cat <- factor(Overall$BH_cat, levels = c("> 0.05", "≤ 0.05", "≤ 0.01", "≤ 0.001", "≤ 0.0001"))

#Order by decreasing estimates
Overall <- Overall[order(-Overall$estimate),]

#remove suffixes from lipid subclasses' name
Overall$term=gsub("_log2_std","",Overall$term)

#Start_rain_plot ------------------------------------------------------
#extract the 16 lipid species displaying the biggest estimates
Overall_top <- subset(Overall, term == "PC 16:0_18:0" | term == "SM 33:1;2" | term == "SM 34:2;2" | term == "LPE O-18:1" | term == "PE O-20:1_20:4" | term == "Cer 18:1;2/16:0" | term == "PE O-18:1_22:5" | term == "PE O-18:1_18:1" | term == "LPC 22:5" | term == "SM 34:1;2" | term == "PC 17:0_18:1" | term == "LPC 17:0" | term == "PC 16:0_22:5" | term == "LPI 18:0" | term == "PC O-37:5" | term == "PC O-38:6" , 
                      select=c(term, estimate, std.error, statistic, p.value, response, BH, BH_cat))

#Import ------------------------------------------------------
plot_data <- Overall_top

# Theme+ Palette ------------------------------------------------------

## Palette

palette <-
  # Blue
  c("#053061",
    "#313695",
    "#4575b4",
    "#74add1",
    "#abd9e9",
    "#e0f3f8",
    "#fee090",
    "#fdae61",
    "#f46d43",
    "#d73027",
    "#a50026",
    '#67001f')
# Red

# Calculate symmetric limits based on most extreme value
max_abs_estimate <- max(abs(plot_data$estimate))

max_lim <- max_abs_estimate
min_lim = -1 * max_lim

## theme

thm <-
  # Good starting theme + set text size
  theme_light(base_size = 7) +
  theme(
    # Remove axis ticks and titles
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    # Remove gridlines and boxes
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    legend.key = element_blank(),
    
    # White backgrounds
    panel.background = element_rect(fill = 'white'),
    plot.background = element_rect(fill = 'white'),
    legend.background = element_rect(fill = 'white'),
    
    # Angle text
    axis.text.x.top = element_text(angle = 45, hjust = 0)
  )

# Bare-bones rainplot ------------------------------------------------------

rainplot <-
  ggplot(plot_data) +
  geom_point(aes(x = response, y = term, colour = estimate, size = BH_cat))

print(rainplot)

# Basic Rainplot ------------------------------------------------------

rainplot <-
  ggplot(plot_data) +
  geom_point(aes(x = response, y = term, colour = estimate, size = BH_cat))  +
  scale_x_discrete(position = 'top') +
  scale_size_manual(name = expression("BH p-value"), values = c(2, 4, 6, 8, 10), drop = FALSE) +
  #scale_size_manual(name = expression(italic(BH)*" p-value"), values = c(2, 4, 6, 8, 10), drop = FALSE) +
  #scale_size_area(expression(paste(-log[10]('BH p-value'))), max_size = 7) +
  #scale_size_area(expression(paste(-log[10]('Q-value'))), max_size = 7) +
  scale_color_gradientn(
    'Effect Size (β coefficient)',
    colors = palette,
    limits = c(min_lim, max_lim),
    breaks = c(min_lim,  min_lim / 2, 0 , max_lim/2, max_lim)
  ) +
  guides(colour = guide_colourbar(order = 1),
         size = guide_legend(order = 2)) +
  thm

print(rainplot)

# Ordering by Cluster ------------------------------------------------------

# Convert to matrix and reshape for clustering.
cluster_data <-
  plot_data %>%
  select(response, term, estimate) %>%
  spread(response, estimate)

rnms <-
  cluster_data$term

cluster_data <-
  cluster_data %>%
  select(-term) %>%
  as.matrix()

rownames(cluster_data) <- rnms

# cluster dependent variable terms
clust <- hclust(dist(cluster_data), method = 'ward.D2')

# `clust$order` orders `term` into clusters
term_order <-
  clust$labels[clust$order]

# Convert term to a factor, ordered by `term_order`
plot_data_clo <-
  plot_data %>%
  mutate(term = factor(term, levels = term_order))


rainplot <-
  # Use cluter ordered data
  ggplot(plot_data_clo) +
  geom_point(aes(x = response, y = term, colour = estimate, size = BH_cat)) +
  scale_x_discrete(position = 'top') +
  scale_size_manual(name = expression("BH p-value"), values = c(2, 4, 6, 8, 10), drop = FALSE) +
  #scale_size_manual(name = expression(italic(BH)*" p-value"), values = c(2, 4, 6, 8, 10), drop = FALSE) +
  #scale_size_area(expression(paste(-log[10]('BH p-value'))), max_size = 7) +
  #scale_size_area(expression(paste(-log[10]('Q-value'))), max_size = 7) +
  scale_color_gradientn(
    'Effect Size\n(β coefficient)',
    colors = palette,
    limits = c(min_lim, max_lim),
    breaks = c(min_lim, min_lim / 2, 0 , max_lim / 2, max_lim)
  ) +
  guides(colour = guide_colourbar(order = 1),
         size = guide_legend(order = 2)) +
  thm

print(rainplot)


# Adding dendrograms ------------------------------------------------------

dendro_dat <- segment(dendro_data(clust))


dendro <-
  # Empty ggplot with same y-scale as rainplot
  ggplot() +
  geom_blank(aes(y = term), data = plot_data) +
  theme_dendro() +
  # 'expand' controls whitespace around the dendrogram. The non-zero argument
  # may need to be increasesed if the line thickness of the dendrogram is
  # increased to make sure the entire dendrogram is plotted
  scale_x_discrete(position = 'top', expand = c(0, 0.03, 0, 0)) +
  # Draw dendrogram
  geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend),
               colour = 'black',
               data = dendro_dat)


p <- ggarrange(dendro, rainplot, ncol = 2, widths = c(1, 5))

p

#export rain plot as a .png
ggsave(paste(graphics_path, paste("rain_plot_top16_lipid_species.png", sep = ""), sep = "/"), p, width = 4.3*1, height = 6*1, units = "in", dpi = 300)

#rename columns
names(Overall)[1] <- "Dependent_variable"
names(Overall)[2] <- "β coefficient"
names(Overall)[3] <- "standard error"
names(Overall)[5] <- "p-value"
names(Overall)[6] <- "Independent variables"
names(Overall)[7] <- "BH p-value"
names(Overall)[8] <- "Categorical BH p-value"

#export excel file
write_xlsx(Overall,"/Volumes/nuxyde32/Uni-Basel/DSBG/Forschung/COmPLETE Health/Metabolomics/Pilotstudie/Book chapter/R/output/text//assocations.xlsx")

#------------------------------------------------------------------------------------------
# Preparation for ANCOVAs
#------------------------------------------------------------------------------------------

# Define variables to plot
(vars_to_plot <- names(dat)[c(471:684)])

# Reshape to long
dat_long <- reshape2::melt(
  dat
  , id.vars = c("Decade", "Sex", "Statins", "PBF_std", "Total_PA_std", "HbA1c_std", "Sampling_time_cat","Fasting_std")
  , measure.vars = vars_to_plot
)

#rename columns
names(dat_long)[9] <- "Lipid_species"

# Recode factors
dat_long$Sex <- factor(dat_long$Sex)
dat_long$Decade <- factor(dat_long$Decade)
dat_long$Sampling_time_cat <- factor(dat_long$Sampling_time_cat)

str(dat_long)
head(dat_long)

#------------------------------------------------------------------------------------------
# ANCOVAs
#------------------------------------------------------------------------------------------

p_val_frame <- data.frame(
  Lipid_species = character()
  , p_val_sex = numeric()
  , p_val_decade = numeric()
  , p_val_interaction = numeric()
  , mean_diff_decade_f = numeric()
  , mean_diff_decade_m = numeric()
  , mean_diff_sex_y = numeric()
  , mean_diff_sex_o = numeric()
  , p_mean_diff_decade_f = numeric() 
  , p_mean_diff_decade_m = numeric() 
  , p_mean_diff_sex_y = numeric() 
  , p_mean_diff_sex_o = numeric()
)


combinations <- unique(dat_long$Lipid_species)

for (i in seq_along(combinations)) {
  
  dat_tmp <- subset(dat_long, Lipid_species %in% combinations[i])

  mod <- lm(value ~ Decade*Sex + Statins + PBF_std + Total_PA_std + HbA1c_std + Sampling_time_cat + Fasting_std, data = dat_tmp)
  
  ANOVA_res <- Anova(mod, type = "II")
  
  mean_diff_decade <- as.data.frame(emmeans::contrast(emmeans(mod, "Decade", by = "Sex"), "revpairwise"))
  mean_diff_sex <- as.data.frame(emmeans::contrast(emmeans(mod, "Sex", by = "Decade"), "pairwise"))
  
  tmp_frame <- data.frame(
    Lipid_species = combinations[i]
    , p_val_sex = ANOVA_res$`Pr(>F)`[which(rownames(ANOVA_res) %in% "Sex")]
    , p_val_decade = ANOVA_res$`Pr(>F)`[which(rownames(ANOVA_res) %in% "Decade")]
    , p_val_interaction = ANOVA_res$`Pr(>F)`[which(rownames(ANOVA_res) %in% "Decade:Sex")]
    , mean_diff_decade_f = mean_diff_decade$estimate[mean_diff_decade$Sex %in% "0"]
    , mean_diff_decade_m = mean_diff_decade$estimate[mean_diff_decade$Sex %in% "1"]
    , mean_diff_sex_y = mean_diff_sex$estimate[mean_diff_sex$Decade %in% "0"]
    , mean_diff_sex_o = mean_diff_sex$estimate[mean_diff_sex$Decade %in% "1"]
    , p_mean_diff_decade_f = mean_diff_decade$p.value[mean_diff_decade$Sex %in% "0"]
    , p_mean_diff_decade_m = mean_diff_decade$p.value[mean_diff_decade$Sex %in% "1"]
    , p_mean_diff_sex_y = mean_diff_sex$p.value[mean_diff_sex$Decade %in% "0"]
    , p_mean_diff_sex_o = mean_diff_sex$p.value[mean_diff_sex$Decade %in% "1"]
  )
  
  p_val_frame <- rbind(p_val_frame, tmp_frame)
  
  rm(dat_tmp, mod, ANOVA_res, tmp_frame, mean_diff_decade, mean_diff_sex)
  
}

p_vals_adjusted_sex_decade <- p.adjust(c(
  p_val_frame$p_val_sex
  , p_val_frame$p_val_decade
  , p_val_frame$p_val_interaction
), method = "BH")

p_val_frame$p_val_sex_adj <- p_vals_adjusted_sex_decade[seq(1, dim(p_val_frame)[1])]
p_val_frame$p_val_decade_adj <- p_vals_adjusted_sex_decade[seq(1, dim(p_val_frame)[1]) + dim(p_val_frame)[1]]
p_val_frame$p_val_interaction_adj <- p_vals_adjusted_sex_decade[seq(1, dim(p_val_frame)[1]) + 2*dim(p_val_frame)[1]]


p_vals_adjusted_fmyo <- p.adjust(c(
  p_val_frame$p_mean_diff_decade_f
  , p_val_frame$p_mean_diff_decade_m
  , p_val_frame$p_mean_diff_sex_y
  , p_val_frame$p_mean_diff_sex_o
), method = "BH")

p_val_frame$p_val_decade_f_adj <- p_vals_adjusted_fmyo[seq(1, dim(p_val_frame)[1]) + 0*dim(p_val_frame)[1]]
p_val_frame$p_val_decade_m_adj <- p_vals_adjusted_fmyo[seq(1, dim(p_val_frame)[1]) + 1*dim(p_val_frame)[1]]
p_val_frame$p_val_sex_y_adj <-    p_vals_adjusted_fmyo[seq(1, dim(p_val_frame)[1]) + 2*dim(p_val_frame)[1]]
p_val_frame$p_val_sex_o_adj <-    p_vals_adjusted_fmyo[seq(1, dim(p_val_frame)[1]) + 3*dim(p_val_frame)[1]]


p_val_frame$p_val_sex_adj_cat <- NA
p_val_frame$p_val_decade_adj_cat <- NA
p_val_frame$p_val_interaction_adj_cat <- NA

p_val_frame$p_val_sex_adj_cat[p_val_frame$p_val_sex_adj > 0.05] <- "> 0.05"
p_val_frame$p_val_sex_adj_cat[p_val_frame$p_val_sex_adj <= 0.05 & p_val_frame$p_val_sex_adj > 0.01] <- "≤ 0.05"
p_val_frame$p_val_sex_adj_cat[p_val_frame$p_val_sex_adj <= 0.01 & p_val_frame$p_val_sex_adj > 0.001] <- "≤ 0.01"
p_val_frame$p_val_sex_adj_cat[p_val_frame$p_val_sex_adj <= 0.001 & p_val_frame$p_val_sex_adj > 0.0001] <- "≤ 0.001"
p_val_frame$p_val_sex_adj_cat[p_val_frame$p_val_sex_adj <= 0.0001] <- "≤ 0.0001"

p_val_frame$p_val_decade_adj_cat[p_val_frame$p_val_decade_adj > 0.05] <- "> 0.05"
p_val_frame$p_val_decade_adj_cat[p_val_frame$p_val_decade_adj <= 0.05 & p_val_frame$p_val_decade_adj > 0.01] <- "≤ 0.05"
p_val_frame$p_val_decade_adj_cat[p_val_frame$p_val_decade_adj <= 0.01 & p_val_frame$p_val_decade_adj > 0.001] <- "≤ 0.01"
p_val_frame$p_val_decade_adj_cat[p_val_frame$p_val_decade_adj <= 0.001 & p_val_frame$p_val_decade_adj > 0.0001] <- "≤ 0.001"
p_val_frame$p_val_decade_adj_cat[p_val_frame$p_val_decade_adj <= 0.0001] <- "≤ 0.0001"

p_val_frame$p_val_interaction_adj_cat[p_val_frame$p_val_interaction_adj > 0.05] <- "> 0.05"
p_val_frame$p_val_interaction_adj_cat[p_val_frame$p_val_interaction_adj <= 0.05 & p_val_frame$p_val_interaction_adj > 0.01] <- "≤ 0.05"
p_val_frame$p_val_interaction_adj_cat[p_val_frame$p_val_interaction_adj <= 0.01 & p_val_frame$p_val_interaction_adj > 0.001] <- "≤ 0.01"
p_val_frame$p_val_interaction_adj_cat[p_val_frame$p_val_interaction_adj <= 0.001 & p_val_frame$p_val_interaction_adj > 0.0001] <- "≤ 0.001"
p_val_frame$p_val_interaction_adj_cat[p_val_frame$p_val_interaction_adj <= 0.0001] <- "≤ 0.0001"

p_val_frame$p_val_sex_adj_cat <- factor(p_val_frame$p_val_sex_adj_cat, levels = c("> 0.05", "≤ 0.05", "≤ 0.01", "≤ 0.001", "≤ 0.0001"))
p_val_frame$p_val_decade_adj_cat <- factor(p_val_frame$p_val_decade_adj_cat, levels = c("> 0.05", "≤ 0.05", "≤ 0.01", "≤ 0.001", "≤ 0.0001"))
p_val_frame$p_val_interaction_adj_cat <- factor(p_val_frame$p_val_interaction_adj_cat, levels = c("> 0.05", "≤ 0.05", "≤ 0.01", "≤ 0.001", "≤ 0.0001"))

# Difference in age within gender
p_val_frame$p_val_decade_f_adj_cat <- NA
p_val_frame$p_val_decade_m_adj_cat <- NA
p_val_frame$p_val_sex_y_adj_cat <- NA
p_val_frame$p_val_sex_o_adj_cat <- NA

p_val_frame$p_val_decade_f_adj_cat[p_val_frame$p_val_decade_f_adj > 0.05] <- "> 0.05"
p_val_frame$p_val_decade_f_adj_cat[p_val_frame$p_val_decade_f_adj <= 0.05 & p_val_frame$p_val_decade_f_adj > 0.01] <- "≤ 0.05"
p_val_frame$p_val_decade_f_adj_cat[p_val_frame$p_val_decade_f_adj <= 0.01 & p_val_frame$p_val_decade_f_adj > 0.001] <- "≤ 0.01"
p_val_frame$p_val_decade_f_adj_cat[p_val_frame$p_val_decade_f_adj <= 0.001 & p_val_frame$p_val_decade_f_adj > 0.0001] <- "≤ 0.001"
p_val_frame$p_val_decade_f_adj_cat[p_val_frame$p_val_decade_f_adj <= 0.0001] <- "≤ 0.0001"

p_val_frame$p_val_decade_m_adj_cat[p_val_frame$p_val_decade_m_adj > 0.05] <- "> 0.05"
p_val_frame$p_val_decade_m_adj_cat[p_val_frame$p_val_decade_m_adj <= 0.05 & p_val_frame$p_val_decade_m_adj > 0.01] <- "≤ 0.05"
p_val_frame$p_val_decade_m_adj_cat[p_val_frame$p_val_decade_m_adj <= 0.01 & p_val_frame$p_val_decade_m_adj > 0.001] <- "≤ 0.01"
p_val_frame$p_val_decade_m_adj_cat[p_val_frame$p_val_decade_m_adj <= 0.001 & p_val_frame$p_val_decade_m_adj > 0.0001] <- "≤ 0.001"
p_val_frame$p_val_decade_m_adj_cat[p_val_frame$p_val_decade_m_adj <= 0.0001] <- "≤ 0.0001"

p_val_frame$p_val_sex_y_adj_cat[p_val_frame$p_val_sex_y_adj > 0.05] <- "> 0.05"
p_val_frame$p_val_sex_y_adj_cat[p_val_frame$p_val_sex_y_adj <= 0.05 & p_val_frame$p_val_sex_y_adj > 0.01] <- "≤ 0.05"
p_val_frame$p_val_sex_y_adj_cat[p_val_frame$p_val_sex_y_adj <= 0.01 & p_val_frame$p_val_sex_y_adj > 0.001] <- "≤ 0.01"
p_val_frame$p_val_sex_y_adj_cat[p_val_frame$p_val_sex_y_adj <= 0.001 & p_val_frame$p_val_sex_y_adj > 0.0001] <- "≤ 0.001"
p_val_frame$p_val_sex_y_adj_cat[p_val_frame$p_val_sex_y_adj <= 0.0001] <- "≤ 0.0001"

p_val_frame$p_val_sex_o_adj_cat[p_val_frame$p_val_sex_o_adj > 0.05] <- "> 0.05"
p_val_frame$p_val_sex_o_adj_cat[p_val_frame$p_val_sex_o_adj <= 0.05 & p_val_frame$p_val_sex_o_adj > 0.01] <- "≤ 0.05"
p_val_frame$p_val_sex_o_adj_cat[p_val_frame$p_val_sex_o_adj <= 0.01 & p_val_frame$p_val_sex_o_adj > 0.001] <- "≤ 0.01"
p_val_frame$p_val_sex_o_adj_cat[p_val_frame$p_val_sex_o_adj <= 0.001 & p_val_frame$p_val_sex_o_adj > 0.0001] <- "≤ 0.001"
p_val_frame$p_val_sex_o_adj_cat[p_val_frame$p_val_sex_o_adj <= 0.0001] <- "≤ 0.0001"

p_val_frame$p_val_decade_f_adj_cat <- factor(p_val_frame$p_val_decade_f_adj_cat, levels = c("> 0.05", "≤ 0.05", "≤ 0.01", "≤ 0.001", "≤ 0.0001"))
p_val_frame$p_val_decade_m_adj_cat <- factor(p_val_frame$p_val_decade_m_adj_cat, levels = c("> 0.05", "≤ 0.05", "≤ 0.01", "≤ 0.001", "≤ 0.0001"))
p_val_frame$p_val_sex_y_adj_cat <- factor(p_val_frame$p_val_sex_y_adj_cat, levels = c("> 0.05", "≤ 0.05", "≤ 0.01", "≤ 0.001", "≤ 0.0001"))
p_val_frame$p_val_sex_o_adj_cat <- factor(p_val_frame$p_val_sex_o_adj_cat, levels = c("> 0.05", "≤ 0.05", "≤ 0.01", "≤ 0.001", "≤ 0.0001"))

head(p_val_frame)

#add percentage of difference aged vs. young and male vs. female
p_val_frame$Percentage_difference_aged_vs_young_female = 100*((2^p_val_frame$mean_diff_decade_f)-1)
p_val_frame$Percentage_difference_aged_vs_young_male = 100*((2^p_val_frame$mean_diff_decade_m)-1)
p_val_frame$Percentage_difference_female_vs_male_young = 100*((2^p_val_frame$mean_diff_sex_y)-1)
p_val_frame$Percentage_difference_female_vs_male_aged = 100*((2^p_val_frame$mean_diff_sex_o)-1)

#remove suffixes from Lipids name
p_val_frame$Lipid_species=gsub("_std","",p_val_frame$Lipid_species)
p_val_frame$Lipid_species=gsub("_log2","",p_val_frame$Lipid_species)

#export excel file
write_xlsx(p_val_frame,"/Volumes/nuxyde32/Uni-Basel/DSBG/Forschung/COmPLETE Health/Metabolomics/Pilotstudie/Book chapter/R/output/text//ANCOVA.xlsx")

#------------------
#Circular bar plots
#------------------

# Important file containing subclass classification
dat_header <- as.data.frame(read_xlsx(paste0(data_path, "/", "Lipid_subclasses.xlsx")
                                      , sheet = 1
                                      , skip = 0
                                      , col_names = FALSE
                                      , col_types = "text"
                                      , trim_ws = TRUE
))

dat_header <- dat_header[1:2, ]

head(dat_header)

# Reshape the variable names

dat_header_t <- as.data.frame(t(dat_header))

names(dat_header_t) <- dat_header_t[1, ]

dat_header <- dat_header_t[-1, ]

dim(dat_header)

#### AGE DIFFERENCES

# Create dataset
dat_cbp_1 <- p_val_frame[, names(p_val_frame)[c(1, 23, 24, 27, 28)]]
dat_cbp_2 <- dat_header[, names(dat_header)[c(1, 2)]]
dat_cbp <- merge(dat_cbp_1, dat_cbp_2, by = "Lipid_species")

dat_cbp$color <- factor(dat_cbp$Lipid_subclass, levels=c("CE", "Cer", "GSL", "LPC", "LPE", "LPE-O", "LPI", "PC", "PC-O", "PE", "PE-O", "PE-P", "PI", "SM", "TG"), labels=rainbow(15))

### FEMALES
# create subset for females with p< 0.05 only
dat_cbp_f <- filter(dat_cbp, p_val_decade_f_adj_cat == "≤ 0.05"|p_val_decade_f_adj_cat =="≤ 0.01"|p_val_decade_f_adj_cat =="≤ 0.001"|p_val_decade_f_adj_cat =="≤ 0.0001")

#Define Lipid_subclass as factor
dat_cbp_f$Lipid_subclass <- factor(dat_cbp_f$Lipid_subclass)
dat_cbp_f$p_val_decade_f_adj_cat <- factor(dat_cbp_f$p_val_decade_f_adj_cat)
levels(dat_cbp_f$p_val_decade_f_adj_cat) <- c("≤ 0.05", "≤ 0.01", "≤ 0.001", "≤ 0.0001")

# Order data:
dat_cbp_f = dat_cbp_f %>% arrange(Lipid_subclass, Percentage_difference_aged_vs_young_female)

# Set a number of 'empty bar' to add at the end of each Lipid_subclass
empty_bar <- 4
to_add <- data.frame(matrix(NA, empty_bar*nlevels(dat_cbp_f$Lipid_subclass), ncol(dat_cbp_f)) )
colnames(to_add) <- colnames(dat_cbp_f)
to_add$Lipid_subclass <- rep(levels(dat_cbp_f$Lipid_subclass), each=empty_bar)
dat_cbp_f <- rbind(dat_cbp_f, to_add)
dat_cbp_f<- dat_cbp_f %>% arrange(Lipid_subclass)
dat_cbp_f$id <- seq(1, nrow(dat_cbp_f))

# prepare a data frame for base lines
base_data <- dat_cbp_f %>% 
  group_by(Lipid_subclass) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(dat_cbp_f, aes(x=as.factor(id), y=Percentage_difference_aged_vs_young_female, fill=color, alpha = factor(p_val_decade_f_adj_cat))) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=Percentage_difference_aged_vs_young_female, fill=color, alpha = factor(p_val_decade_f_adj_cat)), stat="identity") +
  
  # Add a val=60/40/20/0 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 150, xend = start, yend = 150), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(dat_cbp_f$id),4), y = c(0, 20, 40, 60), label = c("0%", "20%", "40%", "60%") , color="grey", size=2 , angle=0, fontface="bold", hjust=1) +
  annotate("text", x = rep(max(dat_cbp_f$id),5), y = c(0, 50, 100, 150, 200), label = c("0%", "50%", "100%", "150%", "200%") , color="black", size=2 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=Percentage_difference_aged_vs_young_female, fill=color, alpha = factor(p_val_decade_f_adj_cat)), stat="identity") +
  scale_alpha_discrete(name = expression("BH p-value"), range=c(0.10, 0.40, 0.70, 1)) +
  ylim(-100,230) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  theme(
    legend.position = c(0.95, 0.5),
    legend.text=element_text(size=6),
    legend.title = element_text(size = 8),) +
  scale_fill_manual(values = c("#FF0000", "#FF6600", "#FFCC00", "#CCFF00", "#66FF00", "#00FF00", "#00FF66", "#00FFCC", "#00CCFF", "#0066FF", "#0000FF", "#6600FF", "#CC00FF", "#FF00CC", "#FF0066")) +
  guides(fill = "none") +
  #geom_text(data=label_data, aes(x=id, y=Percentage_difference_aged_vs_young_female+10, label=Lipid_ID, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=1.8, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=Lipid_subclass), hjust=c(0.55, 0.55, 0.55, 0.55, 0.55, 0.7, 0.55, 0.6, 0.55, 0.55, 0.4, 0.4, 0.55, 0.55), colour = "black", alpha=0.8, size=3, fontface="bold", inherit.aes = FALSE)
p

ggsave(paste(graphics_path, paste("circular_bar_plot_aged_vs_young_female.png", sep = ""), sep = "/"), p, width = 8, height = 7, units = "in", dpi = 300)

### MALES
# create subset for males with p< 0.05 only
dat_cbp_m <- filter(dat_cbp, p_val_decade_m_adj_cat == "≤ 0.05"|p_val_decade_m_adj_cat =="≤ 0.01"|p_val_decade_m_adj_cat =="≤ 0.001"|p_val_decade_m_adj_cat =="≤ 0.0001")

#Define Lipid_subclass as factor
dat_cbp_m$Lipid_subclass <- factor(dat_cbp_m$Lipid_subclass)
dat_cbp_m$p_val_decade_m_adj_cat <- factor(dat_cbp_m$p_val_decade_m_adj_cat)
levels(dat_cbp_m$p_val_decade_m_adj_cat) <- c("≤ 0.05", "≤ 0.01", "≤ 0.001", "≤ 0.0001")

# Order data:
dat_cbp_m = dat_cbp_m %>% arrange(Lipid_subclass, Percentage_difference_aged_vs_young_male)

# Set a number of 'empty bar' to add at the end of each Lipid_subclass
empty_bar <- 4
to_add <- data.frame(matrix(NA, empty_bar*nlevels(dat_cbp_m$Lipid_subclass), ncol(dat_cbp_m)) )
colnames(to_add) <- colnames(dat_cbp_m)
to_add$Lipid_subclass <- rep(levels(dat_cbp_m$Lipid_subclass), each=empty_bar)
dat_cbp_m <- rbind(dat_cbp_m, to_add)
dat_cbp_m<- dat_cbp_m %>% arrange(Lipid_subclass)
dat_cbp_m$id <- seq(1, nrow(dat_cbp_m))

# prepare a data frame for base lines
base_data <- dat_cbp_m %>% 
  group_by(Lipid_subclass) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(dat_cbp_m, aes(x=as.factor(id), y=Percentage_difference_aged_vs_young_male, fill=color, alpha = factor(p_val_decade_m_adj_cat))) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=Percentage_difference_aged_vs_young_male, fill=color, alpha = factor(p_val_decade_m_adj_cat)), stat="identity") +
  
  # Add a val=60/40/20/0 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 150, xend = start, yend = 150), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(dat_cbp_m$id),4), y = c(0, 20, 40, 60), label = c("0%", "20%", "40%", "60%") , color="grey", size=2 , angle=0, fontface="bold", hjust=1) +
  annotate("text", x = rep(max(dat_cbp_m$id),5), y = c(0, 50, 100, 150, 200), label = c("0%", "50%", "100%", "150%", "200%") , color="black", size=2 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=Percentage_difference_aged_vs_young_male, fill=color, alpha = factor(p_val_decade_m_adj_cat)), stat="identity") +
  scale_alpha_discrete(name = expression("BH p-value"), range=c(0.10, 0.40, 0.70, 1)) +
  ylim(-100,200) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  theme(
    legend.position = c(0.95, 0.5),
    legend.text=element_text(size=6),
    legend.title = element_text(size = 8),) +
  scale_fill_manual(values = c("#FF0000", "#FF6600", "#FFCC00", "#CCFF00", "#66FF00", "#00FF66", "#00FFCC", "#00CCFF", "#0066FF", "#0000FF", "#CC00FF", "#FF00CC", "#FF0066")) +
  guides( fill = "none") +
  #geom_text(data=label_data, aes(x=id, y=Percentage_difference_aged_vs_young_male+10, label=Lipid_ID, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=1.8, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=Lipid_subclass), hjust=c(0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.7, 0.55, 0.35, 0.55, 0.4, 0.55), colour = "black", alpha=0.8, size=3, fontface="bold", inherit.aes = FALSE)
p

ggsave(paste(graphics_path, paste("circular_bar_plot_aged_vs_young_male.png", sep = ""), sep = "/"), p, width = 8, height = 7, units = "in", dpi = 300)

#### SEX DIFFERENCES

# Create dataset
dat_cbp_3 <- p_val_frame[, names(p_val_frame)[c(1, 25, 26, 29, 30)]]
dat_cbp_4 <- dat_header[, names(dat_header)[c(1, 2)]]
dat_cbp_sex <- merge(dat_cbp_3, dat_cbp_4, by = "Lipid_species")

dat_cbp_sex$color <- factor(dat_cbp$Lipid_subclass, levels=c("CE", "Cer", "DG", "GSL", "LPC", "LPE", "LPE-O", "LPI", "PC", "PC-O", "PE", "PE-O", "PE-P", "PI", "SM", "TG"), labels=rainbow(16))

### YOUNG
# create subset for young with p< 0.05 only
dat_cbp_y <- filter(dat_cbp_sex, p_val_sex_y_adj_cat == "≤ 0.05"|p_val_sex_y_adj_cat =="≤ 0.01"|p_val_sex_y_adj_cat =="≤ 0.001"|p_val_sex_y_adj_cat =="≤ 0.0001")

#Define Lipid_subclass as factor
dat_cbp_y$Lipid_subclass <- factor(dat_cbp_y$Lipid_subclass)
dat_cbp_y$p_val_sex_y_adj_cat <- factor(dat_cbp_y$p_val_sex_y_adj_cat)
levels(dat_cbp_y$p_val_sex_y_adj_cat) <- c("≤ 0.05", "≤ 0.01", "≤ 0.001", "≤ 0.0001")

# Order data:
dat_cbp_y = dat_cbp_y %>% arrange(Lipid_subclass, Percentage_difference_female_vs_male_young)

# Set a number of 'empty bar' to add at the end of each Lipid_subclass
empty_bar <- 4
to_add <- data.frame(matrix(NA, empty_bar*nlevels(dat_cbp_y$Lipid_subclass), ncol(dat_cbp_y)) )
colnames(to_add) <- colnames(dat_cbp_y)
to_add$Lipid_subclass <- rep(levels(dat_cbp_y$Lipid_subclass), each=empty_bar)
dat_cbp_y <- rbind(dat_cbp_y, to_add)
dat_cbp_y<- dat_cbp_y %>% arrange(Lipid_subclass)
dat_cbp_y$id <- seq(1, nrow(dat_cbp_y))

# Get the name and the y position of each label
label_data <- dat_cbp_y
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- dat_cbp_y %>% 
  group_by(Lipid_subclass) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(dat_cbp_y, aes(x=as.factor(id), y=Percentage_difference_female_vs_male_young, fill=color, alpha = factor(p_val_sex_y_adj_cat))) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=Percentage_difference_female_vs_male_young, fill=color, alpha = factor(p_val_sex_y_adj_cat)), stat="identity") +
  
  # Add a val=60/40/20/0 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 25, xend = start, yend = 25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = -25, xend = start, yend = -25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = -50, xend = start, yend = -50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(dat_cbp_y$id),6), y = c(-50, -25, 0, 25, 50, 100), label = c("-50%", "-25%", "0%", "25%", "50%", "100%") , color="grey", size=2 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=Percentage_difference_female_vs_male_young, fill=color, alpha = factor(p_val_sex_y_adj_cat)), stat="identity") +
  scale_alpha_discrete(name = expression("BH p-value"), range=c(0.10, 0.40, 0.70, 1)) +
  ylim(-80,100) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  theme(
    legend.position = "right",
    legend.text=element_text(size=6),
    legend.title = element_text(size = 8),) +
  scale_fill_manual(values = c("#FF6D00", "#FFDB00", "#B6FF00", "#49FF00", "#00FF24", "#00FF92", "#00FFFF", "#0092FF", "#4900FF", "#B600FF", "#FF00DB")) +
  guides(fill = "none") +
  geom_text(data=label_data, aes(x=id, y=28, label=Lipid_species, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=1.8, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = 0, xend = end, yend = 0), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = 10, label=Lipid_subclass), hjust=c(0.55, 0.55, 0.7, 0.55, 0.55, 0.7, 0.55, 0.55, 0.55, 0.3, 0.55), colour = "black", alpha=0.8, size=3, fontface="bold", inherit.aes = FALSE)
p

ggsave(paste(graphics_path, paste("circular_bar_plot_female_vs_male_young.png", sep = ""), sep = "/"), p, width = 9, height = 7, units = "in", dpi = 300)

### AGED
# create subset for aged with p< 0.05 only
dat_cbp_o <- filter(dat_cbp_sex, p_val_sex_o_adj_cat == "≤ 0.05"|p_val_sex_o_adj_cat =="≤ 0.01"|p_val_sex_o_adj_cat =="≤ 0.001"|p_val_sex_o_adj_cat =="≤ 0.0001")

#Define Lipid_subclass as factor
dat_cbp_o$Lipid_subclass <- factor(dat_cbp_o$Lipid_subclass)
dat_cbp_o$p_val_sex_o_adj_cat <- factor(dat_cbp_o$p_val_sex_o_adj_cat)
levels(dat_cbp_o$p_val_sex_o_adj_cat) <- c("≤ 0.05", "≤ 0.01", "≤ 0.001", "≤ 0.0001")

# Order data:
dat_cbp_o = dat_cbp_o %>% arrange(Lipid_subclass, Percentage_difference_female_vs_male_aged)

# Set a number of 'empty bar' to add at the end of each Lipid_subclass
empty_bar <- 4
to_add <- data.frame(matrix(NA, empty_bar*nlevels(dat_cbp_o$Lipid_subclass), ncol(dat_cbp_o)) )
colnames(to_add) <- colnames(dat_cbp_o)
to_add$Lipid_subclass <- rep(levels(dat_cbp_o$Lipid_subclass), each=empty_bar)
dat_cbp_o <- rbind(dat_cbp_o, to_add)
dat_cbp_o <- dat_cbp_o %>% arrange(Lipid_subclass)
dat_cbp_o$id <- seq(1, nrow(dat_cbp_o))

# Get the name and the y position of each label
label_data <- dat_cbp_o
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- dat_cbp_o %>% 
  group_by(Lipid_subclass) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(dat_cbp_o, aes(x=as.factor(id), y=Percentage_difference_female_vs_male_aged, fill=color, alpha = factor(p_val_sex_o_adj_cat))) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=Percentage_difference_female_vs_male_aged, fill=color, alpha = factor(p_val_sex_o_adj_cat)), stat="identity") +
  
  # Add a val=60/40/20/0 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 75, xend = start, yend = 75), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 25, xend = start, yend = 25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = -25, xend = start, yend = -25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = -50, xend = start, yend = -50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(dat_cbp_o$id),6), y = c(-50, -25, 0, 25, 50, 75), label = c("-50%", "-25%", "0%", "25%", "50%", "75%") , color="grey", size=2 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=Percentage_difference_female_vs_male_aged, fill=color, alpha = factor(p_val_sex_o_adj_cat)), stat="identity") +
  scale_alpha_discrete(name = expression("BH p-value"), range=c(0.10, 0.40, 0.70, 1)) +
  ylim(-50,100) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  theme(
    legend.position = "right",
    legend.text=element_text(size=6),
    legend.title = element_text(size = 8),) +
  scale_fill_manual(values = c("#FF0000", "#FF6D00", "#49FF00", "#00FF24", "#00FFFF", "#B600FF", "#FF00DB", "#FF0066")) +
  guides( fill = "none") +
  geom_text(data=label_data, aes(x=id, y=22, label=Lipid_species, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=1.8, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = 0, xend = end, yend = 0), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = 10, label=Lipid_subclass), hjust=c(0.55, 0.55, 0.55, 0.55, 0.55, 0.7, 0.55, 0.55), colour = "black", alpha=0.8, size=3, fontface="bold", inherit.aes = FALSE)
p

ggsave(paste(graphics_path, paste("circular_bar_plot_female_vs_male_aged.png", sep = ""), sep = "/"), p, width = 9, height = 7, units = "in", dpi = 300)
