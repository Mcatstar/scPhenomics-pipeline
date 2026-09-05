# Lobule analysis for Natalie Porat-Shliom lab
# 07/05/24

# Author: Andy D. Tran, CCR Microscopy Core, LCBG, CCR, NCI

# Libraries and themes---------------------------------------------------------

library(tidyverse)
library(ggsci)

theme <- theme(panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               panel.background=element_blank(),
               axis.line=element_line(color="black", linewidth=1),
               axis.ticks=element_line(color="black", linewidth=1),
               text=element_text(size=18),
               plot.title=element_text(size=24),
               axis.text=element_text(color="black"),
               plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

# Download data----------------------------------------------------------------

input_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [TGIMB] Natalie Porat-Shliom/Lauryn Brown/output/PLIN5 variants in ND'

df <- tibble(condition = character(), image = character())

input_list <- list.files(input_path)

for(condition in input_list){
  condition_path <- file.path(input_path, condition)
  img_list <- list.files(condition_path, pattern = '_output.csv')
  for(img in img_list){
    img_path <- file.path(condition_path, img)
    img_df <- read_csv(img_path) %>% 
      mutate(condition = condition) %>% 
      mutate(image = str_replace(img, '_output.csv', '')) %>% 
      select(!'...1')
    df <- full_join(df, img_df)
  }
}

unique(df$condition)
df$condition <- factor(df$condition, levels = c(
  'Null', 
  'PLIN5 (WT)',
  'PLIN5 (1-424)',
  'PLIN5A',
  'PLIN5A (1-424)',
  'PLIN5E',
  'PLIN5E (1-424)'
))

# Data clearning---------------------------------------------------------------
output_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [TGIMB] Natalie Porat-Shliom/Lauryn Brown/results/PLIN5 variants in ND_02'

cell_df <- df %>% 
  filter(cell_id > 0) %>% 
  filter(!is.na(cell_cv_mean_dist)) %>% 
  mutate(relative_cell_cv_dist = cell_cv_max_dist / pv_cv_dist) %>% 
  mutate(cell_width = cell_cv_max_dist - cell_cv_min_dist) %>% 
  group_by(condition, image, cell_id) %>% 
  summarise(relative_cell_cv_dist = mean(relative_cell_cv_dist),
            cell_pv_min_dist = mean(cell_pv_min_dist),
            cell_cv_mean_dist = mean(cell_cv_mean_dist),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            cell_cv_min_dist = mean(cell_cv_min_dist), 
            cell_width = mean(cell_width),
            cell_area = mean(cell_area),
            overlap_area = mean(overlap_area)) %>% 
  mutate(relative_cell_cv_dist = abs(1 - relative_cell_cv_dist)) %>% 
  replace_na(list(overlap_area = 0)) %>% 
  mutate(percent_overlap = overlap_area / cell_area)

cell_df2 <- cell_df %>% 
  group_by(condition) %>% 
  summarise(cell_width = mean(cell_width))


ggplot(cell_df, aes(x = relative_cell_cv_dist)) +
  geom_freqpoly(aes(y = after_stat(ndensity)), bins = 100) +
  geom_vline(xintercept = c(1/12, 2/12, 3/12, 4/12, 5/12, 6/12,
                            7/12, 8/12, 9/12, 10/12, 11/12),
             col = 'red', linetype = 2) +
  theme +
  facet_grid(rows = vars(condition)) +
  labs(title = 'Cell Max CV distance normalized')

plotname <- 'relative_cell_cv_max_distribution.png'
ggsave(plotname, plot = last_plot(), path = output_path,  height = 11,
       width = 16, dpi = 150, units = c('in'))

tmp_plot <- cell_df %>% 
  filter(condition == 'Null')

ggplot(tmp_plot, aes(x = relative_cell_cv_dist)) +
  geom_freqpoly(aes(y = after_stat(ndensity)), bins = 100) +
  geom_vline(xintercept = c(1/12, 2/12, 3/12, 4/12, 5/12, 6/12,
                            7/12, 8/12, 9/12, 10/12, 11/12),
             col = 'red', linetype = 2) +
  facet_grid(rows = vars(condition)) +
  theme +
  labs(title = 'Cell Max CV distance normalized')

plotname <- 'relative_cell_cv_max_distribution_null.png'
ggsave(plotname, plot = last_plot(), path = output_path,  height = 8,
       width = 16, dpi = 150, units = c('in'))

# tmp_plot <- cell_df %>% 
#   filter(condition == 'Fasting')
# 
# ggplot(tmp_plot, aes(x = relative_cell_cv_dist)) +
#   geom_freqpoly(aes(y = after_stat(ndensity)), bins = 100) +
#   geom_vline(xintercept = c(1/12, 2/12, 3/12, 4/12, 5/12, 6/12,
#                             7/12, 8/12, 9/12, 10/12, 11/12),
#              col = 'red', linetype = 2) +
#   facet_grid(rows = vars(condition)) +
#   theme +
#   labs(title = 'Cell Max CV distance normalized')
# 
# plotname <- 'relative_cell_cv_max_distribution_fasting.png'
# ggsave(plotname, plot = last_plot(), path = output_path,  height = 8,
#        width = 16, dpi = 150, units = c('in'))
# 
# tmp_plot <- cell_df %>% 
#   filter(condition == 'WD')
# 
# ggplot(tmp_plot, aes(x = relative_cell_cv_dist)) +
#   geom_freqpoly(aes(y = after_stat(ndensity)), bins = 100) +
#   geom_vline(xintercept = c(1/12, 2/12, 3/12, 4/12, 5/12, 6/12,
#                             7/12, 8/12, 9/12, 10/12, 11/12),
#              col = 'red', linetype = 2) +
#   facet_grid(rows = vars(condition)) +
#   theme +
#   labs(title = 'Cell Max CV distance normalized')
# 
# plotname <- 'relative_cell_cv_max_distribution_wd.png'
# ggsave(plotname, plot = last_plot(), path = output_path,  height = 8,
#        width = 16, dpi = 150, units = c('in'))

pixel_um <- 0.0694233

df2 <- df %>% 
  filter(cell_id > 0) %>% 
  filter(!is.na(cell_cv_mean_dist)) %>% 
  mutate(relative_cell_cv_dist = cell_cv_max_dist / pv_cv_dist) %>% 
  mutate(relative_cell_cv_dist = abs(1 - relative_cell_cv_dist)) %>% 
  mutate(ring = case_when(
    cell_pv_min_dist == 0 ~ 'R1',
    relative_cell_cv_dist > 11/12 ~ 'R12',
    relative_cell_cv_dist > 10/12 ~ 'R11',
    relative_cell_cv_dist > 9/12 ~ 'R10',
    relative_cell_cv_dist > 8/12 ~ 'R9',
    relative_cell_cv_dist > 7/12 ~ 'R8',
    relative_cell_cv_dist > 6/12 ~ 'R7',
    relative_cell_cv_dist > 5/12 ~ 'R6',
    relative_cell_cv_dist > 4/12 ~ 'R5',
    relative_cell_cv_dist > 3/12 ~ 'R4', 
    relative_cell_cv_dist > 2/12 ~ 'R3',
    relative_cell_cv_dist > 1/12 ~ 'R2',
    TRUE ~ 'R1'
  )) %>% 
  mutate(ring = case_when(
    image == '3920E M2_L1_Cropped PP-PC Axis 4' & cell_id == 86 ~ 'R2',
    image == '3920E M2_L1_Cropped PP-PC Axis 4' & cell_id == 61 ~ 'R2',
    image == '3920E M2_L1_Cropped PP-PC Axis 4' & cell_id == 48 ~ 'R2',
    image == '3920E M2_L1_Cropped PP-PC Axis 4' & cell_id == 46 ~ 'R2',
    image == '3920E M2_L1_Cropped PP-PC Axis 4' & cell_id == 32 ~ 'R2',
    TRUE ~ ring
  )) %>% 
  dplyr::rowwise() %>% 
  mutate(area = area * pixel_um * pixel_um) %>% 
  mutate(perimeter = perimeter * pixel_um) %>% 
  mutate(axis_major_length = axis_major_length *
           pixel_um) %>% 
  mutate(axis_minor_length = axis_minor_length * 
           pixel_um) %>% 
  mutate(feret_diameter_max = feret_diameter_max * 
           pixel_um) %>% 
  mutate(cell_area = cell_area * 
           pixel_um * pixel_um) %>% 
  mutate(organelle_min_dist = organelle_min_dist * 
           pixel_um) %>% 
  mutate(organelle_max_dist = organelle_max_dist *
           pixel_um) %>% 
  mutate(organelle_mean_dist = organelle_mean_dist *
           pixel_um) %>% 
  mutate(overlap_area = overlap_area *
           pixel_um * pixel_um) %>% 
  mutate(circularity = (4*pi*area) / (perimeter^2)) %>% 
  mutate(circularity = ifelse(circularity > 1, 1, circularity))

df2$ring <- factor(df2$ring, levels = c('R1', 'R2', 'R3', 'R4',
                                        'R5', 'R6', 'R7', 'R8',
                                        'R9', 'R10', 'R11', 'R12'))
unique(df2$ring)

cell_df <- cell_df %>% 
  mutate(ring = case_when(
    cell_pv_min_dist == 0 ~ 'R1',
    relative_cell_cv_dist > 11/12 ~ 'R12',
    relative_cell_cv_dist > 10/12 ~ 'R11',
    relative_cell_cv_dist > 9/12 ~ 'R10',
    relative_cell_cv_dist > 8/12 ~ 'R9',
    relative_cell_cv_dist > 7/12 ~ 'R8',
    relative_cell_cv_dist > 6/12 ~ 'R7',
    relative_cell_cv_dist > 5/12 ~ 'R6',
    relative_cell_cv_dist > 4/12 ~ 'R5',
    relative_cell_cv_dist > 3/12 ~ 'R4', 
    relative_cell_cv_dist > 2/12 ~ 'R3',
    relative_cell_cv_dist > 1/12 ~ 'R2',
    TRUE ~ 'R1'
  ))

cell_df <- cell_df %>% 
  mutate(ring = case_when(
    image == '3920E M2_L1_Cropped PP-PC Axis 4' & cell_id == 86 ~ 'R2',
    image == '3920E M2_L1_Cropped PP-PC Axis 4' & cell_id == 61 ~ 'R2',
    image == '3920E M2_L1_Cropped PP-PC Axis 4' & cell_id == 48 ~ 'R2',
    image == '3920E M2_L1_Cropped PP-PC Axis 4' & cell_id == 46 ~ 'R2',
    image == '3920E M2_L1_Cropped PP-PC Axis 4' & cell_id == 32 ~ 'R2',
    TRUE ~ ring
  ))

cell_df$ring <- factor(cell_df$ring, levels = c('R1', 'R2', 'R3', 'R4',
                                                'R5', 'R6', 'R7', 'R8',
                                                'R9', 'R10', 'R11', 'R12'))

unique(cell_df$ring)

test_df <- cell_df %>% 
  group_by(image, condition, ring) %>% 
  summarise()

test_df2 <- test_df %>% 
  filter(grepl('3920E', image))

# Organelle count per cell-----------------------------------------------------
# Mitochondria-----------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'mito') %>% 
  replace_na(list(roi_id = 0)) %>% 
  mutate(count = ifelse(roi_id > 0, 1, 0)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(count = sum(count), cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = count, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  labs(title = 'Mitochondria counts per cell',
       y = 'Count',
       x = 'Region')

plotname <- 'mito_count_region.png'
ggsave(plotname, plot = last_plot(), path = output_path,  height = 8,
       width = 16, dpi = 150, units = c('in'))

# Lipid droplets---------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'lipid') %>% 
  replace_na(list(roi_id = 0)) %>% 
  mutate(count = ifelse(roi_id > 0, 1, 0)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(count = sum(count), cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = count, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  labs(title = 'LD counts per cell',
       y = 'Count',
       x = 'Region')

plotname <- 'lipid_count_region.png'
ggsave(plotname, plot = last_plot(), path = output_path,  height = 8,
       width = 16, dpi = 150, units = c('in'))

# Organelle area---------------------------------------------------------------
# Mitochondria-----------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'mito') %>% 
  filter(!is.na(roi_id)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(area = mean(area), cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = area, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  ylim(0, 3) +
  labs(title = 'Mean mitochondrion area per cell',
       y = 'Area [um^2]',
       x = 'Region')

plotname <- 'mito_area_region.png'
ggsave(plotname, plot = last_plot(), path = output_path,  height = 8,
       width = 16, dpi = 150, units = c('in'))

# Lipid droplets---------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'lipid') %>% 
  filter(!is.na(roi_id)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(area = mean(area), cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = area, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  ylim(0, 10) +
  labs(title = 'Mean LD area per cell',
       y = 'Area [um^2]',
       x = 'Region')

plotname <- 'lipid_area_region.png'
ggsave(plotname, plot = last_plot(), path = output_path,  height = 8,
       width = 16, dpi = 150, units = c('in'))

# Total organelle area---------------------------------------------------------
# Mitochondria-----------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'mito') %>% 
  replace_na(list(area = 0)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(total_area = sum(area), cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = total_area, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  labs(title = 'Total mitochondria area per cell',
       y = 'Area [um^2]',
       x = 'Region')

plotname <- 'mito_total_area_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Lipid droplets---------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'lipid') %>% 
  replace_na(list(area = 0)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(total_area = sum(area), cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = total_area, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  labs(title = 'Total LD area per cell',
       y = 'Area [um^2]',
       x = 'Region')

plotname <- 'lipid_total_area_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Organelle density------------------------------------------------------------
# Mitochondria-----------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'mito') %>% 
  replace_na(list(area = 0)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(total_area = sum(area), cell_area = mean(cell_area),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist)) %>% 
  mutate(density = total_area / cell_area) %>% 
  filter(density <= 1)

ggplot(plot_df, aes(x = ring, y = density, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  #ylim(0, 0.01) +
  labs(title = 'Mitochondria density',
       y = 'Density [Organelle area / Cell area]',
       x = 'Region')

plotname <- 'mito_density_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Lipid droplets---------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'lipid') %>% 
  replace_na(list(area = 0)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(total_area = sum(area), cell_area = mean(cell_area),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist)) %>% 
  mutate(density = total_area / cell_area) %>% 
  filter(density <= 1)

ggplot(plot_df, aes(x = ring, y = density, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  #ylim(0, 0.05) +
  labs(title = 'LD density',
       y = 'Density [Organelle area / Cell area]',
       x = 'Region')

plotname <- 'lipid_density_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Total cellular fluorescence intensity----------------------------------------
# Mitochondria-----------------------------------------------------------------

plot_df <- df2 %>% 
  filter(cell_id > 0) %>% 
  filter(!is.na(cell_cv_mean_dist)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(cell_mito_int = mean(cell_mito_int), 
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = cell_mito_int, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  ylim(0, 150) +
  labs(title = 'Total cell mitochondria intensity',
       y = 'Intensity [AU]',
       x = 'Region')

plotname <- 'mito_total_cell_intensity_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Lipid droplets---------------------------------------------------------------

plot_df <- df2 %>% 
  filter(cell_id > 0) %>% 
  filter(!is.na(cell_cv_mean_dist)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(cell_lipid_int = mean(cell_lipid_int),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = cell_lipid_int, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  ylim(0, 40) +
  labs(title = 'Total cell LD intensity',
       y = 'Intensity [AU]',
       x = 'Region')

plotname <- 'lipid_total_cell_intensity_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Organelle intensity----------------------------------------------------------
# Mitochondria-----------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'mito') %>% 
  filter(!is.na(roi_id)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(mito_int = mean(mito_int),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = mito_int, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  ylim(0, 200) +
  labs(title = 'Mean mitochondria intensity per cell',
       y = 'Intensity [AU]',
       x = 'Region')


plotname <- 'mito_mean_intensity_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Lipid droplets---------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'lipid') %>% 
  filter(!is.na(roi_id)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(lipid_int = mean(lipid_int),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = lipid_int, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  ylim(0, 50) +
  labs(title = 'Mean LD intensity per cell',
       y = 'Intensity [AU]',
       x = 'Region')

plotname <- 'lipid_mean_intensity_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Eccentricity-----------------------------------------------------------------
# Mitochondria-----------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'mito') %>% 
  filter(!is.na(roi_id)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(eccentricity = mean(eccentricity),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = eccentricity, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  labs(title = 'Mean mitochondria eccentricity per cell',
       y = 'Eccentricity',
       x = 'Region')

plotname <- 'mito_eccentricity_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Lipid droplets---------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'lipid') %>% 
  filter(!is.na(roi_id)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(eccentricity = mean(eccentricity),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = eccentricity, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  labs(title = 'Mean LD eccentricity per cell',
       y = 'Eccentricity',
       x = 'Region')

plotname <- 'lipid_eccentricity_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Soliditiy--------------------------------------------------------------------
# Mitochondria-----------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'mito') %>% 
  filter(!is.na(roi_id)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(solidity = mean(solidity),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = solidity, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  ylim(0.8, 1) +
  labs(title = 'Mean mitochondria solidity per cell',
       y = 'Solidity',
       x = 'Region')

plotname <- 'mito_solidity_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Lipid droplets---------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'lipid') %>% 
  filter(!is.na(roi_id)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(solidity = mean(solidity),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = solidity, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  ylim(0.8, 1) +
  labs(title = 'Mean LD solidity per cell',
       y = 'Solidity',
       x = 'Region')

plotname <- 'lipid_solidity_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Feret diameter---------------------------------------------------------------
# Mitochondria-----------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'mito') %>% 
  filter(!is.na(roi_id)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(feret_diameter_max = mean(feret_diameter_max),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = feret_diameter_max, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  ylim(0, 3) +
  labs(title = 'Mean mitochondria Feret diameter per cell',
       y = 'Feret diameter [um]',
       x = 'Region')

plotname <- 'mito_feret_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Lipid droplets---------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'lipid') %>% 
  filter(!is.na(roi_id)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(feret_diameter_max = mean(feret_diameter_max),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = feret_diameter_max, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  ylim(0, 5) +
  labs(title = 'Mean LD Feret diameter per cell',
       y = 'Feret diameter [px]',
       x = 'Region')

plotname <- 'lipid_feret_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Circularity------------------------------------------------------------------
# Mitochondria-----------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'mito') %>% 
  filter(!is.na(roi_id)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(circularity = mean(circularity),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = circularity, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  #ylim(0, 5) +
  labs(title = 'Mean mitochondria circularity per cell',
       y = 'Circularity',
       x = 'Region')

plotname <- 'mito_circularity_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Lipid droplet----------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'lipid') %>% 
  filter(!is.na(roi_id)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(circularity = mean(circularity),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = circularity, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  #ylim(0, 5) +
  labs(title = 'Mean LD circularity per cell',
       y = 'Circularity',
       x = 'Region')

plotname <- 'lipid_circularity_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Distance to corresponding organelle------------------------------------------
# Mitochondria-----------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'mito') %>% 
  filter(!is.na(roi_id)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(organelle_min_dist = mean(organelle_min_dist),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = organelle_min_dist, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  ylim(0, 40) +
  labs(title = 'Mean minimum mitochondria to LD distance per cell',
       y = 'Distance [um]',
       x = 'Region')

plotname <- 'mito_min_dist_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Lipid droplets---------------------------------------------------------------

plot_df <- df2 %>% 
  filter(organelle == 'lipid') %>% 
  filter(!is.na(roi_id)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(organelle_min_dist = mean(organelle_min_dist),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = organelle_min_dist, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  ylim(0, 1) +
  labs(title = 'Mean minimum LD to mitochondria distance per cell',
       y = 'Distance [um]',
       x = 'Region')

plotname <- 'lipid_min_dist_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Overlap area-----------------------------------------------------------------

plot_df <- df2 %>% 
  filter(!is.na(roi_id)) %>% 
  filter(!is.na(overlap_area)) %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(overlap_area = mean(overlap_area),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = overlap_area, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T,
               outlier.size = -1) +
  theme +
  scale_fill_npg() +
  labs(title = 'Organelle overlap area per cell',
       y = 'Area [um^2]',
       x = 'Region')

plotname <- 'overlap_area_region.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Cell area--------------------------------------------------------------------

plot_df <- df2 %>% 
  group_by(condition, image, cell_id, ring) %>% 
  summarise(cell_area = mean(cell_area),
            cell_cv_max_dist = mean(cell_cv_max_dist),
            relative_cell_cv_dist = mean(relative_cell_cv_dist))

ggplot(plot_df, aes(x = ring, y = cell_area, fill = condition)) +
  geom_boxplot(alpha = 1, col = 'black', show.legend = T, 
               outlier.size = -1) +
  theme + 
  scale_fill_npg() +
  labs(title = 'Cell area',
       y = 'Area [um^2]',
       x = 'Region')

plotname <- 'cell_area.png'
ggsave(plotname, plot = last_plot(), path = output_path, height = 8,
       width = 16, dpi = 150, units = c('in'))

# Export spreadsheet-----------------------------------------------------------

csvname <- 'combined_data.csv'
csvpath <- file.path(output_path, csvname)

write_csv(df2, csvpath)

csvname <- 'cell_data.csv'
csvpath <- file.path(output_path, csvname)

write_csv(cell_df, csvpath)

# Format data to FIJI standard-------------------------------------------------

table_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [TGIMB] Natalie Porat-Shliom/Lauryn Brown/results/PLIN5 variants in ND_02/table_data'

df3 <- df2 %>% 
  mutate(count = ifelse(is.na(roi_id), 0, 1)) %>% 
  group_by(condition, image, organelle, ring, cell_id) %>% 
  summarise(count = sum(count), cell_area = mean(cell_area),
            total_area = sum(area), average_area = mean(area),
            average_perimeter = mean(perimeter), 
            average_eccentricity = mean(eccentricity),
            average_solidity = mean(solidity),
            average_feret = mean(feret_diameter_max),
            average_mito_int = mean(mito_int),
            average_lipid_int = mean(lipid_int),
            average_organelle_min_dist = mean(organelle_min_dist),
            cell_area = mean(cell_area),
            average_circularity = mean(circularity)) %>% 
  mutate(density = count / cell_area) %>% 
  filter(count > 0) %>% 
  group_by(condition, image, organelle, ring) %>% 
  summarise(count = sum(count), cell_area = mean(cell_area),
            total_area = sum(total_area), 
            average_area = mean(average_area),
            average_perimeter = mean(average_perimeter),
            average_eccentricity = mean(average_eccentricity),
            average_solidity = mean(average_solidity),
            average_feret = mean(average_feret),
            average_mito_int = mean(average_mito_int),
            average_lipid_int = mean(average_lipid_int),
            average_organelle_min_dist = mean(average_organelle_min_dist),
            average_circularity = mean(average_circularity),
            average_density = mean(density)) %>% 
  arrange(condition, organelle, image, ring)

csvname <- 'tabulated_data.csv'
csvpath <- file.path(table_path, csvname)

write_csv(df3, csvpath)

df3 <- df2 %>% 
  mutate(count = ifelse(is.na(roi_id), 0, 1)) %>% 
  group_by(condition, image, organelle, ring, cell_id) %>% 
  summarise(count = sum(count), cell_area = mean(cell_area),
            total_area = sum(area), average_area = mean(area),
            average_perimeter = mean(perimeter), 
            average_eccentricity = mean(eccentricity),
            average_solidity = mean(solidity),
            average_feret = mean(feret_diameter_max),
            average_mito_int = mean(mito_int),
            average_lipid_int = mean(lipid_int),
            average_organelle_min_dist = mean(organelle_min_dist),
            cell_area = mean(cell_area),
            average_circularity = mean(circularity)) %>% 
  mutate(density = count / cell_area)

ggplot(df3, aes(x = cell_area)) +
  geom_histogram(aes(y = after_stat(..density..)), 
                 fill = 'white', col = 'black') +
  scale_x_log10() +
  theme

df4 <- df3 %>% 
  filter(count > 0) %>% 
  filter(cell_area > 10) %>% 
  group_by(condition, image, organelle, ring) %>% 
  summarise(count = sum(count), cell_area = mean(cell_area),
            total_area = sum(total_area), 
            average_area = mean(average_area),
            average_perimeter = mean(average_perimeter),
            average_eccentricity = mean(average_eccentricity),
            average_solidity = mean(average_solidity),
            average_feret = mean(average_feret),
            average_mito_int = mean(average_mito_int),
            average_lipid_int = mean(average_lipid_int),
            average_organelle_min_dist = mean(average_organelle_min_dist),
            average_circularity = mean(average_circularity),
            average_density = mean(density)) %>% 
  arrange(condition, organelle, image, ring)

csvname <- 'tabulated_data.csv'
csvpath <- file.path(table_path, csvname)

write_csv(df4, csvpath)

df5 <- df3 %>% 
  filter(cell_area > 10) %>% 
  group_by(condition, image, organelle, ring) %>% 
  summarise(cell_area = mean(cell_area),
            average_density = mean(density)) %>% 
  arrange(condition, organelle, image, ring)

csvname <- 'cell_data.csv'
csvpath <- file.path(table_path, csvname)

write_csv(df5, csvpath)

df_test <- df2 %>% 
  mutate(count = ifelse(is.na(roi_id), 0, 1)) %>% 
  group_by(condition, mouse, lobule, image, organelle, ring, cell_id) %>% 
  summarise(count = sum(count), cell_area = mean(cell_area),
            total_area = sum(area), average_area = mean(area),
            average_perimeter = mean(perimeter), 
            average_eccentricity = mean(eccentricity),
            average_solidity = mean(solidity),
            average_feret = mean(feret_diameter_max),
            average_mito_int = mean(mito_int),
            average_lipid_int = mean(lipid_int),
            average_organelle_min_dist = mean(organelle_min_dist),
            cell_area = mean(cell_area),
            average_circularity = mean(circularity)) %>% 
  mutate(density = count / cell_area)

df_test2 <- df_test %>% 
  filter(condition == 'Fasting') %>% 
  filter(lobule == 'Lobule 2') %>% 
  filter(mouse == 'M3')

cell_df3 <- cell_df %>% 
  group_by(condition, image, ring) %>% 
  summarise(percent_overlap = mean(percent_overlap))

csvname <- 'tabulated_overlap.csv'
csvpath <- file.path(table_path, csvname)

write_csv(cell_df3, csvpath)
