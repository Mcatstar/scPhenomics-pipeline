# Data cleanup

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

# Upload data------------------------------------------------------------------

data_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [TGIMB] Natalie Porat-Shliom/Lauryn Brown/results/PLIN5 variants in ND_02/table_data/tabulated_data.csv'

data_path2 <- '/Volumes/LECIMAGE/Analysis/[NCI] [TGIMB] Natalie Porat-Shliom/Lauryn Brown/results/PLIN5 variants in ND_02/table_data/cell_data.csv'

data_path3 <- '/Volumes/LECIMAGE/Analysis/[NCI] [TGIMB] Natalie Porat-Shliom/Lauryn Brown/results/PLIN5 variants in ND_02/table_data/tabulated_overlap.csv'

output_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [TGIMB] Natalie Porat-Shliom/Lauryn Brown/results/PLIN5 variants in ND_02/table_data'

df <- read_csv(data_path)
df2 <- read_csv(data_path2)
df3 <- read_csv(data_path3)

df$condition <- factor(df$condition, levels = c(
  'Null',
  'PLIN5 (WT)',
  'PLIN5 (1-424)',
  'PLIN5A',
  'PLIN5A (1-424)',
  'PLIN5E',
  'PLIN5E (1-424)'
))
unique(df$condition)

df2$condition <- factor(df2$condition, levels = c(
  'Null',
  'PLIN5 (WT)',
  'PLIN5 (1-424)',
  'PLIN5A',
  'PLIN5A (1-424)',
  'PLIN5E',
  'PLIN5E (1-424)'
))
unique(df2$condition)

df3$condition <- factor(df3$condition, levels = c(
  'Null',
  'PLIN5 (WT)',
  'PLIN5 (1-424)',
  'PLIN5A',
  'PLIN5A (1-424)',
  'PLIN5E',
  'PLIN5E (1-424)'
))
unique(df3$condition)

# Count------------------------------------------------------------------------

df_tmp <- df %>% 
  select(image, condition, organelle, ring, count) %>% 
  arrange(organelle, condition, image) %>% 
  pivot_wider(names_from = ring,
              values_from = count)

df_tmp <- df_tmp %>% 
  replace_na(list(R1 = 0, R2 = 0, R3 = 0, R4 = 0,
                  R5 = 0, R6 = 0, R7 = 0, R8 = 0,
                  R9 = 0, R10 = 0, R11 = 0 , R12 = 0))

csvname <- 'count.csv'
csvpath <- file.path(output_path, csvname)

write_csv(df_tmp, csvpath)

# Total area-------------------------------------------------------------------

df_tmp <- df %>% 
  select(image, condition, organelle, ring, total_area) %>% 
  arrange(organelle, condition, image) %>% 
  pivot_wider(names_from = ring,
              values_from = total_area)

df_tmp <- df_tmp %>% 
  replace_na(list(R1 = 0, R2 = 0, R3 = 0, R4 = 0,
                  R5 = 0, R6 = 0, R7 = 0, R8 = 0,
                  R9 = 0, R10 = 0, R11 = 0 , R12 = 0))

csvname <- 'total_area.csv'
csvpath <- file.path(output_path, csvname)

write_csv(df_tmp, csvpath)

# Average area-----------------------------------------------------------------

df_tmp <- df %>% 
  select(image, condition, organelle, ring, average_area) %>% 
  arrange(organelle, condition, image) %>% 
  pivot_wider(names_from = ring,
              values_from = average_area)

# df_tmp <- df_tmp %>% 
#   replace_na(list(R1 = 0, R2 = 0, R3 = 0, R4 = 0,
#                   R5 = 0, R6 = 0, R7 = 0, R8 = 0,
#                   R9 = 0, R10 = 0, R11 = 0 , R12 = 0))

csvname <- 'average_area.csv'
csvpath <- file.path(output_path, csvname)

write_csv(df_tmp, csvpath)

# Average perimeter------------------------------------------------------------

df_tmp <- df %>% 
  select(image, condition, organelle, ring, average_perimeter) %>% 
  arrange(organelle, condition, image) %>% 
  pivot_wider(names_from = ring,
              values_from = average_perimeter)

# df_tmp <- df_tmp %>% 
#   replace_na(list(R1 = 0, R2 = 0, R3 = 0, R4 = 0,
#                   R5 = 0, R6 = 0, R7 = 0, R8 = 0,
#                   R9 = 0, R10 = 0, R11 = 0 , R12 = 0))

csvname <- 'average_perimeter.csv'
csvpath <- file.path(output_path, csvname)

write_csv(df_tmp, csvpath)

# Average eccentricity---------------------------------------------------------

df_tmp <- df %>% 
  select(image, condition, organelle, ring, average_eccentricity) %>% 
  arrange(organelle, condition, image) %>% 
  pivot_wider(names_from = ring,
              values_from = average_eccentricity)

# df_tmp <- df_tmp %>% 
#   replace_na(list(R1 = 0, R2 = 0, R3 = 0, R4 = 0,
#                   R5 = 0, R6 = 0, R7 = 0, R8 = 0,
#                   R9 = 0, R10 = 0, R11 = 0 , R12 = 0))

csvname <- 'average_eccentricity.csv'
csvpath <- file.path(output_path, csvname)

write_csv(df_tmp, csvpath)

# Average solidity-------------------------------------------------------------

df_tmp <- df %>% 
  select(image, condition, organelle, ring, average_solidity) %>% 
  arrange(organelle, condition, image) %>% 
  pivot_wider(names_from = ring,
              values_from = average_solidity)

# df_tmp <- df_tmp %>% 
#   replace_na(list(R1 = 0, R2 = 0, R3 = 0, R4 = 0,
#                   R5 = 0, R6 = 0, R7 = 0, R8 = 0,
#                   R9 = 0, R10 = 0, R11 = 0 , R12 = 0))

csvname <- 'average_solidity.csv'
csvpath <- file.path(output_path, csvname)

write_csv(df_tmp, csvpath)

# Average Feret's--------------------------------------------------------------

df_tmp <- df %>% 
  select(image, condition, organelle, ring, average_feret) %>% 
  arrange(organelle, condition, image) %>% 
  pivot_wider(names_from = ring,
              values_from = average_feret)

# df_tmp <- df_tmp %>% 
#   replace_na(list(R1 = 0, R2 = 0, R3 = 0, R4 = 0,
#                   R5 = 0, R6 = 0, R7 = 0, R8 = 0,
#                   R9 = 0, R10 = 0, R11 = 0 , R12 = 0))

csvname <- 'average_feret.csv'
csvpath <- file.path(output_path, csvname)

write_csv(df_tmp, csvpath)

# Average mito int-------------------------------------------------------------

df_tmp <- df %>% 
  select(image, condition, organelle, ring, average_mito_int) %>% 
  arrange(organelle, condition, image) %>% 
  pivot_wider(names_from = ring,
              values_from = average_mito_int)

csvname <- 'average_mito_int.csv'
csvpath <- file.path(output_path, csvname)

write_csv(df_tmp, csvpath)

# Average lipid int------------------------------------------------------------

df_tmp <- df %>% 
  select(image, condition, organelle, ring, average_lipid_int) %>% 
  arrange(organelle, condition, image) %>% 
  pivot_wider(names_from = ring,
              values_from = average_lipid_int)

csvname <- 'average_lipid_int.csv'
csvpath <- file.path(output_path, csvname)

write_csv(df_tmp, csvpath)

# Average organelle min dist---------------------------------------------------

df_tmp <- df %>% 
  select(image, condition, organelle, ring, 
         average_organelle_min_dist) %>% 
  arrange(organelle, condition, image) %>% 
  pivot_wider(names_from = ring,
              values_from = average_organelle_min_dist)

csvname <- 'average_organelle_min_dist.csv'
csvpath <- file.path(output_path, csvname)

write_csv(df_tmp, csvpath)

# Average circularity----------------------------------------------------------

df_tmp <- df %>% 
  select(image, condition, organelle, ring, average_circularity) %>% 
  arrange(organelle, condition, image) %>% 
  pivot_wider(names_from = ring,
              values_from = average_circularity)

csvname <- 'average_circularity.csv'
csvpath <- file.path(output_path, csvname)

write_csv(df_tmp, csvpath)

# Average density--------------------------------------------------------------

df_tmp <- df2 %>% 
  select(image, condition, organelle, ring, average_density) %>% 
  arrange(organelle, condition, image) %>% 
  pivot_wider(names_from = ring,
              values_from = average_density)

df_tmp <- df_tmp %>% 
  replace_na(list(R1 = 0, R2 = 0, R3 = 0, R4 = 0,
                  R5 = 0, R6 = 0, R7 = 0, R8 = 0,
                  R9 = 0, R10 = 0, R11 = 0 , R12 = 0))


csvname <- 'average_density.csv'
csvpath <- file.path(output_path, csvname)

write_csv(df_tmp, csvpath)

# Cell area--------------------------------------------------------------------

df_tmp <- df2 %>% 
  select(image, condition, organelle, ring, cell_area) %>% 
  arrange(organelle, condition, image) %>% 
  pivot_wider(names_from = ring, 
              values_from = cell_area)

csvname <- 'cell_area.csv'
csvpath <- file.path(output_path, csvname)

write_csv(df_tmp, csvpath)

# Overlap percentage-----------------------------------------------------------

df_tmp <- df3 %>% 
  select(image, condition, ring, percent_overlap) %>% 
  arrange(condition, image) %>% 
  pivot_wider(names_from = ring,
              values_from = percent_overlap)

df_tmp <- df_tmp %>% 
  replace_na(list(R1 = 0, R2 = 0, R3 = 0, R4 = 0,
                  R5 = 0, R6 = 0, R7 = 0, R8 = 0,
                  R9 = 0, R10 = 0, R11 = 0 , R12 = 0))


csvname <- 'percent_overlap.csv'
csvpath <- file.path(output_path, csvname)

write_csv(df_tmp, csvpath)
