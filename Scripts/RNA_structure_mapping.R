#secondary structure map

library(readxl)
library(dplyr)
library(tidyr)


mappings <- read_excel("Tree/wrefs/recombination_5.66/Secondary_structs/structure_mapping/mappings.xlsx")

mappings_expanded <- mappings %>%
  separate(STK, into = c("STK_start", "STK_end"), sep = ";", convert = TRUE)

hcv_1_mapping <- read_csv("Tree/wrefs/recombination_5.66/Secondary_structs/structure_mapping/hcv_1_secondary_structure_mapping.csv")

hcv_2_mapping <- read_csv("Tree/wrefs/recombination_5.66/Secondary_structs/structure_mapping/hcv_2_secondary_structure_mapping.csv")

hcv_3_mapping <- read_csv("Tree/wrefs/recombination_5.66/Secondary_structs/structure_mapping/hcv_3_secondary_structure_mapping.csv")

hcv_4_mapping <- read_csv("Tree/wrefs/recombination_5.66/Secondary_structs/structure_mapping/hcv_4_secondary_structure_mapping.csv")

hcv_5_mapping <- read_csv("Tree/wrefs/recombination_5.66/Secondary_structs/structure_mapping/hcv_5_secondary_structure_mapping.csv")

hcv_6_mapping <- read_csv("Tree/wrefs/recombination_5.66/Secondary_structs/structure_mapping/hcv_6_secondary_structure_mapping.csv")

# Convert both to character
mappings_expanded$STK_start <- as.character(mappings_expanded$STK_start)
mappings_expanded$STK_end   <- as.character(mappings_expanded$STK_end)

hcv_1_mapping$STK_Position <- as.character(hcv_1_mapping$STK_Position)
hcv_2_mapping$STK_Position <- as.character(hcv_2_mapping$STK_Position)
hcv_3_mapping$STK_Position <- as.character(hcv_3_mapping$STK_Position)
hcv_4_mapping$STK_Position <- as.character(hcv_4_mapping$STK_Position)
hcv_5_mapping$STK_Position <- as.character(hcv_5_mapping$STK_Position)
hcv_6_mapping$STK_Position <- as.character(hcv_6_mapping$STK_Position)


library(purrr)

# Store all 6 mapping data frames in a named list
hcv_mappings <- list(
  `HCV-1` = hcv_1_mapping,
  `HCV-2` = hcv_2_mapping,
  `HCV-3` = hcv_3_mapping,
  `HCV-4` = hcv_4_mapping,
  `HCV-5` = hcv_5_mapping,
  `HCV-6` = hcv_6_mapping
)

# Create all columns for HCV-1 to HCV-6
mapped_columns <- map_dfc(names(hcv_mappings), function(hcv_name) {
  mapping_df <- hcv_mappings[[hcv_name]]
  
  start_pos <- mapping_df$RDP_Position[match(mappings_expanded$STK_start, mapping_df$STK_Position)]
  end_pos <- mapping_df$RDP_Position[match(mappings_expanded$STK_end, mapping_df$STK_Position)]
  
  tibble(!!hcv_name := ifelse(!is.na(start_pos) & !is.na(end_pos),
                              paste(start_pos, end_pos, sep = ";"),
                              NA))
})

# Combine with original data
mappings_mapped <- bind_cols(mappings_expanded, mapped_columns)

write.xlsx(mappings_mapped, "Tree/wrefs/recombination_5.66/Secondary_structs/structure_mapping/HCV_structure_STK_map.xlsx")