# Author: Ali Farzam (GitHub: theAliFarzam)
# Correspondence: afarzam@purdue.edu | Dr. Bryon S. Drown (bsdrown@purdue.edu)

library(tidyverse)
library(UpSetR)
library(ComplexUpset)
library(grid)

# concatenate all the different search strategies and mapping/cross-mapping results
all_search_and_cross_map_peptides <- bind_rows(
  downstream_search_upstream_map = downstream_search_upstream_map_classified_peptides,
  upstream_search_downstream_map = upstream_search_downstream_map_classified_peptides,
  upstream_search = upstream_search_classified_peptides,
  downstream_search = downstream_search_classified_peptides,
  .id = "search_strategy"
)

# keep unique peptides and their classification, this ignores the search/mapping strategy
# so that peptides with the same classifications appear once, 
# and only those with different classifications appear more than once
unique_peptides_by_category <- all_search_and_cross_map_peptides |> 
  group_by(category, peptide) |> 
  distinct(peptide)

peptide_binary_matrix_category <- unique_peptides_by_category |> 
  pivot_longer(
    cols = category,
    names_to = "attribute_type",
    values_to = "set_name"
  ) |> 
  distinct(peptide, set_name) |> 
  mutate(present = 1) |> 
  pivot_wider(
    names_from = set_name,
    values_from = present,
    values_fill = list(present = 0)
  ) |> 
  as.data.frame()

category_sets <- setdiff(colnames(peptide_binary_matrix_category), "peptide")

# do as above, but include the search strategy
# so the upset plot shows where the source of diffrent classifications
unique_peptides_by_search_and_category <- all_search_and_cross_map_peptides |> 
  group_by(search_strategy, peptide) |> 
  distinct(peptide, .keep_all = TRUE) |> 
  select(search_strategy, peptide, category)

peptide_binary_matrix_search_and_category <- unique_peptides_by_search_and_category |> 
  pivot_longer(
    cols = c(search_strategy, category),
    names_to = "attribute_type",
    values_to = "set_name"
  ) |> 
  distinct(peptide, set_name) |> 
  mutate(present = 1) |> 
  pivot_wider(
    names_from = set_name,
    values_from = present,
    values_fill = list(present = 0)
  ) |> 
  as.data.frame()

all_sets <- setdiff(colnames(peptide_binary_matrix_search_and_category), "peptide")

# Generate the plot for intersection of categories
UpSetR::upset(
  peptide_binary_matrix_category,
  sets = category_sets,
  order.by = "freq",               # Sort intersections by largest size
  keep.order = FALSE,              # Lets UpSetR sort the sets by size on the left axis
  nsets = length(all_sets),        # Ensures all your generated sets are included
  main.bar.color = "#2c3e50",      
  sets.bar.color = "#e74c3c",      
  text.scale = 1.2                 
)

# generate the plot for search/mapping strategy and classification
UpSetR::upset(
  peptide_binary_matrix_search_and_category,
  sets = all_sets,
  order.by = "freq",               # Sort intersections by largest size
  keep.order = FALSE,              # Lets UpSetR sort the sets by size on the left axis
  nsets = length(all_sets),        # Ensures all your generated sets are included
  main.bar.color = "#2c3e50",      
  sets.bar.color = "#e74c3c",      
  text.scale = 1.2                 
)

# same plots but styling easier to read, maybe
# intersection between peptides and classifications
UpSetR::upset(
  peptide_binary_matrix_category,
  sets = category_sets,
  nsets = length(all_sets),        # Ensures all 12 sets are plotted on the y-axis
  nintersects = 40,                # Caps the plot at the top 40 largest intersections
  order.by = "freq",               # Sorts intersections by frequency (largest first)
  keep.order = FALSE,              # Sorts the 12 sets on the left by their total size
  mb.ratio = c(0.55, 0.45),        # Allocates 45% of the plot height to the matrix to fit 12 rows
  point.size = 2.5,                # Slightly smaller dots to fit the grid
  line.size = 0.8,                 # Slightly thinner lines
  main.bar.color = "#2c3e50",      
  sets.bar.color = "#e74c3c",      
  # text.scale takes a vector: c(intersection size title, intersection size tick labels, 
  # set size title, set size tick labels, set names, numbers above bars)
  text.scale = c(1.3, 1.1, 1.2, 1, 1.1, 1)
)

grid.text("Intersection of \n Peptide Classifications", 
          x = 0.2,          # Centered horizontally (0 is left, 1 is right)
          y = 0.95,         # Placed near the top (0 is bottom, 1 is top)
          gp = gpar(fontsize = 15, fontface = "bold"))

# intersection between peptides, classification, and search/mapping strategy
UpSetR::upset(
  peptide_binary_matrix_search_and_category,
  sets = all_sets,
  nsets = length(all_sets),        # Ensures all 12 sets are plotted on the y-axis
  nintersects = 40,                # Caps the plot at the top 40 largest intersections
  order.by = "freq",               # Sorts intersections by frequency (largest first)
  keep.order = FALSE,              # Sorts the 12 sets on the left by their total size
  mb.ratio = c(0.55, 0.45),        # Allocates 45% of the plot height to the matrix to fit 12 rows
  point.size = 2.5,                # Slightly smaller dots to fit the grid
  line.size = 0.8,                 # Slightly thinner lines
  main.bar.color = "#2c3e50",      
  sets.bar.color = "#e74c3c",      
  # text.scale takes a vector: c(intersection size title, intersection size tick labels, 
  # set size title, set size tick labels, set names, numbers above bars)
  text.scale = c(1.3, 1.1, 1.2, 1, 1.1, 1) 
)

grid.text("Intersection of Search/Mapping \n and Peptide Classifications", 
          x = 0.2,          # Centered horizontally (0 is left, 1 is right)
          y = 0.95,         # Placed near the top (0 is bottom, 1 is top)
          gp = gpar(fontsize = 15, fontface = "bold"))

# Complex Upset Plot
# 1. Transform only the categories into sets
category_upset_data <- unique_peptides_by_search_and_category |>
  # Add a logical flag for presence
  mutate(present = TRUE) |>
  # Pivot only the categories wide
  pivot_wider(
    names_from = category,
    values_from = present,
    values_fill = FALSE # Fills absences with FALSE
  ) |>
  as.data.frame()

# 2. Extract just the category names to pass to the plot
# This gets all column names EXCEPT strategy and peptide
category_sets <- setdiff(colnames(category_upset_data), c("search_strategy", "peptide"))

colored_upset_plot <- upset(
  category_upset_data,
  intersect = category_sets,
  name = "Classifications",
  # Customize the top bar chart
  base_annotations = list(
    'Intersection Size' = intersection_size(
      mapping = aes(fill = search_strategy)
    ) +
      # You can customize colors here if desired using standard ggplot functions
      scale_fill_brewer(palette = "Set2") +
      theme(legend.position = "bottom")
  ),
  # Additional formatting
  sort_intersections = "descending",
  width_ratio = 0.15 # Adjusts the ratio between the left set sizes and main matrix
)

# Display the plot
colored_upset_plot

# make the required data set for bar plots
search_and_category_barplot <-  unique_peptides_by_search_and_category |> 
  ungroup() |> 
  group_by(search_strategy, category) |> 
  summarise(n = n())

# stacked bar plot
ggplot(search_and_category_barplot, 
       aes(x = search_strategy, y = n, fill = category)) +
  geom_col() +
  labs(
    title = "Unique Peptides by Category and Search Strategy",
    x = "Category",
    y = "Count",
    fill = "Mapping Classification"
  ) +
  theme_minimal()

# bar plot with bars split out
ggplot(search_and_category_barplot, aes(x = search_strategy, y = n, fill = category)) +
  geom_col(position = "dodge") +
  labs(
    title = "Unique Peptides by Category and Search Strategy",
    x = "Category",
    y = "Count",
    fill = "Mapping Classification"
  ) +
  theme_minimal()
