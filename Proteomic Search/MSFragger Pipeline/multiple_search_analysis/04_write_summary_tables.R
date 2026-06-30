downstream_search_upstream_map_classified_peptides |> 
  group_by(category) |> 
  summarise(n = n())

upstream_search_downstream_map_classified_peptides |> 
  group_by(category) |> 
  summarise(n = n())

main_search_upstream_map_classified_peptides |> 
  group_by(category) |> 
  summarise(n = n())

upstream_search_classified_peptides  |> 
  group_by(category) |> 
  summarise(n = n())

downstream_search_classified_peptides  |> 
  group_by(category) |> 
  summarise(n = n())

main_search_classified_peptides |> 
  group_by(category) |> 
  summarise(n = n())

# gemini version of turning above summaries to a for loop to write summaries to file:
tibble_names <- c(
  "downstream_search_upstream_map_classified_peptides",
  "upstream_search_downstream_map_classified_peptides",
  "main_search_upstream_map_classified_peptides",
  "upstream_search_classified_peptides",
  "downstream_search_classified_peptides",
  "main_search_classified_peptides"
)

# 2. Set your output file name
output_file <- "peptide_category_summaries.txt"

# 3. Loop through the names, summarize, and write to file
for (name in tibble_names) {
  
  # Fetch the actual tibble using its name
  current_tibble <- get(name)
  
  # Generate the summary
  summary_tibble <- current_tibble |> 
    group_by(category) |> 
    summarise(n = n())
  
  # Write the header to the file
  cat(paste0("========== ", name, " ==========\n"), 
      file = output_file, 
      append = TRUE)
  
  # Write the summarized data (using tab separation for readability)
  write.table(summary_tibble, 
              file = output_file, 
              append = TRUE, 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)
  
  # Add a blank newline for spacing between blocks
  cat("\n", file = output_file, append = TRUE)
}

print(paste("Summaries successfully written to", output_file))

# copilot version:
# list of object names you want to summarize
tbl_names <- c(
  "downstream_search_upstream_map_classified_peptides",
  "upstream_search_downstream_map_classified_peptides",
  "main_search_upstream_map_classified_peptides",
  "upstream_search_classified_peptides",
  "downstream_search_classified_peptides",
  "main_search_classified_peptides"
)

# output file
outfile <- "category_summaries.txt"

# clear file first
write("", file = outfile)

for (nm in tbl_names) {
  
  # retrieve tibble by name
  tbl <- get(nm)
  
  # compute summary
  summary_tbl <- tbl |>
    dplyr::group_by(category) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  # write header
  write(paste0("===== ", nm, " ====="), file = outfile, append = TRUE)
  
  # write summary
  capture.output(print(summary_tbl), file = outfile, append = TRUE)
  
  # blank line between sections
  write("\n", file = outfile, append = TRUE)
}