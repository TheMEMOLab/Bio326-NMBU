library(tidyverse)
library(patchwork)
library(ggpubr)

# Set your working directory in RStudio before running this script.
# Example:
 setwd("/scratch/2026/BIO326")

project_dir <- getwd()

# Define the two NanoStats folders explicitly.
dir_2025 <- file.path(project_dir, "2025", "bio326NanoStats2025.Prok.dir")
dir_2026 <- file.path(project_dir, "2026", "bio326NanoStats2026.Prok.dir")

# Check that the folders exist.
if (!dir.exists(dir_2025)) {
  stop("2025 directory not found: ", dir_2025)
}

if (!dir.exists(dir_2026)) {
  stop("2026 directory not found: ", dir_2026)
}

infer_method <- function(sample_name) {
  case_when(
    str_detect(sample_name, "^FastPrep") ~ "FastPrep",
    str_detect(sample_name, "^Vortex_SRE") ~ "Vortex_SRE",
    str_detect(sample_name, "^Vortex") ~ "Vortex",
    TRUE ~ "Unknown"
  )
}

read_nanostat_file <- function(file, year) {
  lines <- readr::read_lines(file)
  lines <- lines[lines != ""]

  stats <- lines %>%
    str_trim() %>%
    str_replace_all(",", "") %>%
    str_match("^(.*?):\\s+([0-9\\.]+)$") %>%
    as_tibble(.name_repair = ~ c("Match", "Metric", "Value")) %>%
    filter(!is.na(Metric)) %>%
    transmute(
      Metric = Metric,
      Value = as.numeric(Value)
    )

  sample_name <- sub("\\.NanoStats\\.txt$", "", basename(file))

  stats %>%
    mutate(
      Sample = sample_name,
      Year = as.character(year),
      Method = infer_method(sample_name)
    ) %>%
    select(Sample, Metric, Value, Year, Method)
}

# List input files for each year.
files_2025 <- list.files(dir_2025, pattern = "\\.NanoStats\\.txt$", full.names = TRUE)
files_2026 <- list.files(dir_2026, pattern = "\\.NanoStats\\.txt$", full.names = TRUE)

# Parse files into per-year tidy tables.
stats_2025_tidy <- purrr::map_dfr(files_2025, read_nanostat_file, year = "2025")
stats_2026_tidy <- purrr::map_dfr(files_2026, read_nanostat_file, year = "2026")

# Combine both years into one tidy table.
stats_tidy <- bind_rows(stats_2025_tidy, stats_2026_tidy) %>%
  mutate(
    Year = factor(Year, levels = c("2025", "2026")),
    Method = factor(Method, levels = c("FastPrep", "Vortex", "Vortex_SRE"))
  )

# Build the wide table.
stats_wide <- stats_tidy %>%
  pivot_wider(names_from = Metric, values_from = Value)

# Keep sample order explicit for plotting.
sample_levels <- stats_wide %>%
  distinct(Year, Sample) %>%
  arrange(Year, Sample) %>%
  pull(Sample) %>%
  unique()

stats_tidy <- stats_tidy %>%
  mutate(Sample = factor(Sample, levels = sample_levels))

stats_wide <- stats_wide %>%
  mutate(Sample = factor(Sample, levels = sample_levels))

# Metrics from the template.
overview_metrics <- c(
  "Number of reads",
  "Total bases",
  "Median read length",
  "Mean read length",
  "STDEV read length",
  "Read length N50",
  "Mean read quality",
  "Median read quality"
)

qc_metrics <- c(
  "Mean read length",
  "Read length N50",
  "Mean read quality"
)

metric_colors <- setNames(
  c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"),
  overview_metrics
)

create_overview_plots <- function(data_wide, year_value) {
  purrr::map(overview_metrics, function(metric) {
    ggplot(
      data_wide %>% filter(Year == year_value),
      aes(y = forcats::fct_reorder(Sample, .data[[metric]]), x = .data[[metric]])
    ) +
      geom_col(fill = metric_colors[[metric]]) +
      theme_minimal() +
      labs(title = paste(year_value, metric), x = metric, y = "Sample")
  }) %>%
    set_names(gsub(" ", "", overview_metrics))
}

create_qc_plot <- function(data_tidy, metric, grouping = c("Method", "Year", "Method + Year")) {
  grouping <- match.arg(grouping)

  qc_data <- data_tidy %>%
    filter(Metric == metric)

  qc_data <- if (grouping == "Method") {
    mutate(qc_data, Group = as.character(Method))
  } else if (grouping == "Year") {
    mutate(qc_data, Group = as.character(Year))
  } else {
    mutate(qc_data, Group = paste(Year, Method, sep = " | "))
  }

  qc_data$Group <- factor(qc_data$Group, levels = unique(qc_data$Group))

  palette_values <- c(
    "FastPrep" = "#4C78A8",
    "Vortex" = "#F58518",
    "Vortex_SRE" = "#54A24B",
    "2025" = "#B279A2",
    "2026" = "#FF9DA6",
    "2025 | FastPrep" = "#4C78A8",
    "2025 | Vortex" = "#F58518",
    "2025 | Vortex_SRE" = "#54A24B",
    "2026 | FastPrep" = "#9ECAE9",
    "2026 | Vortex" = "#FFBF79",
    "2026 | Vortex_SRE" = "#86BC86"
  )

  ggplot(qc_data, aes(x = Group, y = Value, fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.65) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    scale_fill_manual(values = palette_values[levels(qc_data$Group)]) +
    theme_minimal() +
    labs(title = paste(metric, grouping), x = grouping, y = metric) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
}

# Per-year plot collections.
plots_2025 <- create_overview_plots(stats_wide, "2025")
plots_2026 <- create_overview_plots(stats_wide, "2026")

overview_2025 <- wrap_plots(plotlist = plots_2025, ncol = 2)
overview_2026 <- wrap_plots(plotlist = plots_2026, ncol = 2)

# Final QC comparison plots.
qc_method_plots <- purrr::map(qc_metrics, create_qc_plot, data_tidy = stats_tidy, grouping = "Method") %>%
  set_names(qc_metrics)

qc_year_plots <- purrr::map(qc_metrics, create_qc_plot, data_tidy = stats_tidy, grouping = "Year") %>%
  set_names(qc_metrics)

qc_method_year_plots <- purrr::map(qc_metrics, create_qc_plot, data_tidy = stats_tidy, grouping = "Method + Year") %>%
  set_names(qc_metrics)

qc_method_panel <- wrap_plots(plotlist = qc_method_plots, ncol = 2)
qc_year_panel <- wrap_plots(plotlist = qc_year_plots, ncol = 2)
qc_method_year_panel <- wrap_plots(plotlist = qc_method_year_plots, ncol = 2)


