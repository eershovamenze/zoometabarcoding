# Load prerequisites -------------------------------------------------------

## Required packages -------------------------------------------------------
required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "reshape2", "stringr", "forcats",
  "patchwork", "ggpubr", "ggrepel", "ggh4x",
  "GUniFrac", "vegan", "iNEXT",
  "lulu",  "ComplexUpset",    "ggVennDiagram", "pheatmap", "ggplotify", 
  "colorspace",  "ggpmisc",  "viridis","ade4", "parallel", "pbmcapply", "here")

install.packages(setdiff(required_packages, rownames(installed.packages())))
invisible(lapply(required_packages, library, character.only = TRUE))

## Project paths -----------------------------------------------------------
# All paths are relative to the project root (where this repo was cloned).
proj_root    <- here::here()
data_raw     <- file.path(proj_root, "data/raw")
data_proc    <- file.path(proj_root, "data/processed")
out_figures  <- file.path(proj_root, "output/figures")

## Shared colour palettes --------------------------------------------------
colors_sample_type <- c(
  Fjord_2021    = "#DB3A07FF",
  Mock_2022     = "#FE9B2DFF",
  NorthSea_2024 = "#5DC863FF",
  NorthSea_2021 = "#2C728EFF")

colors_labs <- c(
  "Lab 1 (Illumina)"   = "#27AD81FF",
  "Lab 2 (IonTorrent)" = "#7A0403FF")

colors_preservative <- c(
  Ethanol     = "#F05B12FF",
  "Heat-dried" = "#2C728EFF",
  "Freeze-dried" = "#472D7BFF",
  "TES buffer" = "#21908CFF")


## Load input tables ---------------------------------------------------
motu_wide <- read.csv(file.path(data_raw, "all_reads_MOTU.csv"), check.names = FALSE)
metadata  <- read.csv(file.path(data_raw, "all_reads_metadata.csv"), row.names = 1, check.names = FALSE)
metadata$Sample_ID_derep <- paste(metadata$Sequence_run, metadata$Sample_name, sep="_")
metadata_stations <- metadata %>%
  select(-ID_unique, -Rep) %>%
  distinct()
taxonomy <- read.csv(file.path(data_proc, "merged_taxonomy_final_upd.csv"), check.names = FALSE)
taxonomy_class <- read.csv(file.path(data_proc, "taxonomy_classification.csv"), check.names = FALSE)

## Pivot to long format -----------------------------------------------
reads_long_raw <- motu_wide %>%
  pivot_longer(
    cols      = -MIN,
    names_to  = "ID_unique",
    values_to = "reads") %>%
  left_join(metadata, by = "ID_unique")

# Blank-based decontamination ----------------------------------------
# Applied independently within each Sequence_run, because blanks are only
# meaningful relative to samples processed in the same run.
#
# Logic: for each MIN / Sequence_run, zero out any read count in a real sample
# that does not exceed the maximum count observed across the blanks of that run.

decontaminate_run <- function(df, run_name) {
  blanks  <- df %>% filter(Sample_type == "blank")
  samples <- df %>% filter(Sample_type != "blank")
  
  if (nrow(blanks) == 0) {
    message("  No blanks in run: ", run_name, " — skipping decontamination.")
    return(df)
  }
  
  blank_max <- blanks %>%
    group_by(MIN) %>%
    summarise(blank_threshold = max(reads, na.rm = TRUE), .groups = "drop")
  
  samples <- samples %>%
    left_join(blank_max, by = "MIN") %>%
    mutate(
      blank_threshold = replace_na(blank_threshold, 0),
      reads           = ifelse(reads <= blank_threshold, 0, reads)
    ) %>%
    select(-blank_threshold)
  
  bind_rows(blanks, samples)
}

reads_long_decon <- reads_long_raw %>%
  group_by(Sequence_run) %>%
  group_modify(~ decontaminate_run(.x, .y$Sequence_run)) %>%
  ungroup()

#no motus lost


# Report reads removed
total_before <- sum(reads_long_raw$reads,   na.rm = TRUE)
total_after  <- sum(reads_long_decon$reads, na.rm = TRUE)
message(sprintf("Reads removed by blank decontamination: %s (%.1f%%)",
                format(total_before - total_after, big.mark = ","),
                (total_before - total_after) / total_before * 100))
reads_long_raw %>%
  group_by(Sequence_run) %>%
  summarise(before = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    reads_long_decon %>%
      group_by(Sequence_run) %>%
      summarise(after = sum(reads, na.rm = TRUE), .groups = "drop"),
    by = "Sequence_run") %>%
  mutate(
    removed = before - after,
    pct = removed / before * 100) %>%
  mutate(across(c(before, after, removed), ~format(.x, big.mark = ","))) %>%
  mutate(pct = sprintf("%.1f%%", pct)) %>%
  print()

blank_summary <- reads_long_raw %>%
  filter(Sample_type == "blank") %>%
  group_by(Sequence_run, ID_unique) %>%
  summarise(total_reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  group_by(Sequence_run) %>%
  summarise(
    n_blanks = n_distinct(ID_unique),
    mean_reads_per_blank = mean(total_reads),
    .groups = "drop")

## Remove blanks and zero-read rows for downstream use ----------------
reads_long <- reads_long_decon %>%
  filter(Sample_type != "blank") %>%
  filter(reads > 0)

## Wide version for downstream use ----------------
all_reads_mat <- reads_long %>%
  select(MIN, Sample_ID_derep, reads) %>%
  pivot_wider(names_from = Sample_ID_derep, values_from = reads, values_fill = 0, values_fn = sum) %>%
  as.data.frame()
rownames(all_reads_mat) <- all_reads_mat$MIN

## Quick summary -----------------------------------------------
message("\nSamples loaded: ", n_distinct(reads_long$ID_unique))
message("Experiments: ", paste(unique(reads_long$Experiment), collapse = ", "))
message("MOTUs (>0): ", n_distinct(reads_long$MIN))

reads_long %>%
  group_by(Experiment, Laboratory, Sample_source, Sequence_run) %>%
  summarise(
    n_samples   = n_distinct(ID_unique),
    n_motu      = n_distinct(MIN),
    total_reads = sum(reads)) %>%
  print()


# Rarefaction ------------------------------------------

rare_depth <- 30000  # reads per sample (replicates pooled); change here to update everywhere
rare_depth_rep <- 10000  # reads per sample if per replicate

## Check read depth coverage ------------------------------------------
# Report samples that fall below the rarefaction threshold
depth_check <- reads_long %>%
  group_by(Sample_ID_derep) %>%
  summarise(total_reads = sum(reads))

below_threshold <- depth_check %>% filter(total_reads < rare_depth) %>%
  print()

## Pivot to wide, group replicates, rarefy, pivot back to long --------------------------
# Samples below threshold are excluded automatically by Rarefy().

#Replicates pooled
reads_wide_raw <- reads_long %>%
  select(MIN, Sample_ID_derep, reads) %>%
  pivot_wider(names_from = Sample_ID_derep, values_from = reads, values_fill = 0, values_fn = sum) %>%
  as.data.frame()

rownames(reads_wide_raw) <- reads_wide_raw$MIN

#Replicates individually
reads_wide_raw_rep <- reads_long %>%
  select(MIN, ID_unique, reads) %>%
  pivot_wider(names_from = ID_unique, values_from = reads, values_fill = 0) %>%
  as.data.frame()
rownames(reads_wide_raw_rep) <- reads_wide_raw_rep$MIN

##  Rarefy (GUniFrac::Rarefy expects samples as rows) -----------------------------
set.seed(42)  # for reproducibility
rarefied_mat <- t(Rarefy(t(reads_wide_raw[,-1]), depth = rare_depth)$otu.tab.rff)%>%
  as.data.frame()
rarefied_mat_rep <- t(Rarefy(t(reads_wide_raw_rep[,-1]), depth = rare_depth_rep)$otu.tab.rff)%>%
  as.data.frame()

##  Remove MOTUs with zero reads across all samples after rarefaction -----------------------------
rarefied_mat <- rarefied_mat[rowSums(rarefied_mat) > 0, ]
rarefied_mat$MIN <- rownames(rarefied_mat)
rarefied_mat_rep <- rarefied_mat_rep[rowSums(rarefied_mat_rep) > 0, ]
rarefied_mat_rep$MIN <- rownames(rarefied_mat_rep)


##  Pivot back to long and rejoin metadata -----------------------------
reads_rarefied <- rarefied_mat %>%
  as.data.frame() %>%
  pivot_longer(
    cols      = -MIN,
    names_to  = "Sample_ID_derep",
    values_to = "reads") %>%
  filter(reads > 0) %>%
  left_join(metadata_stations, by = "Sample_ID_derep")
reads_rarefied_rep <- rarefied_mat_rep %>%
  as.data.frame() %>%
  pivot_longer(
    cols      = -MIN,
    names_to  = "ID_unique",
    values_to = "reads") %>%
  filter(reads > 0) %>%
  left_join(metadata, by = "ID_unique")

# make a long table with replicates pooled per stations
reads_long_derep <- reads_long %>%
  group_by(across(-c(ID_unique, Rep, reads))) %>%
  summarise(reads = sum(reads))%>%
  ungroup()


# LULU --------------------------------------------------------------------

lulu_result_file <- file.path(data_proc, "LULU_result.csv")

## Run LULU or load cached result ------------------------------------
if (!file.exists(lulu_result_file)) {
  message("Running LULU curation: this may take several minutes...")
  matchlist <- read.table(
    file.path(data_proc, "match_list.txt"),
    sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  # LULU requires a data frame with MOTUs as row names
  reads_for_lulu <- all_reads_mat[, colnames(all_reads_mat) != "MIN"]
  
  lulu_output <- lulu(
    reads_for_lulu,
    matchlist,
    minimum_match = 84,
    minimum_ratio  = 0.99)
  
  # Format and annotate the OTU map
  lulu_result <- lulu_output$otu_map
  lulu_result$MIN <- rownames(lulu_result)
  lulu_result <- lulu_result %>%
    dplyr::rename(parent_MIN = parent_id,
                  no_samples  = spread,
                  type        = curated) %>%
    mutate(type = ifelse(type == "merged", "pseudogene", "parent"))
  
  write.csv(lulu_result, lulu_result_file, row.names = FALSE)
  message("LULU complete. Results written to ", lulu_result_file)
  
} else {
  lulu_result <- read.csv(lulu_result_file)
  message("Loaded cached LULU result: ", nrow(lulu_result), " MOTUs (",
          sum(lulu_result$type == "parent"), " parents, ",
          sum(lulu_result$type == "pseudogene"), " pseudogenes)")
}

## Filtered (parent-only) read matrices -------------------------------
# Rarefied — used for most diversity analyses
reads_rarefied_lulu <- reads_rarefied[
  reads_rarefied$MIN %in% lulu_result[lulu_result$type == "parent", ]$MIN, ]

# Non-rarefied — used for filtering comparisons
reads_all_lulu <- reads_long[
  reads_long$MIN %in% lulu_result[lulu_result$type == "parent", ]$MIN, ]


# Filtering comparison ----------------------------------------------------

taxonomy_full <- merge(taxonomy, taxonomy_class, by='Merged_final_assignment')
## helper function: per-site MOTU / species counts ----------------------------------
count_per_site <- function(df, type_label) {
  df %>%
    filter(Reads > 0) %>%
    group_by(Sample_ID_derep) %>%
    mutate(No_MOTU = n()) %>%
    group_by(Sample_ID_derep, Merged_final_assignment, Taxa, Merged_rank,
             Plankton_classification, No_MOTU) %>%
    summarise(Reads = sum(Reads), .groups = "drop") %>%
    group_by(Sample_ID_derep, No_MOTU) %>%
    mutate(No_species = n()) %>%
    filter(
      Plankton_classification %in% c("Holoplankton", "Meroplankton"),
      Merged_rank %in% c("species", "genus")
    ) %>%
    mutate(No_zoop = n()) %>%
    distinct(Sample_ID_derep, No_MOTU, No_species, No_zoop) %>%
    mutate(Type = type_label)
}

## Build long-format tables for each filtering strategy ---------------
base_melt <- function(mat) {
  melt(mat, id.vars = "MIN", variable.name = "Sample_ID_derep", value.name = "Reads") %>%
    merge(taxonomy_full, by = "MIN")
}

# Strategy 1: Singletons only (no filtering beyond presence)
reads_unfiltered <- base_melt(all_reads_mat) %>%
  count_per_site("singletons")

# Strategy 2: Global abundance >10 reads
reads_10min <- base_melt(all_reads_mat) %>%
  group_by(MIN) %>%
  mutate(Total_reads = sum(Reads)) %>%
  filter(Total_reads > 10) %>%
  ungroup() %>%
  count_per_site("10_reads")

# Strategy 3: LULU pseudogene removal only
reads_lulu <- base_melt(all_reads_mat) %>%
  filter(MIN %in% lulu_result[lulu_result$type == "parent", ]$MIN) %>%
  count_per_site("lulu")

# Strategy 4: Rarefaction only
reads_rare <- base_melt(rarefied_mat) %>%
  count_per_site("rarefaction")

# Strategy 5: Rarefaction + abundance filter (>0.01% in any sample = max > 2 reads at 30k depth)
reads_rare_filt <- base_melt(rarefied_mat) %>%
  filter(MIN %in% lulu_result[lulu_result$type == "parent", ]$MIN) %>%
  group_by(MIN) %>%
  mutate(Max_abund = max(Reads)) %>%
  filter(Max_abund > 2) %>%
  ungroup() %>%
  count_per_site("abundance")

## Combine and add station metadata -----------------------------------
filtering_all <- bind_rows(
  reads_unfiltered, reads_10min, reads_lulu, reads_rare, reads_rare_filt) %>%
  melt(id.vars = c("Sample_ID_derep", "Type")) %>%
  merge(metadata_stations, by = "Sample_ID_derep")%>%
  mutate(
    Sample_source = ifelse(Sample_source == "North Sea", "North Sea 2021", Sample_source)) %>%
  filter(!Sample_source %in% c("Blank", "Blank_ext"))

# Order filtering levels for plotting
filtering_all$Type <- factor(
  filtering_all$Type,
  levels = c("singletons", "10_reads", "lulu", "rarefaction", "abundance"),
  labels = c("Singletons only", ">10 reads total",
             "Pseudogene filtration", "Rarefaction", ">0.01% in any sample"))

filtering_all$variable <- factor(
  filtering_all$variable,
  labels = c(
    "bold('A')~Number~of~MOTUs",
    "bold('B')~Number~of~unique~taxa",
    "bold('C')~Number~of~unique~holo/meroplankton~species"
  )
)

## Plot
p_filtering <- ggboxplot(
  data       = filtering_all,
  x          = "Type",
  y          = "value",
  add        = "jitter",
  add.params = list(alpha = 0.4, color = "Sample_source", shape = "Laboratory")
) +
  labs(x = "", y = "Value") +
  facet_wrap(~variable, scales = "free_y", nrow = 1, labeller = label_parsed) +
  scale_shape_manual(values = c(20, 3)) +
  scale_colour_manual(values = unname(colors_sample_type)) +
  theme_bw() +
  theme(
    panel.border    = element_rect(colour = "black", linewidth = 0.8),
    axis.text.x     = element_text(angle = 45, vjust = 0.5),
    legend.title    = element_text(face = "bold"),
    strip.background = element_blank()
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))


