# =============================================================================
# taxonomy_formatting.R
# =============================================================================
# Purpose: Process taxonomic assignment outputs from five classifier/database
#          combinations (Wang/MZGdb, Wang/COInr, mkLTG/COInr, mkLTG/MZGdb,
#          BOLDigger), validate species names against WoRMS, and merge into a
#          single consensus taxonomy table.
#
# Inputs
#   OTU_reference.MZGmothur_coi.wang.taxonomy     (Mothur-formatted Wang classifier output)
#   Taxonomy_MOTU_COInr_2026_ltg.tsv              (mkLTG/COInr output)
#   OTU_reference_identification_result_2026.xlsx (BOLDigger output)
#   OTU_reference_bold_result_2026_part_*.xlsx    (BOLDigger hits)
#   bold_tax_final_checked.csv                    (manually curated flags - made manually, see code below)
#   ltg_ranks.csv                                 (table specifying how ltg_rank corresponds to common rank, formatted as below:)
#
#   final_level     ltg_rank mkLTG_rank
#   1           0      no rank    no rank
#   2           1 superkingdom    kingdom
#   3           2      kingdom    kingdom
#   4           3   subkingdom    kingdom
#   5           4       phylum     phylum
#   6           5    subphylum     phylum
#
#   LULU_result.csv                               (LULU results - see main code)
#
# Outputs:
#   data/processed/Merged_taxonomy_final.csv
#   data/processed/WORMS_matched_all.csv   (updated incrementally)
#
# Note: WoRMS lookups require an internet connection.
#       Species already in WORMS_matched_all.csv are
#       skipped to avoid redundant API calls.
# =============================================================================
library(here)
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(reshape2)

source(here("scripts/wm_taxamatch_all.R"))

data_ref  <- here("data", "reference")
data_proc <- here("data", "processed")

source(here("scripts/wm_taxamatch_all.R"))
# WANG / MZGdb taxonomy ---------------------------------------------------------
# Define the taxonomic ranks from the Mothur output file (this will depend on your mothur-formatted taxa file)
tax_ranks <- data.frame(
  final_level = 1:20,
  rank_simple = c(
    "kingdom", "kingdom", "kingdom",
    "phylum",  "phylum",  "phylum",  "phylum",  "phylum",
    "class",   "class",   "class",   "class",
    "order",   "order",   "order",
    "family",  "family",
    "genus",   "genus",
    "species"))
tax_ranks$rank <- paste(tax_ranks$rank_simple, tax_ranks$final_level, sep='_')
# load output taxonomy file
tax_table <- read.delim(file.path(data_ref, "MZGmothur_coi.wang.taxonomy"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(tax_table) <- c("MIN", "taxonomy")

tax_table$taxonomy <- sub(";$", "", tax_table$taxonomy)
tax_table1 <- tax_table %>% 
  separate_wider_delim(taxonomy, ";", names=tax_ranks$rank) %>% 
  melt(id.vars='MIN', value.name='taxa', variable.name='rank') %>% 
  separate_wider_delim(taxa, delim="(", names=c('Taxa', 'Confidence'), too_few='align_start')%>% 
  mutate(Confidence = replace_na(Confidence, '0)'))%>% 
  mutate(Confidence = as.numeric(sub("\\)", "", Confidence)))%>% 
  filter(Confidence >= 95 & !str_detect(Taxa, "unclassified"))%>% 
  merge(tax_ranks, by='rank')%>% 
  group_by(MIN, Taxa)%>%
  summarise(final_level=min(final_level), rank=as.character(rank[final_level==min(final_level)][1]), rank_simple=as.character(rank_simple[final_level==min(final_level)][1]), Final_confidence=Confidence[final_level==min(final_level)][1])%>%
  group_by(MIN)%>%
  mutate(Final_level=max(final_level), Final_rank=as.character(rank_simple[final_level==max(final_level)][1]), Final_confidence=Final_confidence[final_level==max(final_level)][1])%>%
  merge(tax_table[,c(1,2)], by='MIN', all.y=TRUE)%>%
  mutate(rank=replace_na(rank, 'no rank'), Final_level=replace_na(Final_level, 0), Final_confidence=replace_na(Final_confidence, 0), Final_rank=replace_na(Final_rank, 'no rank'))%>%
  mutate(Taxa = gsub("_EXT", "", Taxa))%>%
  mutate(Taxa = gsub("_", " ", Taxa))%>%
  group_by(MIN) %>%
  mutate(Final_assignment = Taxa[Final_level == final_level][1]) %>%
  ungroup()%>%
  mutate(Final_assignment = ifelse(Final_level<13, 'unclassified', Final_assignment), Final_rank = ifelse(Final_level<6, 'no rank', Final_rank))%>%
  dcast(MIN+Final_assignment+Final_rank+Final_confidence~rank, value.var = 'Taxa')%>%
  merge(tax_table, by='MIN', all.y=T)%>%
  mutate(Final_assignment = if_else(Final_assignment == 'Littorinimorpha', "unclassified", Final_assignment))%>%
  mutate(Final_rank = ifelse(Final_assignment=='unclassified', 'no rank', Final_rank)) %>%
  ungroup()

#WORMS taxonomy check
Species_to_check_MZGdb <- tax_table1%>%
  filter(Final_rank == 'species')%>%
  distinct(Final_assignment)
Species_checked_MZGdb <- wm_taxamatch_all(Species_to_check_MZGdb$Final_assignment) %>%
  rename(Final_assignment = query,
         Match.type= match_type,
         ScientificName = scientificname,
         ScientificName_accepted = valid_name,
         AphiaID_accepted = valid_AphiaID) %>%
  select(Final_assignment,
         AphiaID,
         Match.type,
         ScientificName,
         AphiaID_accepted,
         ScientificName_accepted)

#This file is built incrementally (to only look up new hits with each added database)
WORMS_check <- Species_checked_MZGdb

mzgdb_taxonomy_final <- tax_table1 %>%
  left_join(WORMS_check, by = "Final_assignment") %>%
  mutate(Final_assignment = coalesce(ScientificName_accepted, Final_assignment)) %>%
  dplyr::select(
    MIN,
    W_MZGdb_final_assignment = Final_assignment,
    W_MZGdb_rank = Final_rank,
    W_MZGdb_confidence = Final_confidence,
    W_MZGdb_genus = genus_18,
    W_MZGdb_family = family_16,
    W_MZGdb_order = order_13) %>%
  mutate(W_MZGdb_final_assignment = replace_na(W_MZGdb_final_assignment, 'unclassified'),
         W_MZGdb_rank = replace_na(W_MZGdb_rank, 'no rank'),
         W_MZGdb_confidence = replace_na(W_MZGdb_confidence, 0))  %>%
  mutate(MIN = sub(" ", "", MIN))


# mkLTG COInr taxonomy ----------------------------------------------------------
mkltg_tax <- read.delim(file.path(data_ref, "Taxonomy_MOTU_COInr_2026_ltg.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ltg_ranks <- read.csv(file.path(data_ref, 'ltg_ranks.csv'))
#table with matching ltg_rank to common ranks

mkltg_tax <- mkltg_tax %>% mutate(ltg_rank = na_if(ltg_rank, "")) %>%
  mutate(ltg_rank = replace_na(ltg_rank, 'no rank'))%>% 
  merge(ltg_ranks, by='ltg_rank')%>% 
  mutate(ltg_name = gsub(',', '', ltg_name)) 

mkltg_tax$Final_assignment <- mkltg_tax$ltg_name
mkltg_tax[mkltg_tax$final_level<13,]$Final_assignment <- 'unclassified'

Species_to_check_mkltg <- mkltg_tax%>%
  filter(ltg_rank == 'species')%>%
  filter(!Final_assignment %in% WORMS_check$Final_assignment)%>%
  distinct(Final_assignment)

Species_checked_mkltg <- wm_taxamatch_all(Species_to_check_mkltg$Final_assignment) %>%
  rename(Final_assignment = query,
         Match.type= match_type,
         ScientificName = scientificname,
         ScientificName_accepted = valid_name,
         AphiaID_accepted = valid_AphiaID) %>%
  select(Final_assignment,
         AphiaID,
         Match.type,
         ScientificName,
         AphiaID_accepted,
         ScientificName_accepted)

#This file is built incrementally (to only look up new hits with each added database)
WORMS_check <- rbind(WORMS_check, Species_checked_mkltg) %>%
  distinct()


mkltg_coinr_taxonomy_final <- mkltg_tax %>%
  left_join(WORMS_check, by = "Final_assignment") %>%
  mutate(ScientificName_accepted = coalesce(ScientificName_accepted, Final_assignment)) %>%
  dplyr::select(MIN = seqid,
                mkLTG_COInr_final_assignment = ScientificName_accepted,
                mkLTG_COInr_phylum = phylum,
                mkLTG_COInr_order = order,
                mkLTG_COInr_family = family,
                mkLTG_COInr_genus = genus,
                mkLTG_COInr_assigned = ltg_name,
                mkLTG_COInr_rank = mkLTG_rank, 
                mkLTG_COInr_confidence = pid) 


# BOLDigger taxonomy ----------------------------------------------------------
#output file
bold_tax <- read_excel(file.path(data_ref, 'BOLDigger_identification_result_2026.xlsx'))
#metadata output file (depends on how many parts are in results)
bold_hits <- read_excel(file.path(data_ref, 'BOLDigger_result_2026_part_1.xlsx'))

bold_tax <- bold_tax %>% 
  mutate(selected_level = replace_na(selected_level, 'no rank'))
bold_tax$Final_assignment <- apply(bold_tax, 1, function(row) {
  sel <- row["selected_level"]
  if (sel %in% colnames(bold_tax)) {
    return(row[[sel]])
  } else {
    return("unclassified")
  }
})


#Dealing with flagged taxa
bold_tax_flagged <- bold_tax %>%
  filter(str_detect(flags, "\\b2\\b"))

#Species-level double hits
bold_hits_flagged_species <- bold_hits %>%
  filter(id %in% bold_tax_flagged[bold_tax_flagged$selected_level=='species',]$id)%>%
  filter(pct_identity>97)%>%
  filter(!is.na(species)) %>% 
  group_by(id, order, family, genus, species)%>%
  summarise(hits = length(database), pct_identity=mean(pct_identity))%>%
  group_by(id)%>%
  mutate(total_hits=sum(hits))%>%
  mutate(percent_hits=hits/total_hits)%>%
  group_by(id, total_hits) %>%
  summarise(best_id = species[which.max(percent_hits)], percent_hits=max(percent_hits))

#BOLDigger automatically assigns to dominant match
bold_hits_flagged_species$Final_assignment <- bold_hits_flagged_species$best_id
#Flag cases with <10 hits or less than 80% agreement for manual checking
bold_hits_flagged_species[bold_hits_flagged_species$percent_hits < 0.8 | bold_hits_flagged_species$total_hits < 10,]$Final_assignment <- "_Manual_check"
#Automatically downgrade to genus cases with <= 50% agreement
bold_hits_flagged_species[bold_hits_flagged_species$percent_hits <= 0.5 & bold_hits_flagged_species$Final_assignment == "_Manual_check",]$Final_assignment <- "genus"
bold_hits_flagged_species$Level_check <- 'Species'
#Remainder for manual check
bold_hits_flagged_species[bold_hits_flagged_species$Final_assignment == "_Manual_check",]$Level_check <- '_Manual_check'
#Make list of genus double hits (including genus level assignments and "downgrades" from species)
genus_ids <- c(bold_tax_flagged[bold_tax_flagged$selected_level=='genus',]$id, bold_hits_flagged_species[bold_hits_flagged_species$Final_assignment=='genus',]$id)
#exclude genus downgrades in final table
bold_hits_flagged_species <- filter(bold_hits_flagged_species, Final_assignment != 'genus')

#Genus double hits
#note: code below will not produce any hits in sample table, but should work with real output
bold_hits_flagged_genus <- bold_hits %>%
  filter(id %in% genus_ids)%>%
  filter(pct_identity>95)%>%
  filter(!is.na(genus)) %>% 
  group_by(id, genus)%>%
  summarise(hits = length(database), pct_identity=mean(pct_identity))%>%
  group_by(id)%>%
  mutate(total_hits=sum(hits))%>%
  mutate(percent_hits=hits/total_hits)%>%
  group_by(id, total_hits) %>%
  summarise(best_id = genus[which.max(percent_hits)], percent_hits=max(percent_hits))

bold_hits_flagged_genus$Final_assignment <- bold_hits_flagged_genus$best_id
bold_hits_flagged_genus[bold_hits_flagged_genus$percent_hits < 0.8 | bold_hits_flagged_genus$total_hits < 10,]$Final_assignment <- "_Manual_check"
bold_hits_flagged_genus[bold_hits_flagged_genus$percent_hits <= 0.5 & bold_hits_flagged_genus$Final_assignment == "_Manual_check",]$Final_assignment <- "family"
bold_hits_flagged_genus[bold_hits_flagged_genus$percent_hits == 1 | bold_hits_flagged_genus$total_hits == 1,]$Final_assignment <- bold_hits_flagged_genus[bold_hits_flagged_genus$percent_hits == 1,]$best_id
bold_hits_flagged_genus$Level_check <- 'genus'
bold_hits_flagged_genus[bold_hits_flagged_genus$Final_assignment == "_Manual_check",]$Level_check <- '_Manual_check'
family_ids <- c(bold_tax_flagged[bold_tax_flagged$selected_level=='family',]$id, bold_hits_flagged_genus[bold_hits_flagged_genus$Final_assignment=='family',]$id)
bold_hits_flagged_genus <- filter(bold_hits_flagged_genus, Final_assignment != 'family')

#Family double hits
bold_hits_flagged_family <- bold_hits %>%
  filter(id %in% family_ids)%>%
  filter(pct_identity>90)%>%
  filter(!is.na(family)) %>% 
  group_by(id, family)%>%
  summarise(hits = length(database), pct_identity=mean(pct_identity))%>%
  group_by(id)%>%
  mutate(total_hits=sum(hits))%>%
  mutate(percent_hits=hits/total_hits)%>%
  group_by(id, total_hits) %>%
  summarise(best_id = family[which.max(percent_hits)], percent_hits=max(percent_hits))
bold_hits_flagged_family$Final_assignment <- bold_hits_flagged_family$best_id
bold_hits_flagged_family[bold_hits_flagged_family$percent_hits < 0.8 | bold_hits_flagged_family$total_hits < 10,]$Final_assignment <- "_Manual_check"
bold_hits_flagged_family[bold_hits_flagged_family$percent_hits <= 0.5 & bold_hits_flagged_family$Final_assignment == "_Manual_check",]$Final_assignment <- "order"
bold_hits_flagged_family[bold_hits_flagged_family$percent_hits == 1,]$Final_assignment <- bold_hits_flagged_family[bold_hits_flagged_family$percent_hits == 1,]$best_id
bold_hits_flagged_family$Level_check <- 'family'
bold_hits_flagged_family[bold_hits_flagged_family$Final_assignment == "_Manual_check",]$Level_check <- '_Manual_check'
order_ids <- c(bold_tax_flagged[bold_tax_flagged$selected_level=='order',]$id, bold_hits_flagged_family[bold_hits_flagged_family$Final_assignment=='order',]$id)
bold_hits_flagged_family <- filter(bold_hits_flagged_family, Final_assignment != 'order')

#Order double hits
bold_hits_flagged_order <- bold_hits %>%
  filter(id %in% order_ids)%>%
  filter(pct_identity>85)%>%
  filter(!is.na(order)) %>% 
  group_by(id, order)%>%
  summarise(hits = length(database), pct_identity=mean(pct_identity))%>%
  group_by(id)%>%
  mutate(total_hits=sum(hits))%>%
  mutate(percent_hits=hits/total_hits)%>%
  group_by(id, total_hits) %>%
  summarise(best_id = order[which.max(percent_hits)], percent_hits=max(percent_hits))
bold_hits_flagged_order$Final_assignment <- bold_hits_flagged_order$best_id
bold_hits_flagged_order[bold_hits_flagged_order$percent_hits < 0.8,]$Final_assignment <- "unclassified"
bold_hits_flagged_order$Level_check <- 'order'
bold_hits_flagged_order[bold_hits_flagged_order$Final_assignment == "unclassified",]$Level_check <- 'no rank'

flagged_all <- rbind(bold_hits_flagged_species, bold_hits_flagged_genus, bold_hits_flagged_family, bold_hits_flagged_order) %>%
  mutate(Manual_check = ifelse(Final_assignment=='_Manual_check', 'manual', 'automatic')) %>%
  dplyr::rename(Final_assignment_check = Final_assignment)

###Use code below for manual checking
bold_tax_manual <- bold_tax %>% 
  left_join(flagged_all, by=c('id')) %>% 
  mutate(Assignment_checked_BD = NA,
         Level_check_BD = NA,
         Manual_outcome_BD= NA,
         Reason_BD = NA)
colnames(bold_tax_manual)
manual_checks <- bold_hits[bold_hits$id %in% bold_tax_manual[bold_tax_manual$Manual_check == 'manual',]$id,]

# Write tables, sort Final_assignment (use _Manual_check) and manually check and fill in columns: 
# Assignment_checked_BD = Manually checked final assignment
# Level_checked_BD = Manually checked final level
# Manual_outcome_BD = Outcome of check ("Down rank" or "kept rank")
# Reason_BD = reason for outcome

# uncomment to write tables:
# write.csv(bold_tax_manual, file.path(data_proc,'bold_tax_final_to_check.csv'))
# write.csv(manual_checks, file.path(data_proc,'bold_hits_for_manual_check.csv'))

# save new table as bold_tax_final_checked.csv

# Load manually checked table
bold_final_checked <- read.csv(file.path(data_proc, 'bold_tax_final_checked.csv'))

#For test purposes (pretend we checked the )
bold_final_checked <- bold_tax_manual %>%
  mutate(Level_check = ifelse(Level_check == '_Manual_check', selected_level, Level_check)) %>%
  mutate(Final_assignment_check = ifelse(Final_assignment_check == '_Manual_check', Final_assignment, Final_assignment_check))


Species_to_check_BOLD <- bold_final_checked%>%
  filter(selected_level == 'species')%>%
  filter(!Final_assignment %in% WORMS_check$Final_assignment)%>%
  distinct(Final_assignment)

Species_checked_BOLD <- wm_taxamatch_all(Species_to_check_BOLD$Final_assignment) %>%
  rename(Final_assignment = query,
         Match.type= match_type,
         ScientificName = scientificname,
         ScientificName_accepted = valid_name,
         AphiaID_accepted = valid_AphiaID) %>%
  select(colnames(WORMS_check))

WORMS_check <- bind_rows(WORMS_check, Species_to_check_BOLD)

bd_taxonomy_final <- bold_final_checked %>%
  mutate(Final_assignment_check = coalesce(Final_assignment_check, Final_assignment))%>%
  left_join(WORMS_check, by = "Final_assignment")%>%
  mutate(ScientificName_accepted = coalesce(ScientificName_accepted, Final_assignment_check))%>%
  mutate(Level_checked_BD =coalesce(Level_check, selected_level))%>%
  mutate(ScientificName_accepted = ifelse(Level_checked_BD %in% c('class', 'phylum', 'order'), 'unclassified', ScientificName_accepted))%>%
  mutate(Level_checked_BD = ifelse(Level_checked_BD %in% c('class', 'phylum', 'order'), 'no rank', Level_checked_BD)) %>%
  mutate(Manual_check = replace_na(Manual_check, 'none'))%>%
  dplyr::select(MIN = id,
                BD_final_assignment = ScientificName_accepted,
                BD_phylum = phylum,
                BD_class = class,
                BD_order = order,
                BD_family = family,
                BD_genus = genus,
                BD_assigned_raw = Final_assignment,    
                BD_flags = flags,
                BD_rank = Level_checked_BD,
                BD_flag_check = Manual_check,
                BD_manual_check_outcome = Manual_outcome_BD,
                BD_reason_outcome = Reason_BD)

# Merged taxonomy ---------------------------------------------------------
lulu_result_file <- file.path(data_proc, "LULU_result.csv")
lulu_result <- read.csv(lulu_result_file)

#Function to get lowest common rank across taxonomies
get_lowest_common_rank <- function(picked) {
  ranks <- c("genus", "family", "order")
  bd_vals <- picked[1:3]
  mzg_vals <- picked[4:6]
  mkltg_vals <- picked[7:9]
  
  for (i in seq_along(ranks)) {
    vals <- c(bd_vals[[i]], mzg_vals[[i]], mkltg_vals[[i]])
    
    # Remove NA, unclassified, and blank entries
    vals_clean <- vals[!is.na(vals) & vals != "unclassified" & vals != "" & !grepl("^\\s*$", vals)]
    
    # If 2+ values remain and all are equal assign consensus
    if (length(vals_clean) >= 2 && length(unique(vals_clean)) == 1) {
      return(c(vals_clean[[1]], ranks[[i]]))
    }
  }
  
  return(c("unclassified", "no rank"))  # no agreement found
}

merged_taxonomy <- bd_taxonomy_final %>%
  left_join(mzgdb_taxonomy_final, by = "MIN") %>%
  left_join(mkltg_coinr_taxonomy_final, by = "MIN")%>%
  left_join(select(lulu_result, c('MIN', 'type')), by = "MIN")%>%
  mutate(type = replace_na(type, "parent"))%>%
  #Make order level assignments unclassified
  mutate(BD_final_assignment = ifelse(BD_rank == 'order', 'unclassified', BD_final_assignment),
         W_MZGdb_final_assignment = ifelse(W_MZGdb_rank == 'order', 'unclassified', W_MZGdb_final_assignment),
         mkLTG_COInr_final_assignment = ifelse(mkLTG_COInr_rank == 'order', 'unclassified', mkLTG_COInr_final_assignment),
  )%>%
  #Agreement between 3
  mutate(Agreement = ifelse(BD_final_assignment==W_MZGdb_final_assignment & BD_final_assignment==mkLTG_COInr_final_assignment & mkLTG_COInr_final_assignment==W_MZGdb_final_assignment, 'All', 'None'))%>%
  mutate(Agreement = ifelse(Agreement == 'None' & BD_final_assignment==W_MZGdb_final_assignment, 'BD+MZGbd', Agreement))%>%
  mutate(Agreement = ifelse(Agreement == 'None' & mkLTG_COInr_final_assignment==W_MZGdb_final_assignment, 'MZGbd+mkLTG', Agreement))%>%
  mutate(Agreement = ifelse(Agreement == 'None' & mkLTG_COInr_final_assignment==BD_final_assignment, 'BD+mkLTG', Agreement))%>%
  #If all three agree, assign to consensus
  mutate(Merged_final_assignment = ifelse(Agreement=='All', BD_final_assignment, '_check'))%>%
  #If two agree at species level or the third is 'unclassified', assign to consensus. If the consensus is "unclassified", then assign to third
  mutate(Merged_final_assignment = ifelse(Agreement=='BD+MZGbd' & mkLTG_COInr_final_assignment == 'unclassified', BD_final_assignment, Merged_final_assignment))%>%
  mutate(Merged_final_assignment = ifelse(Agreement=='BD+MZGbd' & BD_rank == 'species', BD_final_assignment, Merged_final_assignment))%>%
  mutate(Merged_final_assignment = ifelse(Agreement=='BD+MZGbd' & BD_final_assignment == 'unclassified', mkLTG_COInr_final_assignment, Merged_final_assignment))%>%
  mutate(Merged_final_assignment = ifelse(Agreement=='MZGbd+mkLTG' & BD_final_assignment == 'unclassified', mkLTG_COInr_final_assignment, Merged_final_assignment))%>%
  mutate(Merged_final_assignment = ifelse(Agreement=='MZGbd+mkLTG' & W_MZGdb_rank == 'species', mkLTG_COInr_final_assignment, Merged_final_assignment))%>%
  mutate(Merged_final_assignment = ifelse(Agreement=='MZGbd+mkLTG' & mkLTG_COInr_final_assignment == 'unclassified', BD_final_assignment, Merged_final_assignment))%>%
  mutate(Merged_final_assignment = ifelse(Agreement=='BD+mkLTG' & W_MZGdb_final_assignment == 'unclassified', mkLTG_COInr_final_assignment, Merged_final_assignment))%>%
  mutate(Merged_final_assignment = ifelse(Agreement=='BD+mkLTG' & BD_rank == 'species', mkLTG_COInr_final_assignment, Merged_final_assignment))%>%
  mutate(Merged_final_assignment = ifelse(Agreement=='BD+mkLTG' & mkLTG_COInr_final_assignment == 'unclassified', W_MZGdb_final_assignment, Merged_final_assignment))%>%
  #If ambiguity was manually validated with BOLDigger, assign to BOLDigger
  mutate(Merged_final_assignment = ifelse(BD_flag_check == "manual", BD_final_assignment, Merged_final_assignment))%>%
  #If there is no remaining agreement and ambiguity was automatically validated with BOLDigger, assign to BOLDigger
  mutate(Merged_final_assignment = ifelse(Merged_final_assignment=='_check' & BD_flag_check == "automatic", BD_final_assignment, Merged_final_assignment))%>%
  #If there is no remaining agreement and likely pseudogene, then assign 'unclassified'
  mutate(Merged_final_assignment = ifelse(Merged_final_assignment=='_check' & type == "pseudogene", 'unclassified', Merged_final_assignment))%>%
  mutate(Merged_rank = ifelse(Merged_final_assignment==mkLTG_COInr_final_assignment, mkLTG_COInr_rank, '_check'))%>%
  mutate(Merged_rank = ifelse(Merged_final_assignment==W_MZGdb_final_assignment, W_MZGdb_rank, Merged_rank))%>%
  mutate(Merged_rank = ifelse(Merged_final_assignment==BD_final_assignment, BD_rank, Merged_rank))%>%
  mutate(Merged_rank = ifelse(Merged_final_assignment=='unclassified', 'no rank', Merged_rank))

#For remaining, assign to lowest common rank
merged_taxonomy <- merged_taxonomy %>%
  rowwise() %>%
  mutate(
    lowest_common = list(get_lowest_common_rank(pick(
      BD_genus, BD_family, BD_order,
      W_MZGdb_genus, W_MZGdb_family, W_MZGdb_order,
      mkLTG_COInr_genus, mkLTG_COInr_family, mkLTG_COInr_order
    ))),
    Merged_final_assignment = ifelse(
      Merged_final_assignment == "_check",
      lowest_common[[1]],
      Merged_final_assignment
    ),
    Merged_rank = ifelse(
      Merged_final_assignment == lowest_common[[1]],
      lowest_common[[2]],
      Merged_rank
    )
  ) %>%
  mutate(Merged_final_assignment = replace_na(Merged_final_assignment, BD_final_assignment))%>%
  mutate(Merged_rank = ifelse(Merged_final_assignment == 'unclassified', 'no rank', Merged_rank))%>%
  mutate(Merged_rank = replace_na(Merged_rank, 'no rank'))%>%
  ungroup()%>%
  select(-lowest_common)


write.csv(WORMS_check, file.path(data_proc, "WORMS_matched_all.csv"), row.names = FALSE)
write.csv(merged_taxonomy,
          file.path(data_proc, "Merged_taxonomy_final_test.csv"),
          row.names = FALSE)
