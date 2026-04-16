---
editor_options: 
  markdown: 
    wrap: 72
---

# zoometabarcoding

R code for the analysis presented in:

**Ershova-Menze E., Westgaard J.-I., Hjellnes H., Falkenhaug T.**
*Optimising zooplankton DNA metabarcoding: methodological considerations
for large-scale monitoring.* (2026) Molecular Ecology Resources

------------------------------------------------------------------------

## Scripts

### `data_analysis.R`

Main analysis script. Loads the two global input tables, runs the full
pipeline from data cleaning through to the filtering comparison figure:

-   Package loading and project path setup
-   Loading `all_reads_MOTU.csv` and `all_reads_metadata.csv`
-   Blank-based decontamination (applied per sequencing run)
-   Rarefaction to 30,000 reads (pooled replicates) and 10,000 reads
    (individual replicates) using `GUniFrac::Rarefy`
-   LULU post-clustering curation (run cached from
    `data/processed/LULU_result.csv`)
-   Comparison of five filtering strategies (singletons, global \>10
    reads, LULU only, rarefaction only, rarefaction + \>0.01% abundance
    threshold)

### `taxonomy_formatting.R`

Generates the merged taxonomy table from three classifier/database
combinations:

-   Wang/RDP classifier against MZGdb (MOTHUR pipeline)
-   mkLTG classifier against COInr
-   BOLDigger against BOLD, including resolution of flagged (ambiguous)
    assignments and a manual curation step
-   WoRMS name validation via `wm_taxamatch_all()` (incremental, cached)
-   Consensus assignment based on pairwise agreement between classifiers

Output files are in the default format from the above tools. Short
sample output files (50 entries) are provided.

Output: `data/processed/merged_taxonomy_final.csv`

### `wm_taxamatch_all.R`

Helper function sourced by `taxonomy_formatting.R`. Queries the WoRMS
API in chunks of 50 names using `worrms::wm_records_taxamatch()`, with a
`taxize` fallback for names that fail the primary lookup.

------------------------------------------------------------------------

## Data

Input data are deposited at DOI: 10.21335/NMDC-742608018. Place files in
data/raw
