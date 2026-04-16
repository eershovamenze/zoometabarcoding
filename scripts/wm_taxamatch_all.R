library(worrms)
library(purrr)
library(dplyr)
library(taxize)

wm_taxamatch_all <- function(names, chunk_size = 50) {
  
  # split into chunks of 50 (WoRMS API limit)
  idx_chunks <- split(seq_along(names),
                      ceiling(seq_along(names) / chunk_size))
  
  res_chunks <- map(idx_chunks, function(idx) {
    nm_chunk <- names[idx]
    
    # Try Taxamatch; if it errors, return list(NULL, NULL, ...)
    out <- tryCatch(
      wm_records_taxamatch(nm_chunk),
      error = function(e) {
        # optional: message(e$message)
        vector("list", length(nm_chunk))
      }
    )
    
    # If only one name was given, some versions may return a data.frame instead of list
    if (inherits(out, "data.frame")) {
      out <- list(out)
    }
    
    # Ensure length of out matches length of nm_chunk
    if (length(out) != length(nm_chunk)) {
      out <- rep(list(NULL), length(nm_chunk))
    }
    
    map2(out, nm_chunk, ~ {
      
      ## CASE 1: Taxamatch returned something - use first row
      if (!is.null(.x) && nrow(.x) > 0) {
        return(
          .x %>%
            slice(1) %>% 
            mutate(query = .y, .before = 1)
        )
      }
      
      ## CASE 2: Taxamatch failed -> try taxize::get_wormsid fallback
      aphia_id <- tryCatch(
        taxize::get_wormsid(.y, messages = FALSE, ask=FALSE),
        error = function(e) NA
      )
      
      # coerce to numeric and take the first ID if multiple
      aphia_id_num <- suppressWarnings(as.numeric(aphia_id)[1])
      
      if (!is.na(aphia_id_num)) {
        rec <- tryCatch(
          worrms::wm_record(aphia_id_num),
          error = function(e) NULL
        )
        
        if (!is.null(rec)) {
          rec_df <- tibble::as_tibble(rec)
          
          # ensure match_type exists
          if (!"match_type" %in% names(rec_df)) {
            rec_df$match_type <- "taxize_fallback"
          }
          
          return(
            rec_df %>%
              slice(1) %>%
              mutate(query = .y, .before = 1)
          )
        }
      }
      
      ## CASE 3: nothing worked - mark as not matched
      tibble(
        query = .y,
        status = "not matched",
        match_type = "not matched"
      )
    })
  })
  
  # flatten list and bind into one tibble
  res <- res_chunks %>% flatten() %>% bind_rows()
  
  # make sure status/match_type exist everywhere and are filled
  if (!"status" %in% names(res)) {
    res$status <- NA_character_
  }
  if (!"match_type" %in% names(res)) {
    res$match_type <- NA_character_
  }
  
  res <- res %>%
    mutate(
      status = ifelse(is.na(status), "not matched", status),
      match_type = ifelse(is.na(match_type), "not matched", match_type)
    )
  
  res
}
