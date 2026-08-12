## Project specific functions
##
## fcns_patho_grouping_v1.R
## ==========================================
## Purpose: Define whether each isolated organism is a likely cause of infection
##          (pathogen) or a likely contaminant, then aggregate these classifications
##          up through four hierarchical levels:
##            organism -> specimen -> episode -> person
##
## Reference table required:
##   algorithms/patho_grouping_v1.0_20260809.csv
##   Columns: Org_name (organism name), pathogenicity (high/low/depends/""),
##            acorn_pathogen (definite/possible/contaminant/depends)
##
## Key design principle:
##   - High pathogenicity -> always "likely cause"
##   - Low pathogenicity -> generally "likely contaminant", EXCEPT for
##     Staphylococcus species when the patient has an indwelling medical device
##     (ventilator, peripheral IV, central IV catheter), in which case they
##     are upgraded to "likely cause" (device-associated infection logic).
##
## Data dependencies (must exist in the calling environment):
##   organism_patho  : data frame, the pathogenicity lookup table (loaded from CSV)
##   odat0           : organism-level data, one row per isolated organism
##   edat0           : episode-level data, one row per clinical episode
##   odat2           : odat0 after organism-level processing (now includes
##                     likely_cause and likely_contaminant columns)
##   edat2           : edat0 after episode-level processing
##
## Output columns added by these functions:
##   Organism level:  likely_cause, likely_contaminant
##   Specimen level:  likely_causes_s, likely_contaminants_s,
##                    any_causes_s, any_contaminants_s
##   Episode level:   likely_causes_e, likely_contaminants_e,
##                    any_causes_e, any_contaminants_e
##   Person level:    likely_causes_p, likely_contaminants_p,
##                    any_causes_p, any_contaminants_p
##


#####################################################################
#                                                                   #
###      Define cause and/or contaminant at ORGANISM level        ###
#                                                                   #
# For each isolated organism, classify it as a likely cause of      #
# infection ("Yes") or not ("No"), based on:                        #
#   1. The pathogenicity lookup table (organism_patho)              #
#   2. The presence of indwelling medical devices in the episode    #
#                                                                   #
# Required input columns:                                           #
#   odat0:  episode_id, orgname                                     #
#   edat0:  episode_id, med_device_vent, med_device_pcv,            #
#           med_device_cvc                                           #
#   organism_patho: Org_name, pathogenicity                         #
#                                                                   #
# Output: returns a single classification string per organism       #
#####################################################################

#' Assign likely cause status to a single organism
#'
#' This function creates a new column in the organism-level data frame to
#' define if each reported organism is a likely cause based on pre-specified
#' rules specific to the project's research question.
#'
#' @param orgname    Character. The name of the organism (will be lowercased and trimmed).
#' @param episode_id Character/Integer. The episode identifier, used to look up
#'                   device information from edat0.
#'
#' @return Character. One of:
#'   - "Yes"                            : likely causative pathogen
#'   - "No"                             : likely contaminant (low pathogenicity, no device)
#'   - "No growth"                      : culture showed no growth
#'   - "No organism name / culture result" : missing culture data
#'   - "Ungrouped"                      : organism not found in the pathogenicity lookup table
#'   - "No conclusion"                  : pathogenicity value not recognised
#'
#' @examples
#'   odat0$likely_cause <- mapply(assign_likely_cause,
#'                                 orgname    = odat0$orgname,
#'                                 episode_id = odat0$episode_id)
#'
assign_likely_cause <- function(orgname, episode_id) {

  # ---- Step 1: Clean the organism name ----
  orgname <- tolower(trimws(orgname))

  # ---- Step 2: Look up the episode to get device information ----
  episode_rows <- edat0[which(edat0$episode_id == episode_id), ]

  # Guard: if episode_id is not found in edat0, we cannot check device status.
  # Conservative default: treat as if no device present (low-patho Staph -> No).
  if (nrow(episode_rows) == 0) {
    episode_rows <- data.frame(
      med_device_vent = NA_character_,
      med_device_pcv  = NA_character_,
      med_device_cvc  = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  # ---- Step 3: Handle special / missing values ----
  ## [YW-2025-12-11: There's a need to monitor unique values in orgname column
  ## in the raw data, to ensure new values can be handled properly in case the
  ## data entering practice changes]
  if (orgname == "no growth") {
    return("No growth")
  } else if (orgname == "no organism name / culture result") {
    return("No organism name / culture result")
  }

  # ---- Step 4: Match organism against the pathogenicity lookup table ----
  # Primary: resolve orgname to a mo code via the AMR package, then match on mo_code.
  # Fallback: if AMR cannot resolve the name (UNKNOWN), fall back to direct string
  # matching on the cleaned Org_name column (backward compatibility).
  org_mo <- suppressMessages(as.character(as.mo(orgname, keep_synonyms = FALSE)))

  if (org_mo != "UNKNOWN") {
    matched_row <- organism_patho[!is.na(organism_patho$mo_code) &
                                   organism_patho$mo_code == org_mo, ]
    # mo-code match found — done
  } else {
    # AMR could not resolve — fall back to cleaned string matching
    organism_patho$Org_name <- tolower(trimws(organism_patho$Org_name))
    matched_row <- organism_patho[organism_patho$Org_name == orgname, ]
  }

  # If the organism is not in the lookup table, we cannot classify it
  if (nrow(matched_row) == 0) {
    return("Ungrouped")
  }

  # ---- Step 5: Decision logic based on pathogenicity level ----
  pathogenicity <- tolower(trimws(matched_row$pathogenicity[1]))

  # --- 5a. High pathogenicity: always a likely cause ---
  if (pathogenicity == "high") {
    return("Yes")
  }

  # --- 5b. Low pathogenicity (or empty, treated as low): depends on device ---
  if (pathogenicity == "low" || pathogenicity == "") {

    # For Staphylococcus species specifically, check device status.
    # Rationale: coagulase-negative staphylococci (CoNS) are common skin
    # commensals but become true pathogens when intravascular devices or
    # endotracheal tubes provide a route for invasion.
    if (grepl("staph", orgname)) {

      # Ventilator / intubation -> respiratory device-associated infection
      if (!is.na(episode_rows$med_device_vent) &&
          episode_rows$med_device_vent == "Intubation / Mechanical ventilation") {
        return("Yes")

      # Peripheral IV catheter -> line-associated infection
      } else if (!is.na(episode_rows$med_device_pcv) &&
                 episode_rows$med_device_pcv == "Peripheral IV catheter") {
        return("Yes")

      # Central IV catheter -> central line-associated bloodstream infection (CLABSI)
      } else if (!is.na(episode_rows$med_device_cvc) &&
                 episode_rows$med_device_cvc == "Central IV catheter") {
        return("Yes")

      # Staphylococcus with low pathogenicity AND no indwelling device -> contaminant
      } else {
        return("No")
      }

    # Non-staphylococcus with low pathogenicity: always a contaminant
    } else {
      return("No")
    }
  }

  # --- 5c. Any other pathogenicity value: cannot reach a conclusion ---
  return("No conclusion")
}


#' Assign likely contaminant status (inverse of likely cause)
#'
#' For organisms that are NOT a likely cause, they are classified as a likely
#' contaminant, and vice versa.  Special sentinel values (No growth, Ungrouped,
#' etc.) are passed through unchanged so they can be handled downstream.
#'
#' @param likely_cause Character. Output from assign_likely_cause().
#'
#' @return Character. "Yes" if the organism is a likely contaminant,
#'   "No" if it is a likely cause, and sentinel values passed through as-is.
#'
#' @examples
#'   odat0$likely_contaminant <- sapply(odat0$likely_cause,
#'                                       assign_likely_contaminant)
#'
assign_likely_contaminant <- function(likely_cause) {

  likely_cause <- trimws(likely_cause)

  # Sentinel values: pass through unchanged (neither cause nor contaminant)
  if (likely_cause == "No growth") {
    return("No growth")
  }

  if (likely_cause == "No organism name / culture result") {
    return("No organism name / culture result")
  }

  if (likely_cause == "Ungrouped") {
    return("Ungrouped")
  }

  # Inverse logic: if it is NOT a cause, it IS a contaminant
  if (likely_cause == "No") {
    return("Yes")
  }

  if (likely_cause == "Yes") {
    return("No")
  }

  # Catch-all for any unexpected values (e.g., "No conclusion")
  return("No conclusion")
}


#####################################################################
#                                                                   #
###      Define cause and/or contaminant at SPECIMEN level       ###
#                                                                   #
# Aggregate organism-level classifications up to the specimen.      #
# A single specimen (e.g., blood culture bottle, urine sample) may  #
# yield multiple organisms. These functions summarise which         #
# organisms are considered causes / contaminants within that        #
# specimen.                                                         #
#                                                                   #
# Required input columns:                                           #
#   odat2: specimen_id, orgname, likely_cause, likely_contaminant   #
#                                                                   #
#####################################################################

#' Collapse organism-level likely causes to specimen level
#'
#' For a given specimen, collect all organism names that were classified as
#' likely cause ("Yes"), and return them as a single summary string.
#'
#' @param specimen_id Character/Integer. The specimen identifier.
#'
#' @return Character string:
#'   - "Not reported" if no organisms were marked "Yes"
#'   - The organism name if exactly one
#'   - "Multiple: 'orgA','orgB',..." if more than one
#'
assign_likely_causes_s <- function(specimen_id) {

  # Subset to this specimen
  org_rows <- odat2[which(odat2$specimen_id == specimen_id), ]
  org_rows$likely_cause <- trimws(org_rows$likely_cause)

  # Extract the names of all organisms flagged as likely cause
  likely_causes <- trimws(unique(
    org_rows[which(org_rows$likely_cause == "Yes"), "orgname", drop = TRUE]
  ))

  # Return appropriate summary
  if (length(likely_causes) == 0) {
    "Not reported"
  } else if (length(likely_causes) == 1) {
    likely_causes[1]
  } else {
    paste0("Multiple: '", paste0(likely_causes, collapse = "','"), "'")
  }
}


#' Collapse organism-level likely contaminants to specimen level
#'
#' Mirror of assign_likely_causes_s for contaminants.
#'
assign_likely_contaminants_s <- function(specimen_id) {

  org_rows <- odat2[which(odat2$specimen_id == specimen_id), ]
  org_rows$likely_contaminant <- trimws(org_rows$likely_contaminant)

  likely_contaminants <- trimws(unique(
    org_rows[which(org_rows$likely_contaminant == "Yes"), "orgname", drop = TRUE]
  ))

  if (length(likely_contaminants) == 0) {
    "Not reported"
  } else if (length(likely_contaminants) == 1) {
    likely_contaminants[1]
  } else {
    paste0("Multiple: '", paste0(likely_contaminants, collapse = "','"), "'")
  }
}


#' Convert specimen-level cause string to binary Yes/No
#'
#' Useful for tabulation and filtering: turns the descriptive specimen-level
#' string into a simple Yes/No indicator.
#'
assign_any_causes_s <- function(likely_causes_s) {

  likely_causes_s <- trimws(likely_causes_s)

  if (likely_causes_s == "Not reported") {
    return("No")
  }
  if (likely_causes_s == "No") {
    return("No")
  }
  if (likely_causes_s == "Ungrouped") {
    return("Ungrouped")
  }

  # Any actual organism name (or "Multiple: ...") counts as Yes
  return("Yes")
}


#' Convert specimen-level contaminant string to binary Yes/No
#'
assign_any_contaminants_s <- function(likely_contaminants_s) {

  likely_contaminants_s <- trimws(likely_contaminants_s)

  # All sentinel / negative values -> "No"
  if (likely_contaminants_s == "Not reported") {
    return("No")
  }
  if (likely_contaminants_s == "No organism name / culture result") {
    return("No organism name / culture result")
  }
  if (likely_contaminants_s == "No") {
    return("No")
  }
  if (likely_contaminants_s == "No growth") {
    return("No")
  }
  if (likely_contaminants_s == "Ungrouped") {
    return("Ungrouped")
  }

  # Any actual contaminant organism -> "Yes"
  return("Yes")
}


#####################################################################
#                                                                   #
###      Define cause and/or contaminant at EPISODE level        ###
#                                                                   #
# Aggregate organism-level classifications to the episode.          #
# An episode (e.g., a hospital admission for suspected sepsis) may  #
# include multiple specimens, each potentially with multiple         #
# organisms.                                                        #
#                                                                   #
# Special rule: When multiple likely causes are present and one of  #
# them is the generic "Staphylococcus sp", it is removed from the   #
# cause list (the more specific Staphylococcus species takes        #
# priority) and may be reclassified as a contaminant if another     #
# cause exists.                                                     #
#                                                                   #
#####################################################################

#' Collapse organism-level likely causes to episode level
#'
#' Collects all unique organism names flagged "Yes" across all specimens in
#' an episode.  When multiple causes exist and the non-specific
#' "Staphylococcus sp" is among them, it is dropped in favour of more
#' precisely identified Staphylococcus species.
#'
#' @param episode_id Character/Integer. The episode identifier.
#'
#' @return Character string:
#'   - "No growth" if no organism was marked "Yes"
#'   - A single organism name, or "Multiple: 'orgA','orgB',..."
#'
assign_likely_causes_e <- function(episode_id) {

  # Subset to this episode
  org_rows <- odat2[odat2$episode_id == episode_id, ]

  # Get unique organism names flagged as likely cause
  likely_causes <- unique(
    org_rows[trimws(org_rows$likely_cause) == "Yes", "orgname", drop = TRUE]
  )
  likely_causes <- trimws(likely_causes)
  likely_causes <- tolower(likely_causes)   # lowercase for safe comparison

  # No causes identified
  if (length(likely_causes) == 0) {
    return("No growth")
  }

  # Deduplication rule: if there are multiple causes and one is the generic
  # "Staphylococcus sp", remove it.  The more specific identification
  # (e.g., "Staphylococcus aureus") takes precedence.
  if (length(likely_causes) > 1 && "staphylococcus sp" %in% likely_causes) {
    likely_causes <- setdiff(likely_causes, "staphylococcus sp")
  }

  # Return summary
  if (length(likely_causes) == 1) {
    return(likely_causes)
  } else {
    return(paste0("Multiple: '", paste0(likely_causes, collapse = "','"), "'"))
  }
}


#' Collapse organism-level likely contaminants to episode level
#'
#' Mirror of assign_likely_causes_e, with the complementary logic:
#' if "Staphylococcus sp" was removed from the cause list (because another
#' Staphylococcus species was a confirmed cause), it is added to the
#' contaminant list.
#'
assign_likely_contaminants_e <- function(episode_id) {

  org_rows <- odat2[odat2$episode_id == episode_id, ]

  # Clean columns
  org_rows$likely_cause       <- trimws(org_rows$likely_cause)
  org_rows$likely_contaminant <- trimws(org_rows$likely_contaminant)

  # Get lists of causes and contaminants (lowercase for case-insensitive comparison)
  likely_causes       <- tolower(trimws(unique(
    org_rows[org_rows$likely_cause == "Yes", "orgname", drop = TRUE]
  )))
  likely_contaminants <- tolower(trimws(unique(
    org_rows[org_rows$likely_contaminant == "Yes", "orgname", drop = TRUE]
  )))

  # Mirror of the episode cause deduplication rule:
  # If "Staphylococcus sp" is among multiple causes, it is downgraded to a
  # contaminant (because a more specific Staph species is the true cause).
  if (length(likely_causes) > 1 && "staphylococcus sp" %in% likely_causes) {
    if (!("staphylococcus sp" %in% likely_contaminants)) {
      likely_contaminants <- c(likely_contaminants, "staphylococcus sp")
    }
  }

  # Return summary
  if (length(likely_contaminants) == 0) {
    return("No growth")
  } else if (length(likely_contaminants) == 1) {
    return(likely_contaminants)
  } else {
    return(paste0("Multiple: '", paste0(likely_contaminants, collapse = "','"), "'"))
  }
}


#' Convert episode-level cause string to binary Yes/No
#'
assign_any_causes_e <- function(likely_causes_e) {

  if (likely_causes_e == "Not reported") {
    return("No")
  }
  if (likely_causes_e == "No growth") {
    return("No")
  }
  if (likely_causes_e == "Ungrouped") {
    return("Ungrouped")
  }

  # Any other value (organism name or "Multiple: ...") -> Yes
  return("Yes")
}


#' Convert episode-level contaminant string to binary Yes/No
#'
assign_any_contaminants_e <- function(likely_contaminants_e) {

  # All sentinel / negative values -> "No"
  if (likely_contaminants_e == "Not reported") {
    return("No")
  }
  if (likely_contaminants_e == "No organism name / culture result") {
    return("No organism name / culture result")
  }
  if (likely_contaminants_e == "No") {
    return("No")
  }
  if (likely_contaminants_e == "No growth") {
    return("No")
  }
  if (likely_contaminants_e == "Ungrouped") {
    return("Ungrouped")
  }

  # Any actual contaminant organism -> "Yes"
  return("Yes")
}


#####################################################################
#                                                                   #
###      Define cause and/or contaminant at PERSON level         ###
#                                                                   #
# Aggregate episode-level classifications to the person (patient).  #
# A single patient may have multiple episodes over time. These      #
# functions summarise all organisms ever identified as causes or     #
# contaminants across all of the patient's episodes.                 #
#                                                                   #
# Required input columns:                                           #
#   edat2: patient_id, likely_causes_e, likely_contaminants_e       #
#                                                                   #
#####################################################################

#' Collapse episode-level likely causes to person level
#'
#' Across all episodes for a given patient, collect the unique likely cause
#' organisms, excluding episodes with "No growth".
#'
#' @param patient_id Character/Integer. The patient identifier.
#'
#' @return Character string:
#'   - "No growth" if no organisms were found across any episode
#'   - A single organism name, or "Multiple: '...','...'"
#'
assign_likely_causes_p <- function(patient_id) {

  # Get all episodes for this patient
  episodes_of_person <- edat2[edat2$patient_id == patient_id, ]
  episodes_of_person$likely_causes_e <- trimws(episodes_of_person$likely_causes_e)

  # Collect unique causes, excluding "No growth" episodes
  likely_causes <- trimws(unique(
    episodes_of_person[
      episodes_of_person$likely_causes_e != "No growth",
      "likely_causes_e",
      drop = TRUE
    ]
  ))

  # Return summary
  if (length(likely_causes) == 0) {
    "No growth"
  } else if (length(likely_causes) == 1) {
    likely_causes
  } else {
    paste0("Multiple: '", paste0(likely_causes, collapse = "','"), "'")
  }
}


#' Collapse episode-level likely contaminants to person level
#'
assign_likely_contaminants_p <- function(patient_id) {

  episodes_of_person <- edat2[edat2$patient_id == patient_id, ]
  episodes_of_person$likely_contaminants_e <- trimws(
    episodes_of_person$likely_contaminants_e
  )

  # Collect unique contaminants, excluding "No growth" episodes
  likely_contaminants <- trimws(unique(
    episodes_of_person[
      which(episodes_of_person$likely_contaminants_e != "No growth"),
      "likely_contaminants_e",
      drop = TRUE
    ]
  ))

  if (length(likely_contaminants) == 0) {
    "Not reported"
  } else if (length(likely_contaminants) == 1) {
    likely_contaminants[1]
  } else {
    paste0("Multiple: '", paste0(likely_contaminants, collapse = "','"), "'")
  }
}


#' Convert person-level cause string to binary Yes/No
#'
assign_any_causes_p <- function(likely_causes_p) {

  if (likely_causes_p == "Not reported") {
    return("No")
  }
  if (likely_causes_p == "No growth") {
    return("No")
  }
  if (likely_causes_p == "Ungrouped") {
    return("Ungrouped")
  }
  return("Yes")
}


#' Convert person-level contaminant string to binary Yes/No
#'
assign_any_contaminants_p <- function(likely_contaminants_p) {

  # All sentinel / negative values
  if (likely_contaminants_p == "Not reported") {
    return("No")
  }
  if (likely_contaminants_p == "No organism name / culture result") {
    return("No organism name / culture result")
  }
  if (likely_contaminants_p == "No") {
    return("No")
  }
  if (likely_contaminants_p == "No growth") {
    return("No")
  }
  if (likely_contaminants_p == "Ungrouped") {
    return("Ungrouped")
  }
  return("Yes")
}
