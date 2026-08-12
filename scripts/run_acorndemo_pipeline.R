# =========================================================
# ACORN Demo Data Processing Pipeline
# =========================================================

# Package loader
dynamic_require <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Required packages
  pkgs <- c("dplyr", "readxl", "AMR")
invisible(lapply(pkgs, dynamic_require))

# =========================================================
# Raw data and output location
# =========================================================
raw_data_location_all = "./data/raw/acorndemo_data.xlsx"

output_data_location_persons = "./data/output/acorndemo_persons.csv"
output_data_location_episodes = "./data/output/acorndemo_episodes.csv"
output_data_location_organisms = "./data/output/acorndemo_organisms.csv"
output_data_location_specimens = "./data/output/acorndemo_specimens.csv"

# =========================================================
# Load Organism Pathogenicity Grouping Data
# =========================================================
organism_patho <- read.csv("algorithms/patho_grouping_v1.0_20260809.csv")

# --- AMR taxonomy integration ---
# Resolve each curated organism name to a stable mo code.
# The mo code becomes the primary join key — never match on free-text.
# Unresolved entries (UNKNOWN) keep their original Org_name for the fallback path.
organism_patho <- organism_patho %>%
  mutate(
    mo_code = suppressMessages(as.character(as.mo(Org_name, keep_synonyms = FALSE))),
    mo_current_name = suppressMessages(mo_name(mo_code, keep_synonyms = FALSE))
  )

# Load project specific functions
source("scripts/fcns_patho_grouping_v1.R")

# =========================================================
# Initial ETL
# =========================================================

source("./scripts/initial_acorndemo_data_extract.R")   # Select relevant records and columns
source("./scripts/initial_acorndemo_data_transform.R") # Apply transformations
source("./scripts/initial_acorndemo_data_load.R")      # Write processed data (CSV)

