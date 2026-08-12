# =========================================================
# Load Raw Tables
# =========================================================

# Full data from ACORN's REDCap CRF
rawdat.all   <- read_excel(raw_data_location_all, sheet = "redcap_f01f05_dta", col_types = "text")

# ACORN dashboard-linked episodes
rawdat.linked <- read_excel(raw_data_location_all, sheet = "acorn_dta", col_types = "text")


# =========================================================
# Define Relevant Variable Sets
# =========================================================
  
relevantOrganisms <- c(
  'patient_id'
  ,'episode_id'
  ,'specid'
  ,'orgname'
)

relevantEpisodes <- c(
  'patient_id'
  ,'episode_id'
  ,'med_device_vent'
  ,'med_device_pcv'
  ,'med_device_cvc'
)

relevantPersons <- c(
  'patient_id'
)

# =========================================================
# Create Initial Data Frames
# =========================================================
odat0      <- rawdat.linked[, relevantOrganisms]
edat0      <- rawdat.all[, relevantEpisodes]
pdat0      <- rawdat.all[, relevantPersons] %>% distinct(patient_id, .keep_all = TRUE)

