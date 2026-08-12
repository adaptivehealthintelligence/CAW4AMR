#==================================
##### Run steps
#- Apply assign_likely_cause and  assign_likely_contaminant to odat0
#- Apply assign_likely_causes_s and  assign_likely_contaminants_s to sdat0
#- Apply assign_likely_causes_e and  assign_likely_contaminants_e to edat0
#- Apply assign_likely_causes_p and  assign_likely_contaminants_p to pdat0
#==================================

#==================================
##### Output:
##odat2: new columns created: likely_cause and likely_contaminant
##sdat2: new columns created: likely_causes_s and likely_contaminants_s
##edat2: new columns created: likely_causes_e and likely_contaminants_e
##pdat2: new columns created: likely_causes_p and likely_contaminants_p
#==================================

#==================================
##assign_likely_cause: To run this function, we need both edat0 and odat0, and we go through each organism row in odat0 and apply this function to assign ‘likely
## cause’ based one both the pathogenicity grouping and medical device status at the point of episode management. 
## The output will be odat1 table, which will have the ‘likely_cause’ column.
#==================================

# specimen_id variable creation 
odat0$specimen_id <- paste0(odat0$patient_id, "_", odat0$specid) 
odat0$specid <- NULL

# specimen table creation
sdat0 <- data.frame(specimen_id = unique(odat0$specimen_id))


odat1 <- odat0 %>%   
  rowwise() %>% 
  mutate(likely_cause = assign_likely_cause(orgname, episode_id)) %>%  
  ungroup() 

#==================================
##### assign_likely_contaminant: To run this function, we need odat1 
## and we go through each organism row depending on previously assigned likely cause status. 
## The output will be odat2 table, which will have the ‘likely_cause’ column.
#==================================

odat2 <- odat1 %>%   
  rowwise() %>% 
  mutate(likely_contaminant = assign_likely_contaminant(likely_cause)) %>%  
  ungroup() 

#==================================
## assign_likely_causes_s: To run this function, we need sdat0 (the specimen table) and odat1 (the organism table). 
# For each specimen_id in the organism table, this function identifies organisms marked as a likely cause.  
# The output will be sdat1 table, which will have the ‘likely_causes_s’ column. 
#==================================

sdat1 <- sdat0 %>%   
  rowwise() %>% 
  mutate(likely_causes_s = assign_likely_causes_s(specimen_id)) %>%  
  ungroup() 

#==================================
# assign_likely_contaminants_s: 
#To run this function, we need sdat1 and odat1. 
# For each specimen_id in the organism table, this function summarizes all organisms marked as likely contaminants. 
# The output will be sdat2 table, which will have the ‘likely_contaminants_s’ column. 
#==================================

sdat2 <- sdat1 %>%   
  rowwise() %>% 
  mutate(likely_contaminants_s = assign_likely_contaminants_s(specimen_id)) %>%  
  ungroup() 

#==================================
# assign_likely_causes_e: 
# To run this function, we need edat0 and odat1. 
# For each episode_id in the organism table, this function summarizes all organisms marked as likely causes. 
# The output will be edat1 table, which will have the ‘likely_causes_e’ column. 
#==================================

edat1 <- edat0 %>% 
  rowwise() %>% 
  mutate( 
    likely_causes_e = assign_likely_causes_e(episode_id) 
  ) %>% 
  ungroup() 

#==================================
#assign_likely_contaminants_e: 
#To run this function, we need edat1 and odat1. 
#For each episode_id in the organism table, this function summarizes all organisms marked as likely contaminants. 
#The output will be edat2 table, which will have the ‘likely_contaminants_e’ column.
#==================================

edat2 <- edat1 %>% 
  rowwise() %>% 
  mutate( 
    likely_contaminants_e = assign_likely_contaminants_e(episode_id) 
  ) %>% 
  ungroup() 

#==================================
#assign_likely_causes_p: 
# To run this function, we need pdat1 and edat1. 
#For each patient_id in the organism table, this function summarizes all organisms marked as likely contaminants. 
# The output will be pdat2 table, which will have the ‘likely_contaminants_p’ column.
#==================================

pdat1 <- pdat0 %>%   
  rowwise() %>% 
  mutate(likely_causes_p = assign_likely_causes_p(patient_id)) %>%  
  ungroup() 

#==================================
# assign_likely_contaminants_p: 
# To run thi sfunction, we need pdat0 and edat1. 
# For each patient_id in the organism table, this function summarizes all organisms marked as likely causes. 
# The output will be pdat1 table, which will have the ‘likely_causes_p’ column.
#==================================

pdat2 <- pdat1 %>%   
  rowwise() %>% 
  mutate(likely_contaminants_p = assign_likely_contaminants_p(patient_id)) %>%  
  ungroup() 
