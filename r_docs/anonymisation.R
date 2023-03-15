if (!require("stringi")) install.packages("stringi")
library(tidyverse)
library(here)
library(qualtRics)

# creating 1517 pseudonyms to replace 1517 unique Prolific IDs in existing datasets
#stri_rand_strings(1517, 10) %>% write.csv(., "data/pseudonyms.csv")

export_T1 <- qualtRics::read_survey(here("confidential/export_T1.csv"))
export_T2 <- qualtRics::read_survey(here("confidential/export_T2.csv"))
pseudo    <- read_csv(here("confidential/participant_key.csv"))

export_T1$PROLIFIC_ID <- pseudo$pseudonym[match(export_T1$PROLIFIC_ID, pseudo$PROLIFIC_ID)]
export_T2$PROLIFIC_ID <- pseudo$pseudonym[match(export_T2$PROLIFIC_ID, pseudo$PROLIFIC_ID)]

saveRDS(export_T1, file = "data/deidentified_T1.rds")
saveRDS(export_T2, file = "data/deidentified_T2.rds")
