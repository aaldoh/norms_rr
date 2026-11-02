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

# Remove the original Prolific ID. Qualtrics also stores it from the survey URL
# in a separate column, PROLIFIC_PID, which must not be shared.
export_T1 <- select(export_T1, -any_of("PROLIFIC_PID"))
export_T2 <- select(export_T2, -any_of("PROLIFIC_PID"))

# Safety check: stop if any real Prolific IDs (24-character hex strings) remain
has_real_ids <- function(d) any(vapply(d, function(x) any(grepl("^[0-9a-f]{24}$", x)), logical(1)))
if (has_real_ids(export_T1) || has_real_ids(export_T2)) {
  stop("Real Prolific IDs still present - do not save or share these data.")
}

saveRDS(export_T1, file = "data/deidentified_T1.rds")
saveRDS(export_T2, file = "data/deidentified_T2.rds")
