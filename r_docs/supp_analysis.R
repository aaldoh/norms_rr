# author: Alaa Aldoh
# contact: az.aldoh@gmail.com/a.aldoh@sussex.ac.uk

# Set up ---------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, emmeans, papaja, bfrr, broom)

no_out <- readRDS("data/no_outliers.rds")
max_FFQ_meat <- round(mad(no_out$FFQ_T1_meat)*3 + median(no_out$FFQ_T1_meat)) # max threshold for univariate outliers on FFQ

complete_cases <- no_out %>%
  mutate(FFQ_T2_meat       = rowSums(across(FFQ_T2_2:FFQ_T2_7)),
         interest_diff     = interest_T2 - interest_T1,
         attitude_diff     = attitude_T2 - attitude_T1,
         intention_diff    = intentionuni_T2 - intentionuni_T1,
         intention_bi_diff = intentionbi_T2 - intentionbi_T1,
         FFQ_meat_diff     = FFQ_T2_meat - FFQ_T1_meat) %>%
  filter(!is.na(no_out$FFQ_T2_7))

complete_no_mad <- filter(complete_cases, FFQ_T1_meat <= max_FFQ_meat)

no_fail <- no_out %>%
  filter(norm_info == "dynamic" & man_check == "More people are eating less meat"| norm_info == "static" & man_check == "People are eating about the same amount of meat"| norm_info == "none") %>%
  mutate(FFQ_T2_meat       = rowSums(across(FFQ_T2_2:FFQ_T2_7)),
         interest_diff     = interest_T2 - interest_T1,
         attitude_diff     = attitude_T2 - attitude_T1,
         intention_diff    = intentionuni_T2 - intentionuni_T1,
         intention_bi_diff = intentionbi_T2 - intentionbi_T1,
         FFQ_meat_diff     = FFQ_T2_meat - FFQ_T1_meat)

no_fail_mad <-filter(no_fail, FFQ_T1_meat <= max_FFQ_meat)
  


# No manipulation check fails ---------------------------
## H1 ---------------------------
contrasts_norm <- list(none_vs_dy = c(-1,0,1),
                       stat_vs_dy = c(0,-1,1))

outcomes_T1 <- c("attitude_T1", "interest_T1", "intentionuni_T1") %>% set_names(.)
h1check.models <- map_df(outcomes_T1, ~ lm(substitute(i ~ norm_info, list(i = as.name(.))), data = no_fail) %>%
                           emmeans(., "norm_info", contr = contrasts_norm, infer = T) %>%
                           .$contrasts %>%
                           summary(), .id = "outcome")

## H2 ---------------------------
h2check.df <- summary(lm(attitude_T1 ~ norm_info*format, data = no_fail)) %>% .$df

h2check.models <- map_df(outcomes_T1, ~ 
                           lm(formula(paste0(.x, "~", "norm_info*format")), data = no_fail) %>%
                           tidy(conf.int=T), .id = "outcome") %>%
  filter(term == "norm_infostatic:formattext") %>%
  cbind(., h2check.df[2]) %>%
  select(1:4, 9, 7:8, 5:6)

## H3 ---------------------------
contrast_interact = list(dytext_vs_dyvis = c(1, -1,  0,  0,  0)) #text vs visual for dynamic only

h3check.models <- map_df(outcomes_T1, ~ lm(formula(paste0(.x, "~", "condition")), data = no_fail) %>%
                           emmeans(., "condition", contr = contrast_interact, infer = T) %>%
                           .$contrasts %>%
                           summary(), .id = "outcome")

## T1 check output ---------------------------
h.tests <- c("stat_vs_dy", "none_vs_dy", "norm_infostatic:formattext", "dytext_vs_dyvis")
T1check_outcomes <- bind_rows(h1check.models,setNames(h2check.models,names(h1check.models)), h3check.models)
T1check_bf <- sapply(1:12, function(x) 
  bfrr(T1check_outcomes$estimate[x], 
       T1check_outcomes$SE[x], 
       sample_df = T1check_outcomes$df[x], 
       model = "normal",
       mean = 0, 
       sd = 5, 
       tail = 1, 
       criterion = 3, 
       rr_interval = list(mean = c(-10, 10), sd = c(0, 10)), 
       precision = 0.05)[-14])
T1check_rr <- sapply(1:12, function(x) paste0("[", toString(T1check_bf[,x]$RR$sd), "]"))
T1check_table <- cbind(T1check_outcomes, unlist(T1check_bf[3,]), T1check_rr, unlist(T1check_bf[5,])) %>% 
  mutate(p.value = printp(p.value)) %>%
  modify_if(., ~is.numeric(.), ~round(., 2)) %>%
  unite(CL, c("lower.CL", "upper.CL"), sep = ", ") %>%
  arrange(match(contrast, h.tests))

write_csv(T1check_table,"outputs/supplementary/T1check_table.csv")

## H4 ---------------------------
outcomes_T2 <- c("attitude_diff", "interest_diff", "intention_diff", "FFQ_meat_diff") %>% set_names(.)
h4_check.models <- map_df(outcomes_T2[1:3], ~lm(substitute(i ~ norm_info, list(i = as.name(.))), data = no_fail) %>%
                           emmeans(., "norm_info", contr = contrasts_norm, infer = T) %>%
                           .$contrasts %>%
                           summary(), .id = "outcome") %>%
  bind_rows(map_df(outcomes_T2[4], ~lm(FFQ_meat_diff ~ norm_info, data = no_fail_mad) %>%
                     emmeans(., "norm_info", contr = contrasts_norm, infer = T) %>%
                     .$contrasts %>%
                     summary(), .id = "outcome"))

## H5 ---------------------------
h5_check.df <- summary(lm(attitude_diff ~ norm_info*format, data = no_fail)) %>% .$df
h5_check.df[4:6] <- summary(lm(attitude_diff ~ norm_info*format, data = no_fail_mad)) %>% .$df

h5_check.models <- map_df(outcomes_T2[1:3], ~ 
                           lm(formula(paste0(.x, "~", "norm_info*format")), data = no_fail) %>%
                           tidy(conf.int=T), .id = "outcome") %>%
  filter(term == "norm_infostatic:formattext") %>%
  cbind(., h5_check.df[2]) %>%
  bind_rows(map_df(outcomes_T2[4], ~ 
                     lm(formula(paste0(.x, "~", "norm_info*format")), data = no_fail_mad) %>%
                     tidy(conf.int=T), .id = "outcome") %>%
              filter(term == "norm_infostatic:formattext") %>%
              cbind(., h5_check.df[5])) %>%
  mutate(`h5_check.df[2]` = coalesce(`h5_check.df[2]`, `h5_check.df[5]`)) %>%
  select(1:4, 9, 7:8, 5:6)


## H6 ---------------------------
h6check.models <- map_df(outcomes_T2[1:3], ~ lm(formula(paste0(.x, "~", "condition")), data = no_fail) %>%
                          emmeans(., "condition", contr = contrast_interact, infer = T) %>%
                          .$contrasts %>%
                          summary(), .id = "outcome") %>%
  bind_rows(map_df(outcomes_T2[4], ~ lm(FFQ_meat_diff ~ condition, data = no_fail_mad) %>%
                     emmeans(., "condition", contr = contrast_interact, infer = T) %>%
                     .$contrasts %>%
                     summary(), .id = "outcome"))

## T2 check output ---------------------------
T2check_outcomes <- bind_rows(h4_check.models,setNames(h5_check.models,names(h4_check.models)), h6check.models) %>%
  arrange(match(outcome, outcomes_T2))

T2check_meat_bf <- sapply(13:16, function(x) 
  bfrr(T2check_outcomes$estimate[x]*-1, 
       T2check_outcomes$SE[x], 
       sample_df = T2check_outcomes$df[x], 
       model = "normal",
       mean = 0, 
       sd = 1, 
       tail = 1, 
       criterion = 3, 
       rr_interval = list(mean = c(-5, 5), sd = c(0, 5)), 
       precision = 0.05)[-14])

T2check_bf <- sapply(1:12, function(x) 
  bfrr(T2check_outcomes$estimate[x], 
       T2check_outcomes$SE[x], 
       sample_df = T2check_outcomes$df[x], 
       model = "normal",
       mean = 0, 
       sd = 5, 
       tail = 1, 
       criterion = 3, 
       rr_interval = list(mean = c(-10, 10), sd = c(0, 10)), 
       precision = 0.05)[-14]) %>%
  cbind(., T2check_meat_bf)

T2check_rr <- sapply(1:16, function(x) paste0("[", toString(T2check_bf[,x]$RR$sd), "]"))
T2check_table <- cbind(T2check_outcomes, unlist(T2check_bf[3,]), T2check_rr, unlist(T2check_bf[5,])) %>% 
  mutate(p.value = printp(p.value)) %>%
  modify_if(., ~is.numeric(.), ~round(., 2)) %>%
  unite(CL, c("lower.CL", "upper.CL"), sep = ", ") %>%
  arrange(match(contrast, h.tests))

write_csv(T2check_table,"outputs/supplementary/T2check_table.csv")

# Complete case analysis ---------------------------
## H4 ---------------------------
h4_comp.models <- map_df(outcomes_T2[1:3], ~lm(substitute(i ~ norm_info, list(i = as.name(.))), data = complete_cases) %>%
                      emmeans(., "norm_info", contr = contrasts_norm, infer = T) %>%
                      .$contrasts %>%
                      summary(), .id = "outcome") %>%
  bind_rows(map_df(outcomes_T2[4], ~lm(FFQ_meat_diff ~ norm_info, data = complete_no_mad) %>%
                     emmeans(., "norm_info", contr = contrasts_norm, infer = T) %>%
                     .$contrasts %>%
                     summary(), .id = "outcome"))

## H5 ---------------------------
h5_comp.df <- summary(lm(attitude_diff ~ norm_info*format, data = complete_cases)) %>% .$df
h5_comp.df[4:6] <- summary(lm(attitude_diff ~ norm_info*format, data = complete_no_mad)) %>% .$df

h5_comp.models <- map_df(outcomes_T2[1:3], ~ 
                      lm(formula(paste0(.x, "~", "norm_info*format")), data = complete_cases) %>%
                      tidy(conf.int=T), .id = "outcome") %>%
  filter(term == "norm_infostatic:formattext") %>%
  cbind(., h5_comp.df[2]) %>%
  bind_rows(map_df(outcomes_T2[4], ~ 
                     lm(formula(paste0(.x, "~", "norm_info*format")), data = complete_no_mad) %>%
                     tidy(conf.int=T), .id = "outcome") %>%
              filter(term == "norm_infostatic:formattext") %>%
              cbind(., h5_comp.df[5])) %>%
  mutate(`h5_comp.df[2]` = coalesce(`h5_comp.df[2]`, `h5_comp.df[5]`)) %>%
              select(1:4, 9, 7:8, 5:6)


## H6 ---------------------------
h6comp.models <- map_df(outcomes_T2[1:3], ~ lm(formula(paste0(.x, "~", "condition")), data = complete_cases) %>%
                      emmeans(., "condition", contr = contrast_interact, infer = T) %>%
                      .$contrasts %>%
                      summary(), .id = "outcome") %>%
  bind_rows(map_df(outcomes_T2[4], ~ lm(FFQ_meat_diff ~ condition, data = complete_no_mad) %>%
                     emmeans(., "condition", contr = contrast_interact, infer = T) %>%
                     .$contrasts %>%
                     summary(), .id = "outcome"))
## complete case output ---------------------------
T2comp_outcomes <- bind_rows(h4_comp.models,setNames(h5_comp.models,names(h4_comp.models)), h6comp.models) %>%
  arrange(match(outcome, outcomes_T2))

T2comp_meat_bf <- sapply(13:16, function(x) 
  bfrr(T2comp_outcomes$estimate[x]*-1, 
       T2comp_outcomes$SE[x], 
       sample_df = T2comp_outcomes$df[x], 
       model = "normal",
       mean = 0, 
       sd = 1, 
       tail = 1, 
       criterion = 3, 
       rr_interval = list(mean = c(-5, 5), sd = c(0, 5)), 
       precision = 0.05)[-14])

T2comp_bf <- sapply(1:12, function(x) 
  bfrr(T2comp_outcomes$estimate[x], 
       T2comp_outcomes$SE[x], 
       sample_df = T2comp_outcomes$df[x], 
       model = "normal",
       mean = 0, 
       sd = 5, 
       tail = 1, 
       criterion = 3, 
       rr_interval = list(mean = c(-10, 10), sd = c(0, 10)), 
       precision = 0.05)[-14]) %>%
  cbind(., T2comp_meat_bf)

T2comp_rr <- sapply(1:16, function(x) paste0("[", toString(T2comp_bf[,x]$RR$sd), "]"))
T2comp_table <- cbind(T2_outcomes, unlist(T2comp_bf[3,]), T2comp_rr, unlist(T2comp_bf[5,])) %>% 
  mutate(p.value = printp(p.value)) %>%
  modify_if(., ~is.numeric(.), ~round(., 2)) %>%
  unite(CL, c("lower.CL", "upper.CL"), sep = ", ") %>%
  arrange(match(contrast, h.tests))

write_csv(T2comp_table,"outputs/supplementary/T2complete_table.csv")

