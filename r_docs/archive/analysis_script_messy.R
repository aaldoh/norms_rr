# author: Alaa Aldoh
# contact: az.aldoh@gmail.com/a.aldoh@sussex.ac.uk

# Set up ---------------------------
if (!require("devtools")) install.packages("devtools")
if (!require("papaja")) devtools::install_github("crsh/papaja")
if (!require("bfrr")) devtools::install_github("debruine/bfrr")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(qualtRics, here, tidyverse, emmeans, papaja, bfrr, ggfortify, mice, kableExtra)

set.seed(1234)

export_T1 <- qualtRics::read_survey(here("data/export_T1.csv"))

# Data cleaning ---------------------------
## structuring data
clean_T1 <- export_T1 %>%
  dplyr::rename_with(., ~gsub("\\s|_1$", "", .x), .cols = interest_T1_1:politics_1) %>%
  dplyr::mutate(across(starts_with("FFQ"), 
                       ~recode(., "Never" = 0, "Once (during those 7 days)" = 1, 
                               "Twice" = 2, "Three times" = 3, "Four times" = 4, 
                               "Five times" = 5, "Six times" = 6, "Seven times" = 7,
                               "2 times per day" = 14, "3 times per day" = 21, 
                               "4 or more times per day" = 28))) %>%
  dplyr::mutate(cons_proj   = rowMeans(dplyr::across(cons_now_perc:cons_six_perc), na.rm = T),
                cons_proj_c = cons_proj - mean(cons_proj, na.rm = T),
                FFQ_T1_meat = rowSums(across(FFQ_2:FFQ_7)),
                genderbi    = na_if(gender, "Other"),
                age_c       = age - mean(age, na.rm = T),
                prob_time   = ifelse(`norm_time_Page Submit` < 5 & condition != "none", 1, 0) %>% as_factor,
                dplyr::across(where(is.character), as_factor),
                condition   = fct_relevel(condition, "dy_graph", "dy_text", "st_graph", "st_text", "none"),
                norm_info   = fct_relevel(norm_info, "none", "static", "dynamic"),
                format      = fct_relevel(format, "none", "text", "visual")) %>%
  droplevels()

## exclusions
complete <- filter(clean_T1, !is.na(veg)) #1) filter out incomplete responses
no_veg <- filter(clean_T1, veg != "Yes") # 2) filter out vegetarians
no_care <- filter(no_veg, prob_time == 0) # 3) filter out careless responding

library(MASS) # 4) removing multivariate outliers
mcd     <- cov.mcd(no_care[c(35:37, 51)], quantile.used = nrow(no_care[c(35:37, 51)])*.75)
mcd_mah <- mahalanobis(no_care[c(35:37, 51)], mcd$center,mcd$cov)
cutoff  <- qchisq(p = 0.99, df = ncol(no_care[c(35:37, 51)]))
no_out  <- no_care[mcd_mah <= cutoff, ]
detach("package:MASS", unload=TRUE)

max_FFQ_meat <- round(mad(no_out$FFQ_T1_meat)*3 + median(no_out$FFQ_T1_meat)) # max threshold for univariate outliers on FFQ

saveRDS(no_out, file = "data/no_outliers.rds")
write_csv(complete,"data/complete.csv")


# Data inspection ---------------------------
## randomization check
random_check <- list(age_stat = apa_print(aov(age ~ norm_info, no_out)), ## age
                     pol_stat = apa_print(aov(politics ~ norm_info, no_out)), ## political position
                     gender_stat = chisq.test(no_out$norm_info, no_out$gender), n = nrow(no_out), ## gender
                     nation_stat = chisq.test(no_out$norm_info, no_out$country), n = nrow(no_out), ## nation
                     meat_T1_check = apa_print(aov(FFQ_T1_meat ~ norm_info, no_out))) ## baseline meat consumption

## manipulation checks
man_checks <- list(man_desc  = round(100*prop.table(table(no_out$norm_info, no_out$man_check), margin = 1),2),
                   man.chitest = chisq.test(no_out$norm_info, no_out$man_check), n = nrow(no_out))

# Data overview ---------------------------
participants <- c(export = nrow(export_T1),
                  complete = nrow(complete),
                  no_veg = nrow(no_veg),
                  no_care = nrow(no_care),
                  no_out = nrow(no_out))
age_desc      <- sapply(no_out[,"age"], function(x) c(min = min(x), max = max(x), avg = mean(x), sd = sd(x))) ## age
gender_freq   <- 100 * prop.table(table(no_out$gender)) ## gender
baseline_meat <-  group_by(no_out, condition) %>% summarise(m = mean(FFQ_T1_meat), sd = sd(FFQ_T1_meat)) ## meat consumption at T1

# T1 overview ---------------------------
## visualization
T1_violin <- no_out %>% ## VIOLIN PLOTS + MEAN/CI
  pivot_longer(cols = c("attitude_T1", "interest_T1", "intentionuni_T1"), names_to = "variable", values_to = "value") %>%
  ggplot2::ggplot(aes(norm_info, value, colour = format)) +
  facet_wrap(~ variable, labeller = as_labeller(c('attitude_T1'="Attitude", 'interest_T1'="Interest", 'intentionuni_T1'="Intention"))) +
  geom_violin(trim = T,  position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_normal", position = position_dodge(width = 0.9)) +
  coord_cartesian(ylim = c(0,100)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  labs(y = "Value (%)", x = "Norm type", colour = "Format", title = "Distribution of outcomes by condition at T1") + 
  scale_x_discrete(labels= c("None", "Static", "Dynamic")) +
  scale_color_discrete(labels= c("None", "Text", "Visual")) +
  papaja::theme_apa()

T1_jitter <- no_out %>% #interaction and simple effects WITH JITTER
  filter(condition != "none") %>%
  pivot_longer(cols = c("attitude_T1", "interest_T1", "intentionuni_T1"), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = norm_info, color = format, group = format, y = value)) +
  facet_wrap(~ variable, labeller = as_labeller(c('attitude_T1'="Attitude", 'interest_T1'="Interest", 'intentionuni_T1'="Intention"))) +
  geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.3)+
  stat_summary(fun = mean, geom = "line", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_se,
               geom = "pointrange", position = position_dodge(width = 0.9)) +
  labs(y = "Value (%)", x = "Norm type", colour = "Format")+
  scale_x_discrete(labels= c("Static", "Dynamic")) +
  scale_color_discrete(labels= c("Text", "Visual")) +
  papaja::theme_apa()

T1_interact <- no_out %>% #interaction and simple effects PLOT WITHOUT JITTER
  filter(condition != "none") %>%
  select(c("norm_info", "format", "attitude_T1", "interest_T1", "intentionuni_T1")) %>%
  pivot_longer(cols = ends_with("_T1"), names_to = "outcome", values_to = "value") %>%
  group_by(norm_info, format, outcome) %>%
  summarise(n= n(),
            m = mean(value),
            se = sd(value)/sqrt(n())) %>%
  ggplot(aes(x = norm_info, color = format, group = format, y = m, ymin = m-se, ymax = m+se)) +
  facet_wrap(~ outcome, labeller = as_labeller(c('attitude_T1'="Attitude", 'interest_T1'="Interest", 'intentionuni_T1'="Intention"))) +
  geom_line(aes(linetype = format, group = format), position = position_dodge(width = 0.9), show.legend = FALSE) +
  geom_point(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = m-se, ymax = m+se, group = format), width = 0.2, position = position_dodge(width = 0.9)) +
  coord_cartesian(ylim = c(0,100)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  labs(y = "Value (%)", x = "Norm type", colour = "Format")+
  scale_x_discrete(labels= c("Static", "Dynamic")) +
  scale_color_discrete(labels= c("Text", "Visual")) +
  papaja::theme_apa()

## descriptive statistics
outcomes_T1_desc <- group_by(no_out, condition) %>% ## outcomes by each condition
  summarise(n = n(),
            across(c("attitude_T1", "interest_T1", "intentionuni_T1"), 
                   list(mean = mean, sd = sd, se = ~sd(.)/sqrt(n())), 
                   .names = "{.fn}_{.col}"), .groups = "rowwise")

## H1 ---------------------------
contrasts_norm <- list(none_vs_dy = c(-1,0,1),
                       stat_vs_dy = c(0,-1,1))

outcomes_T1 <- c("attitude_T1", "interest_T1", "intentionuni_T1") %>% set_names(.)
h1.models <- map_df(outcomes_T1, ~ lm(substitute(i ~ norm_info, list(i = as.name(.))), data = no_out) %>%
                      emmeans(., "norm_info", contr = contrasts_norm, infer = T) %>%
                      .$contrasts %>%
                      summary(), .id = "outcome")

## H2 ---------------------------
h2.df <- summary(lm(attitude_T1 ~ norm_info*format, data = no_out)) %>% .$df

h2.models <- map_df(outcomes_T1, ~ 
                      lm(formula(paste0(.x, "~", "norm_info*format")), data = no_out) %>%
                      tidy(conf.int=T), .id = "outcome") %>%
  filter(term == "norm_infostatic:formattext") %>%
  cbind(., h2.df[2]) %>%
  select(1:4, 9, 7:8, 5:6)

## H3 ---------------------------
contrast_interact = list(dytext_vs_dyvis = c(1, -1,  0,  0,  0)) #text vs visual for dynamic only

h3.models <- map_df(outcomes_T1, ~ lm(formula(paste0(.x, "~", "condition")), data = no_out) %>%
                      emmeans(., "condition", contr = contrast_interact, infer = T) %>%
                      .$contrasts %>%
                      summary(), .id = "outcome")

## T1 output ---------------------------
T1_outcomes <- bind_rows(h1.models,setNames(h2.models,names(h1.models)), h3.models)
T1_bf <- sapply(1:12, function(x) 
  bfrr(T1_outcomes$estimate[x], 
       T1_outcomes$SE[x], 
       sample_df = T1_outcomes$df[x], 
       model = "normal",
       mean = 0, 
       sd = 5, 
       tail = 1, 
       criterion = 3, 
       rr_interval = list(mean = c(-10, 10), sd = c(0, 10)), 
       precision = 0.05)[-14])
T1_rr <- sapply(1:12, function(x) paste0("HN[", toString(T1_bf[,x]$RR$sd), "]"))
T1_table <- cbind(T1_outcomes, unlist(T1_bf[3,]), T1_rr, unlist(T1_bf[5,])) %>% 
  mutate(p.value = printp(p.value)) %>%
  modify_if(., ~is.numeric(.), ~round(., 2))

write_csv(T1_table,"outputs/T1_table.csv")
save(participants, T1_table, T1_violin, T1_jitter, T1_interact, file =  "outputs/T1_out.RData")

saveRDS(imputed_data, "imputed_data.rds")

# Data merge ---------------------------
export_T2 <- qualtRics::read_survey(here("data/export_T2.csv"))

## cleaning combined dataset
clean_T2 <- export_T2 %>%
  dplyr::rename_with(., ~gsub("\\s|_1$", "", .x), .cols = interest_T2_1:intentionbi_T2_1)

merged_tib <- merge(no_out, clean_T2, by= "PROLIFIC_ID", all.x = TRUE, all.y = FALSE) %>%
  select(-ends_with(c(".x", ".y")), -starts_with("norm_time")) %>%
  mutate(across(starts_with("FFQ_T2"), 
                ~recode(., "Never" = 0, "Once (during those 7 days)" = 1, 
                        "Twice" = 2, "Three times" = 3, "Four times" = 4, 
                        "Five times" = 5, "Six times" = 6, "Seven times" = 7,
                        "2 times per day" = 14, "3 times per day" = 21, 
                        "4 or more times per day" = 28))) 

apply(merged_tib, 2, function(col)sum(is.na(col))/length(col)) #check % missing from all columns

merged_tib$complete <- as.numeric(complete.cases(merged_tib[, c("FFQ_T2_7")])) #create variable indicating if case in incomplete

## Multiple imputation ---------------------------
### variables predicting missingness
xtabs(~condition+complete, data = merged_tib) # condition x missingness
xtabs(~gender+complete, data = merged_tib) # gender x missingness

aux_logit <- glm(complete ~ condition+age+genderbi+politics, data = merged_tib, family = "binomial") #check if auxiliary variables predict missingness
summary(aux_logit)


T2_vars_imp <- c("interest_T2", "attitude_T2", "intentionuni_T2", "intentionbi_T2", 
                 "FFQ_T2_1", "FFQ_T2_2", "FFQ_T2_3", "FFQ_T2_4", "FFQ_T2_5", "FFQ_T2_6", 
                 "FFQ_T2_7")

exclude_pred <- c("nationality", "norm_info", "condition", "PROLIFIC_ID", "Consent", "veg", "prob_time", "think_reduce", 
                  "why_think", "active_reduce", "cons_proj", "cons_proj_c", "genderbi", "age_c", "FFQ_T1_meat", T2_vars_imp,
                  "FFQ_T2_8", "FFQ_T2_9",	"FFQ_T2_10",	"FFQ_T2_11",	"FFQ_T2_12",	"FFQ_T2_13",	"FFQ_T2_14",	"FFQ_T2_15")

exclude_outcome <- c("PROLIFIC_ID", "Consent","nationality","country","resident","FFQ_1","FFQ_2","FFQ_3","FFQ_4","FFQ_5","FFQ_6","FFQ_7",          
                     "FFQ_8","FFQ_9","FFQ_10","FFQ_11","FFQ_12","FFQ_13","FFQ_14","FFQ_15","interest_T1","attitude_T1","intentionuni_T1", "intentionbi_T1", 
                     "cons_now_perc",   "cons_next_perc",  "cons_six_perc",   "man_check","politics","age","gender","veg","condition","format","norm_info","cons_proj",      
                     "cons_proj_c","FFQ_T1_meat","genderbi","age_c","prob_time")

init <- mice(merged_tib, maxit = 0)
meth <- init$meth
pred <- init$pred
meth[!names(meth) %in% T2_vars_imp] <- ""
pred[, !names(merged_tib) %in% names(clean_T2)] <- 1
pred[, exclude_pred] <- 0

###### split by group
imputed <- mice(merged_tib, m = 2, maxit = 5, ridge=0.01, meth = meth, group = "norm_info")
imputed_long <- mice::complete(imputed, "long", include = T)

imputed_calc <- imputed_long %>%
  mutate(FFQ_T2_meat = rowSums(across(FFQ_T2_2:FFQ_T2_7)),
         interest_diff = interest_T2 - interest_T1,
         attitude_diff = attitude_T2 - attitude_T1,
         intention_diff = intentionuni_T2 - intentionuni_T1,
         intention_bi_diff = intentionbi_T2 - intentionbi_T1,
         FFQ_meat_diff = FFQ_T2_meat - FFQ_T1_meat)

imputed_dat <- as.mids(imputed_calc)


# T2 overview ---------------------------
## visualization
final_tib <- merged_tib %>%
  mutate(FFQ_T2_meat = rowSums(across(FFQ_T2_2:FFQ_T2_7)),
         interest_diff = interest_T2 - interest_T1,
         attitude_diff = attitude_T2 - attitude_T1,
         intention_diff = intentionuni_T2 - intentionuni_T1,
         intention_bi_diff = intentionbi_T2 - intentionbi_T1,
         FFQ_meat_diff = FFQ_T2_meat - FFQ_T1_meat)

T2_violin <- final_tib %>%
  pivot_longer(cols = c("attitude_diff", "interest_diff", "intention_diff"), names_to = "variable", values_to = "value") %>%
  ggplot2::ggplot(aes(norm_info, value, colour = format)) +
  facet_wrap(~ variable, labeller = as_labeller(c('attitude_diff'="Attitude", 'interest_diff'="Interest", 'intention_diff'="Intention"))) +
  geom_violin(trim = T,  position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_normal", position = position_dodge(width = 0.9)) +
  labs(y = "Value (%)", x = "Norm type", colour = "Format", title = "Distribution of change in outcomes by condition at T2") + 
  scale_x_discrete(labels= c("None", "Static", "Dynamic")) +
  scale_color_discrete(labels= c("None", "Text", "Visual")) +
  papaja::theme_apa()

T2_interact <- final_tib %>% #interaction and simple effects PLOT WITHOUT JITTER
  filter(condition != "none") %>%
  select(c("norm_info", "format", "attitude_diff", "interest_diff", "intention_diff")) %>%
  pivot_longer(cols = ends_with("_diff"), names_to = "outcome", values_to = "value") %>%
  group_by(norm_info, format, outcome) %>%
  summarise(n= n(),
            m = mean(value, na.rm = T),
            se = sd(value, na.rm =T)/sqrt(n())) %>%
  ggplot(aes(x = norm_info, color = format, group = format, y = m, ymin = m-se, ymax = m+se)) +
  facet_wrap(~ outcome, labeller = as_labeller(c('attitude_diff'="Attitude", 'interest_diff'="Interest", 'intention_diff'="Intention"))) +
  geom_line(aes(linetype = format, group = format), position = position_dodge(width = 0.9), show.legend = FALSE) +
  geom_point(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = m-se, ymax = m+se, group = format), width = 0.2, position = position_dodge(width = 0.9)) +
  coord_cartesian(ylim = c(-5,5)) +
  #scale_y_continuous(breaks = seq(0, 100, 10)) +
  labs(y = "Value (%)", x = "Norm type", colour = "Format")+
  scale_x_discrete(labels= c("Static", "Dynamic")) +
  scale_color_discrete(labels= c("Text", "Visual")) +
  papaja::theme_apa()


outcomes_change <- final_tib %>% 
  select(c("condition", ends_with("_T1"), ends_with("_T2"))) %>%
  pivot_longer(cols = ends_with(c("_T1", "_T2")), 
               names_to = c("outcome", "time"), 
               values_to = "value",
               names_sep = "_") %>%
  mutate(outcome = factor(outcome, levels=c('attitude','interest','intentionuni','intentionbi'))) %>%
  group_by(condition, time, outcome) %>%
  summarise(n= n(),
            m = mean(value, na.rm = T),
            se = sd(value, na.rm =T)/sqrt(n())) %>%
  ggplot(aes(x = time, color = condition, group = condition, y = m, ymin = m-se, ymax = m+se)) +
  facet_wrap(~ outcome, labeller = as_labeller(c('attitude'="Attitude", 'interest'="Interest", 'intentionuni'="Intention (unipolar)", 'intentionbi' ="Intention (bipolar/exploratory)"))) +
  geom_line(aes(group = condition), show.legend = T) +
  scale_color_discrete(labels= c("Visual dynamic", "Text dynamic", "Visual static", "Text static", "None")) +
  labs(y = "Value (%)", x = "Time", colour = "Condition")+
  papaja::theme_apa()

FFQ_change <- final_tib %>%
  select(c("condition", "FFQ_T1_meat", "FFQ_T2_meat")) %>%
  pivot_longer(cols = contains(c("_T1", "_T2")), 
               names_to = c("outcome", "time"), 
               values_to = "value",
               names_sep = "_") %>%
  group_by(condition, time) %>%
  summarise(n= n(),
            m = mean(value, na.rm = T),
            se = sd(value, na.rm =T)/sqrt(n())) %>%
  ggplot(aes(x = time, color = condition, group = condition, y = m, ymin = m-se, ymax = m+se)) +
  geom_line(aes(group = condition), show.legend = T) +
  scale_color_discrete(labels= c("Visual dynamic", "Text dynamic", "Visual static", "Text static", "None")) +
  labs(y = "Value (%)", x = "Time", colour = "Condition")+
  papaja::theme_apa()


## descriptive statistics
change_by_cond <- group_by(final_tib, condition) %>%
  summarise(n = n(),
            across(c("attitude_diff", "interest_diff","intention_diff", "FFQ_meat_diff"), 
                   list(mean = ~mean(., na.rm = TRUE), 
                        sd = ~sd(., na.rm = TRUE), 
                        se = ~sd(., na.rm = TRUE)/sqrt(n())), 
                   .names = "{.fn}_{.col}"), .groups = "rowwise")

diffs_by_cond <- group_by(final_tib, condition) %>%
  summarise(n = n(),
            across(c("attitude_T2", "interest_T2","intentionuni_T2"), 
                   list(mean = ~mean(., na.rm = TRUE), 
                        sd = ~sd(., na.rm = TRUE), 
                        se = ~sd(., na.rm = TRUE)/sqrt(n())), 
                   .names = "{.fn}_{.col}"), .groups = "rowwise")

## H4 ---------------------------
outcomes_T2 <- c("attitude_diff", "interest_diff", "intention_diff", "FFQ_meat_diff") %>% set_names(.)
h4.models <- map_df(outcomes_T2, ~ with(imputed_dat , lm(formula(paste0(.x, "~", "norm_info"))))  %>%
                     emmeans(., "norm_info", contr = contrasts_norm, infer = T) %>%
                     .$contrasts %>%
                     summary(), .id = "outcome")



## H5 ---------------------------
h5.df <- summary(pool(with(imputed_dat, (lm(attitude_diff ~ norm_info*format))))) %>% .$df

h5.models <- map_df(outcomes_T2, ~ 
                      with(imputed_dat, lm(formula(paste0(.x, "~", "norm_info*format")), data = no_out)) %>%
                      pool() %>%
                      tidy(conf.int=T), .id = "outcome") %>%
  filter(term == "norm_infostatic:formattext") %>%
  cbind(., h5.df[2]) %>%
  select(1:4, 11, 7:8, 5:6)


## H6 ---------------------------
h6.models <- map_df(outcomes_T2, ~ 
                      with(imputed_dat, lm(formula(paste0(.x, "~", "condition")))) %>%
                      emmeans(., "condition", contr = contrast_interact, infer = T) %>%
                      .$contrasts %>%
                      summary(), .id = "outcome")


## T2 output ---------------------------
h.tests <- c("stat_vs_dy", "none_vs_dy", "norm_infostatic:formattext", "dytext_vs_dyvis")
T2_outcomes <- bind_rows(h4.models,setNames(h5.models,names(h4.models)), h6.models) %>%
  arrange(match(outcome, outcomes_T2), match(contrast, h.tests))

T2_bf <- sapply(1:12, function(x) 
  bfrr(T2_outcomes$estimate[x], 
       T2_outcomes$SE[x], 
       sample_df = T2_outcomes$df[x], 
       model = "normal",
       mean = 0, 
       sd = 5, 
       tail = 1, 
       criterion = 3, 
       rr_interval = list(mean = c(-10, 10), sd = c(0, 10)), 
       precision = 0.05)[-14])

T2_meat_bf <- sapply(13:16, function(x) 
  bfrr(T2_outcomes$estimate[x]*-1, 
       T2_outcomes$SE[x], 
       sample_df = T2_outcomes$df[x], 
       model = "normal",
       mean = 0, 
       sd = 1, 
       tail = 1, 
       criterion = 3, 
       rr_interval = list(mean = c(-5, 5), sd = c(0, 5)), 
       precision = 0.05)[-14])

T2_rr <- sapply(1:12, function(x) paste0("HN[", toString(T2_bf[,x]$RR$sd), "]"))
T2_table <- cbind(T2_outcomes, unlist(T2_bf[3,]), T2_rr, unlist(T2_bf[5,])) %>% 
  mutate(p.value = printp(p.value)) %>%
  modify_if(., ~is.numeric(.), ~round(., 2))

write_csv(T2_table,"outputs/T2_table.csv")
save(T2_table, file =  "outputs/T2_out.RData")


# Correct estimates ---------------------------
desc <- no_out %>%
  group_by(norm_info) %>%
  summarise(n = n(),
            m_att = mean(attitude_T1),
            sd_att = sd(attitude_T1),
            se_att = sd(attitude_T1)/sqrt(n())) 

desc2 <- no_out %>%
  group_by(format) %>%
  summarise(n = n(),
            m_att = mean(attitude_T1),
            sd_att = sd(attitude_T1),
            se_att = sd(attitude_T1)/sqrt(n())) 

desc3 <- no_out %>%
  group_by(condition) %>%
  summarise(m_att = mean(attitude_T1),
            se_att = sd(attitude_T1)/sqrt(n()))

correct_est <- c(stat_vs_dy = desc$m_att[1] - desc$m_att[2],
                 text_vs_vis = desc2$m_att[3] - desc2$m_att[1],
                 interact = (desc3[1,2] - desc3[2,2]) - (desc3[3,2] - desc3[4,2]),
                 simple_eff = desc3[1,2] - desc3[2,2])

correct_se <- c(stat_vs_dy = sqrt((desc$se_att[2]^2)+(desc$se_att[1]^2)),
                text_vs_vis = sqrt((desc2$se_att[3]^2)+(desc2$se_att[1]^2)),
                interact = sqrt((desc3$se_att[1] - desc3$se_att[2]) - (desc3$se_att[3] - desc3$se_att[4]))) ### this is wrong
