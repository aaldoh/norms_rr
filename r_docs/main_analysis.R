# author: Alaa Aldoh
# contact: az.aldoh@gmail.com/a.aldoh@uva.nl

# Set up ---------------------------
if (!require("devtools")) install.packages("devtools")
if (!require("papaja")) devtools::install_github("crsh/papaja")
if (!require("bfrr")) devtools::install_github("debruine/bfrr")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(qualtRics, here, tidyverse, emmeans, papaja, bfrr, ggfortify, mice, kableExtra)

set.seed(1234)

deidentified_T1 <- readRDS("data/deidentified_T1.rds")
deidentified_T2 <- readRDS("data/deidentified_T2.rds")

combined_tib <- merge(deidentified_T1, deidentified_T2, by= "PROLIFIC_ID", all.x = TRUE, all.y = FALSE) 

# Data cleaning ---------------------------
## structuring data
clean <- combined_tib %>%
  dplyr::rename_with(., ~gsub("\\s|_1$", "", .x), .cols = c(interest_T1_1:politics_1, interest_T2_1:intentionbi_T2_1)) %>%
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
                read_time   = `norm_time_Page Submit`,
                dplyr::across(where(is.character), as_factor),
                condition   = fct_relevel(condition, "dy_graph", "dy_text", "st_graph", "st_text", "none"),
                norm_info   = fct_relevel(norm_info, "none", "static", "dynamic"),
                format      = fct_relevel(format, "none", "text", "visual"),
                complete    = as.numeric(complete.cases(.[, c("FFQ_T2_7")]))) %>%
  select(-ends_with(c(".x", ".y")), -starts_with("norm_time")) %>%
  droplevels()

## exclusions
no_dup <- ungroup(filter(group_by(clean, PROLIFIC_ID), n() < 2)) # 1) filter out participants who completed study > once
complete <- filter(no_dup, !is.na(veg)) # 2) filter out incomplete responses
no_veg <- filter(complete, veg != "Yes") %>% select(-veg) # 3) filter out vegetarians
no_care <- filter(no_veg, prob_time == 0) # 4) filter out careless responding

library(MASS) # 4) removing multivariate outliers
mcd     <- cov.mcd(no_care[c(21:23, 57)], quantile.used = nrow(no_care[c(21:23, 57)])*.75)
mcd_mah <- mahalanobis(no_care[c(21:23, 57)], mcd$center,mcd$cov)
cutoff  <- qchisq(p = 0.99, df = ncol(no_care[c(21:23, 57)]))
no_out  <- no_care[mcd_mah <= cutoff, ]
detach("package:MASS", unload=TRUE)

max_FFQ_meat <- round(mad(no_out$FFQ_T1_meat)*3 + median(no_out$FFQ_T1_meat)) # max threshold for univariate outliers on FFQ

saveRDS(clean, file = "data/clean.rds")
saveRDS(no_out, file = "data/no_outliers.rds")
write_csv(combined_tib,"data/01_merged.csv")
write_csv(clean,"data/02_clean.csv")
write_csv(no_dup,"data/03_unique.csv")
write_csv(complete,"data/04_complete.csv")

# Data inspection ---------------------------
## randomization check
random_check <- list(age_stat = apa_print(aov(age ~ norm_info, no_out)), ## age
                     pol_stat = apa_print(aov(politics ~ norm_info, no_out)), ## political position
                     gender_stat = chisq.test(no_out$norm_info, no_out$gender), ## gender
                     nation_stat = chisq.test(no_out$norm_info, no_out$country), ## nation
                     meat_T1_check = apa_print(aov(FFQ_T1_meat ~ norm_info, no_out))) ## baseline meat consumption

## manipulation checks
man_checks <- list(man_desc  = round(100*prop.table(table(no_out$norm_info, no_out$man_check), margin = 1),2),
                   man.chitest = chisq.test(no_out$norm_info, no_out$man_check), n = nrow(no_out))

# Data overview ---------------------------
participants <- c(export = nrow(clean),
                  unique = nrow(no_dup),
                  complete = nrow(complete),
                  no_veg = nrow(no_veg),
                  no_care = nrow(no_care),
                  no_out = nrow(no_out),
                  nomiss_full = sum(!is.na(no_dup$FFQ_T2_7)),
                  nomiss_out = sum(!is.na(no_out$FFQ_T2_7)))
age_desc      <- sapply(no_out[,"age"], function(x) c(min = min(x), max = max(x), avg = mean(x), sd = sd(x))) ## age
gender_freq   <- 100 * prop.table(table(no_out$gender)) ## gender
baseline_meat <-  group_by(no_out, condition) %>% summarise(m = mean(FFQ_T1_meat), sd = sd(FFQ_T1_meat)) ## meat consumption at T1

library(apaTables)
no_out %>% 
  mutate(FFQ_T2_meat       = rowSums(across(FFQ_T2_2:FFQ_T2_7)),
         FFQ_meat_diff     = FFQ_T2_meat - FFQ_T1_meat) %>%
  select(ends_with("T1"), ends_with("T2"), "cons_proj", "FFQ_meat_diff") %>%
  apa.cor.table(filename = "outputs/corrtab.doc", show.conf.interval = FALSE)

no_out %>%
  mutate(FFQ_T2_meat       = rowSums(across(FFQ_T2_2:FFQ_T2_7)),
         FFQ_meat_diff     = FFQ_T2_meat - FFQ_T1_meat) %>%
  select(ends_with("T1"), ends_with("T2"), "cons_proj", "FFQ_meat_diff") %>%
  summarise_all(funs(sum(!is.na(.))))

library(psych)
no_out %>%
  select("cons_now_perc", "cons_next_perc",  "cons_six_perc") %>%
  psych::alpha()

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

T1_interact_overlaid2 <- 
  no_out %>% #interaction and simple effects PLOT WITHOUT OVERLAY
    filter(condition != "none") %>%
    select(c("norm_info", "format", "attitude_T1", "interest_T1", "intentionuni_T1")) %>%
    pivot_longer(cols = ends_with("_T1"), names_to = "outcome", values_to = "value") %>%
    group_by(norm_info, format, outcome) %>%
    summarise(n= n(),
              m = mean(value),
              se = sd(value)/sqrt(n())) %>%
    ggplot(aes(x = norm_info, color = format, group = format, shape = format, y = m, ymin = m-se, ymax = m+se)) +
    facet_wrap(~ outcome, labeller = as_labeller(c('attitude_T1'="Attitude", 'interest_T1'="Interest", 'intentionuni_T1'="Intention"))) +
    geom_line(aes(linetype = format, group = format), linewidth = 0.5) +
    geom_point() +
    geom_errorbar(aes(ymin = m-se, ymax = m+se, group = format), width = 0.2) +
    coord_cartesian(ylim = c(0,100)) +
    scale_y_continuous(breaks = seq(0, 100, 10)) +
    labs(y = "Value of Outcome (%)", x = "Norm Type", colour = "Format", shape = "Format", linetype = "Format")+
    scale_x_discrete(labels= c("Static", "Dynamic")) +
    scale_color_manual(values = c("black", "azure4"), labels= c("Text", "Text +\nVisual Cue")) +
    scale_shape_discrete(labels= c("Text", "Text +\nVisual Cue")) +
    scale_linetype_discrete(labels= c("Text", "Text +\nVisual Cue")) +
    papaja::theme_apa() +
    theme(legend.box.background = element_rect(color = "black", linewidth = 0.6), legend.position = c(.85,.75), strip.text = element_text(face = "bold"), axis.title = element_text(face="bold"), legend.title = element_text(face="bold"))

T1_interact_overlaid3 <- 
  no_out %>% #interaction and simple effects PLOT WITHOUT OVERLAY
  filter(condition != "none") %>%
  select(c("norm_info", "format", "attitude_T1", "interest_T1", "intentionuni_T1")) %>%
  pivot_longer(cols = ends_with("_T1"), names_to = "outcome", values_to = "value") %>%
  group_by(norm_info, format, outcome) %>%
  summarise(n= n(),
            m = mean(value),
            se = sd(value)/sqrt(n())) %>%
  ggplot(aes(x = norm_info, color = format, group = format, shape = format, y = m, ymin = m-se, ymax = m+se)) +
  facet_wrap(~ outcome, labeller = as_labeller(c('attitude_T1'="Attitude", 'interest_T1'="Interest", 'intentionuni_T1'="Intention"))) +
  geom_line(aes(linetype = format, group = format), linewidth = 0.5) +
  geom_point() +
  geom_errorbar(aes(ymin = m-se, ymax = m+se, group = format), width = 0.2) +
  coord_cartesian(ylim = c(0,100)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  labs(y = "Value of Outcome (%)", x = "Norm Type", colour = "Format", shape = "Format", linetype = "Format")+
  scale_x_discrete(labels= c("Static", "Dynamic")) +
  scale_color_manual(values = c("black", "azure4"), labels= c("Text + Visual Cue Absent", "Text + Visual Cue Present")) +
  scale_shape_discrete(labels= c("Text + Visual Cue Absent", "Text + Visual Cue Present")) +
  scale_linetype_discrete(labels= c("Text + Visual Cue Absent", "Text + Visual Cue Present")) +
  papaja::theme_apa() +
  theme(legend.box.background = element_rect(color = "black", linewidth = 0.6), legend.position = "top", strip.text = element_text(face = "bold"), axis.title = element_text(face="bold"), legend.title = element_blank())


no_out %>% #plot for SPSP
  filter(condition != "none") %>%
  select(c("norm_info", "format", "attitude_T1", "interest_T1", "intentionuni_T1")) %>%
  pivot_longer(cols = ends_with("_T1"), names_to = "outcome", values_to = "value") %>%
  group_by(norm_info, format, outcome) %>%
  summarise(n= n(),
            m = mean(value),
            se = sd(value)/sqrt(n())) %>%
  ggplot(aes(x = norm_info, color = format, group = format, y = m, ymin = m-se, ymax = m+se)) +
  facet_wrap(~ outcome, labeller = as_labeller(c('attitude_T1'="Attitude", 'interest_T1'="Interest", 'intentionuni_T1'="Intention"))) +
  geom_line(aes(linetype = format, group = format), position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = m-se, ymax = m+se, group = format), width = 0.2, position = position_dodge(width = 0.3)) +
  coord_cartesian(ylim = c(30,70)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  labs(y = "Value (%)", x = "Norm Type", colour = "Format")+
  scale_x_discrete(labels= c("Static", "Dynamic")) +
  scale_color_manual(values = c("#E69F00", "#0072B2"), labels= c("Text", "Text + Visual")) +
  papaja::theme_apa(base_size = 16)

## descriptive statistics
outcomes_T1_desc <- group_by(no_out, condition) %>% ## outcomes by each condition
  summarise(n = n(),
            across(c("attitude_T1", "interest_T1", "intentionuni_T1"), 
                   list(mean = mean, sd = sd, se = ~sd(.)/sqrt(n())), 
                   .names = "{.fn}_{.col}"), .groups = "rowwise")

write_csv(outcomes_T1_desc,"outputs/t1_desc.csv")

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

contrast_posthoc = list(stext_vs_dtext = c(0, 1, 0, -1,0),
                        stvis_vs_dyvis = c(1, 0, -1, 0, 0))

h2posthoc.models <- map_df(outcomes_T1, ~ lm(formula(paste0(.x, "~", "condition")), data = no_out) %>%
                      emmeans(., "condition", contr = contrast_posthoc, infer = T) %>%
                      .$contrasts %>%
                      summary(), .id = "outcome")

h2posthoc_bf <- sapply(1:6, function(x) 
  bfrr(h2posthoc.models$estimate[x], 
       h2posthoc.models$SE[x], 
       sample_df = h2posthoc.models$df[x], 
       model = "normal",
       mean = 0, 
       sd = 5, 
       tail = 1, 
       criterion = 3, 
       rr_interval = list(mean = c(-10, 10), sd = c(0, 10)), 
       precision = 0.05)[-14])
h2posthoc_rr <- sapply(1:6, function(x) paste0("[", toString(h2posthoc_bf[,x]$RR$sd), "]"))
h2posthoc_table <- cbind(h2posthoc.models, unlist(h2posthoc_bf[3,]), h2posthoc_rr, unlist(h2posthoc_bf[5,])) %>% 
  mutate(p.value = printp(p.value)) %>%
  modify_if(., ~is.numeric(.), ~round(., 2)) %>%
  unite(CL, c("lower.CL", "upper.CL"), sep = ", ") 

h2posthoc_bf_invert <- sapply(1:6, function(x) 
  bfrr(h2posthoc.models$estimate[x]*-1, 
       h2posthoc.models$SE[x], 
       sample_df = h2posthoc.models$df[x], 
       model = "normal",
       mean = 0, 
       sd = 5, 
       tail = 1, 
       criterion = 3, 
       rr_interval = list(mean = c(-10, 10), sd = c(0, 10)), 
       precision = 0.05)[-14])
h2posthoc_rr_invert <- sapply(1:6, function(x) paste0("[", toString(h2posthoc_bf_invert[,x]$RR$sd), "]"))
h2posthoc_table_invert <- cbind(h2posthoc_table, unlist(h2posthoc_bf_invert[3,]), h2posthoc_rr_invert, unlist(h2posthoc_bf_invert[5,])) %>% 
  modify_if(., ~is.numeric(.), ~round(., 2)) %>%
  arrange(contrast)

write_csv(h2posthoc_table_invert,"outputs/h2posthoc.csv")

## T1 output ---------------------------
h.tests <- c("stat_vs_dy", "none_vs_dy", "norm_infostatic:formattext", "dytext_vs_dyvis")
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
T1_rr <- sapply(1:12, function(x) paste0("[", toString(T1_bf[,x]$RR$sd), "]"))

T1_table <- cbind(T1_outcomes, unlist(T1_bf[3,]), T1_rr, unlist(T1_bf[5,])) %>% 
  mutate(p.value = printp(p.value)) %>%
  modify_if(., ~is.numeric(.), ~round(., 2)) %>%
  unite(CL, c("lower.CL", "upper.CL"), sep = ", ") %>%
  #mutate_all(~ str_replace(., "-", "−")) %>%
  arrange(match(contrast, h.tests))

write_csv(T1_table,"outputs/T1_table.csv")
save(participants, outcomes_T1_desc, T1_table, T1_violin, T1_interact_overlaid2, T1_interact_overlaid3, file =  "outputs/T1_out.RData")

## Multiple imputation ---------------------------
### variables predicting missingness
prop_missing <- apply(no_out, 2, function(col)sum(is.na(col))/length(col)) #check % missing from all columns
prop_missing_T2 <- prop_missing[match('interest_T2', names(prop_missing)):match("FFQ_T2_15", names(prop_missing))]*100
complete_Vars <- names(prop_missing[prop_missing==0])

xtabs(~condition+complete, data = no_out) # condition x missingness
xtabs(~gender+complete, data = no_out) # gender x missingness

aux_logit <- glm(complete ~ condition+age+genderbi+politics, data = no_out, family = "binomial") #check if auxiliary variables predict missingness
summary(aux_logit)

### setting up predictor matrix
init <- mice(no_out, maxit = 0)
meth <- init$meth
pred <- init$pred
pred[, names(no_out)[c(1:5,32, 54:64)]] <- 0
pred[names(no_out)[c(1:5,32, 54:64)],] <- 0
meth[names(meth) %in% names(no_out)[c(1:5,32, 54:64)]] <- ""

###### split by group
imputed <- no_out %>%
  split(.$norm_info) %>%
  map(., ~ .x %>% droplevels) %>%
  map(~ mice(.x, m = 100, maxit = 10, meth = meth, pred = pred))

imputed_rbind1 <- rbind(imputed$dynamic, imputed$static)
imputed_rbindfull <- rbind(imputed_rbind1, imputed$none)
imputed_long <- mice::complete(imputed_rbindfull, "long", include = T)

imputed_calc <- imputed_long %>%
  mutate(FFQ_T2_meat       = rowSums(across(FFQ_T2_2:FFQ_T2_7)),
         interest_diff     = interest_T2 - interest_T1,
         attitude_diff     = attitude_T2 - attitude_T1,
         intention_diff    = intentionuni_T2 - intentionuni_T1,
         intention_bi_diff = intentionbi_T2 - intentionbi_T1,
         FFQ_meat_diff     = FFQ_T2_meat - FFQ_T1_meat,
         format            = fct_relevel(format, "none", "text", "visual"),
         condition         = fct_relevel(condition, "dy_graph", "dy_text", "st_graph", "st_text", "none"),
         norm_info         = fct_relevel(norm_info, "none", "static", "dynamic"))

imputed_dat <- as.mids(imputed_calc)
imputed_filtered <- filter(imputed_dat, FFQ_T1_meat <= max_FFQ_meat) #excludes MAD outliers on the FFQ

saveRDS(imputed, file = "data/imputed.rds")

# T2 overview ---------------------------
## visualization
T2_violin <- filter(imputed_calc, .imp == 50) %>%
  pivot_longer(cols = c("attitude_diff", "interest_diff", "intention_diff"), names_to = "variable", values_to = "value") %>%
  ggplot2::ggplot(aes(norm_info, value, colour = format)) +
  facet_wrap(~ variable, labeller = as_labeller(c('attitude_diff'="Attitude", 'interest_diff'="Interest", 'intention_diff'="Intention"))) +
  geom_violin(trim = T,  position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_normal", position = position_dodge(width = 0.9)) +
  labs(y = "Value (%)", x = "Norm type", colour = "Format", title = "Distribution of change in outcomes by condition at T2") + 
  scale_x_discrete(labels= c("None", "Static", "Dynamic")) +
  scale_color_discrete(labels= c("None", "Text", "Visual")) +
  papaja::theme_apa()

T2_interact <- filter(imputed_calc, .imp == 50) %>% #interaction and simple effects PLOT WITHOUT JITTER
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

filter(imputed_calc, .imp == 50) %>% #interaction and simple effects PLOT WITHOUT JITTER
  filter(condition != "none") %>%
  select(c("norm_info", "format", "attitude_T2", "interest_T2", "intentionuni_T2")) %>%
  pivot_longer(cols = ends_with("_T2"), names_to = "outcome", values_to = "value") %>%
  group_by(norm_info, format, outcome) %>%
  summarise(n= n(),
            m = mean(value),
            se = sd(value)/sqrt(n())) %>%
  ggplot(aes(x = norm_info, color = format, group = format, y = m, ymin = m-se, ymax = m+se)) +
  facet_wrap(~ outcome, labeller = as_labeller(c('attitude_T2'="Attitude", 'interest_T2'="Interest", 'intentionuni_T2'="Intention"))) +
  geom_line(aes(linetype = format, group = format), position = position_dodge(width = 0.9), show.legend = FALSE) +
  geom_point(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = m-se, ymax = m+se, group = format), width = 0.2, position = position_dodge(width = 0.9)) +
  coord_cartesian(ylim = c(0,100)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  labs(y = "Value (%)", x = "Norm type", colour = "Format")+
  scale_x_discrete(labels= c("Static", "Dynamic")) +
  scale_color_discrete(labels= c("Text", "Text + Visual")) +
  papaja::theme_apa()

FFQ_change <- filter(imputed_calc, .imp == 50 & FFQ_T1_meat <= max_FFQ_meat) %>%
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
  geom_line(aes(group = condition), show.legend = T, position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = m-se, ymax = m+se), width = 0.2, position = position_dodge(width = 0.9)) +
  scale_color_discrete(labels= c("Visual dynamic", "Text dynamic", "Visual static", "Text static", "None")) +
  labs(y = "Meat consumed (Servings)", x = "Time", colour = "Condition")+
  papaja::theme_apa()


## descriptive statistics
change_by_cond <- filter(imputed_calc, .imp == 50) %>%
  group_by(condition) %>%
  summarise(n = n(),
            across(c("attitude_diff", "interest_diff","intention_diff", "FFQ_meat_diff"), 
                   list(mean = ~mean(., na.rm = TRUE), 
                        sd = ~sd(., na.rm = TRUE), 
                        se = ~sd(., na.rm = TRUE)/sqrt(n())), 
                   .names = "{.fn}_{.col}"), .groups = "rowwise")

diffs_by_cond <- filter(imputed_calc, .imp == 50) %>%
  group_by(condition) %>%
  summarise(n = n(),
            across(c("attitude_T2", "interest_T2","intentionuni_T2"), 
                   list(mean = ~mean(., na.rm = TRUE), 
                        sd = ~sd(., na.rm = TRUE), 
                        se = ~sd(., na.rm = TRUE)/sqrt(n())), 
                   .names = "{.fn}_{.col}"), .groups = "rowwise")

meat_by_cond <- filter(imputed_calc, .imp == 50 & FFQ_T1_meat <= max_FFQ_meat) %>%
  group_by(condition) %>%
  summarise(n = n(),
            across(c("FFQ_T2_meat"), 
                   list(mean = ~mean(., na.rm = TRUE), 
                        sd = ~sd(., na.rm = TRUE), 
                        se = ~sd(., na.rm = TRUE)/sqrt(n())), 
                   .names = "{.fn}_{.col}"), .groups = "rowwise")

write_csv(diffs_by_cond,"outputs/t2_desc.csv")
write_csv(meat_by_cond,"outputs/t2meat_desc.csv")

## H4 ---------------------------
outcomes_T2 <- c("attitude_diff", "interest_diff", "intention_diff", "FFQ_meat_diff") %>% set_names(.)
h4.models <- map_df(outcomes_T2[1:3], ~ with(imputed_dat , lm(formula(paste0(.x, "~", "norm_info"))))  %>%
                     emmeans(., "norm_info", contr = contrasts_norm, infer = T) %>%
                     .$contrasts %>%
                     summary(), .id = "outcome") %>%
  bind_rows(map_df(outcomes_T2[4], ~ with(imputed_filtered , lm(formula(paste0(.x, "~", "norm_info"))))  %>%
                     emmeans(., "norm_info", contr = contrasts_norm, infer = T) %>%
                     .$contrasts %>%
                     summary(), .id = "outcome"))

## H5 ---------------------------
h5.df <- summary(pool(with(imputed_dat, (lm(attitude_diff ~ norm_info*format))))) %>% .$df
h5d.df <- summary(pool(with(imputed_filtered, (lm(attitude_diff ~ norm_info*format))))) %>% .$df

h5.models <- map_df(outcomes_T2[1:3], ~ 
                      with(imputed_dat, lm(formula(paste0(.x, "~", "norm_info*format")))) %>%
                      pool() %>%
                      tidy(conf.int=T), .id = "outcome") %>%
  filter(term == "norm_infostatic:formattext") %>%
  cbind(., h5.df[2]) %>%
  select(1:4, 11, 7:8, 5:6) %>%
  bind_rows(map_df(outcomes_T2[4], ~ 
                     with(imputed_filtered, lm(formula(paste0(.x, "~", "norm_info*format")))) %>%
                     pool() %>%
                     tidy(conf.int=T), .id = "outcome") %>%
              filter(term == "norm_infostatic:formattext") %>%
              cbind(., h5d.df[2]) %>%
              select(1:4, 11, 7:8, 5:6))


## H6 ---------------------------
h6.models <- map_df(outcomes_T2[1:3], ~ 
                      with(imputed_dat, lm(formula(paste0(.x, "~", "condition")))) %>%
                      emmeans(., "condition", contr = contrast_interact, infer = T) %>%
                      .$contrasts %>%
                      summary(), .id = "outcome") %>%
  bind_rows(map_df(outcomes_T2[4], ~ 
                     with(imputed_filtered, lm(formula(paste0(.x, "~", "condition")))) %>%
                     emmeans(., "condition", contr = contrast_interact, infer = T) %>%
                     .$contrasts %>%
                     summary(), .id = "outcome") )


## T2 output ---------------------------
T2_outcomes <- bind_rows(h4.models,setNames(h5.models,names(h4.models)), h6.models) %>%
  arrange(match(outcome, outcomes_T2))

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
       precision = 0.05)[-14]) %>%
  cbind(., T2_meat_bf)

T2_rr <- sapply(1:16, function(x) paste0("[", toString(T2_bf[,x]$RR$sd), "]"))
T2_table <- cbind(T2_outcomes, unlist(T2_bf[3,]), T2_rr, unlist(T2_bf[5,])) %>% 
  mutate(p.value = printp(p.value)) %>%
  modify_if(., ~is.numeric(.), ~round(., 2)) %>%
  unite(CL, c("lower.CL", "upper.CL"), sep = ", ") %>%
  arrange(match(contrast, h.tests))

write_csv(T2_table,"outputs/T2_table.csv")

save(prop_missing_T2, T2_violin, T2_interact, FFQ_change, change_by_cond, diffs_by_cond, T2_table, file =  "outputs/T2_out.RData")

hypotheses_table <- bind_rows(T1_table,setNames(T2_table,names(T1_table)))
write_csv(hypotheses_table,"outputs/hypotheses_table.csv")

## secondary analyses ---------------------------
## manipulation checks
stdy_only <- filter(no_out, norm_info != "none") %>% droplevels()
think_reduceT2 <- list(think_desc  = round(100*prop.table(table(stdy_only$condition, stdy_only$think_reduce), margin = 1),2),
                   think.chitest = chisq.test(stdy_only$condition, stdy_only$think_reduce), n = nrow(stdy_only))

active_reduceT2 <- list(active_desc  = round(100*prop.table(table(stdy_only$condition, stdy_only$active_reduce), margin = 1),2),
                       active.chitest = chisq.test(stdy_only$condition, stdy_only$active_reduce), n = nrow(stdy_only))

measures_T2 <- c("attitude_T2", "interest_T2", "intentionuni_T2") %>% set_names(.)
sec_norm.models <- map_df(measures_T2, ~ with(imputed_dat , lm(formula(paste0(.x, "~", "norm_info"))))  %>%
                      emmeans(., "norm_info", contr = contrasts_norm, infer = T) %>%
                      .$contrasts %>%
                      summary(), .id = "outcome")
sec_interact.df <- summary(pool(with(imputed_dat, (lm(attitude_T2 ~ norm_info*format))))) %>% .$df

sec_interact.models <- map_df(measures_T2, ~ 
                      with(imputed_dat, lm(formula(paste0(.x, "~", "norm_info*format")))) %>%
                      pool() %>%
                      tidy(conf.int=T), .id = "outcome") %>%
  filter(term == "norm_infostatic:formattext") %>%
  cbind(., sec_interact.df[2]) %>%
  select(1:4, 11, 7:8, 5:6)

sec_simple.models <- map_df(measures_T2, ~ 
                      with(imputed_dat, lm(formula(paste0(.x, "~", "condition")))) %>%
                      emmeans(., "condition", contr = contrast_interact, infer = T) %>%
                      .$contrasts %>%
                      summary(), .id = "outcome")

secposthoc.models <- map_df(measures_T2, ~ 
                              with(imputed_dat, lm(formula(paste0(.x, "~", "condition")))) %>%
                             emmeans(., "condition", contr = contrast_posthoc, infer = T) %>%
                             .$contrasts %>%
                             summary(), .id = "outcome")

secposthoc_bf <- sapply(1:6, function(x) 
  bfrr(secposthoc.models$estimate[x], 
       secposthoc.models$SE[x], 
       sample_df = secposthoc.models$df[x], 
       model = "normal",
       mean = 0, 
       sd = 5, 
       tail = 1, 
       criterion = 3, 
       rr_interval = list(mean = c(-10, 10), sd = c(0, 10)), 
       precision = 0.05)[-14])
secposthoc_rr <- sapply(1:6, function(x) paste0("[", toString(secposthoc_bf[,x]$RR$sd), "]"))
secposthoc_table <- cbind(secposthoc.models, unlist(secposthoc_bf[3,]), secposthoc_rr, unlist(secposthoc_bf[5,])) %>% 
  mutate(p.value = printp(p.value)) %>%
  modify_if(., ~is.numeric(.), ~round(., 2)) %>%
  unite(CL, c("lower.CL", "upper.CL"), sep = ", ") %>%

secposthoc_bf_invert <- sapply(1:6, function(x) 
  bfrr(secposthoc.models$estimate[x]*-1, 
       secposthoc.models$SE[x], 
       sample_df = secposthoc.models$df[x], 
       model = "normal",
       mean = 0, 
       sd = 5, 
       tail = 1, 
       criterion = 3, 
       rr_interval = list(mean = c(-10, 10), sd = c(0, 10)), 
       precision = 0.05)[-14])
secposthoc_rr_invert <- sapply(1:6, function(x) paste0("[", toString(secposthoc_bf_invert[,x]$RR$sd), "]"))
secposthoc_table_invert <- cbind(secposthoc_table, unlist(secposthoc_bf_invert[3,]), secposthoc_rr_invert, unlist(secposthoc_bf_invert[5,])) %>% 
  modify_if(., ~is.numeric(.), ~round(., 2)) %>%
  unite(CL, c("lower.CL", "upper.CL"), sep = ", ") %>%
  arrange(contrast)

write_csv(secposthoc_table_invert,"outputs/posthoc_sec_table.csv")

## T2 secondary output ---------------------------
T2_secondary <- bind_rows(sec_norm.models,setNames(sec_interact.models,names(sec_norm.models)), sec_simple.models) %>%
  arrange(match(outcome, measures_T2))

T2_sec_bf <- sapply(1:12, function(x) 
  bfrr(T2_secondary$estimate[x], 
       T2_secondary$SE[x], 
       sample_df = T2_secondary$df[x], 
       model = "normal",
       mean = 0, 
       sd = 5, 
       tail = 1, 
       criterion = 3, 
       rr_interval = list(mean = c(-10, 10), sd = c(0, 10)), 
       precision = 0.05)[-14])

sec_rr <- sapply(1:12, function(x) paste0("[", toString(T2_sec_bf[,x]$RR$sd), "]"))
sec_table <- cbind(T2_secondary, unlist(T2_sec_bf[3,]), sec_rr, unlist(T2_sec_bf[5,])) %>% 
  mutate(p.value = printp(p.value)) %>%
  modify_if(., ~is.numeric(.), ~round(., 2)) %>%
  unite(CL, c("lower.CL", "upper.CL"), sep = ", ") %>%
  arrange(match(contrast, h.tests))

write_csv(sec_table,"outputs/sec_table.csv")


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
