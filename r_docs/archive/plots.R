
T1_jitter <- no_out %>% #interaction and simple effects WITH JITTER
  filter(condition != "none") %>%
  pivot_longer(cols = c("attitude_T1", "interest_T1", "intentionuni_T1"), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = norm_info, color = format, group = format, y = value)) +
  facet_wrap(~ variable, labeller = as_labeller(c('attitude_T1'="Attitude", 'interest_T1'="Interest", 'intentionuni_T1'="Intention"))) +
  geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.3)+
  stat_summary(fun = mean, geom = "line", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_se,
               geom = "pointrange", position = position_dodge(width = 0.9)) +
  labs(y = "Value (%)", x = "Norm Type", colour = "Format")+
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
  ggplot(aes(x = norm_info, color = format, group = format, shape = format, y = m, ymin = m-se, ymax = m+se)) +
  facet_wrap(~ outcome, labeller = as_labeller(c('attitude_T1'="Attitude", 'interest_T1'="Interest", 'intentionuni_T1'="Intention"))) +
  geom_line(aes(linetype = format, group = format), position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = m-se, ymax = m+se, group = format), width = 0.2, position = position_dodge(width = 0.3)) +
  coord_cartesian(ylim = c(0,100)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  labs(y = "Value (%)", x = "Norm Type", colour = "Format", shape = "Format")+
  scale_x_discrete(labels= c("Static", "Dynamic")) +
  scale_color_discrete(labels= c("Text", "Text +\n visual cue")) +
  scale_shape_discrete(labels= c("Text", "Text +\n visual cue")) +
  papaja::theme_apa()

T1_interact_overlaid1 <- 
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
  scale_color_manual(values = c("black", "azure4"), labels= c("Text", "Text + Visual Cue")) +
  scale_shape_discrete(labels= c("Text", "Text + Visual Cue")) +
  scale_linetype_discrete(labels= c("Text", "Text + Visual Cue")) +
  papaja::theme_apa() +
  theme(legend.box.background = element_rect(color = "black", linewidth = 0.6), legend.position = "top", strip.text = element_text(face = "bold"), axis.title = element_text(face="bold"), legend.title = element_text(face="bold"))

outcomes_change <- 
  filter(imputed_calc, .imp == 50) %>% 
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
  geom_line(aes(group = condition), show.legend = T, position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = m-se, ymax = m+se), width = 0.2, position = position_dodge(width = 0.9)) +
  scale_color_discrete(labels= c("Visual dynamic", "Text dynamic", "Visual static", "Text static", "None")) +
  labs(y = "Value (%)", x = "Time", colour = "Condition")+
  papaja::theme_apa()