## Associations between VAF and blood counts and blood biochemistry

### Here we calculate if/how changes in VAF are affected by changes in blood count and blood chemistry data.

```{r load data and merge}
stat_sign_individuals <- statistics_data %>%
  group_by(individual) %>%
  summarise(sig = sum(sum_stat < 0.7) > 0) %>%
  subset(sig == F)

blood_count_data <- load_blood_count_data()
fuller_data <- merge(load_data(),blood_count_data,by = c("SardID","Phase")) %>%
  subset(!(SardID %in% c(load_excluded_individuals(),load_excluded_individuals_lymph()))) %>%
  subset(!(SardID %in% stat_sign_individuals$individual))
```

```{r data preprocessing}
complete_blood_counts <- fuller_data %>%
  group_by(SardID) %>%
  mutate(relative_timepoint = Phase - min(Phase) + 1) %>% 
  select(SardID,Phase,VAF,Age,amino_acid_change,relative_timepoint,
         RBC,HB,MCV,MCH,WBC,NEUT_PERC,LYMPH_PERC,MONO_PERC,EOS_PERC,BASO_PERC,
         PLT,Glucose,SerumInsulin,BLoodUreaNitrogen,Creatinine,ALT = ALT.y,AST,GGT,Fibrinogen,
         TotalCholesterol,HDLcholesterol,Triglycerides,Iron,Transferrin,UricAcid,ESR,Phase,Age) %>%
  group_by(SardID) %>%
  mutate(complete = sum(is.na(RBC)) == 0) %>%
  ungroup() %>%
  subset(complete == T) %>%
  mutate(
    NEUT = WBC * NEUT_PERC/100,
    MONO = WBC * MONO_PERC/100,
    LYMPH = WBC * LYMPH_PERC/100,
    BASO = WBC * BASO_PERC/100,
    EOS = WBC * EOS_PERC/100
  ) %>%
  mutate(GRAN = NEUT + BASO + EOS)
```

Latent class mixture models allow us to define different classes for time series. For our case, the time series we consider here are those described by blood count and chemistry parameters through time, allowing us to stratify individuals based on their growth. This was inspired by [this work](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6819541/pdf/12889_2019_Article_7752.pdf).

```{r prepare_lcmm}
complete_blood_counts_lcmm <- complete_blood_counts %>% 
  subset(!(SardID %in% c(load_excluded_individuals(),load_excluded_individuals_lymph()))) %>%
  merge(distinct(select(full_data,SardID,Gender)),by = "SardID",all = F) %>%
  group_by(SardID) %>%
  mutate(MinAge = min(Age),
         Age = Age - min(Age),
         subject = as.numeric(SardID)) %>%
  distinct()

LINK <- "linear"
MIXTURE_FORMULA <- formula(" ~ Age + Age:factor(Gender) + factor(Gender)")

LCMM_MODELS <- list()
```

```{r lcmm_esr, fig.height=10,fig.width=5}
LCMM_MODELS$ESR <- list()

lcmm_fit <- lcmm(
  fixed = ESR ~ Age +  Age:factor(Gender) + factor(Gender),
  random = ~ 1,
  subject = "subject",
  mixture = MIXTURE_FORMULA,
  data = as.data.frame(complete_blood_counts_lcmm),
  link = LINK,
  idiag = T,
  ng = 2)

LCMM_MODELS$ESR$MODEL_FULL <- lcmm(
  fixed = ESR ~ Age +  Age:factor(Gender) + factor(Gender),
  random = ~ 1,
  subject = "subject",
  data = as.data.frame(complete_blood_counts_lcmm),
  link = LINK,
  idiag = 1)

summary(lcmm_fit)

tmp <- which.max(table(lcmm_fit$pprob$class))

if (tmp == 1) {
  lcmm_fit$pprob <- data.frame(
    subject = lcmm_fit$pprob$subject,
    class = lcmm_fit$pprob$class,
    prob1 = lcmm_fit$pprob$prob1,
    prob2 = lcmm_fit$pprob$prob2
  )
} else {
  lcmm_fit$pprob <- data.frame(
    subject = lcmm_fit$pprob$subject,
    class = 3 - lcmm_fit$pprob$class,
    prob1 = lcmm_fit$pprob$prob2,
    prob2 = lcmm_fit$pprob$prob1
  )
}

LCMM_MODELS$ESR$MODEL <- lcmm_fit
LCMM_MODELS$ESR$TRAJECTORIES <- complete_blood_counts_lcmm %>%
  merge(lcmm_fit$pprob,by = "subject") %>%
  mutate(class = factor(ifelse(class == 1,'A','B'),levels=c("A","B"))) %>%
  ggplot(aes(x = Age,y = ESR)) + 
  geom_line(aes(group = SardID),alpha = 0.5,colour = 'grey') +
  facet_wrap( ~ class) +
  geom_smooth(method = 'loess',formula = y ~ x,span = 0.8,method.args = list(degree = 1),
              colour = 'black',fill = 'grey20') + 
  theme_gerstung(base_size = 15) + 
  theme(legend.position = 'bottom') + 
  xlab("Years since study entry")

LCMM_MODELS$ESR$BOXPLOT <- r_values %>%
  select(individual,coefficient,b_clone) %>%
  group_by(individual) %>%
  filter((b_clone) == max(b_clone)) %>%
  distinct() %>%
  merge(lcmm_fit$pprob,by.x = "individual",by.y = "subject",all.y = T) %>%
  mutate(class = factor(ifelse(class == 1,'A','B'),levels=c("A","B"))) %>%
  mutate(b_clone = ifelse(is.na(b_clone),0,b_clone)) %>% 
  ggplot(aes(x = as.factor(class),y = b_clone)) + 
  geom_boxplot(alpha = 0.8) + 
  geom_signif(comparisons = list(c("A","B")),
              map_signif_level = function(x) return(paste("p-value =",round(x,4))),
              test = wilcox.test) + 
  stat_summary(geom = 'text',
               fontface = 'italic',
               fun.data = function(x) return(list(y = 0.22,label = sprintf('n = %s',length(x))))) +
  xlab("Trajectory type") + 
  ylab("Growth effect") + 
  theme_gerstung(base_size = 15) 
```

```{r lcmm_wbc, fig.height=10,fig.width=5}
LCMM_MODELS$WBC <- list()

lcmm_fit <- lcmm(
  fixed = WBC ~ Age + Age:factor(Gender) + factor(Gender),
  random = ~ 1,
  subject = "subject",
  mixture = MIXTURE_FORMULA,
  data = as.data.frame(complete_blood_counts_lcmm),
  link = LINK,
  idiag = T,
  ng = 2)

LCMM_MODELS$WBC$MODEL_FULL <- lcmm(
  fixed = WBC ~ Age +  Age:factor(Gender) + factor(Gender),
  random = ~ 1,
  subject = "subject",
  data = as.data.frame(complete_blood_counts_lcmm),
  link = LINK,
  idiag = 1)

summary(lcmm_fit)

tmp <- which.max(table(lcmm_fit$pprob$class))

if (tmp == 1) {
  lcmm_fit$pprob <- data.frame(
    subject = lcmm_fit$pprob$subject,
    class = lcmm_fit$pprob$class,
    prob1 = lcmm_fit$pprob$prob1,
    prob2 = lcmm_fit$pprob$prob2
  )
} else {
  lcmm_fit$pprob <- data.frame(
    subject = lcmm_fit$pprob$subject,
    class = 3 - lcmm_fit$pprob$class,
    prob1 = lcmm_fit$pprob$prob2,
    prob2 = lcmm_fit$pprob$prob1
  )
}

LCMM_MODELS$WBC$MODEL <- lcmm_fit
LCMM_MODELS$WBC$TRAJECTORIES <- complete_blood_counts_lcmm %>%
  merge(lcmm_fit$pprob,by = "subject") %>%
  mutate(class = factor(ifelse(class == 1,'A','B'),levels=c("A","B"))) %>%
  ggplot(aes(x = Age,y = WBC)) + 
  geom_line(aes(group = SardID),alpha = 0.5,colour = 'grey') +
  facet_wrap( ~ class) +
  geom_smooth(method = 'loess',formula = y ~ x,span = 0.8,method.args = list(degree = 1),
              colour = 'black',fill = 'grey20') + 
  theme_gerstung(base_size = 15) + 
  theme(legend.position = 'bottom') + 
  xlab("Years since study entry")

LCMM_MODELS$WBC$BOXPLOT <- r_values %>%
  select(individual,coefficient,b_clone) %>%
  group_by(individual) %>%
  filter((b_clone) == max(b_clone)) %>%
  distinct() %>%
  merge(lcmm_fit$pprob,by.x = "individual",by.y = "subject",all.y = T) %>%
  mutate(class = factor(ifelse(class == 1,'A','B'),levels=c("A","B"))) %>%
  mutate(b_clone = ifelse(is.na(b_clone),0,b_clone)) %>% 
  ggplot(aes(x = as.factor(class),y = b_clone)) + 
  geom_boxplot(alpha = 0.8) + 
  geom_signif(comparisons = list(c("A","B")),
              map_signif_level = function(x) return(paste("p-value =",round(x,4))),
              test = wilcox.test) + 
  stat_summary(geom = 'text',
               fontface = 'italic',
               fun.data = function(x) return(list(y = 0.22,label = sprintf('n = %s',length(x))))) +
  xlab("Trajectory type") + 
  ylab("Growth effect") + 
  theme_gerstung(base_size = 15) 
```

```{r lcmm_hb, fig.height=10,fig.width=5}
LCMM_MODELS$HB <- list()

lcmm_fit <- lcmm(
  fixed = HB ~ Age + Age:factor(Gender) + factor(Gender),
  random = ~ 1,
  subject = "subject",
  mixture = MIXTURE_FORMULA,
  data = as.data.frame(complete_blood_counts_lcmm),
  link = LINK,
  idiag = T,
  ng = 2)

LCMM_MODELS$HB$MODEL_FULL <- lcmm(
  fixed = HB ~ Age +  Age:factor(Gender) + factor(Gender),
  random = ~ 1,
  subject = "subject",
  data = as.data.frame(complete_blood_counts_lcmm),
  link = LINK,
  idiag = T)

summary(lcmm_fit)

tmp <- which.max(table(lcmm_fit$pprob$class))

if (tmp == 1) {
  lcmm_fit$pprob <- data.frame(
    subject = lcmm_fit$pprob$subject,
    class = lcmm_fit$pprob$class,
    prob1 = lcmm_fit$pprob$prob1,
    prob2 = lcmm_fit$pprob$prob2
  )
} else {
  lcmm_fit$pprob <- data.frame(
    subject = lcmm_fit$pprob$subject,
    class = 3 - lcmm_fit$pprob$class,
    prob1 = lcmm_fit$pprob$prob2,
    prob2 = lcmm_fit$pprob$prob1
  )
}

LCMM_MODELS$HB$MODEL <- lcmm_fit
LCMM_MODELS$HB$TRAJECTORIES <- complete_blood_counts_lcmm %>%
  merge(lcmm_fit$pprob,by = "subject") %>%
  mutate(class = factor(ifelse(class == 1,'A','B'),levels=c("A","B"))) %>%
  ggplot(aes(x = Age,y = HB)) + 
  geom_line(aes(group = SardID),alpha = 0.5,colour = 'grey') +
  facet_wrap( ~ class) +
  geom_smooth(method = 'loess',formula = y ~ x,span = 0.8,method.args = list(degree = 1),
              colour = 'black',fill = 'grey20') + 
  theme_gerstung(base_size = 15) + 
  theme(legend.position = 'bottom') + 
  xlab("Years since study entry")

LCMM_MODELS$HB$BOXPLOT <- r_values %>%
  select(individual,coefficient,b_clone) %>%
  group_by(individual) %>%
  filter((b_clone) == max(b_clone)) %>%
  distinct() %>%
  merge(lcmm_fit$pprob,by.x = "individual",by.y = "subject",all.y = T) %>%
  mutate(class = factor(ifelse(class == 1,'A','B'),levels=c("A","B"))) %>%
  mutate(b_clone = ifelse(is.na(b_clone),0,b_clone)) %>% 
  ggplot(aes(x = as.factor(class),y = b_clone)) + 
  geom_boxplot(alpha = 0.8) + 
  geom_signif(comparisons = list(c("A","B")),
              map_signif_level = function(x) return(paste("p-value =",round(x,4))),
              test = wilcox.test) + 
  stat_summary(geom = 'text',
               fontface = 'italic',
               fun.data = function(x) return(list(y = 0.22,label = sprintf('n = %s',length(x))))) +
  xlab("Trajectory type") + 
  ylab("Growth effect") + 
  theme_gerstung(base_size = 15) 
```

```{r lcmm_plt, fig.height=10,fig.width=5}
LCMM_MODELS$PLT <- list()

lcmm_fit <- lcmm(
  fixed = PLT ~ Age + Age:factor(Gender) + factor(Gender),
  random = ~ 1,
  subject = "subject",
  mixture = MIXTURE_FORMULA,
  data = as.data.frame(complete_blood_counts_lcmm),
  link = LINK,
  idiag = T,
  ng = 2)

LCMM_MODELS$PLT$MODEL_FULL <- lcmm(
  fixed = PLT ~ Age +  Age:factor(Gender) + factor(Gender),
  random = ~ 1,
  subject = "subject",
  data = as.data.frame(complete_blood_counts_lcmm),
  link = LINK,
  idiag = 1)

summary(lcmm_fit)

tmp <- which.max(table(lcmm_fit$pprob$class))

if (tmp == 1) {
  lcmm_fit$pprob <- data.frame(
    subject = lcmm_fit$pprob$subject,
    class = lcmm_fit$pprob$class,
    prob1 = lcmm_fit$pprob$prob1,
    prob2 = lcmm_fit$pprob$prob2
  )
} else {
  lcmm_fit$pprob <- data.frame(
    subject = lcmm_fit$pprob$subject,
    class = 3 - lcmm_fit$pprob$class,
    prob1 = lcmm_fit$pprob$prob2,
    prob2 = lcmm_fit$pprob$prob1
  )
}

LCMM_MODELS$PLT$MODEL <- lcmm_fit
LCMM_MODELS$PLT$TRAJECTORIES <- complete_blood_counts_lcmm %>%
  merge(lcmm_fit$pprob,by = "subject") %>%
  mutate(class = factor(ifelse(class == 1,'A','B'),levels=c("A","B"))) %>%
  ggplot(aes(x = Age,y = PLT)) + 
  geom_line(aes(group = SardID),alpha = 0.5,colour = 'grey') +
  facet_wrap( ~ class) +
  geom_smooth(method = 'loess',formula = y ~ x,span = 0.8,method.args = list(degree = 1),
              colour = 'black',fill = 'grey20') + 
  theme_gerstung(base_size = 15) + 
  theme(legend.position = 'bottom') + 
  xlab("Years since study entry")

LCMM_MODELS$PLT$BOXPLOT <- r_values %>%
  select(individual,coefficient,b_clone) %>%
  group_by(individual) %>%
  filter((b_clone) == max(b_clone)) %>%
  distinct() %>%
  merge(lcmm_fit$pprob,by.x = "individual",by.y = "subject",all.y = T) %>%
  mutate(class = factor(ifelse(class == 1,'A','B'),levels=c("A","B"))) %>%
  mutate(b_clone = ifelse(is.na(b_clone),0,b_clone)) %>% 
  ggplot(aes(x = as.factor(class),y = b_clone)) + 
  geom_boxplot(alpha = 0.8) + 
  geom_signif(comparisons = list(c("A","B")),
              map_signif_level = function(x) return(paste("p-value =",round(x,4))),
              test = wilcox.test) + 
  stat_summary(geom = 'text',
               fontface = 'italic',
               fun.data = function(x) return(list(y = 0.22,label = sprintf('n = %s',length(x))))) +
  xlab("Trajectory type") + 
  ylab("Growth effect") + 
  theme_gerstung(base_size = 15) 
```

```{r visualise_all_lcmm,fig.height = 6.5,fig.width = 14}
plot_grid(
  plot_grid(LCMM_MODELS$WBC$TRAJECTORIES + 
              theme(plot.title = element_text(hjust = 0.5)) + 
              xlab("") + 
              ylab(bquote('WBC'~(10^3/uL))),
            LCMM_MODELS$WBC$BOXPLOT + 
              xlab("") + 
              ylab("Unknown cause effect") + 
              theme(panel.grid.major.x = element_blank()),
            ncol = 1,
            rel_heights = c(0.9,0.8)),
  plot_grid(LCMM_MODELS$HB$TRAJECTORIES + 
              theme(plot.title = element_text(hjust = 0.5)) + 
              ylab("Hb (g/dL)"),
            LCMM_MODELS$HB$BOXPLOT + 
              ylab("") + 
              theme(panel.grid.major.x = element_blank()),
            ncol = 1,
            rel_heights = c(0.9,0.8)),
  plot_grid(LCMM_MODELS$PLT$TRAJECTORIES + 
              theme(plot.title = element_text(hjust = 0.5)) + 
              xlab("") +
              ylab(bquote('Plt'~(10^3/uL))),
            LCMM_MODELS$PLT$BOXPLOT + 
              xlab("") + 
              ylab("") + 
              theme(panel.grid.major.x = element_blank()),
            ncol = 1,
            rel_heights = c(0.9,0.8)),
  plot_grid(LCMM_MODELS$ESR$TRAJECTORIES + 
              theme(plot.title = element_text(hjust = 0.5)) + 
              xlab("") + 
              ylab("ESR (mm/h)"),
            LCMM_MODELS$ESR$BOXPLOT + 
              xlab("") + 
              ylab("") + 
              theme(panel.grid.major.x = element_blank()),
            ncol = 1,
            rel_heights = c(0.9,0.8)),
  nrow = 1
) + 
  ggsave(useDingbats=FALSE,filename = sprintf("figures/%s/phenotype_associations/FigureBloodCounts.pdf",model_id),
         height = 6.5,width = 14) +
  ggsave(filename = sprintf("figures/%s/phenotype_associations/FigureBloodCounts.svg",model_id),
         height = 6.5,width = 14)
```

We can also analyse this data if we consider each blood count/chemistry parameter trajectory independently and calculate its slope and average value. These values can then be used in a simple model as predictors of unknown growth effects. We also use the LCMM model trajectory probabilities in a second linear model to assess their utility in these scenarios. 

```{r linear_model_blood_counts,fig.width=10,fig.height=4} 
model_probs <- list(select(LCMM_MODELS$ESR$MODEL$pprob,SardID=subject,ESR_traj_B=prob2), 
                    select(LCMM_MODELS$WBC$MODEL$pprob,SardID=subject,WBC_traj_B=prob2),
                    select(LCMM_MODELS$HB$MODEL$pprob,SardID=subject,HB_traj_B=prob2),
                    select(LCMM_MODELS$PLT$MODEL$pprob,SardID=subject,PLT_traj_B=prob2)
) %>% 
  reduce(left_join, by = "SardID") 

slope <- function(x,y) {
  s_x <- sum(x)
  s_x2 <- sum(x^2)
  s_y <- sum(y)
  s_xy <- sum(x*y)
  return((s_xy - s_x * s_y) / (s_x2 - s_x^2))
}                             

data_linear_model <- merge(load_data(),blood_count_data,by = c("SardID","Phase"),all.y = T) %>% 
  group_by(SardID) %>%
  mutate(MeanAge = mean(Age)) %>% 
  group_by(SardID,Gender,MeanAge) %>%
  summarise(
    WBC_mean = mean(WBC),
    Hb_mean = mean(HB),
    Plt_mean = mean(PLT),
    ESR_mean = mean(ESR),
    WBC_slope = slope(Age,WBC),
    Hb_slope = slope(Age,HB),
    Plt_slope = slope(Age,PLT),
    ESR_slope = slope(Age,ESR)
  ) %>% 
  merge(distinct(select(r_values,"individual","b_clone","coefficient","site")),by.x = "SardID",by.y = "individual",all.x = T)  %>% 
  merge(model_probs,by = "SardID",all=T)

data_linear_model <- data_linear_model[
  apply(data_linear_model[,4:11],1,function(x) !any(is.na(x))),
  ] 
data_linear_model[,c(3:11,15:18)] <- apply(data_linear_model[,c(3:11,15:18)],2,scale)
data_linear_model <- data_linear_model %>% 
  mutate(b_clone = ifelse(is.na(b_clone),0,b_clone),
         coefficient = exp(ifelse(is.na(coefficient),0,coefficient))) %>%
  group_by(SardID) %>%
  filter(abs(b_clone) == max(abs(b_clone)))

linear_model_slopes <- lm(
  b_clone ~ WBC_mean + Hb_mean + Plt_mean + ESR_mean + WBC_slope + Hb_slope + Plt_slope + ESR_slope + MeanAge + Gender,
  data = data_linear_model) 

linear_model_traj_probs <- lm(
  b_clone ~ WBC_mean + Hb_mean + Plt_mean + ESR_mean + ESR_traj_B + WBC_traj_B + HB_traj_B + PLT_traj_B + MeanAge + Gender,
  data = data_linear_model
)

lm_res_plot <- data.frame(
  b_clone = data_linear_model$b_clone,
  Slopes = predict(linear_model_slopes,data_linear_model),
  Trajectories = predict(linear_model_traj_probs,data_linear_model)
) %>% 
  gather("key","value",Slopes,Trajectories) %>%
  ggplot(aes(x = b_clone,y = b_clone)) +
  geom_abline(slope = 1,size = 1.5) +
  geom_segment(aes(xend = b_clone,yend = value,colour = key),alpha = 0.2) + 
  geom_point(aes(y = value,colour = key),alpha=0.5) + 
  scale_color_aaas(name = "Model") + 
  theme_gerstung(base_size = 15) + 
  theme(legend.position = "bottom") + 
  xlab("Unknown effect") + 
  ylab("Unknown effect")

lmss <- summary(linear_model_slopes)
lmtps <- summary(linear_model_traj_probs)
lm_coef_plot <- data.frame(coefficient = c(lmss$coefficients[,1],lmtps$coefficients[,1]),
                           std_error = c(lmss$coefficients[,2],lmtps$coefficients[,2]),
                           p.val = c(lmss$coefficients[,4],lmtps$coefficients[,4]),
                           names = c(rownames(lmss$coefficients),rownames(lmtps$coefficients)),
                           model = c(rep(sprintf("Slopes\n(Adj. R sq. = %.4f)",lmss$adj.r.squared),nrow(lmss$coefficients)),
                                     rep(sprintf("Trajectories\n(Adj. R sq. = %.4f)",lmtps$adj.r.squared),nrow(lmtps$coefficients)))
) %>% 
  mutate(signif = cut(p.val,breaks = c(0,0.001,0.01,0.05,0.1,1),labels = c("***","**","*",".","n.s."))) %>%
  mutate(namesGroup = ifelse(
    grepl("_mean",names),
    "Mean",
    ifelse(
      grepl("_slope",names),
      "Slope",
      ifelse(
        grepl("_traj_",names),
        "Trajectories",
        "Other"
      )
    )
  )) %>% 
  mutate(namesGroup = factor(namesGroup,levels = c("Mean","Slope","Trajectories","Other"))) %>% 
  mutate(names = str_match(names,"[A-Za-z()]+")) %>% 
  ggplot(aes(x = names,y = coefficient,colour = model)) + 
  geom_hline(yintercept = 0,alpha = 0.7) +
  geom_point() + 
  geom_errorbar(aes(ymin = coefficient - std_error,ymax = coefficient + std_error),
                width = 0.2) + 
  geom_text(aes(label = signif,y = coefficient + std_error + 0.02),size=3) +
  theme_gerstung(base_size = 15) + 
  facet_wrap(~ model + namesGroup,scales = "free_x",nrow = 1) + 
  rotate_x_text() + 
  xlab("Parameters") + 
  ylab("Coefficient") + 
  scale_color_aaas() + 
  theme(legend.position = "none",
        strip.text.x = element_text(margin = margin(0,0,0,0, "cm")))

plot_grid(lm_res_plot + theme(legend.position = "none"),lm_coef_plot,rel_widths = c(2,3.8)) %>%
  plot_grid(get_legend(lm_res_plot),ncol=1,rel_heights = c(1,0.1))
```

From this analysis we can observe that using the slopes in this context does allows us to make relatively little conclusions between clonal growth and other factors (average parameter value ,trajectory membership and slopes). Apart from a positive association with Plt slope and a negative association with Plt mean in the slopes model, little else can be inferred from this.

```{r compare_trajectories}
trajectory_class_data <- rbind(
  data.frame(LCMM_MODELS$ESR$MODEL$pprob[,c(1,2)],
             Parameter = "ESR"),
  data.frame(LCMM_MODELS$WBC$MODEL$pprob[,c(1,2)],
             Parameter = "WBC"),
  data.frame(LCMM_MODELS$HB$MODEL$pprob[,c(1,2)],
             Parameter = "Hb"),
  data.frame(LCMM_MODELS$PLT$MODEL$pprob[,c(1,2)],
             Parameter = "Plt")
) %>%
  group_by(subject) %>%
  mutate(N = sum(class == 2)) %>%
  ungroup() %>%
  arrange(N) %>%
  mutate(subject = factor(subject,
                          levels = rev(unique(subject)),
                          ordered = T)) 

N_trajectory_plot <- trajectory_class_data %>%
  ggplot(aes(y = "N",x = subject,fill = N)) + 
  geom_tile() +
  scale_fill_material(palette = "blue") + 
  theme_gerstung() +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  theme(axis.text.x = element_blank(),
        legend.position = 'none') + 
  ylab("") + 
  xlab("")

class_trajectory_plot <- trajectory_class_data %>%
  ggplot(aes(x = subject,y = Parameter,fill = as.factor(class))) + 
  geom_tile(color = 'black') + 
  theme_gerstung() + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_blank(),
        legend.position = 'bottom') +
  scale_fill_lancet(name="Trajectory") + 
  xlab("Individual")

plot_grid(
  N_trajectory_plot,
  class_trajectory_plot,
  ncol=1,
  rel_heights = c(0.15,1.0),align = T
)
```

```{r table_trajectory_coefficients}
trajectory_coefficient_data <- rbind(
  data.frame(coeff = c(coef(LCMM_MODELS$WBC$MODEL)[c(3,4)],
                       coef(LCMM_MODELS$WBC$MODEL_FULL)[1]),
             Parameter = "WBC",
             class = c(1,2,"Full")),
  data.frame(coeff = c(coef(LCMM_MODELS$HB$MODEL)[c(3,4)],
                       coef(LCMM_MODELS$HB$MODEL_FULL)[1]),
             Parameter = "Hb",
             class = c(1,2,"Full")),
  data.frame(coeff = c(coef(LCMM_MODELS$PLT$MODEL)[c(3,4)],
                       coef(LCMM_MODELS$PLT$MODEL_FULL)[1]),
             Parameter = "Plt",
             class = c(1,2,"Full")),
  data.frame(coeff = c(coef(LCMM_MODELS$ESR$MODEL)[c(3,4)],
                       coef(LCMM_MODELS$ESR$MODEL_FULL)[1]),
             Parameter = "ESR",
             class = c(1,2,"Full"))
)

trajectory_coefficient_data %>% 
  spread(key = "class",value = "coeff") %>%
  mutate_if(.predicate = is.numeric,function(x) round(x,4))
```

To better assess the relevance of these data we need to compare them with reference values. A few studies ([here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3446744/) and [here](https://www.jstage.jst.go.jp/article/geriatrics1964/28/4/28_4_509/_pdf/-char/ja)) have assessed this, so we can use them to approximate a reference value. 

```{r reference_values}
reference_values_women_HD <- data.frame(
  Decade = c("50-59","60-69","70-79","80-89","+90"),
  Decade_numeric = c(5,6,7,8,9),
  Hb = c(12.1,10.5,10.2,9.3,8.8),
  Plt = c(442.40,385.91,443.72,388.04,360.02),
  WBC = c(4.6,3.2,4.5,4.2,4.0),
  N = c(122,100,121,111,49),
  Sex = "F",
  Place = "HD"
)
reference_values_men_HD <- data.frame(
  Decade = c("50-59","60-69","70-79","80-89","+90"),
  Decade_numeric = c(5,6,7,8,9),
  Hb = c(12.1,11.7,9.9,9.6,9.7),
  Plt = c(136.5,123.9,82.8,104.8,137.4),
  WBC = c(4.3,3.9,4.1,4.2,NA),
  N = c(99,102,105,113,20),
  Sex = "M",
  Place = "HD"
)

reference_values_women_JP_1 <- data.frame(
  Decade = c("50-59","60-69","70-79","80-89","+90"),
  Decade_numeric = c(5,6,7,8,9),
  Hb = c(13.6,13.3,NA,NA,NA),
  Plt = c(235,244,NA,NA,NA),
  WBC = c(5.4,5.5,NA,NA,NA),
  N = c(152,33,NA,NA,NA),
  Sex = "F",
  Place = "JP"
)
reference_values_men_JP_1 <- data.frame(
  Decade = c("50-59","60-69","70-79","80-89","+90"),
  Decade_numeric = c(5,6,7,8,9),
  Hb = c(15.3,15.0,14.6,NA,NA),
  Plt = c(231,232,237,NA,NA),
  WBC = c(6.3,6.1,5.9,NA,NA),
  N = c(300,110,6,NA,NA),
  Sex = "M",
  Place = "JP"
)

reference_values_women_JP_2 <- data.frame(
  Decade = c("50-59","60-69","70-79","80-89","+90"),
  Decade_numeric = c(5,6,7,8,9),
  Hb = c(NA,13.8,13.8,13.4,13.5),
  Plt = c(NA,279,268,269,205),
  WBC = c(NA,5.9,5.8,5.8,5.3),
  N = c(NA,81,208,109,4),
  Sex = "F",
  Place = "JP"
)
reference_values_men_JP_2 <- data.frame(
  Decade = c("50-59","60-69","70-79","80-89","+90"),
  Decade_numeric = c(5,6,7,8,9),
  Hb = c(NA,13.8,13.8,13.4,13.5),
  Plt = c(NA,269,261,239,269),
  WBC = c(NA,6.7,6.3,6.6,6.2),
  N = c(NA,72,106,61,4),
  Sex = "M",
  Place = "JP"
)

reference_values <- rbind(
  reference_values_women_HD,
  reference_values_women_JP_1,
  reference_values_women_JP_2,
  reference_values_men_HD,
  reference_values_men_JP_1,
  reference_values_men_JP_2
)

glm(Hb ~ Decade_numeric + Place,data = reference_values,weights = reference_values$N) %>%
  summary
glm(Plt ~ Decade_numeric + Place,data = reference_values,weights = reference_values$N) %>%
  summary 
glm(WBC ~ Decade_numeric + Place,data = reference_values,weights = reference_values$N) %>%
  summary