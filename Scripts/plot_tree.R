# setup -------------------------------------------------------------------
tm <- 0.352777778
pop_size <- 1e5
source("Scripts/vaf_dynamics_functions.R")
library(phylodyn)
INLA:::inla.dynload.workaround()
model_file_name <- "models/model_ch.RDS"
set.seed(42)
model_id <- 'model_ch'
load(file = sprintf('data_output/vaf_modelling_coefficients_%s.Rdata',model_id))
values_model <- readRDS(file = model_file_name)
c('trees') %>%
  sapply(
    function(x) dir.create(sprintf("figures/%s/%s",model_id,x),showWarnings = F)
  )
source("Scripts/prepare_data.R")

clone_assignment <- list(
  `2259` = list(PD41276k_25 = "SF3B1-K666N"),
  `500` = list(PD41305k_130 = "Unknown driver 1",
               PD41305k_109 = "Unknown driver 3",
               PD41305k_100 = "Unknown driver 2"),
  `3877` = list(PD34493k_45 = "Unknown driver 2",
                PD34493k_32 = "Unknown driver 1",
                PD34493k_13 = "SF3B1-K666N",
                PD34493k_79 = "U2AF1-Q157R")
)
trees <- load_trees()
tree_details <- load_tree_details()
ages <- list(`500`=73.4 + 0.75, 
             `2259`=79.4 + 0.75, 
             `3877`=83.8 + 0.75)
ages_vector <- do.call(c,ages)

change_point <- function(x, b0, m1, m2, delta) { 
  b0 + (x*m1) + (sapply(x-delta, function (t) max(0, t)) * m2)
}

site_samples <- values_model$draws$`11`[,grep('site',colnames(values_model$draws$`11`))]
gene_samples <- values_model$draws$`11`[,grep('gene',colnames(values_model$draws$`11`))]
clone_samples <- values_model$draws$`11`[,grep('clone',colnames(values_model$draws$`11`))]
offset_samples <- values_model$draws$`11`[,grep('u',colnames(values_model$draws$`11`))]
beta_samples <- values_model$draws$`11`[,grep('beta',colnames(values_model$draws$`11`))] %>%
  unlist()

sampling_idxs <- r_values_full %>%
  subset(b_site_095 + b_gene_095 + b_clone_095 > 0) %>%
  group_by(
    individual,site,
    site_numeric,
    gene_numeric,
    clone_numeric
  ) %>% 
  summarise(
    minAge = min(age[which(round(true/coverage,4) >= 0.005)]),
    maxAge = max(age),
    maxVAF = max(true/coverage)
  ) %>%
  arrange(site)

c_names <- colnames(values_model$draws$`11`)

u_c_names <- grep('^u',c_names,perl = T)
b_gene_c_names <- grep('^b_gene',c_names,perl = T) 
b_site_c_names <- grep('^b_site',c_names,perl = T)
b_clone_c_names <- grep('^b_clone',c_names,perl = T)

maximum_probability_value <- function(x) {
  d <- density(x)
  return(d$x[which.max(d$y)])
}

grab_samples <- function(site_numeric,gene_numeric,clone_numeric,sample_size=1000) {
  if (!is.na(site_numeric)) {
    b_site <- sample(site_samples[,site_numeric],sample_size,replace = F)
  } else {
    b_site <- 0
  }
  b_gene <- sample(gene_samples[,gene_numeric],sample_size,replace = F)
  b_clone <- sample(clone_samples[,clone_numeric],sample_size,replace = F)
  u <- sample(offset_samples[,clone_numeric],sample_size,replace = F)
  return(list(b_site,b_clone,b_gene,u))
}

bnpr_at_coalescence <- function(bnpr_estimate) {
  time_at_coal <- bnpr_estimate$coal_times
  X <- bnpr_estimate$summary$time
  Y <- bnpr_estimate$summary$quant0.5
  Y025 <- bnpr_estimate$summary$quant0.025
  Y975 <- bnpr_estimate$summary$quant0.975
  Y_log <- log(na.omit(approx(X,Y,xout = time_at_coal,rule = 2)$y))
  values_at_coalescence <- data.frame(
    X = time_at_coal,
    Y = Y_log,
    Y_exp = approx(X,Y,xout = time_at_coal)$y,
    Y025 = approx(X,Y025,xout = time_at_coal)$y,
    Y975 = approx(X,Y975,xout = time_at_coal)$y
  ) %>%
    na.omit()
  return(values_at_coalescence)
}

bnpr_at_all <- function(bnpr_estimate) {
  X <- bnpr_estimate$summary$time
  Y <- bnpr_estimate$summary$quant0.5
  Y025 <- bnpr_estimate$summary$quant0.025
  Y975 <- bnpr_estimate$summary$quant0.975
  Y_log <- log(Y)
  values_at_coalescence <- data.frame(
    X = X,
    Y = Y_log,
    Y_exp = Y,
    Y025 = Y025,
    Y975 = Y975
  ) %>%
    na.omit()
  return(values_at_coalescence)
}

sample_ages <- function(sampling_idxs,pop_size,gen_time,sample_size=1000) {
  sampling_idxs %>% 
    apply(1,function(x) {
      site_numeric <- as.numeric(x[3])
      gene_numeric <- as.numeric(x[4])
      clone_numeric <- as.numeric(x[5])
      
      S <- grab_samples(site_numeric,gene_numeric,clone_numeric,sample_size)
      b_site <- S[[1]]
      b_clone <- S[[2]]
      b_gene <- S[[3]]
      u <- S[[4]]
      
      t1 <- as.numeric(x[6])
      age_distribution <- t0_adjusted(u,b_site+b_gene+b_clone,gen_time,pop_size*2) + values_model$min_age
      
      tmp <- data.frame(
        age_distribution = age_distribution,
        effect_distribution = b_clone + b_site + b_gene,
        offset_distribution = u,
        site = x[2],
        individual = gsub(" ","",x[1]),
        gene = str_match(x[2],'[A-Z0-9]+'),
        t1 = t1
      ) 
      colnames(tmp) <- c("age_distribution","effect_distribution","offset_distribution","site","individual","gene","t1")
      return(tmp)
    }) %>%
    return
}

characterise_age_samples <- function(age_samples) {
  age_samples %>% 
    lapply(function(x) {
      A <- x$age_distribution
      t1 <- x$t1[1]
      A <- ifelse(A < 0,0, A)
      A <- ifelse(A > t1, t1, A)
      A <- A[!is.na(A)]
      M <- mean(A)
      Qs <- quantile(A,c(0.05,0.10,0.16,0.5,0.84,0.90,0.95))
      D <- density(A)
      MAP <- D$x[which.max(D$y)]
      return(
        c(Mean = M,Qs,MAP = MAP)
      )
    }) %>% 
    do.call(what = rbind) %>% 
    as_tibble() %>%
    return
}

# calculate bnpr estimates ------------------------------------------------
bnpr_estimates <- trees %>% 
  lapply(
    function(x) {
      
      if (x$ID == "2259") {
        clades <- cut_tree(x$tree_ultra, 0.5)
      } else {
        clades <- cut_tree(x$tree_ultra, 0.1)
      }
      ID <- x$ID
      tree_chron <- x$tree_ultra
      tree_chron$edge.length <- tree_chron$edge.length * ages[[x$ID]]
      
      if(x$ID=="3877"){
        late_clades <- cut_tree(x$tree_ultra, 0.7)
        late_clades <- late_clades[sapply(late_clades, length) > 5]
        clades <- c(lapply(clades, setdiff, unlist(late_clades)), late_clades)
      }
      
      x$clade_dynamics <- list()
      i <- 0
      for (clade in clades) { 
        if (length(clade) >= 5) {
          i <- i + 1
          tree_subset <- keep.tip(tree_chron,clade)
          tree_subset <- multi2di(tree_subset)
          
          bnpr_estimate = BNPR(data = tree_subset)
          x$clade_dynamics[[as.character(i)]] <- list(
            tree_subset = tree_subset,
            bnpr_estimate = bnpr_estimate,
            clade = clade
          )
        } 
      }
      return(x)
    }
  )

# calculate linear fits ---------------------------------------------------
par(mfrow = c(4,4),mar=c(1,1,1,1))
for (item_name in names(bnpr_estimates)) {
  tree_size <- length(trees[[item_name]]$tree_ultra$tip.label)
  for (clade_name in names(bnpr_estimates[[item_name]][["clade_dynamics"]])) {
    clade <- bnpr_estimates[[item_name]][["clade_dynamics"]][[clade_name]]
    sard_id <- item_name
    age_range <- full_data %>% 
      subset(SardID == sard_id) %>%
      summarise((min(Age)-ages[[item_name]]),(max(Age)-ages[[item_name]])) %>%
      unlist %>%
      c
    
    bnpr_estimate <- clade$bnpr_estimate
    vaf <- length(clade$clade) / tree_size
    Y <- log(bnpr_estimate$summary$quant0.5)
    X <- -bnpr_estimate$summary$time
    W <- (log(bnpr_estimate$summary$quant0.975) - log(bnpr_estimate$summary$quant0.025))^2/16
    linear_estimate <- lm(Y ~ X,w = 1/W)
    mp <- mean(X[Y > max(Y)/2])
    Y_exp <- exp(Y)
    non_linear_estimate <- nls(exp(Y) ~ SSlogis(X,theta1,theta2,theta3),
                               weights = 1/W,
                               control = nls.control(warnOnly = T,
                                                     maxiter = 1000,
                                                     minFactor = 1 / (1024^2)),
                               start = c(theta1 = max(exp(Y)),
                                         theta2 = mp,
                                         theta3 = 1))
    non_linear_estimate_vaf <- nls(Y_exp / Y_exp[1] * vaf ~ SSlogis(X,1,theta2,theta3),
                                   weights = 1/W,
                                   control = nls.control(warnOnly = T,
                                                         maxiter = 1000,
                                                         minFactor = 1 / (1024^16)),
                                   start = c(theta2 = max(X),
                                             theta3 = 5))

    opt_fn <- function(par) {
      se <- (Y - change_point(X,par[1],par[2],par[3],par[4]))^2
      return(mean(se / W))
    }
    change_point_regression <- optim(
      par = c(linear_estimate$coefficients[1],
              linear_estimate$coefficients[2],
              1,
              mean(X)),
      method = "L-BFGS-B",
      lower = c(NA,0,NA,min(X) + diff(range(X)) * 0.25),
      upper = c(NA,NA,NA,max(X) - diff(range(X)) * 0.25),
      fn = opt_fn)
    
    b <- -log(200e3 - 1)
    expected_vaf <- inv.logit(change_point_regression$par[2] * diff(range(X)) + b)
    pp <- change_point_regression$par
    cp <- pp[4]
    ff <- change_point(X,pp[1],pp[2],pp[3],pp[4])
    se <- sqrt(sum((Y[X>cp]-ff[X>cp])^2)/(sum(X>cp)-2)/sqrt(sum((X[X>cp]-mean(X[X>cp]))^2)))
    ci_95 <- c(pp[2]+pp[3]-1.96*se,pp[2] + pp[3]+1.96*se)
    
    relevant_growth_df <- data.frame(
      X = X[X > age_range[1] & X < age_range[2]],
      Y = Y[X > age_range[1] & X < age_range[2]],
      W = 1/W[X > age_range[1] & X < age_range[2]])
    relevant_growth <- lm(Y ~ X,data = relevant_growth_df,
                          weights = relevant_growth_df$W)
    plot(X,Y)
    points(relevant_growth_df$X,relevant_growth_df$Y,col="red")
    lines(relevant_growth_df$X,predict(relevant_growth),col="blue")
    
    bnpr_estimates[[item_name]][["clade_dynamics"]][[clade_name]]$linear_estimate <- linear_estimate
    bnpr_estimates[[item_name]][["clade_dynamics"]][[clade_name]]$non_linear_estimate <- non_linear_estimate
    bnpr_estimates[[item_name]][["clade_dynamics"]][[clade_name]]$change_point_regression <- change_point_regression
    bnpr_estimates[[item_name]][["clade_dynamics"]][[clade_name]]$non_linear_estimate_vaf <- non_linear_estimate_vaf
    bnpr_estimates[[item_name]][["clade_dynamics"]][[clade_name]]$late_ci <- ci_95
    bnpr_estimates[[item_name]][["clade_dynamics"]][[clade_name]]$vaf <- vaf
    bnpr_estimates[[item_name]][["clade_dynamics"]][[clade_name]]$expected_vaf <- expected_vaf
    bnpr_estimates[[item_name]][["clade_dynamics"]][[clade_name]]$relevant_growth <- relevant_growth
  }
}

bnpr_estimates_coalescence <- bnpr_estimates
for (item_name in names(bnpr_estimates_coalescence)) {
  tree_size <- length(trees[[item_name]]$tree_ultra$tip.label)
  for (clade_name in names(bnpr_estimates_coalescence[[item_name]][["clade_dynamics"]])) {
    clade <- bnpr_estimates_coalescence[[item_name]][["clade_dynamics"]][[clade_name]]
    bnpr_estimate <- clade$bnpr_estimate
    vaf <- length(clade$clade) / tree_size
    tmp_df <- bnpr_at_all(bnpr_estimate)
    Y <- tmp_df$Y
    X <- -tmp_df$X
    W <- (log(tmp_df$Y975) - log(tmp_df$Y025))^2/16
    linear_estimate <- lm(Y ~ X,w = 1/W)
    mp <- mean(X[Y > max(Y)/2])
    Y_exp <- exp(Y)
    non_linear_estimate <- nls(exp(Y) ~ SSlogis(X,theta1,theta2,theta3),
                               weights = 1/W,
                               control = nls.control(warnOnly = T,
                                                     maxiter = 1000,
                                                     minFactor = 1 / (1024^2)),
                               start = c(theta1 = max(exp(Y)),
                                         theta2 = mp,
                                         theta3 = 1))
    non_linear_estimate_vaf <- nls(Y_exp / Y_exp[1] * vaf ~ SSlogis(X,1,theta2,theta3),
                                   weights = 1/W,
                                   control = nls.control(warnOnly = T,
                                                         maxiter = 1000,
                                                         minFactor = 1 / (1024^16)),
                                   start = c(theta2 = max(X),
                                             theta3 = 5))
    
    opt_fn <- function(par) {
      se <- (Y - change_point(X,par[1],par[2],par[3],par[4]))^2
      return(mean(se / W))
    }
    change_point_regression <- optim(
      par = c(linear_estimate$coefficients[1],
              linear_estimate$coefficients[2],
              1,
              mean(X)),
      method = "L-BFGS-B",
      lower = c(NA,0,NA,min(X) + diff(range(X)) * 0.25),
      upper = c(NA,NA,NA,max(X) - diff(range(X)) * 0.25),
      fn = opt_fn)
    bnpr_estimates_coalescence[[item_name]][["clade_dynamics"]][[clade_name]]$linear_estimate <- linear_estimate
    bnpr_estimates_coalescence[[item_name]][["clade_dynamics"]][[clade_name]]$non_linear_estimate <- non_linear_estimate
    bnpr_estimates_coalescence[[item_name]][["clade_dynamics"]][[clade_name]]$change_point_regression <- change_point_regression
    bnpr_estimates_coalescence[[item_name]][["clade_dynamics"]][[clade_name]]$non_linear_estimate_vaf <- non_linear_estimate_vaf
  }
}

# make plots --------------------------------------------------------------

legend_position <- c(0.33,1.1)
plot_list <- list()
for (curr_id in c("2259","500","3877")) {
  O <- bnpr_estimates[[curr_id]]
  max_age <- ages[[curr_id]]
  if (curr_id == "2259") {
    clade_dynamics <- O$clade_dynamics
  } else if (curr_id == "500") {
    clade_dynamics <- O$clade_dynamics
  } else {
    clade_dynamics <- O$clade_dynamics
  }
  
  plot_list[[curr_id]] <- list()
  
  clades_plot <- list()
  mutation_timings <- list()
  mrca <- c()
  mrca_names <- c()
  tree_ultra <- O$tree_ultra
  C <- 0
  for (clade in clade_dynamics) {
    C <- C + 1
    driver_id <- clone_assignment[[curr_id]][names(clone_assignment[[curr_id]]) %in% clade$tree_subset$tip.label] %>%
      unlist 
    clades_plot[[driver_id]] <- clade$clade
  
    mrca[C] <- MRCA(tree_ultra,clade$clade)
    mrca_names[C] <- driver_id
    mutation_timings[[driver_id]] <- time_mutation(tree_ultra,mrca[C]) * max_age
  }
  tree_detail <- tree_details[[curr_id]]$details
  tree_detail <- tree_detail[tree_detail$coding_change_CHgeneInModel != 'no',]
  tree_groups <- groupOTU(tree_ultra,clades_plot)
  
  edge_labels <- data.frame(tree_ultra$edge)
  edge_labels$LL <- tree_ultra$edge %>%
    apply(1,function(x) {
      ix <- which(mrca %in% x[2])
      if (length(ix) > 0) {
        return(mrca_names[ix])
      } else {
        return(NA)
      }
    }) 
  colnames(edge_labels)=c("parent", "node", "edge_num")
  edge_labels$top <- ifelse(
    grepl('-',edge_labels$edge_num),
    sprintf('italic("%s")\n',str_match(edge_labels$edge_num,'[A-Z0-9]+')),
    NA)
  edge_labels$top <- ifelse(
    grepl('Unknown driver',edge_labels$edge_num),
    'Unknown\n',
    edge_labels$top)
  edge_labels$bottom <- ifelse(
    grepl('-',edge_labels$edge_num),
    sprintf('\n%s',gsub('-','',str_match(edge_labels$edge_num,'-[A-Z0-9]+'))),
    NA)
  edge_labels$bottom <- ifelse(
    grepl('Unknown driver',edge_labels$edge_num),
    sprintf('\ndriver %s',str_match(edge_labels$edge_num,'[0-9]')),
    edge_labels$bottom)
  edge_labels$edge_num <- ifelse(
    !is.na(edge_labels$top),
    sprintf(
      'atop(atop(textstyle(%s),textstyle("%s")),NA)',
      edge_labels$top,
      edge_labels$bottom),
    NA
  )
  tree_detail <- tree_detail %>%
    mutate(
      Gene = ifelse(
        !(variant_ID %in% c("delY","del20q")),
        Gene,
        variant_ID)
    )
  other_mutation_annotation <- data.frame(
    gene = rep(NA,tree_ultra$Nnode*2),
    site = rep(NA,tree_ultra$Nnode*2),
    stringsAsFactors = F)
  other_mutation_annotation$annotation <- ""
  for (N in 1:length(tree_detail$node)) {
    nn <- tree_detail$node[N]
    site <- gsub('p.',"",tree_detail$Protein[N])
    site <- ifelse(is.na(site),"",site)
    other_mutation_annotation[edge_labels[,2] == nn,1] <- tree_detail$Gene[N]
    other_mutation_annotation[edge_labels[,2] == nn,2] <- site
    if (site != "") { 
      L <- sprintf(
        'atop(atop(textstyle(italic("%s")),textstyle("%s")),NA)',
        tree_detail$Gene[N],
        site)  
    } else {
      L <- sprintf(
        '"%s"',
        tree_detail$Gene[N])
    }
    other_mutation_annotation[edge_labels[,2] == nn,3] <- L
  }
  other_mutation_annotation <- cbind(edge_labels[,c(1,2)],
                                     other_mutation_annotation) 
  other_mutation_annotation <- other_mutation_annotation %>%
    mutate(annotation = ifelse(!(paste(gene,site,sep = '-') %in% mrca_names),
                               annotation,
                               NA)) %>%
    mutate(plot = ifelse(annotation != "","Y",NA)) %>%
    mutate(annotation = ifelse(annotation == '',NA,annotation))
  
  tree_plot <- ggtree(tree_groups,
                      aes(colour = str_match(group,'[A-Z0-9]+'),
                          size = ifelse(group == 0,"N","Y")),
                      ladderize=T) + 
    scale_size_manual(values = c(N = 0.2,Y = 0.4),
                      guide = F) + 
    scale_colour_manual(values = c(`0` = 'grey80',U = 'grey30',gene_colours,Y = "red4"),
                        guide=F) + 
    theme_gerstung(base_size = 6) +
    xlab("Age") + 
    ggtitle(sprintf("Patient ID=%s",curr_id)) +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 6),
          plot.title = element_text(margin = margin(b = -5),size = 6)) + 
    coord_cartesian(xlim = c(0,max_age) / max_age) +
    scale_x_continuous(expand = c(0.001,0),
                       breaks = (c(0,20,40,60,80)+0.75) / max_age,
                       labels = (c(0,20,40,60,80))) +
    scale_y_continuous(trans = "reverse")
  
  tree_plot <- tree_plot %<+% other_mutation_annotation + 
    geom_linerange(aes(x=branch,xmin = branch - branch.length/2,
                       xmax = branch + branch.length/2,
                       shape = plot),
                   size = 0.4,
                   colour = "red4") + 
    geom_label_repel(aes(x=branch, label = annotation),
                     alpha = 1,
                     size = 6 * tm,
                     hjust = 0,
                     vjust = 0.5,
                     segment.size = 0.2,
                     segment.color = "red",
                     force = 0.5,
                     min.segment.length = 0.05,
                     label.padding = unit(0.02,"cm"),
                     label.r = unit(0,"cm"),
                     label.size = unit(0,"cm"),
                     alpha = 0.7,
                     parse = T,
                     colour = "red4")
  
  tree_plot <- tree_plot %<+% edge_labels + 
    geom_label(aes(x=branch, label=edge_num),
               parse = T,
               alpha = 0.9,
               size = 6 * tm,
               vjust = 0.6,
               label.padding = unit(0,"cm"),
               label.r = unit(0,"cm"),
               label.size = unit(0,"cm"),
               colour = "black")

  for (clade in clade_dynamics) {
    driver_id <- clone_assignment[[curr_id]][names(clone_assignment[[curr_id]]) %in% clade$tree_subset$tip.label] %>%
      unlist
    linear_estimate <- clade$linear_estimate
    non_linear_estimate <- clade$linear_estimate
    tree_ultra <- O$tree_ultra
    tree_ultra$edge.length <- tree_ultra$edge.length * max_age
    
    bnpr_x <- -clade$bnpr_estimate$summary$time + max_age
    bnpr_y <- clade$bnpr_estimate$summary$mean
    bnpr_y_low <- clade$bnpr_estimate$summary$quant0.025
    bnpr_y_high <- clade$bnpr_estimate$summary$quant0.975
    W <- (log(clade$bnpr_estimate$summary$quant0.975) - log(clade$bnpr_estimate$summary$quant0.025))^2/16
    cp <- clade$change_point_regression$par
    observed_vaf <- clade$vaf
    expected_vaf <- clade$expected_vaf
    clade_frequency <- length(clade$tree_subset$tip.label) / length(trees[[curr_id]]$tree_ultra$tip.label)
    if (!grepl("Unknown driver",driver_id)) {
      tmp <- Parameters_Age %>% 
        subset(individual == as.numeric(curr_id) & site == driver_id)
      tmp_data <- full_data %>%
        subset(SardID == as.numeric(curr_id) & amino_acid_change == driver_id)
      
      curr_sampling_idxs <- sampling_idxs %>%
        subset(individual == curr_id & site == driver_id)
      beta_value <- values_model$beta_values[[1]][1,1]
      samples <- grab_samples(
        curr_sampling_idxs$site_numeric,
        curr_sampling_idxs$gene_numeric,
        curr_sampling_idxs$clone_numeric,
        sample_size = 2500)
      
      age_estimates <- t0_adjusted(
        samples[[4]],samples[[1]]+samples[[2]]+samples[[3]],
        2,pop_size) + values_model$min_age
      plausible_idx <- which(age_estimates > tmp$Q05 & !is.na(age_estimates) & age_estimates < tmp$Q95)
      X_pre <- data.frame(
        u = samples[[4]],
        b = samples[[1]] + samples[[2]] + samples[[3]])[plausible_idx,] %>%
        apply(1,function(x) {
          do.call(data.frame,trajectory_from_parameters_2(
            x[2],x[1],
            seq(-values_model$min_age,-values_model$min_age+max_age,by = 0.1),
            N = pop_size,g = 2)) 
          }) %>%
        do.call(what = rbind)

      X <- X_pre %>%
        mutate(y = ifelse(is.na(y),0,y)) %>% 
        group_by(x) %>% 
        summarise(y_ = median(y,na.rm=T),
                  y05 = quantile(y,0.05,na.rm=T),
                  y95 = quantile(y,0.95,na.rm=T)) %>% 
        select(x,y_,y05,y95)
      
      u <- mean(samples[[4]])
      b <- mean(samples[[1]] + samples[[2]] + samples[[3]])
      b_05 <- quantile(samples[[1]] + samples[[2]] + samples[[3]],0.05)
      u_05 <- quantile(samples[[4]],0.05)
      b_95 <- quantile(samples[[1]] + samples[[2]] + samples[[3]],0.95)
      u_95 <- quantile(samples[[4]],0.95)
      
      gen_per_year <- 2
      fitness <- b / gen_per_year
      pop_size <- 1e5 * 2
      drift_threshold <- 1 / (fitness * pop_size)
      time_at_onset <- t0_adjusted(u,b,2,pop_size) + values_model$min_age

      y <- inv.logit(b * (bnpr_x - values_model$min_age) + u) / 2
      y_05 <- inv.logit(b_05 * (bnpr_x - values_model$min_age) + u_05) / 2
      y_95 <- inv.logit(b_95 * (bnpr_x - values_model$min_age) + u_95) / 2
      
      vaf_trajectory <- data.frame(
        x = X$x + values_model$min_age,
        y = X$y_,
        y_05 = X$y05,
        y_95 = X$y95) %>%
        mutate(y = ifelse(y < 1 / pop_size,NA,y))
      
      bnpr_inferred_trajectory <- data.frame(
        x = bnpr_x,
        y = exp(predict(linear_estimate,x = bnpr_x))
      )
      
      bnpr_x_at_coalescence <- -clade$bnpr_estimate$coal_times+max_age
      bnpr_y_at_coalescence <- approx(
        y = bnpr_y,
        x = bnpr_x,
        xout = bnpr_x_at_coalescence
      )$y

      bnpr_at_coalescence <- data.frame(
        x = bnpr_x_at_coalescence,
        y = bnpr_y_at_coalescence
      ) %>% 
        na.omit()
      
      # bnpr_x <- Filter(f = function(x) {
      #   tmp_data$Age
      # })
      y_ <- approx(X$x+values_model$min_age,X$y_,
                   xout = bnpr_x)$y

      optim_func <- function(params) {
        scaling_factor <- params[1]
        return(mean(1/W * (y_ * scaling_factor - bnpr_y)^2))
      }
      solution <- optim(par = c(scaling_factor = 1e3),optim_func,method = "Brent",
                        lower = 100, upper = 1e10)
      
      scaling_factor <- solution$par
      scaling_factor <- bnpr_y[1] / clade_frequency
      print(c(bnpr_y[1] / clade_frequency,
              bnpr_y_low[1] / clade_frequency,
              bnpr_y_high[1] / clade_frequency,
              tmp_data$SardID[1],
              tmp_data$amino_acid_change[1]))
      vaf_trajectory$y <- vaf_trajectory$y * scaling_factor
      vaf_trajectory$y_05 <- vaf_trajectory$y_05 * scaling_factor
      vaf_trajectory$y_95 <- vaf_trajectory$y_95 * scaling_factor
    
      tmp_data <- tmp_data %>%
        mutate(VAF = scaling_factor * VAF) %>%
        mutate(VAF005 = qbeta(0.05,MUTcount_Xadj,TOTALcount - MUTcount_Xadj) * scaling_factor,
               VAF095 = qbeta(0.95,MUTcount_Xadj,TOTALcount - MUTcount_Xadj) * scaling_factor)
      y_lims <- c(
        min(vaf_trajectory$y[vaf_trajectory$y>0],na.rm = T),
        max(c(y * scaling_factor,bnpr_y,tmp_data$VAF)) * 5
      )
      y_lims <- c(1,1e6)
      
      x_axis_at <- unique(floor(bnpr_x / 10) * 10) %>%
        Filter(f = function(x) x > min(bnpr_x))
      if (length(x_axis_at) < 4) {
        x_axis_at <- unique(floor(bnpr_x / 5) * 5) %>%
          Filter(f = function(x) x > min(bnpr_x))
      }
      y_axis_at <- unique(10^floor(log(bnpr_y))) %>%
        Filter(f = function(x) x > y_lims[1])
      y_axis_at <- c(1,y_axis_at)
      y_axis_at_vaf <- c(10^c(-3:-1),0.5)
      
      colour_for_trajectory <- ifelse(
        grepl('-',driver_id),
        gene_colours[str_match(driver_id,'[0-9A-Z]+')],
        'grey60'
      )
      
      age_midpoint <- max_age / 2
      if (time_at_onset < age_midpoint) {
        annotate_hjust = 1
        annotate_x = max_age
      } else {
        annotate_hjust = 0
        annotate_x = 0
      }
      
      mutation_timing <- data.frame(
        x = c(mutation_timings[[driver_id]]),
        y = 0
      )
      
      driver_id_split <- str_split(driver_id,pattern = '-')[[1]]
      title_label <- substitute(
        expression(paste(
          italic(gene),'-', site,sep=''
        )),
        list(
          gene = driver_id_split[1],
          site = driver_id_split[2]
        )
      ) %>%
        eval()

      trajectory_plot_ <- data.frame(
        y = y * scaling_factor,
        y_05 = y_05 * scaling_factor,
        y_95 = y_95 * scaling_factor,
        bnpr_x = bnpr_x,
        bnpr_y = bnpr_y,
        bnpr_y025 = clade$bnpr_estimate$summary$quant0.025,
        bnpr_y975 = clade$bnpr_estimate$summary$quant0.975
      ) %>%
        ggplot() +
        geom_line(data = mutation_timing,
                  inherit.aes = F,
                  size = 3,
                  colour = colorspace::lighten(colour_for_trajectory,0.7),
                  aes(x = x,y = y)) + 
        geom_hline(yintercept = 2e5,colour = "black",linetype = 2) +
        geom_ribbon(aes(x = bnpr_x,
                        ymin = bnpr_y025,
                        ymax = bnpr_y975,
                        fill = "EPS (fit and 95% CI)",
                        shape = "EPS (fit and 95% CI)",
                        linetype = "EPS (fit and 95% CI)")) +
        geom_line(aes(x = bnpr_x,y = bnpr_y,
                      fill = "EPS (fit and 95% CI)",
                      shape = "EPS (fit and 95% CI)",
                      colour = "EPS (fit and 95% CI)",
                      linetype = "EPS (fit and 95% CI)",
                      size = "EPS (fit and 95% CI)")) +
        geom_line(data = vaf_trajectory,
                  aes(x = x,y = y,
                      fill = "VAF (inferred)",
                      colour = "VAF (inferred)",
                      shape = "VAF (inferred)",
                      linetype = "VAF (inferred)",
                      size = "VAF (inferred)")) + 
        geom_ribbon(data = vaf_trajectory,
                    aes(x = x,ymin = y_05,ymax = y_95,
                        fill = "VAF (interval)",
                        colour = "VAF (interval)",
                        shape = "VAF (interval)",
                        linetype = "VAF (interval)",
                        size = "VAF (interval)"),
                    alpha = 0.4) + 
        geom_point(inherit.aes = F,
                   data = tmp_data,
                   aes(x = Age,
                       y = VAF,
                       size = "Data",
                       shape = "Data",
                       fill = "Data",
                       colour = "Data",
                       linetype = "Data")) +
        geom_linerange(inherit.aes = F,
                       data = tmp_data,
                       linetype = 1,size = 0.25,
                       aes(x = Age,
                           ymin = VAF005,
                           ymax = VAF095,
                           shape = "Data",
                           fill = "Data",
                           colour = "Data")) +
        geom_hline(yintercept = 2e5,colour = "black",linetype = 2) +
        theme_gerstung(base_size = 6) + 
        scale_y_continuous(
          trans = 'log10',
          sec.axis = sec_axis(trans = ~ .,
                              name = "Variant allele frequency",
                              breaks = y_axis_at_vaf * scaling_factor,
                              labels = scientific(y_axis_at_vaf)),
          breaks = y_axis_at,
          labels = scientific) + 
        scale_x_continuous(limits = c(0,max_age),
                           breaks = seq(0,100,by = 20) + 0.75,
                           labels = seq(0,100,by = 20),
                           expand = c(0,0,0,0)) +
        ylab("Effective population size") + 
        xlab("Age") + 
        coord_cartesian(ylim = y_lims) + 
        scale_colour_manual(values = c(`EPS (fit and 95% CI)` = colour_for_trajectory,
                                       `VAF (inferred)` = "black",
                                       `VAF (interval)` = NA,
                                       `Data` = "black",
                                       `EPS (at coalescence)` = colour_for_trajectory),name = NULL,
                            breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)",
                                       "VAF (interval)",
                                       "VAF (inferred)","Data")) + 
        scale_fill_manual(values = c(`EPS (fit and 95% CI)` = colorspace::lighten(colour_for_trajectory,0.8),
                                     `VAF (inferred)` = NA,
                                     `VAF (interval)` = "grey",
                                     `Data` = NA,
                                     `EPS (at coalescence)` = NA),name = NULL,
                          breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)",
                                     "VAF (interval)",
                                     "VAF (inferred)","Data")) + 
        scale_shape_manual(values = c(`EPS (fit and 95% CI)` = NA,
                                      `VAF (inferred)` = NA,
                                      `VAF (interval)` = NA,
                                      `Data` = 16,
                                      `EPS (at coalescence)` = 16),
                           name = NULL,
                           breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)",
                                      "VAF (interval)",
                                      "VAF (inferred)","Data")) + 
        scale_size_manual(values = c(`EPS (fit and 95% CI)` = 0.5,
                                     `VAF (inferred)` = 0.5,
                                     `VAF (interval)` = 0,
                                     `Data` = 1,
                                     `EPS (at coalescence)` = 1),
                          name = NULL,
                          breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)",
                                     "VAF (interval)",
                                     "VAF (inferred)","Data")) + 
        scale_linetype_manual(values = c(`EPS (fit and 95% CI)` = 1,
                                         `VAF (inferred)` = 3,
                                         `VAF (interval)` = 0,
                                         `Data` = 0,
                                         `EPS (at coalescence)` = 0),
                              name = NULL,
                              breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)",
                                         "VAF (interval)",
                                         "VAF (inferred)","Data")) + 
        labs(title = title_label) +
        theme(legend.position = legend_position,
              legend.direction = "vertical",
              legend.justification = c(0.5,1),
              legend.text.align = 0,
              legend.key.height = unit(0.35,"cm"),
              legend.text = element_text(size = 6),
              axis.title = element_text(size = 6),
              plot.title = element_text(size = 6,margin = margin(b = 0)),
              axis.ticks.y.left = element_line(colour = colour_for_trajectory),
              axis.text.y.left = element_text(colour = colour_for_trajectory),
              axis.line.y.left = element_line(colour = colour_for_trajectory),
              axis.title.y.left = element_text(colour = colour_for_trajectory))
      
      trajectory_plot <- trajectory_plot_ +
        annotate("text",
                 label = c("Growth/yr\n\n",
                           sprintf("\n%.2f (phylogeny)",linear_estimate$coefficients[2]),
                           sprintf("\n\n%.2f (time series)",exp(b)-1)),
                 x = annotate_x,
                 y = y_lims[1],
                 colour = c(colour_for_trajectory,colour_for_trajectory,"black"), 
                 hjust = annotate_hjust,
                 size = 6 * tm,
                 vjust = 0)
      
      trajectory_plot <- trajectory_plot_ +
        annotate("text",
                 label = c("Growth/yr\n\n",
                           sprintf("\n%.2f(phylogeny)",linear_estimate$coefficients[2]),
                           sprintf("\n\n%.2f (time series)",exp(b)-1)),
                 max_age * 0.02,
                 y = y_lims[2],
                 colour = c(colour_for_trajectory,colour_for_trajectory,"black"), 
                 hjust = 0,
                 size = 6 * tm,
                 vjust = 0.9)
      
      C <- as.data.frame(coefficients(summary(linear_estimate)))
      rg <- summary(clade$relevant_growth)$coefficients
      coefficient_df <- data.frame(
        b_mean = c(exp(b)-1,C$Estimate[2]),
        b_q05 = c(exp(b_05)-1,C$Estimate[2] - C$`Std. Error`[2]),
        b_q95 = c(exp(b_95)-1,C$Estimate[2] + C$`Std. Error`[2]),
        source = c("VAF","EPS"),
        expected_late = c(NA,1/clade$non_linear_estimate_vaf$m$getPars()[2]),
        b_early = c(NA,cp[2]),
        b_late = c(NA,cp[3]),
        ci_late_low = clade$late_ci[1],
        ci_late_high = clade$late_ci[2],
        midpoint = c(NA,cp[4]),
        relevant_growth = rg[2,1],
        relevant_growth_se = rg[2,2],
        observed_vaf = observed_vaf,
        expected_vaf = expected_vaf,
        age_mean = c(tmp$Mean,mean(mutation_timing$x)),
        age_q05 = c(tmp$Q05,mutation_timing$x[1]),
        age_q95 = c(tmp$Q95,mutation_timing$x[2]),
        gene = driver_id_split[1],
        site = driver_id_split[2],
        individual = curr_id
      ) 
      ages_plot <- data.frame(
        b_mean = c(b,C$Estimate[2]),
        b_q05 = c(b_05,C$Estimate[2] - C$`Std. Error`[2]),
        b_q95 = c(b_95,C$Estimate[2] + C$`Std. Error`[2]),
        source = c("VAF","EPS")
      )
      
    } else {
      bnpr_inferred_trajectory <- data.frame(
        x = bnpr_x,
        y = exp(predict(linear_estimate,x = bnpr_x))
      )
      
      bnpr_x_at_coalescence <- -clade$bnpr_estimate$coal_times+max_age
      bnpr_y_at_coalescence <- approx(
        y = bnpr_y,
        x = bnpr_x,
        xout = bnpr_x_at_coalescence
      )$y
      
      bnpr_at_coalescence <- data.frame(
        x = bnpr_x_at_coalescence,
        y = bnpr_y_at_coalescence
      ) %>% 
        na.omit()

      y_lims <- c(
        1,
        max(bnpr_y) * 5
      )
      y_lims <- c(1,1e6)
      x_axis_at <- unique(floor(bnpr_x / 10) * 10) %>%
        Filter(f = function(x) x > min(bnpr_x))
      if (length(x_axis_at) < 4) {
        x_axis_at <- unique(floor(bnpr_x / 5) * 5) %>%
          Filter(f = function(x) x > min(bnpr_x))
      }
      y_axis_at <- unique(10^floor(log(bnpr_y))) %>%
        Filter(f = function(x) x > y_lims[1])
      y_axis_at <- c(1,y_axis_at)
      
      age_midpoint <- max_age / 2
      if (min(bnpr_x) < age_midpoint) {
        annotate_hjust = 1
        annotate_x = max_age
      } else {
        annotate_hjust = 0
        annotate_x = 0
      }
      
      mutation_timing <- data.frame(
        x = c(mutation_timings[[driver_id]]),
        y = 0
      )
      
      title_label <- sprintf("%s",driver_id,curr_id)
      
      trajectory_plot_ <- data.frame(bnpr_x = bnpr_x,bnpr_y = bnpr_y,
                                    bnpr_y025 = clade$bnpr_estimate$summary$quant0.025,
                                    bnpr_y975 = clade$bnpr_estimate$summary$quant0.975
      ) %>%
        ggplot() +
        geom_line(data = mutation_timing,
                  inherit.aes = F,
                  size = 3,
                  colour = colorspace::lighten("black",0.7),
                  aes(x = x,y = y)) + 
        geom_ribbon(aes(x = bnpr_x,ymin = bnpr_y025,ymax = bnpr_y975,
                        fill = "EPS (fit and 95% CI)",
                        shape = "EPS (fit and 95% CI)")) +
        geom_line(aes(x = bnpr_x,y = bnpr_y,
                      colour = "EPS (fit and 95% CI)",
                      size = "EPS (fit and 95% CI)",
                      linetype = "EPS (fit and 95% CI)")) +
        theme_gerstung(base_size = 6) + 
        scale_y_continuous(
          breaks = y_axis_at,
          labels = scientific,
          trans = 'log10') +
        scale_x_continuous(limits = c(0,max_age),
                           breaks = seq(0,100,by = 20) + 0.75,
                           labels = seq(0,100,by = 20),
                           expand = c(0,0,0,0)) +
        ylab("Effective population size") + 
        xlab("Age") + 
        coord_cartesian(ylim = y_lims) + 
        geom_hline(yintercept = 2e5,colour = "black",linetype = 2) +
        scale_colour_manual(values = c(`EPS (fit and 95% CI)` = "grey60",
                                       `EPS (at coalescence)` = "grey60"),
                            name = NULL,
                            breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)")) + 
        scale_fill_manual(values = c(`EPS (fit and 95% CI)` = "grey90",
                                     `EPS (at coalescence)` = NA),
                          name = NULL,
                          breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)")) + 
        scale_shape_manual(values = c(`EPS (fit and 95% CI)` = NA,
                                      `EPS (at coalescence)` = 16),
                           name = NULL,
                           breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)")) + 
        scale_size_manual(values = c(`EPS (fit and 95% CI)` = 0.5,
                                     `EPS (at coalescence)` = 1),
                          name = NULL,
                          breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)")) + 
        scale_linetype_manual(values = c(`EPS (fit and 95% CI)` = 1,
                                         `EPS (at coalescence)` = 0),
                              name = NULL,
                              breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)")) + 
        labs(title = title_label) +
        theme(legend.position = legend_position,
              legend.direction = "vertical",
              legend.justification = c(0.5,1),
              legend.text.align = 0,
              legend.key.height = unit(0.35,"cm"),
              axis.title = element_text(size = 6),
              legend.text = element_text(size = 6),
              plot.title = element_text(size = 6,margin = margin(b = 0)))
      
      trajectory_plot <- trajectory_plot_ + 
        annotate("text",
                 label = c("\nGrowth/yr\n",
                           sprintf("\n\n%.2f (phylogeny)",linear_estimate$coefficients[2])),
                 x = annotate_x,
                 y = y_lims[1],
                 colour = c("black"), 
                 hjust = annotate_hjust,
                 size = 6 * tm,
                 vjust = 0) 
      
      trajectory_plot <- trajectory_plot_ + 
        annotate("text",
                 label = c("Growth/yr\n\n",
                           sprintf("\n%.2f (phylogeny)",linear_estimate$coefficients[2]),
                           '\n\n0'),
                 x = max_age * 0.02,
                 y = y_lims[2],
                 colour = c("grey60","grey60","white"), 
                 hjust = 0,
                 size = 6 * tm,
                 vjust = 0.9)
      
      C <- summary(linear_estimate) %>% 
        coefficients %>%
        as.data.frame()
      rg <- summary(clade$relevant_growth)$coefficients
      coefficient_df <- data.frame(
        b_mean = c(C$Estimate[2]),
        b_q05 = c(C$Estimate[2] - C$`Std. Error`[2]),
        b_q95 = c(C$Estimate[2] + C$`Std. Error`[2]),
        source = c("EPS"),
        expected_late = 1/c(clade$non_linear_estimate_vaf$m$getPars()[2]),
        b_early = c(cp[2]),
        b_late = c(cp[3]),
        ci_late_low = clade$late_ci[1],
        ci_late_high = clade$late_ci[2],
        midpoint = c(cp[4]),
        relevant_growth = rg[2,1],
        relevant_growth_se = rg[2,2],
        observed_vaf = observed_vaf,
        expected_vaf = expected_vaf,
        age_mean = c(mean(mutation_timing$x)),
        age_q05 = c(mutation_timing$x[1]),
        age_q95 = c(mutation_timing$x[2]),
        gene = "",
        site = driver_id,
        individual = curr_id) 
    }
    plot_list[[curr_id]][[sprintf('%s',driver_id)]] <- trajectory_plot
    plot_list[[curr_id]][[sprintf('%s_age_coef',driver_id)]] <- coefficient_df
  }
  plot_list[[curr_id]][["tree"]] <- tree_plot
}


# finish and save plots ---------------------------------------------------

all_trajectories <- list()

for (individual in names(bnpr_estimates)) {
  bnpr_estimate <- bnpr_estimates[[individual]]
  for (clade in bnpr_estimate$clade_dynamics) {
    tip_labels <- bnpr_estimate$tree_ultra$tip.label[clade$clade]
    clone_idx <- which(names(clone_assignment[[individual]]) %in% tip_labels)
    driver <- clone_assignment[[individual]][[clone_idx]]
    all_trajectories[[length(all_trajectories) + 1]] <- data.frame(
      X = ages[[individual]] - clade$bnpr_estimate$summary$time,
      Y = clade$bnpr_estimate$summary$mean,
      W = (log(clade$bnpr_estimate$summary$quant0.975) - log(clade$bnpr_estimate$summary$quant0.025))^2/16,
      individual = individual,
      driver = driver
    )
  }
}

all_trajectories_df <- all_trajectories %>% do.call(what = rbind) %>%
  mutate(alpha_W = W <= 5) 

all_trajectories_df <- all_trajectories_df %>% 
  group_by(driver,individual) %>% 
  mutate(X = X - min(X)) %>% 
  mutate(gene = str_match(driver,"[A-Z0-9]+"))

all_trajectories_df %>%
  subset(alpha_W == T) %>% 
  ggplot(aes(x = X,y = Y,group = driver,colour = gene)) + 
  geom_rect(aes(xmin = 0,xmax = max(all_trajectories_df$X),
                ymin = 50e3,ymax = 200e3),fill = "grey92",colour = NA) +
  geom_line(data = all_trajectories_df,linetype = 3,size = 0.25) + 
  geom_line(size = 0.25) + 
  facet_wrap(~ individual,ncol = 2) + 
  theme_gerstung(base_size = 6) + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(trans = 'log10',expand = c(0,0),breaks = c(1e-2,1,1e2,1e4,1e6)) +
  coord_cartesian(ylim = c(0.5,1e7)) +
  theme(strip.text = element_text(margin = margin(b = 0.5))) +
  xlab("Time since first coal. (years)") +
  ylab("Neff") +
  scale_color_manual(values = c(U = "grey40",gene_colours),
                     guide = F) + 
  ggsave("figures/model_ch/trees/tree_trajectories.pdf",
         height = 1.5,
         width = 1.8)

coefficients_df <- plot_list %>%
  lapply(function(x) {
    grep("*age_coef",names(x)) %>%
      as.vector() %>%
      lapply(function(y) {
        return(x[[y]])
      }) %>%
      do.call(what = rbind)
  }) %>%
  do.call(what = rbind) %>% 
  as.data.frame() %>% 
  mutate(plot_label = ifelse(gene != "",sprintf('%s-%s',gene,site),as.character(site))) %>%
  mutate(plot_label = gsub("Unknown driver ","UD",plot_label)) %>% 
  mutate(colour_code = ifelse(gene != "",names(gene_colours[as.character(gene)]),"Unknown")) %>%
  mutate(source = ifelse(source == "EPS","Colonies","Longitudinal")) %>%
  mutate(age_q05 = age_q05 - ifelse(source == "Colonies",0.75,0),
         age_q95 = age_q95 - ifelse(source == "Colonies",0.75,0)) %>%
  mutate(late_growth = b_early + b_late) %>%
  mutate(b = ifelse(source == "Longitudinal",b_mean,relevant_growth),
         b05 = ifelse(source == "Longitudinal",b_q05,relevant_growth-relevant_growth_se*1.96),
         b95 = ifelse(source == "Longitudinal",b_q95,relevant_growth+relevant_growth_se*1.96))

site_labels <- sapply(
  strsplit(rev(c("SF3B1-K666N","U2AF1-Q157R","UD1","UD2","UD3")), "-"), 
  function(x) ifelse(
    length(x) > 1,
    parse(text = sprintf("atop(atop(textstyle(italic('%s')),textstyle('%s')),'')",x[1],x[2])),
    x[1]
  )
)

coefficients_df_long <- rbind(
  transmute(
    coefficients_df,
    mean = b,
    q05 = b05,
    q95 = b95,
    source,gene,site,individual,plot_label,colour_code,
    id = "Growth per year"
  ),
  transmute(
    coefficients_df,
    mean = age_mean,
    q05 = age_q05,
    q95 = age_q95,
    source,gene,site,individual,plot_label,colour_code,
    id = "Age at onset"
  )
) %>%
  mutate(id = factor(id,levels = c("Growth per year","Age at onset"))) %>%
  mutate(individual = sprintf("id = %s",individual)) 

coefficient_age_comparison_plot <- coefficients_df_long %>% 
  mutate(plot_label = factor(
    plot_label,
    levels = rev(c("SF3B1-K666N","U2AF1-Q157R","UD1","UD2","UD3")))) %>%
  ggplot(aes(x = as.numeric(as.factor(plot_label)),
             y = mean,
             ymin = q05,
             ymax = q95,
             group = paste(individual,source,site),
             colour = colour_code)) + 
  geom_point(aes(shape = source),position = position_dodge(width = 0.8),
             size = 1.5) + 
  geom_linerange(position = position_dodge(width = 0.8)) + 
  theme_gerstung(base_size = 6) + 
  xlab("Driver") + 
  ylab("Growth per year") + 
  scale_shape_discrete(name = NULL) +
  scale_colour_manual(values = c(gene_colours,Unknown = "black"),
                      guide = F) + 
  scale_x_continuous(labels = site_labels,
                     breaks = c(1,2,3,4,5),
                     trans = "reverse",
                     expand = c(0,0)) + 
  facet_grid(id ~ individual,
             scales = "free",switch = "y",space = "free_x") + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 6),
        strip.placement = "outside",
        legend.position = "bottom",
        legend.spacing = unit(0,"cm"),
        legend.box.margin = margin(),
        legend.margin = margin(),
        legend.background = element_blank(),
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 6,margin = margin(b = 0))) + 
  xlab("") + 
  ylab("") +
  theme(strip.text.y = element_text(angle = 0))

linear_coefficients_only_plot <- coefficients_df %>% 
  mutate(plot_label = factor(plot_label,
    levels = rev(c("SF3B1-K666N","U2AF1-Q157R","UD1","UD2","UD3")))) %>%
  mutate(individual = factor(individual,levels = c(2259,3877,500))) %>% 
  ggplot(aes(x = as.numeric(as.factor(plot_label)),
             y = b_mean,
             ymin = b_q05,
             ymax = b_q95,
             group = paste(individual,source,site),
             colour = colour_code)) + 
  geom_point(aes(shape = source),position = position_dodge(width = 0.8),
             size = 1.5) + 
  geom_linerange(position = position_dodge(width = 0.8)) + 
  theme_gerstung(base_size = 6) + 
  xlab("Driver") + 
  ylab("Growth per year") + 
  scale_shape_discrete(name = NULL) +
  scale_colour_manual(values = c(gene_colours,Unknown = "black"),
                      guide = F) + 
  scale_x_continuous(labels = site_labels,
                     breaks = c(1,2,3,4,5),
                     trans = "reverse",
                     expand = c(0,0)) + 
  facet_grid( ~ individual,
             scales = "free",switch = "y",space = "free_x") + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 6),
        strip.placement = "outside",
        legend.position = "bottom",
        legend.spacing = unit(0,"cm"),
        legend.box.margin = margin(),
        legend.margin = margin(),
        legend.background = element_blank(),
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 6,margin = margin(b = 0))) + 
  xlab("") + 
  ylab("") +
  theme(strip.text.y = element_text(angle = 0))

load_data_6th_tp() %>%
  subset(SardID %in% c(names(trees))) %>% 
  mutate(site = as.character(str_match(amino_acid_change,"[A-Z0-9]+$")),
         individual = as.character(SardID)) %>% 
  merge(mutate(coefficients_df,individual = as.character(individual),
               site = as.character(site)),
        by = c("individual","site"),all = F) %>% 
  select(gene,site,individual,observed_vaf,VAF) %>%
  distinct %>%
  mutate(observed_vaf = observed_vaf / 2) %>%
  rowwise() %>%
  mutate(tree_size = length(trees[[as.character(individual)]]$tree_ultra$tip.label)) %>% 
  transmute(
    site = site,
    gene = gene,
    individual = individual,
    tree = observed_vaf,
    tree_05 = qbeta(0.05,observed_vaf*tree_size,tree_size-observed_vaf*tree_size),
    tree_95 = qbeta(0.95,observed_vaf*tree_size,tree_size-observed_vaf*tree_size),
    long = VAF,
    long_05 = qbeta(0.05,VAF*1e3,1e3-VAF*1e3),
    long_95 = qbeta(0.95,VAF*1e3,1e3-VAF*1e3)
  ) %>% 
  rbind(
    data.frame(
      site = "P95H",
      gene = "SRSF2",
      individual = "500",
      tree = 1 / 92,
      tree_05 = qbeta(0.05,1,92-1),
      tree_95 = qbeta(0.95,1,92-1),
      long = 0.07355,
      long_05 = qbeta(0.05,0.07355*1e3,1e3-0.07355*1e3),
      long_95 = qbeta(0.95,0.07355*1e3,1e3-0.07355*1e3)
    )
  ) %>% 
  ggplot(aes(x = long,y = tree)) + 
  geom_abline(slope = 1) +
  geom_point(size = 0.5) +
  geom_linerange(aes(ymin = tree_05,ymax = tree_95)) +
  geom_errorbarh(aes(xmin = long_05,xmax = long_95),height = 0) +
  geom_text_repel(aes(label = sprintf("%s-%s (%s)",gene,site,individual)),
                  size = 2.7) + 
  coord_cartesian(xlim = c(0.001,0.5),ylim = c(0.001,0.5)) + 
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') + 
  xlab("Clone size (sequencing)") + 
  ylab("Clone size (fraction of tips in clade)") + 
  theme_gerstung(base_size = 6) +
  ggsave(filename = sprintf("figures/%s/trees/compare_vaf.pdf",model_id),
         width = 2,height = 2,
         useDingbats = F)

output_width = 12
output_height = 4.5

plot_3877 <- plot_grid(
  plot_list$`3877`$tree,
  plot_grid(
    plot_list$`3877`$`SF3B1-K666N`,
    plot_list$`3877`$`U2AF1-Q157R`,
    ncol = 1,
    align='hv'
  ),
  nrow = 1,
  rel_widths = c(1.0,0.7)
) 
plot_2259 <- plot_grid(
  plot_list$`2259`$tree,
  plot_list$`2259`$`SF3B1-K666N` + 
    theme(legend.position = c(legend_position[1],1.05)),
  rel_widths = c(1.0,0.7),
  nrow = 1)

plot_500 <- plot_grid(
  plot_list$`500`$tree,
  plot_grid(
    plot_list$`500`$`Unknown driver 1`,
    plot_list$`500`$`Unknown driver 2`,
    plot_list$`500`$`Unknown driver 3`,
    ncol = 1,
    align='hv'
  ),
  nrow = 1,
  rel_widths = c(1.0,0.7)
)

pos_adjust = "none"
ts <- 6
tree_plot_reordered <- plot_grid(
  plot_grid(
    plot_list$`2259`$tree + theme(plot.title = element_text(size = ts),
                                  axis.text = element_text(size = ts),
                                  plot.margin = margin(0.3,0.3,0.3,0.3)),
    plot_list$`2259`$`SF3B1-K666N` + theme(legend.position = pos_adjust,
                                           axis.text = element_text(size = ts),
                                           axis.title = element_text(size = ts)),
    get_legend(plot_list$`2259`$`SF3B1-K666N` + theme(legend.key.height = unit(0.6,"cm"))),
    ncol = 1,
    align = "v",
    axis = "tblr",
    rel_heights = c(2,1,1)
  ),
  ggplot() + theme_nothing(),
  plot_grid(
    plot_list$`3877`$tree + theme(plot.title = element_text(size = ts),
                                  axis.text = element_text(size = ts),
                                  plot.margin = margin(0.3,0.3,0.3,0.3)),
    plot_list$`3877`$`SF3B1-K666N` + theme(legend.position = pos_adjust,
                                           axis.text = element_text(size = ts),
                                           axis.title = element_text(size = ts)),
    plot_list$`3877`$`U2AF1-Q157R` + theme(legend.position = pos_adjust,
                                           axis.text = element_text(size = ts),
                                           axis.title = element_text(size = ts)),
    ncol = 1,
    align = "v",
    axis = "tblr",
    rel_heights = c(2,1,1)
  ),
  ggplot() + theme_nothing(),
  plot_grid(
    plot_list$`500`$tree + theme(plot.title = element_text(size = ts),
                                 axis.text = element_text(size = ts),
                                 plot.margin = margin(0.3,0.3,0.3,0.3)),
    plot_list$`500`$`Unknown driver 1` + theme(legend.position = pos_adjust,
                                               axis.text = element_text(size = ts),
                                               axis.title = element_text(size = ts)),
    plot_list$`500`$`Unknown driver 3` + theme(legend.position = pos_adjust,
                                               axis.text = element_text(size = ts),
                                               axis.title = element_text(size = ts)),
    ncol = 1,
    align = "v",
    axis = "tblr",
    rel_heights = c(2,1,1)
  ),
  ncol = 5,
  rel_widths = c(1,0.04,1,0.04,0.84)
) 

tree_plot_reordered %>% 
  ggsave(filename = sprintf("figures/%s/trees/all_trees_and_dynamics_full.pdf",model_id),
         height = 8,
         width = 8,
         useDingbats = F)
tree_plot_reordered %>% 
  ggsave(filename = sprintf("figures/%s/trees/all_trees_and_dynamics_full_redux.pdf",model_id),
         height = 5,
         width = 5,
         useDingbats = F)

ggsave(plot = coefficient_age_comparison_plot,
       filename = sprintf("figures/%s/trees/coef_age_plot.pdf",model_id),
       width = 8.6 / 3,height = 2.7,
       useDingbats = F)

ggsave(plot = linear_coefficients_only_plot + ylab("Growth per year"),
       filename = sprintf("figures/%s/trees/linear_coef_only_plot.pdf",model_id),
       width = 2.7,height = 2,
       useDingbats = F)

write.csv(coefficients_df,"data_output/tree_bnpr_coefficients.csv",row.names = T)
