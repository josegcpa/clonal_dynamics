legend_position <- c(0.33,1.11)

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
  tree_groups <- groupOTU(tree_ultra,clades_plot)
  tree_plot <- ggtree(tree_groups,aes(colour=str_match(group,'[A-Z0-9]+')),size = 0.7,
                      ladderize=T) + 
    scale_colour_manual(values = c(`0` = 'grey80',U = 'black',gene_colours),
                        guide=F) + 
    theme_gerstung(base_size = 15) +
    xlab("Age") + 
    ggtitle(sprintf("SardID %s",curr_id)) +
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 11)) + 
    coord_flip() + 
    scale_x_continuous(trans = "reverse",
                       expand = c(0.001,0),
                       limits = c(max_age,0) / max_age,
                       breaks = c(0,20,40,60,80) / max_age,
                       labels = c(0,20,40,60,80))
  
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
      'atop(%s,"%s")',
      edge_labels$top,
      edge_labels$bottom),
    NA
  )
  tree_plot <- tree_plot %<+% edge_labels + 
    geom_label(aes(x=branch, label=edge_num),
               parse = T,
               label.r = unit(0,"cm"),label.size = unit(0.01,"cm"),
               colour = "black")
  
  for (clade in clade_dynamics) {
    driver_id <- clone_assignment[[curr_id]][names(clone_assignment[[curr_id]]) %in% clade$tree_subset$tip.label] %>%
      unlist
    linear_estimate <- clade$linear_estimate
    tree_ultra <- O$tree_ultra
    tree_ultra$edge.length <- tree_ultra$edge.length * max_age
    
    bnpr_x <- -clade$bnpr_estimate$summary$time + max_age
    bnpr_y <- clade$bnpr_estimate$summary$mean
    if (!grepl("Unknown driver",driver_id)) {
      tmp <- Parameters_Age %>% 
        subset(individual == as.numeric(curr_id) & site == driver_id)
      tmp_data <- full_data %>%
        subset(SardID == as.numeric(curr_id) & amino_acid_change == driver_id)

      u <- (tmp$u_values_095 + tmp$u_values_005)/2
      u_005 <- tmp$u_values_005
      u_095 <- tmp$u_values_095
      b <- tmp$b_site_mean + tmp$b_gene_mean + tmp$b_clone_mean
      b_05 <- tmp$b_site_005 + tmp$b_gene_005 + tmp$b_clone_005
      b_95 <- tmp$b_site_095 + tmp$b_gene_095 + tmp$b_clone_095
      
      gen_per_year <- 13
      fitness <- b / gen_per_year
      pop_size <- 2e5 * 2
      drift_threshold <- 1 / (fitness * pop_size)
      time_at_drift <- (logit(drift_threshold) - u)/b + values_model$min_age
      time_at_onset <- t0_adjusted(u,b,13,pop_size) + values_model$min_age
      
      y <- inv.logit(b * (bnpr_x - values_model$min_age) + u) / 2
      y_05 <- inv.logit(b_05 * (bnpr_x - values_model$min_age) + u_005) / 2
      y_95 <- inv.logit(b_95 * (bnpr_x - values_model$min_age) + u_095) / 2
      vaf_trajectory <- data.frame(
        x = c(
          seq(time_at_onset,time_at_drift,length.out = 100),
          seq(time_at_drift,max_age,length.out = 100)),
        y = 0.5 * c(
          seq(1 / 2e5,drift_threshold,length.out = 100),
          inv.logit(b * (seq(time_at_drift,max_age,length.out = 100)-values_model$min_age) + u)))
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
      
      y_ <- inv.logit(b * (bnpr_at_coalescence$x - values_model$min_age) + u) / 2
      
      optim_func <- function(params) {
        scaling_factor <- params[1]
        return(mean((y_ * scaling_factor - bnpr_at_coalescence$y)^2))
      }
      solution <- optim(par = c(scaling_factor = 1e3),optim_func,method = "Brent",
                        lower = 100, upper = 1e10)
      
      vaf_trajectory$y <- vaf_trajectory$y * solution$par
    
      tmp_data <- tmp_data %>%
        mutate(VAF = solution$par * VAF)
      y_lims <- c(
        min(vaf_trajectory$y),#max(1,min(c(y * solution$par,bnpr_y,tmp_data$VAF))),
        max(c(y * solution$par,bnpr_y,tmp_data$VAF)) * 5
      )
      
      x_axis_at <- unique(floor(bnpr_x / 10) * 10) %>%
        Filter(f = function(x) x > min(bnpr_x))
      if (length(x_axis_at) < 4) {
        x_axis_at <- unique(floor(bnpr_x / 5) * 5) %>%
          Filter(f = function(x) x > min(bnpr_x))
      }
      y_axis_at <- unique(10^floor(log(bnpr_y))) %>%
        Filter(f = function(x) x > y_lims[1])
      y_axis_at <- c(1,y_axis_at)
      y_axis_at_vaf <- c(10^c(-6:-1),0.5)
      
      colour_for_trajectory <- ifelse(
        grepl('-',driver_id),
        gene_colours[str_match(driver_id,'[0-9A-Z]+')],
        'black'
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
          italic(gene),'-', site
        )),
        list(
          gene = driver_id_split[1],
          site = driver_id_split[2]
        )
      ) %>%
        eval()
      expression(paste("Sepal length and sepal width of various ",italic("Iris")," species"))

      trajectory_plot <- data.frame(
        y = y * solution$par,
        bnpr_x = bnpr_x,bnpr_y = bnpr_y,
        bnpr_y025 = clade$bnpr_estimate$summary$quant0.025,
        bnpr_y975 = clade$bnpr_estimate$summary$quant0.975
      ) %>%
        ggplot() +
        geom_line(data = mutation_timing,
                  inherit.aes = F,
                  size = 5,
                  colour = colorspace::lighten(colour_for_trajectory,0.7),
                  aes(x = x,y = y)) + 
        geom_ribbon(aes(x = bnpr_x,
                        ymin = bnpr_y025,
                        ymax = bnpr_y975,
                        fill = "EPS (fit and 95% CI)",
                        shape = "EPS (fit and 95% CI)",
                        linetype = "EPS (fit and 95% CI)")) +
        # geom_line(aes(x = bnpr_x,y = bnpr_y,colour = "EPS (fit and 95% CI)",
        #               shape = "EPS (fit and 95% CI)",
        #               linetype = "EPS (fit and 95% CI)"),
        #           alpha = 0.3,
        #           size = 1) + 
        geom_line(aes(x = x,y = y,
                      fill = "EPS (fit and 95% CI)",
                      shape = "EPS (fit and 95% CI)",
                      colour = "EPS (fit and 95% CI)",
                      linetype = "EPS (fit and 95% CI)",
                      size = "EPS (fit and 95% CI)"),
                  data = bnpr_inferred_trajectory) +
        geom_point(aes(x = x,y = y,
                       colour = "EPS (at coalescence)",
                       shape = "EPS (at coalescence)",
                       linetype = "EPS (at coalescence)",
                       size = "EPS (at coalescence)",
                       fill = "EPS (at coalescence)"),
                   data = bnpr_at_coalescence) +
        geom_line(data = vaf_trajectory,
                  aes(x = x,y = y,
                      fill = "VAF (inferred)",
                      colour = "VAF (inferred)",
                      shape = "VAF (inferred)",
                      linetype = "VAF (inferred)",
                      size = "VAF (inferred)")) + 
        geom_point(inherit.aes = F,
                   data = tmp_data,
                   aes(x = Age,
                       y = VAF,
                       size = "Data",
                       shape = "Data",
                       fill = "Data",
                       colour = "Data",
                       linetype = "Data")) +
        annotate("text",
                 label = c(sprintf("Growth/year (VAF) = %.2f\n",b),
                           sprintf("\nGrowth/year (EPS) = %.2f",linear_estimate$coefficients[2])),
                 x = annotate_x,
                 y = y_lims[1],
                 colour = c("black",colour_for_trajectory), 
                 hjust = annotate_hjust,
                 size = 3.5,
                 fontface = 'bold',
                 vjust = 0) +  
        theme_gerstung(base_size = 15) + 
        scale_y_continuous(
          trans = 'log10',
          sec.axis = sec_axis(trans = ~ .,
                              name = "Variant allele frequency\n(VAF)",
                              breaks = y_axis_at_vaf * solution$par,
                              labels = y_axis_at_vaf),
          breaks = y_axis_at,
          labels = y_axis_at) + 
        scale_x_continuous(limits = c(0,max_age),expand = c(0.02,0,0.02,0)) +
        ylab("Effective  pop. size (EPS)") + 
        xlab("Age") + 
        coord_cartesian(ylim = y_lims) + 
        scale_colour_manual(values = c(`EPS (fit and 95% CI)` = colour_for_trajectory,
                                       `VAF (inferred)` = "black",
                                       `Data` = "black",
                                       `EPS (at coalescence)` = colour_for_trajectory),name = NULL,
                            breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)","VAF (inferred)","Data")) + 
        scale_fill_manual(values = c(`EPS (fit and 95% CI)` = colorspace::lighten(colour_for_trajectory,0.8),
                                     `VAF (inferred)` = NA,
                                     `Data` = NA,
                                     `EPS (at coalescence)` = NA),name = NULL,
                          breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)","VAF (inferred)","Data")) + 
        scale_shape_manual(values = c(`EPS (fit and 95% CI)` = NA,
                                      `VAF (inferred)` = NA,
                                      `Data` = 1,
                                      `EPS (at coalescence)` = 16),
                           name = NULL,
                           breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)","VAF (inferred)","Data")) + 
        scale_size_manual(values = c(`EPS (fit and 95% CI)` = 1,
                                     `VAF (inferred)` = 1,
                                     `Data` = 3,
                                     `EPS (at coalescence)` = 1.5),
                          name = NULL,
                          breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)","VAF (inferred)","Data")) + 
        scale_linetype_manual(values = c(`EPS (fit and 95% CI)` = 1,
                                         `VAF (inferred)` = 3,
                                         `Data` = 0,
                                         `EPS (at coalescence)` = 0),
                              name = NULL,
                              breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)","VAF (inferred)","Data")) + 
        labs(title = title_label) +
        theme(legend.position = legend_position,
              legend.direction = "vertical",
              legend.justification = c(0.5,1),
              legend.text.align = 0,
              plot.title = element_text(size = 12),
              axis.ticks.y.left = element_line(colour = colour_for_trajectory),
              axis.text.y.left = element_text(colour = colour_for_trajectory),
              axis.line.y.left = element_line(colour = colour_for_trajectory),
              axis.title.y.left = element_text(colour = colour_for_trajectory))
      
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
        1,#min(bnpr_y),
        max(bnpr_y) * 5
      )
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
      
      trajectory_plot <- data.frame(bnpr_x = bnpr_x,bnpr_y = bnpr_y,
                                    bnpr_y025 = clade$bnpr_estimate$summary$quant0.025,
                                    bnpr_y975 = clade$bnpr_estimate$summary$quant0.975
      ) %>%
        ggplot() +
        geom_line(data = mutation_timing,
                  inherit.aes = F,
                  size = 5,
                  colour = colorspace::lighten("black",0.7),
                  aes(x = x,y = y)) + 
        geom_ribbon(aes(x = bnpr_x,ymin = bnpr_y025,ymax = bnpr_y975,
                        fill = "EPS (fit and 95% CI)",
                        shape = "EPS (fit and 95% CI)")) +
        # geom_line(aes(x = bnpr_x,y = bnpr_y,
        #               colour = "EPS (fit and 95% CI)"),
        #           size = 1,
        #           alpha = 0.3) + 
        geom_point(aes(x = x,y = y,
                       size = "EPS (at coalescence)",
                       shape = "EPS (at coalescence)",
                       linetype = "EPS (at coalescence)",
                       fill = "EPS (at coalescence)",
                       colour = "EPS (at coalescence)"),
                   data = bnpr_at_coalescence) +
        geom_line(aes(x = x,y = y,
                      colour = "EPS (fit and 95% CI)",
                      size = "EPS (fit and 95% CI)",
                      linetype = "EPS (fit and 95% CI)"),
                  data = bnpr_inferred_trajectory) +
        annotate("text",
                 label = sprintf("\nGrowth/year (EPS) = %.2f",linear_estimate$coefficients[2]),
                 x = annotate_x,
                 y = y_lims[1],
                 colour = c("black"), 
                 hjust = annotate_hjust,
                 size = 3.5,
                 fontface = 'bold',
                 vjust = 0) +  
        theme_gerstung(base_size = 15) + 
        scale_y_continuous(
          breaks = y_axis_at,
          labels = y_axis_at,
          trans = 'log10') +
        scale_x_continuous(limits = c(0,max_age),
                           expand = c(0.02,0,0.02,0)) +
        ylab("Effective  pop. size (EPS)") + 
        xlab("Age") + 
        coord_cartesian(ylim = y_lims) + 
        scale_colour_manual(values = c(`EPS (fit and 95% CI)` = "black",
                                       `EPS (at coalescence)` = "black"),
                            name = NULL,
                            breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)")) + 
        scale_fill_manual(values = c(`EPS (fit and 95% CI)` = "grey80",
                                     `EPS (at coalescence)` = NA),
                          name = NULL,
                          breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)")) + 
        scale_shape_manual(values = c(`EPS (fit and 95% CI)` = NA,
                                      `EPS (at coalescence)` = 16),
                           name = NULL,
                           breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)")) + 
        scale_size_manual(values = c(`EPS (fit and 95% CI)` = 1,
                                     `EPS (at coalescence)` = 1.5),
                          name = NULL,
                          breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)")) + 
        scale_linetype_manual(values = c(`EPS (fit and 95% CI)` = 1,
                                         `EPS (at coalescence)` = 0),
                              name = NULL,
                              breaks = c("EPS (fit and 95% CI)","EPS (at coalescence)")) + 
        ggtitle(driver_id) +
        theme(legend.position = legend_position,
              legend.direction = "vertical",
              legend.justification = c(0.5,1),
              legend.text.align = 0,
              plot.title = element_text(size = 12)) 
      }
    plot_list[[curr_id]][[driver_id]] <- trajectory_plot
  }
  plot_list[[curr_id]][["tree"]] <- tree_plot
}

output_width = 13.3
output_height = 8.0

plot_3877 <- plot_grid(
  plot_list$`3877`$tree,
  plot_grid(
    plot_list$`3877`$`SF3B1-K666N`,
    plot_list$`3877`$`U2AF1-Q157R`,
    plot_list$`3877`$`Unknown driver 1`,
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

plot_3877 %>%
  ggsave(filename = sprintf("figures/%s/trees/all_trees_and_dynamics_subset_ggplot_3877.pdf",model_id),
         height = output_height,
         width = output_width,
         useDingbats = F)

plot_500 %>%
  ggsave(filename = sprintf("figures/%s/trees/all_trees_and_dynamics_subset_ggplot_500.pdf",model_id),
         height = output_height,
         width = output_width,
         useDingbats = F)

plot_2259 %>%
  ggsave(filename = sprintf("figures/%s/trees/all_trees_and_dynamics_subset_ggplot_2259.pdf",model_id),
         height = 4,
         width = output_width,
         useDingbats = F)

plot_grid(
  plot_3877,
  plot_2259,
  rel_heights = c(output_height,4),
  align = "v",
  ncol = 1
) %>% 
  ggsave(filename = sprintf("figures/%s/trees/all_trees_and_dynamics_subset_ggplot_3877_2259.pdf",model_id),
         height = output_height + 4,
         width = output_width,
         useDingbats = F)
