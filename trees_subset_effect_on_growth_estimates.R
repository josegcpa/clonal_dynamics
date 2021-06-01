library(tidyverse)
library(ape)
library(dendextend)
library(cowplot)
library(phylodyn)
library(ggsci)
library(phytools)
INLA:::inla.dynload.workaround()

run_clonex <- function(d=1000, p=10000, o=0, mu=1e-9, nu=1e-5, s=0.01,
                       t=0.01, n=1, N=200000, a=2, g=1000, G=100, r=42, 
                       X = 0,Y = 0,Z = 0.0001,L = 0, R = 1,
                       path_to_clonex=".",folder="tmp"){
  dir.create(folder,showWarnings = F)
  x <- sprintf("%s/clonex -d %i -p %i -o %i -u %e -v %e -s %e -t %e -n %i -N %i -g %i -G %i -r %i -f %s -X %i -Y %i -Z %f -L %i -R %i",
               path_to_clonex,d,p,o,mu,nu,s,t,n,N,g,G,r,folder,X,Y,Z,L,R)
  print(x)
  system(x)
}

read_clonex <- function(path,d = 1000,pop_size_calc = T) {
  x <- read_tsv(path,
                col_names = c("Gen","NClones","CloneID","MutID"),
                col_types = c(col_integer(),col_integer(),col_integer(),col_integer()),
                progress = T) %>%
    mutate(Driver = MutID <= d)
  if (pop_size_calc == T) {
    pop_size <- x %>%
      group_by(Gen,CloneID) %>%
      summarise(CloneSize = max(NClones),.groups = "drop") %>% 
      group_by(Gen) %>%
      summarise(PopulationSize = sum(CloneSize),.groups = "drop")
    x <- merge(x,pop_size,by = "Gen")
  }
  return(x)
}

get_allele_frequencies <- function(clonex_data) {
  clonex_data %>%
    group_by(Gen,MutID,Driver) %>%
    summarise(CloneSize = sum(NClones),
              PopulationSize = PopulationSize[1],
              .groups = "drop") %>% 
    mutate(AF = CloneSize / PopulationSize) %>% 
    return
}

theme_elegant <- function(...) {
  theme_minimal(...) + 
    theme(
      axis.ticks = element_line(),
      axis.line = element_line(),
      panel.grid = element_blank()
    )
}

build_tree <- function(sub_data,subsample_size=100,
                       detection_threshold=0) {
  MaxClone <- max(sub_data$CloneID)
  MaxMut <- max(sub_data$MutID)
  sub_af <- sub_data %>% 
    group_by(CloneID) %>% 
    summarise(N = max(NClones),.groups = "drop") %>%
    arrange(CloneID)
  sub_af_sparse <- rep(0,MaxClone)
  sub_af_sparse[sub_af$CloneID] <- sub_af$N 
  
  Presence <- matrix(0,nrow = MaxClone,ncol = MaxMut)
  for (i in 1:nrow(sub_data)) {
    Presence[sub_data$CloneID[i],sub_data$MutID[i]] <- 1
  }
  Presence[is.na(Presence)] <- 0
  
  if (sum(as.logical(sub_af_sparse)) < subsample_size) {
    sub_af_sparse <- rep(1,length(sub_af_sparse))
  }
  
  S <- sample(seq(1,MaxClone),subsample_size,replace = F,prob = sub_af_sparse)
  Presence <- Presence[S,]
  Presence <- Presence[,colSums(Presence) > 0]
  Presence <- rbind(Presence,wt=0)
  
  dst <- as.matrix(dist(Presence > 0,method = "manhattan"))
  tree <- root(njs(dst),
               outgroup = "wt",
               resolve.root = T,
               edgelabel = T) 
  tree <- drop.tip(tree,"wt")
  tree$tip.label <- S
  return(list(tree = tree,S = S))
}

get_root <- function(tree){
  e <- tree$edge
  return(setdiff(e[,1], e[,2]))
}
cut_tree <- function(tree, depth){
  t <- length(tree$tip.label)
  d <- node.depth.edgelength(tree)
  e <- tree$edge
  a <- d
  a[e[,2]] <- d[e[,1]]
  w <- which(d>depth & a<=depth)
  c <- w[w>t]
  p <- prop.part(tree) 
  clades <- c(as.list(w[w<=t]), p[c-t])
  return(clades)
}
calibrate_tree <- function(tree){
  t <- length(tree$tip.label)
  d <- node.depth.edgelength(tree)
  m <- mean(d[1:t])
  tree$edge.length <- tree$edge.length/m
  return(tree)
}

prune_tips <- function(tree, depth){
  w <- tree$edge[,2] %in% 1:length(tree$tip.label)
  tree$edge.length[w] <- tree$edge.length[w] - depth
  return(tree)
}

run_clonex(path_to_clonex = "$HOME/clonex",
           folder = sprintf("clonex_runs/%s",0.01),
           n = 2e5,N = 2e5,p = 10e3,d = 1e3,g = 1000,t = 0,
           mu = 1e-9,nu=1e-5,G = 50,s = 0.005,R = 3)

FILE <- "clonex_runs/0.005/r0005.csv"
DATA <- read_clonex(FILE,pop_size_calc = F)
DATA$PopulationSize <- 200e3
af <- get_allele_frequencies(DATA)
af %>% 
  subset(MutID != 0) %>%
  subset(Driver == T) %>% 
  ggplot(aes(x = Gen,y = AF,colour = as.factor(MutID))) + 
  geom_point() + 
  geom_line() + 
  scale_colour_aaas(guide = F) + 
  theme_elegant(base_size = 15)

sub_data <- DATA %>%
  subset(Gen == 1000)

clones_w_driver <- sub_data %>%
  subset(Driver == T) %>% 
  select(CloneID,MutID)  %>%
  distinct

table(clones_w_driver$MutID)

ESTIMATES <- list()
for (i in seq(1,10)) {
  print(i)
  TR <- build_tree(sub_data) 
  clades <- cut_tree(TR$tree, 0.1)
  ESTIMATES[[i]] <- list()
  CC <- 0
  for (c in clades) {
    CC <- CC + 1
    if (length(c) >= 5) {
      tree_subset <- keep.tip(TR$tree,c) 
      DRIVER <- clones_w_driver$MutID[clones_w_driver$CloneID %in% tree_subset$tip.label] %>%
        unique 
      DRIVER <- paste(DRIVER,collapse = ',')
      if (length(DRIVER) > 0) {
        tree_subset <- multi2di(tree_subset)
        bnpr_estimate <- BNPR(tree_subset)
        ESTIMATES[[i]][[as.character(CC)]] <- list(
          bnpr_estimate = bnpr_estimate,
          tree = TR$tree,
          driver = DRIVER,
          ID = paste(i,CC,sep = ','),
          clade = c)
      }
    }
  }
}

pdf("tmp.pdf",8,5,pointsize = 10)

M_af <- af %>% 
  subset(MutID != 0) %>%
  subset(Driver == T) %>% 
  group_by(MutID) %>% 
  filter(length(unique(Gen)) > 3) %>% 
  summarise(MA = max(AF),
            MG = max(Gen))

af %>% 
  subset(MutID != 0) %>%
  subset(Driver == T) %>% 
  ggplot(aes(x = Gen,y = AF,colour = as.factor(MutID))) + 
  geom_point() + 
  geom_line() + 
  scale_colour_aaas(guide = F) + 
  theme_elegant(base_size = 15) +
  geom_text(
    inherit.aes = F,
    data = M_af,
    aes(x = MG,MA,label = MutID),
    hjust = 1.5,
    vjust = -0.5,
    size = 4
  )

DR <- c()
EST <- c()
for (I in 1:length(ESTIMATES)) {
  par(mfrow=c(length(ESTIMATES[[I]]),3),mar = c(4,4,1,1))
  for (element in ESTIMATES[[I]]) {
    colour <- ifelse(c(1:100) %in% element$clade,
                     "red",
                     "grey40")
    plot(element$tree,tip.color = colour)
    title(sprintf("Driver = %s",element$driver))
    plot_BNPR(element$bnpr_estimate)
    tmp <- data.frame(x = -element$bnpr_estimate$summary$time,
                      y = element$bnpr_estimate$summary$mean,
                      ID = element$ID)
    w <- (log(element$bnpr_estimate$summary$quant0.975) - log(element$bnpr_estimate$summary$quant0.025))^2
    E <- lm(log(y) ~ x,data = tmp,weights = w)
    x <- tmp$x
    y <- tmp$y
    m <- E$coefficients[2]
    b <- E$coefficients[1]
    EST <- c(EST,m)
    DR <- c(DR,element$driver)
    L <- exp(x * m + b)
    plot(x=x,y=y,log = "y")
    lines(x = x,y = L)

  }
}

data.frame(
  Driver = DR,
  Estimates = EST
) %>% 
  ggplot(aes(x = as.character(Driver),y = Estimates)) + 
  geom_jitter(height = 0,width = 0.2) + 
  geom_boxplot(fill = NA,size = 0.3) + 
  theme_elegant(base_size = 15) + 
  xlab("Driver") + 
  ylab("Coefficient estimates") + 
  stat_summary(geom = "point",
               fun = mean,
               alpha = 0.5,
               size = 5)

dev.off()