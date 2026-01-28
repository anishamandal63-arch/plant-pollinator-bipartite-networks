# Load packages
library(tidyverse)
library(bipartite)
library(vegan)
library(igraph)

#read data
pp.data <- read.csv(file= "pp_data.csv", na.strings = c("NA", ""))

#matrix for IDH (invasive dominated habits) only nectaring for interactions
dat.idh <- droplevels(subset(pp.data, 
                             Activity == "Nectaring" & 
                               FT == "IDH"))

plants <- sort(unique(dat.idh$Plants))
polls <- sort(unique(dat.idh$Scientificn))

mat <- matrix(NA, nrow = length(plants), ncol = length(polls))
dimnames(mat) <- list(plants, polls)

for(i in 1:length(plants)){
  for(j in 1:length(polls)){
    temp <- subset(dat.idh, Plants == plants[i] &
                     Scientificn == polls[j])
    mat[i,j] <- sum(temp$PN)}}

mat
t(mat)
x <- mat

#matrix for PF (primary forest)
dat.primary <- droplevels(subset(pp.data, 
                                 Activity == "Nectaring" & 
                                   FT == "PF" ))


plants1 <- sort(unique(dat.primary$Plants))
polls1 <- sort(unique(dat.primary$Scientificn))

mat1 <- matrix(NA, nrow = length(plants1), ncol = length(polls1))
dimnames(mat1) <- list(plants1, polls1)

for(K in 1:length(plants1)){
  for(L in 1:length(polls1)){
    temp1 <- subset(dat.primary, Plants == plants1[K] &
                      Scientificn == polls1[L])
    mat1[K,L] <- sum(temp1$PN)}}

mat1
mat[is.na(mat)] <- 0   # Replace NA with 0
mat[mat < 0] <- 0      # Replace negative values with 0
mat1[is.na(mat1)] <- 0
#############################################################interactions_plot
plotweb (mat)
plotweb (mat1)

# Set global font to italic (font = 3)
par(font = 3)

plotweb(mat,  method = "cca", labsize = 0.7, arrow="both", y.width.high=0.06, y.width.low=0.06, text.rot=90, col.high="blue", col.low="green")
plotweb(mat1,  method = "cca", labsize = 0.7, arrow="both", y.width.high=0.06, y.width.low=0.06, text.rot=90, col.high="blue", col.low="green")

###################################################______________visweb____________________###########

plot_matrix_discrete <- function(M, main = "Interactions (discrete colors)") {
  M <- as.matrix(M)
  M[is.na(M)] <- 0
  
  # Map counts to category codes
  cats <- matrix(1L, nrow = nrow(M), ncol = ncol(M))
  cats[M == 1] <- 2L
  cats[M == 2] <- 3L
  cats[M >= 3 & M <= 5] <- 4L
  cats[M > 5] <- 5L
  
  pal <- c("white", "green", "blue", "pink", "black")
  nr <- nrow(cats); nc <- ncol(cats)
  cats_flip <- cats[nr:1, , drop = FALSE]
  
  op <- par(mar = c(8, 8, 4, 2))
  on.exit(par(op), add = TRUE)
  
  image(x = 1:nc, y = 1:nr, z = t(cats_flip),
        col = pal, zlim = c(1, 5),
        xaxt = "n", yaxt = "n", xlab = "", ylab = "",
        main = main, font.main = 1)  # keep (A), (B) normal
  
  # Italic axis labels for species names
  axis(1, at = 1:nc, labels = colnames(M), las = 2, cex.axis = 0.8, font = 3)
  axis(2, at = 1:nr, labels = rev(rownames(M)), las = 2, cex.axis = 0.8, font = 3)
  
  abline(h = 0.5 + 1:nr, v = 0.5 + 1:nc, col = "grey90")
}

# Plot both side-by-side
par(mfrow = c(1, 2), mar = c(8, 5, 4, 2))

plot_matrix_discrete(mat1, main = "(A)")
plot_matrix_discrete(mat,  main = "(B)")

########################################################################
##             modularity and nestedness
###################################################################

mod <- computeModules(mat)
plotModuleWeb(mod)

mod1 <- computeModules(mat1)
plotModuleWeb(mod1)

mod@likelihood  # modularity score for IDH
mod1@likelihood   # modularity score for PF
#################################################plot_modules
# Compute modules
set.seed(123)                 # for reproducibility
mod1 <- computeModules(mat1)  # PF
mod  <- computeModules(mat)   # IDH

# Split plotting window: 1 row, 2 columns
op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 1))

# Plot PF (left)
plotModuleWeb(mod1,
              displayAlabels = TRUE,
              displayBlabels = TRUE,
              labsize = 0.7)


# Plot IDH (right)
plotModuleWeb(mod,
              displayAlabels = TRUE,
              displayBlabels = TRUE,
              labsize = 0.7)

#####################################################
# after the above code Restore par
par(op)
###############################################################
# Function to build matrix from data
###############################################################
make_matrix <- function(data) {
  plants <- sort(unique(na.omit(data$Plants)))
  polls  <- sort(unique(na.omit(data$Scientificn)))
  mat <- matrix(0, nrow = length(plants), ncol = length(polls),
                dimnames = list(plants, polls))
  for (i in seq_along(plants)) {
    for (j in seq_along(polls)) {
      temp <- subset(data, Plants == plants[i] & Scientificn == polls[j])
      mat[i, j] <- sum(temp$PN, na.rm = TRUE)
    }
  }
  mat[is.na(mat)] <- 0
  return(mat)
}
# Filter only nectaring interactions
nectar_data <- droplevels(subset(pp.data, Activity == "Nectaring"))

###############################################################
# Bootstrap NODF and Modularity (1000 runs) for PF and IDH
###############################################################
bootstrap_network <- function(data, n_iter = 1000) {
  nodf_vals <- numeric(n_iter)
  mod_vals  <- numeric(n_iter)
  
  for (i in 1:n_iter) {
    # Resample interactions with replacement
    samp <- data[sample(1:nrow(data), nrow(data), replace = TRUE), ]
    
    #Build matrix
    
    mat <- make_matrix(samp)
    
    # Calculate metrics
    mat_bin <- (mat > 0) * 1
    nodf_vals[i] <- tryCatch(networklevel(mat_bin, index = "NODF"), error = function(e) NA)
    mod_vals[i]  <- tryCatch(computeModules(mat)@likelihood, error = function(e) NA)
    
    if (i %% 100 == 0) cat("Iteration:", i, "completed\n")
  }
  
  # Return summary results
  list(
    NODF_mean   = mean(nodf_vals, na.rm = TRUE),
    NODF_CI     = quantile(nodf_vals, probs = c(0.025, 0.975), na.rm = TRUE),
    NODF_all    = nodf_vals,
    Mod_mean    = mean(mod_vals, na.rm = TRUE),
    Mod_CI      = quantile(mod_vals, probs = c(0.025, 0.975), na.rm = TRUE),
    Mod_all     = mod_vals
  )
}
results_PF  <- bootstrap_network(subset(nectar_data, FT == "PF"), n_iter = 1000)
results_IDH <- bootstrap_network(subset(nectar_data, FT == "IDH"), n_iter = 1000)

###############################################################
# Display summary results
###############################################################
cat("\n--- Primary Forest (PF) ---\n")
cat("NODF mean ± 95% CI:", results_PF$NODF_mean, "(", 
    round(results_PF$NODF_CI[1], 2), "-", round(results_PF$NODF_CI[2], 2), ")\n")
cat("Modularity mean ± 95% CI:", results_PF$Mod_mean, "(", 
    round(results_PF$Mod_CI[1], 3), "-", round(results_PF$Mod_CI[2], 3), ")\n")

cat("\n--- Invasive-Dominated Habitat (IDH) ---\n")
cat("NODF mean ± 95% CI:", results_IDH$NODF_mean, "(", 
    round(results_IDH$NODF_CI[1], 2), "-", round(results_IDH$NODF_CI[2], 2), ")\n")
cat("Modularity mean ± 95% CI:", results_IDH$Mod_mean, "(", 
    round(results_IDH$Mod_CI[1], 3), "-", round(results_IDH$Mod_CI[2], 3), ")\n")


###############################################################
# Save results as CSV

results_df <- data.frame(
  Habitat = c("PF", "IDH"),
  NODF_mean = c(results_PF$NODF_mean, results_IDH$NODF_mean),
  NODF_lower = c(results_PF$NODF_CI[1], results_IDH$NODF_CI[1]),
  NODF_upper = c(results_PF$NODF_CI[2], results_IDH$NODF_CI[2]),
  Modularity_mean = c(results_PF$Mod_mean, results_IDH$Mod_mean),
  Modularity_lower = c(results_PF$Mod_CI[1], results_IDH$Mod_CI[1]),
  Modularity_upper = c(results_PF$Mod_CI[2], results_IDH$Mod_CI[2])
)

write.csv(results_df, "Bootstrap_1000_PFvsIDH.csv", row.names = FALSE)



# ########################################################
# NMDS + PERMANOVA + BETADISPER Summary (Nectaring only)
# ###############################################################

# --- Load packages ---
library(dplyr)
library(tidyr)
library(stringr)
library(vegan)

# --- Load and clean data ---
pp.data <- read.csv("pp_data.csv", na.strings = c("NA", "")) %>%
  rename_with(~ gsub("[^A-Za-z0-9_]", "_", .x))

# --- Filter only nectaring activity ---
df <- pp.data %>%
  mutate(
    Season = case_when(
      grepl("pre", Season, ignore.case = TRUE)  ~ "Pre-monsoon",
      grepl("post", Season, ignore.case = TRUE) ~ "Post-monsoon",
      TRUE ~ Season
    ),
    FT = case_when(
      grepl("IDH", FT, ignore.case = TRUE) ~ "IDH",
      grepl("PF", FT, ignore.case = TRUE)  ~ "PF",
      TRUE ~ FT
    ),
    PN = as.numeric(PN),
    Plants = na_if(str_squish(Plants), ""),
    Pollinators = na_if(str_squish(Pollinators), "")
  ) %>%
  filter(!is.na(Activity) & grepl("^nectar", Activity, ignore.case = TRUE),
         FT %in% c("PF","IDH"),
         Season %in% c("Pre-monsoon","Post-monsoon"))

# Choose sampling unit 

df <- df %>%
  mutate(SamplingUnit = as.character(Transect))

# --- Helper functions ---
build_comm_meta <- function(data, species_col) {
  sp <- rlang::ensym(species_col)
  comm_wide <- data %>%
    filter(!is.na(!!sp)) %>%
    group_by(SamplingUnit, FT, Season, !!sp) %>%
    summarise(abund = sum(PN, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = !!sp, values_from = abund, values_fill = 0) %>%
    ungroup()
  
  meta <- comm_wide %>% select(SamplingUnit, FT, Season)
  comm <- comm_wide %>% select(-SamplingUnit, -FT, -Season)
  comm <- comm[, colSums(comm, na.rm = TRUE) > 0, drop = FALSE]
  keep_rows <- rowSums(comm, na.rm = TRUE) > 0
  meta <- meta[keep_rows, , drop = FALSE]
  comm <- comm[keep_rows, , drop = FALSE]
  list(comm = comm, meta = meta)
}

run_community_analysis <- function(comm, meta, label = "Community") {
  set.seed(123)
  dist <- vegdist(comm, method = "bray")
  nmds <- metaMDS(comm, distance = "bray", k = 2, trymax = 100, autotransform = FALSE, trace = FALSE)
  stress <- nmds$stress
  ad <- adonis2(dist ~ Season * FT, data = meta, permutations = 999)
  bd_season <- betadisper(dist, meta$Season)
  bd_ft <- betadisper(dist, meta$FT)
  list(stress = stress,
       adonis = ad,
       betadisper_season = anova(bd_season),
       betadisper_ft = anova(bd_ft))
}

# --- Run for PLANTS ---
plants_obj <- build_comm_meta(df, "Plants")
plants_res <- run_community_analysis(plants_obj$comm, plants_obj$meta, "Plants")

# --- Run for POLLINATORS ---
polls_obj <- build_comm_meta(df, "Pollinators")
polls_res <- run_community_analysis(polls_obj$comm, polls_obj$meta, "Pollinators")

# --- Extract key values into compact summary table ---
extract_summary <- function(label, res) {
  perm <- as.data.frame(res$adonis)
  data.frame(
    Group = label,
    Stress = round(res$stress, 3),
    Season_R2 = round(perm["Season", "R2"], 3),
    Season_F = round(perm["Season", "F"], 2),
    Season_p = perm["Season", "Pr(>F)"],
    FT_R2 = round(perm["FT", "R2"], 3),
    FT_F = round(perm["FT", "F"], 2),
    FT_p = perm["FT", "Pr(>F)"],
    Interaction_R2 = round(perm["Season:FT", "R2"], 3),
    Interaction_F = round(perm["Season:FT", "F"], 2),
    Interaction_p = perm["Season:FT", "Pr(>F)"],
    Disp_Season_F = round(res$betadisper_season$`F value`[1], 2),
    Disp_Season_p = res$betadisper_season$`Pr(>F)`[1],
    Disp_FT_F = round(res$betadisper_ft$`F value`[1], 2),
    Disp_FT_p = res$betadisper_ft$`Pr(>F)`[1]
  )
}

summary_table <- bind_rows(
  extract_summary("Plants (nectaring)", plants_res),
  extract_summary("Pollinators (nectaring)", polls_res)
)

cat("\n\n============== SUMMARY TABLE ==============\n")
print(summary_table)
cat("===========================================\n")



# Print full PERMANOVA table with df
plants_res$adonis
polls_res$adonis

##############################################################
#############_____________nmds____plot
#####################################################
library(vegan)
library(ggplot2)
library(dplyr)

set.seed(123)
nmds_plants <- metaMDS(
  plants_obj$comm,
  distance = "bray",
  k = 2,
  trymax = 100,
  autotransform = FALSE,
  trace = FALSE
)
set.seed(123)
nmds_polls <- metaMDS(
  polls_obj$comm,
  distance = "bray",
  k = 2,
  trymax = 100,
  autotransform = FALSE,
  trace = FALSE
)
nmds_plants$stress
nmds_polls$stress

nmds_df <- function(nmds, meta){
  scores_df <- as.data.frame(scores(nmds, display = "sites"))
  df <- cbind(scores_df, meta)
  df$Season <- factor(df$Season, levels = c("Pre-monsoon", "Post-monsoon"))
  df$FT     <- factor(df$FT, levels = c("PF", "IDH"))
  df
}
plant_meta <- plants_obj$meta
poll_meta  <- polls_obj$meta

plants_df <- nmds_df(nmds_plants, plant_meta)
polls_df  <- nmds_df(nmds_polls, poll_meta)
# Groups with >=3 samples (for ellipse)
groups_ok <- function(df) {
  df %>% count(FT, Season) %>% filter(n >= 3) %>% mutate(g = interaction(FT, Season)) %>% pull(g)
}
plants_ok <- groups_ok(plants_df)
polls_ok  <- groups_ok(polls_df)

# shapes
shape_vals <- c(PF = 16, IDH = 17)

# helper: safe compute corner coords 
corner_ll_safe <- function(df, xcol = "NMDS1", ycol = "NMDS2", xpad = 0.05, ypad = 0.05) {
  if (is.null(df) || !all(c(xcol, ycol) %in% names(df))) return(NULL)
  if (nrow(df) < 1) return(NULL)
  xr <- range(df[[xcol]], na.rm = TRUE)
  yr <- range(df[[ycol]], na.rm = TRUE)
  if (any(is.infinite(xr)) || any(is.infinite(yr))) return(NULL)
  x <- xr[1] + (xr[2] - xr[1]) * xpad
  y <- yr[1] + (yr[2] - yr[1]) * ypad
  list(x = x, y = y)
}

# safe annotate 
safe_annotate_text <- function(p, coords, label_text, size = 4, ...) {
  if (is.null(coords) || is.null(label_text) || !nzchar(as.character(label_text))) return(p)
  # coords must be numeric
  if (!is.finite(coords$x) || !is.finite(coords$y)) return(p)
  p + annotate("text", x = coords$x, y = coords$y, label = label_text, hjust = 0, vjust = 0, size = size, ...)
}

# Build plants plot 
if (!is.null(plants_df) && nrow(plants_df) > 0) {
  llp <- corner_ll_safe(plants_df)
  stress_label_plants <- if (!is.null(nmds_plants) && !is.null(nmds_plants$stress)) paste0("stress = ", sprintf("%.3f", nmds_plants$stress)) else NA_character_
  p_plants <- ggplot(plants_df, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(color = Season, shape = FT), size = 3.5, alpha = 0.9) +
    stat_ellipse(
      data = plants_df %>% filter(interaction(FT, Season) %in% plants_ok),
      aes(group = interaction(FT, Season), color = Season),
      type = "norm", level = 0.68, linewidth = 1
    ) +
    scale_shape_manual(values = shape_vals) +
    coord_equal() +
    labs(x = "NMDS1", y = "NMDS2", color = "Season", shape = "Forest type") +
    theme_bw(base_size = 14) +
    theme(panel.grid = element_blank(), plot.title = element_blank(),
          legend.position = "bottom", legend.title = element_text(size = 12),
          legend.text = element_text(size = 11), axis.title = element_text(size = 13),
          axis.text = element_text(size = 11))
  p_plants <- safe_annotate_text(p_plants, llp, stress_label_plants, size = 4)
} else {
  p_plants <- ggplot() + theme_void() + ggtitle("Insufficient plant data for NMDS")
}

# Build pollinators plot 
if (!is.null(polls_df) && nrow(polls_df) > 0) {
  llp2 <- corner_ll_safe(polls_df)
  stress_label_polls <- if (!is.null(nmds_polls) && !is.null(nmds_polls$stress)) paste0("stress = ", sprintf("%.3f", nmds_polls$stress)) else NA_character_
  p_polls <- ggplot(polls_df, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(color = Season, shape = FT), size = 3.5, alpha = 0.9) +
    stat_ellipse(
      data = polls_df %>% filter(interaction(FT, Season) %in% polls_ok),
      aes(group = interaction(FT, Season), color = Season),
      type = "norm", level = 0.68, linewidth = 1
    ) +
    scale_shape_manual(values = shape_vals) +
    coord_equal() +
    labs(x = "NMDS1", y = "NMDS2", color = "Season", shape = "Forest type") +
    theme_bw(base_size = 14) +
    theme(panel.grid = element_blank(), plot.title = element_blank(),
          legend.position = "bottom", legend.title = element_text(size = 12),
          legend.text = element_text(size = 11), axis.title = element_text(size = 13),
          axis.text = element_text(size = 11))
  p_polls <- safe_annotate_text(p_polls, llp2, stress_label_polls, size = 4)
} else {
  p_polls <- ggplot() + theme_void() + ggtitle("Insufficient pollinator data for NMDS")
}

# Combine side-by-side with shared legend 
combine_and_print_safe <- function(p1, p2) {
  if ("patchwork" %in% rownames(installed.packages())) {
    library(patchwork)
    combined <- (p1 + p2) + plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "bottom")
    print(combined)
    return(invisible(TRUE))
  }
  if ("cowplot" %in% rownames(installed.packages())) {
    library(cowplot)
    # attempt to get legend from p1, if p1 has none, try p2; catch errors
    legend <- tryCatch(get_legend(p1 + theme(legend.position = "bottom")), error = function(e) {
      tryCatch(get_legend(p2 + theme(legend.position = "bottom")), error = function(e2) NULL)
    })
    p1_noleg <- p1 + theme(legend.position = "none")
    p2_noleg <- p2 + theme(legend.position = "none")
    if (!is.null(legend)) {
      combined <- plot_grid(p1_noleg, p2_noleg, ncol = 2, align = "hv")
      final <- plot_grid(combined, legend, ncol = 1, rel_heights = c(1, 0.12))
      print(final)
      return(invisible(TRUE))
    } else {
      print(plot_grid(p1_noleg, p2_noleg, ncol = 2, align = "hv"))
      return(invisible(TRUE))
    }
  }
  if ("gridExtra" %in% rownames(installed.packages())) {
    library(gridExtra)
    grid.arrange(p1 + theme(legend.position = "bottom"), p2 + theme(legend.position = "none"), ncol = 2)
    return(invisible(TRUE))
  }
  par(mfrow = c(1,2)); print(p1); print(p2); par(mfrow = c(1,1))
}

# call combine
combine_and_print_safe(p_plants, p_polls)





#####################################################################
# Shannon diversity (H') and Evenness (Pielou J) for Plants & Pollinators
#########################################################################3
# clean names
pp.data <- pp.data %>% rename_with(~ gsub("[^A-Za-z0-9_]", "_", .x))

# prepare and FILTER to nectaring ONLY
df <- pp.data %>%
  mutate(
    Season = case_when(
      grepl("pre", Season, ignore.case = TRUE)  ~ "Pre-monsoon",
      grepl("post", Season, ignore.case = TRUE) ~ "Post-monsoon",
      TRUE ~ Season
    ),
    FT = case_when(
      grepl("IDH", FT, ignore.case = TRUE) ~ "IDH",
      grepl("PF",  FT, ignore.case = TRUE) ~ "PF",
      TRUE ~ FT
    ),
    PN = as.numeric(PN),
    Plants = na_if(str_squish(Plants), ""),
    Pollinators = na_if(str_squish(Pollinators), "")
  ) %>%
  # <<< IMPORTANT: keep only nectaring interactions >>>
  filter(!is.na(Activity), grepl("^nectar", Activity, ignore.case = TRUE)) %>%
  filter(FT %in% c("PF", "IDH"))

compute_H_J <- function(abund_vector) {
  abund_vector <- as.numeric(abund_vector)
  abund_vector <- abund_vector[abund_vector > 0]
  S <- length(abund_vector)
  if (S == 0) {
    return(tibble(richness = 0, H = NA_real_, J = NA_real_))
  }
  N <- sum(abund_vector)
  p <- abund_vector / N
  H <- -sum(p * log(p))        # natural log
  J <- if (S <= 1) NA_real_ else H / log(S)
  tibble(richness = S, H = H, J = J)
}

# POOLED per FT (Plants) 
plants_pooled <- df %>%
  filter(!is.na(Plants)) %>%
  group_by(FT, Plants) %>%
  summarise(abund = sum(PN, na.rm = TRUE), .groups = "drop") %>%
  group_by(FT) %>%
  summarise(metrics = list(compute_H_J(abund)), .groups = "drop") %>%
  unnest(cols = c(metrics))

# POOLED per FT (Pollinators) 
polls_pooled <- df %>%
  filter(!is.na(Pollinators)) %>%
  group_by(FT, Pollinators) %>%
  summarise(abund = sum(PN, na.rm = TRUE), .groups = "drop") %>%
  group_by(FT) %>%
  summarise(metrics = list(compute_H_J(abund)), .groups = "drop") %>%
  unnest(cols = c(metrics))

cat("\n--- Plants (pooled per FT, nectaring only) ---\n"); print(plants_pooled)
cat("\n--- Pollinators (pooled per FT, nectaring only) ---\n"); print(polls_pooled)

#POOLED per FT x SEASON (nectaring only)
plants_ft_season <- df %>%
  filter(!is.na(Plants)) %>%
  group_by(FT, Season, Plants) %>%
  summarise(abund = sum(PN, na.rm = TRUE), .groups = "drop") %>%
  group_by(FT, Season) %>%
  summarise(metrics = list(compute_H_J(abund)), .groups = "drop") %>%
  unnest(cols = c(metrics))

polls_ft_season <- df %>%
  filter(!is.na(Pollinators)) %>%
  group_by(FT, Season, Pollinators) %>%
  summarise(abund = sum(PN, na.rm = TRUE), .groups = "drop") %>%
  group_by(FT, Season) %>%
  summarise(metrics = list(compute_H_J(abund)), .groups = "drop") %>%
  unnest(cols = c(metrics))

cat("\n--- Plants (pooled per FT x Season, nectaring only) ---\n"); print(plants_ft_season)
cat("\n--- Pollinators (pooled per FT x Season, nectaring only) ---\n"); print(polls_ft_season)

###################bootstrap graphs###########in supplementary_material

library(ggplot2)

df <- results_df
df$Modularity_obs <- df$Modularity_mean
df$NODF_obs <- df$NODF_mean

required_cols <- c("Habitat","NODF_obs","NODF_lower","NODF_upper","Modularity_obs","Modularity_lower","Modularity_upper")
missing_cols <- setdiff(required_cols, colnames(df))
if(length(missing_cols)>0) stop("Missing columns in result data frame: ", paste(missing_cols, collapse=", "))



####Observed modularity from original networks (triangles) ----
observed_Q_PF <- if (exists("mod1") && !is.null(mod1)) mod1@likelihood else NA
observed_Q_IDH <- if (exists("mod") && !is.null(mod)) mod@likelihood else NA

# Create plotting data
df_mod <- data.frame(
  Habitat = as.character(df$Habitat),
  Mean = as.numeric(df$Modularity_obs),
  Lower = as.numeric(df$Modularity_lower),
  Upper = as.numeric(df$Modularity_upper),
  stringsAsFactors = FALSE
)

df_nodf <- data.frame(
  Habitat = as.character(df$Habitat),
  Mean = as.numeric(df$NODF_obs),
  Lower = as.numeric(df$NODF_lower),
  Upper = as.numeric(df$NODF_upper),
  stringsAsFactors = FALSE
)
# Modularity plot
p_mod <- ggplot(df_mod, aes(x = Habitat, y = Mean, color = Habitat)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.12, size = 0.9) +
  
  geom_point(data = df_mod, aes(x = Habitat, y = ifelse(Habitat == "PF", observed_Q_PF, observed_Q_IDH)),
             shape = 17, size = 3, color = "black", inherit.aes = FALSE) +
  labs(title = "Modularity (Q) with 95% CI", y = "Modularity (Q)", x = "") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

print(p_mod)  

# Save to files
ggsave("Modularity_Q_CI.png", p_mod, width = 6, height = 4, dpi = 300)
ggsave("Modularity_Q_CI.pdf", p_mod, width = 6, height = 4)

# ---- Nestedness (NODF) plot ----
p_nodf <- ggplot(df_nodf, aes(x = Habitat, y = Mean, color = Habitat)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.12, size = 0.9) +
  labs(title = "Nestedness (NODF) with 95% CI", y = "NODF", x = "") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

print(p_nodf)

# Save to files
ggsave("Nestedness_NODF_CI.png", p_nodf, width = 6, height = 4, dpi = 300)
ggsave("Nestedness_NODF_CI.pdf", p_nodf, width = 6, height = 4)



