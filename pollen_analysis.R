# Pollen experiment analyses
# Created by D. Denney Nov. 2025
# Last updated - Sept. 2026

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(car)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(terra)
library(tigris)
library(ggrepel)
library(patchwork)
library(elevatr)
library(sf)
library(ggspatial)
library(maps)
library(cowplot)
library(osmdata)
library(gt)
library(writexl)

options(tigris_use_cache = TRUE)

#############################
## Read in datasets
##########################

p <- read.table("summary_by_flower.txt",header = T)
p$Treatment <- as.factor(p$Treatment)
p$Pop <- as.factor(p$Pop)
p$Elevation_Cat <- as.factor(p$Elevation_Cat)
p$s_elev <- scale(p$Elevation,center=TRUE, scale=TRUE)
p$Treat_Block <- paste0(p$Treatment,"_",p$Block)
p$Treat_Block <- as.factor(p$Treat_Block)

str(p)
summary(p)

pollen1 <- read.csv("Pollen1_DD.csv", header = T)
pollen1$Treatment <- as.factor(pollen1$Treatment)
pollen1$Population <- as.factor(pollen1$Population)
pollen1$Replicate <- as.numeric(pollen1$Replicate)
# Updated tilingii names based on taxonomic revisions. See Gabbie's paper (2021) and Nesom 2014 & 2012
pollen1$Updated_species <- as.factor(pollen1$Updated_species)
pollen1$s_elev <- scale(pollen1$Elevation,center=TRUE, scale=TRUE)
pollen1$PlantID <- paste0(pollen1$Population, "_",pollen1$Replicate)

pollen2 <- read.csv("Pollen2_DD.csv", header = T)
pollen2$Treatment <- as.factor(pollen2$Treatment)
pollen2$Population <- as.factor(pollen2$Population)
pollen2$Replicate <- as.numeric(pollen2$Replicate)
pollen2$PlantID <- with(pollen2, paste(Population, Replicate, sep = "_"))
pollen2$Group <- as.factor(pollen2$Group)
pollen2$s_elev <- scale(pollen2$Elevation,center=TRUE, scale=TRUE)
# Updated taxonomy of tilingii. See above note
pollen2$Species <- as.factor(pollen2$Species)

# Standardize treatment labels 
pollen1 <- pollen1 %>%
  mutate(Treatment_short = as.character(Treatment),
    Treatment_label = recode(as.character(Treatment), "0" = "Control", "1" = "Heat", "C" = "Control", "H" = "Heat",.default = as.character(Treatment)),
    Treatment_label = factor(Treatment_label,levels = c("Control", "Heat")))

pollen2 <- pollen2 %>%
  mutate(Treatment_short = as.character(Treatment),
    Treatment_label = recode(as.character(Treatment), "0" = "Control", "1" = "Heat", "C" = "Control", "H" = "Heat",.default = as.character(Treatment)),
    Treatment_label = factor(
      Treatment_label,
      levels = c("Control", "Heat") ) )

# Pollen 1 grouping
pollen1 <- pollen1 %>%
  mutate(elev_group = factor(if_else(standardize_species(Updated_species) == "guttatus",  "widespread",  "restricted" ),
      levels = c("widespread", "restricted") ) )

# Pollen 2 grouping
pollen2 <- pollen2 %>%  mutate(type = factor(if_else(standardize_species(Species) == "guttatus", "guttatus",  "tilingii"),
      levels = c("guttatus", "tilingii")))


####################################
### Sierra Nevada experiment 
#################################
## =========================
## Map of populations from Sierra Nevada expt
## =========================
pollen_gps <- p %>%
  group_by(Pop) %>%
  summarise(
    Pop = dplyr::first(as.character(Pop)),           
    Latitude = mean(Latitude, na.rm = TRUE),
    Longitude = mean(Longitude, na.rm = TRUE),
    Elevation_Cat = dplyr::first(Elevation_Cat),
    .groups = "drop"
  )

## Convert to sf points (WGS84 lon/lat)
pts_sf <- st_as_sf(pollen_gps, coords = c("Longitude", "Latitude"), crs = 4326)

## Make a plain data.frame with X/Y columns for ggplot geoms + ggrepel
pts_xy <- cbind(st_drop_geometry(pts_sf), st_coordinates(pts_sf))

## US, CA, NV polygons
states_sf <- tigris::states(cb = TRUE) %>% st_transform(4326)
ca <- states_sf %>% filter(STUSPS == "CA")
nv <- states_sf %>% filter(STUSPS == "NV")


## Define map extent from points
bb <- st_bbox(pts_sf)
pad <- 0.3

x_lim <- c(as.numeric(bb[["xmin"]]) - pad, as.numeric(bb[["xmax"]]) + pad)
y_lim <- c(as.numeric(bb[["ymin"]]) - pad, as.numeric(bb[["ymax"]]) + pad)

## Polygon for the inset map
win_poly <- st_as_sfc(
  st_bbox(
    c(
      xmin = x_lim[1], xmax = x_lim[2],
      ymin = y_lim[1], ymax = y_lim[2]
    ),
    crs = 4326
  )
)

## Elevation raster
dem_raster <- elevatr::get_elev_raster(locations = ca, z = 7, clip = "locations")

## Convert raster to df for ggplot
dem_df <- as.data.frame(dem_raster, xy = TRUE, na.rm = TRUE)
names(dem_df) <- c("x", "y", "elev_m")


## Generate map of inset, main CA locations
us <- states_sf

p_inset <- ggplot() +
  geom_sf(data = us, fill = "grey95", color = "grey40", linewidth = 0.25) +
  geom_sf(data = ca, fill = "grey85", color = "grey30", linewidth = 0.35) +
  geom_sf(data = win_poly, fill = NA, color = "red3", linewidth = 0.9) +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = "grey10", linewidth = 0.4)
  )

# California plot with scale bar!
p_main <- ggplot() +
  geom_sf(data = nv, fill = "grey92", color = NA) +
  geom_raster(data = dem_df, aes(x = x, y = y, fill = elev_m)) +
  scale_fill_viridis_c(name = "Elevation (m)", option = "D") +
  geom_sf(data = ca, fill = NA, color = "grey20", linewidth = 0.4) +
  geom_point(data = pts_xy, aes(x = X, y = Y), size = 3) +
  ggrepel::geom_text_repel(data = pts_xy, aes(x = X, y = Y, label = Pop), size = 3) +
  ggspatial::annotation_scale(
    location = "bl",
    width_hint = 0.28,
    text_cex = 0.7,
    line_width = 0.6,
    unit_category = "metric"
  ) +
  coord_sf(xlim = x_lim, ylim = y_lim, expand = FALSE) +
  theme_void()

sierra_plot  <- p_main +
  patchwork::inset_element(
    p_inset,
    left = 0.67, bottom = 0.75, right = 0.98, top = 0.995
  )

ggsave(filename = "map_sierra.svg", plot = sierra_plot, width = 8, height = 5.5, units = "in")
##=======================
## Filter pollen data for models
##======================

# Only consider plants that produced flowers
dat <- p %>% filter(Flowered == 1)

treatment_labels <- c("C" = "Control", "H" = "Heat")

dat <- dat %>%
  mutate(Population = reorder(Pop, Elevation, FUN = mean),
    Elev_num = as.numeric(as.character(Elevation)),   
    Elev_round = round(Elev_num),
    Elev_label = paste0(Elev_round),
    Population_ord = reorder(Elev_label, Elev_num, FUN = mean))

high <- dat %>% filter(Elevation_Cat == "High")
low <- dat %>% filter(Elevation_Cat == "Low")

### ======================
## Control vs Heat by population boxplots
### ====================

m_pop_fix <- glmmTMB(
  cbind(Viable.Pollen, Inviable.Pollen) ~ Treatment * Population + (1 | Treat_Block),
  data = dat,
  family = betabinomial())

emm_pop <- emmeans(m_pop_fix, ~ Treatment | Population, type = "response")

con_pop <- contrast(
  emm_pop,
  method = "revpairwise",   # gives H - C if Treatment levels are C, H
  by = "Population")

con_pop_sum <- summary(con_pop, infer = c(TRUE, TRUE), adjust = "holm")
con_pop_sum

# convert p values to * for annotations on box plots
p_to_stars <- function(p) {
  case_when(
    is.na(p)      ~ "",
    p < 0.001     ~ "***",
    p < 0.01      ~ "**",
    p < 0.05      ~ "*",
    TRUE          ~ "ns"
  )
}

ann <- as.data.frame(con_pop_sum) %>%
  mutate(sig = p_to_stars(p.value))

# Labels should appear at y = 1
ypos <- dat %>%
  group_by(Population) %>%
  summarise(
    y = pmin(1, max(Prop.Viable, na.rm = TRUE) + 0.06),
    .groups = "drop"
  )

ann <- ann %>% left_join(ypos, by = "Population")

# Boxplots by population
p_viab <- ggplot(dat, aes(x = Population, y = Prop.Viable, color = Treatment, fill = Treatment)) +
   geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
  labs(
    title = "Viable Pollen by Sierran Population",
    x = "Population",
    y = "Proportion Viable"
  ) + scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("C" = "#6c95ad", "H" = "#8a404e"))+
  scale_fill_manual(values = c("C" = "#6c95ad", "H" = "#8a404e"))

p_viab


# Zoom in to high vs low categories
ggplot(high, aes(x = Population_ord, y = Prop.Viable, color = Treatment, fill = Treatment)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
  labs(title = "High elevation populations",
    x = "Elevation (m)",
    y = "Proportion Viable Pollen",
    color = "Treatment", fill = "Treatment"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(
    values = c("C" = "#6c95ad", "H" = "#8a404e"),
    labels = c("C" = "Control", "H" = "Heat")
  ) +
  scale_fill_manual(
    values = c("C" = "#6c95ad", "H" = "#8a404e"),
    labels = c("C" = "Control", "H" = "Heat")
  ) +
  theme_bw()

ggplot(low, aes(x = Population_ord, y = Prop.Viable, color = Treatment, fill = Treatment)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
  labs( title = "Low elevation populations",
    x = "Elevation (m)",
    y = "Proportion Viable Pollen",
    color = "Treatment", fill = "Treatment"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(
    values = c("C" = "#6c95ad", "H" = "#8a404e"),
    labels = c("C" = "Control", "H" = "Heat")
  ) +
  scale_fill_manual(
    values = c("C" = "#6c95ad", "H" = "#8a404e"),
    labels = c("C" = "Control", "H" = "Heat")
  ) +
  theme_bw()

# Total pollen production
# Including significance values for the plot
pd <- position_dodge(width = 0.6)
m_total <- glmmTMB(
  Total.Pollen ~ Population * Treatment + (1 | Treat_Block),
  data = dat,
  family = nbinom2()
)

emm_total <- emmeans(m_total, ~ Treatment | Population, type = "link")
con_total <- contrast(emm_total, method = "revpairwise", by = "Population")  # H - C
con_total_sum <- summary(con_total, infer = c(TRUE, TRUE), adjust = "holm")

# Use the p to stars function from above
ann_total <- as.data.frame(con_total_sum) %>%
  mutate(sig = p_to_stars(p.value))

# Variable y positions
ypos_total <- dat %>%
  group_by(Population) %>%
  summarise(
    y = max(Total.Pollen, na.rm = TRUE) * 1.08,
    .groups = "drop"
  )

ann_total <- ann_total %>%
  left_join(ypos_total, by = "Population")

p_total <- ggplot(dat, aes(x = Population, y = Total.Pollen, color = Treatment, fill = Treatment)) +
  #geom_violin(trim = FALSE, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
  #facet_wrap(~ Treatment, labeller = labeller(Treatment = treatment_labels)) +
  labs(
    title = "Total Pollen Counts by Sierran Population",
    x = "Population",
    y = "Total Pollen Counted"
  ) +
  theme_bw() + scale_color_manual(
    values = c("C" = "#6c95ad", "H" = "#8a404e"),
    labels = c("C" = "Control", "H" = "Heat")
  ) +
  scale_fill_manual(
    values = c("C" = "#6c95ad", "H" = "#8a404e"),
    labels = c("C" = "Control", "H" = "Heat")
  ) 

p_total



## ====================
## Pollen viability models
## ====================

m_beta <- glmmTMB(
  cbind(Viable.Pollen, Inviable.Pollen) ~ s_elev * Treatment + (1 | Treat_Block) + (1| Population),
  data = dat,
  family = betabinomial()
)
summary(m_beta)
Anova(m_beta, type = "III")
simulationOutput <- simulateResiduals(fittedModel = m_beta, plot = T, re.form = NULL)
emmip(m_beta, ~ Treatment, type = "response", CIs = T) +theme_bw()


# Significant elev x treatment

# Create grid of raw Elevation values
elev_seq <- seq(min(dat$Elevation), max(dat$Elevation), length.out = 100)

# Scale them using the same center/scale as in dat$s_elev
center_val <- attr(scale(dat$Elevation), "scaled:center")
scale_val  <- attr(scale(dat$Elevation), "scaled:scale")

# Predicted data for plotting
pred_grid <- expand.grid(
  Elevation = elev_seq,
  Treatment = levels(dat$Treatment),
  Block = NA,        # placeholder because it needs to be in the data for some reason
  Population = NA    # placeholder 
) %>%
  mutate(s_elev = (Elevation - center_val) / scale_val)

preds <- cbind(
  pred_grid,
  predict(m_beta, newdata = pred_grid, type = "response", se.fit = TRUE, re.form = NA)
) %>%
  rename(fit = fit, se = se.fit) %>%
  mutate(
    lower = fit - 1.96 * se,
    upper = fit + 1.96 * se
  )

#Predicted data trends
ggplot(preds, aes(x = Elevation, y = fit, color = Treatment, fill = Treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("C" = "#6c95ad", "H" = "#8a404e"))+
  scale_fill_manual(values = c("C" = "#6c95ad", "H" = "#8a404e"))+
  labs(
    x = "Elevation (m)",
    y = "Predicted pollen viability",
    title = "Effect of Elevation and Treatment on Pollen Viability"
  ) +  scale_fill_manual(
    values = c("C" = "#6c95ad", "H" = "#8a404e"),
    labels = c("C" = "Control", "H" = "Heat")
  ) +
  theme_bw()

# Raw and predicted data
fig_elev <- ggplot() +
  geom_point(data = dat, aes(x = Elevation, y = Prop.Viable, color = Treatment), alpha = 0.6) +
  geom_line(data = preds, aes(x = Elevation, y = fit, color = Treatment), size = 1.2) +
  geom_ribbon(data = preds, aes(x = Elevation, ymin = lower, ymax = upper, fill = Treatment), alpha = 0.2) +
  scale_color_manual(values = c("C" = "#6c95ad", "H" = "#8a404e"), labels = c("C" = "Control", "H" = "Heat"))+
  scale_fill_manual(values = c("C" = "#6c95ad", "H" = "#8a404e"), labels = c("C" = "Control", "H" = "Heat")) +
  labs(
    x = "Elevation (m)",
    y = "Predicted pollen viability"
  ) +
theme_bw(base_size = 14) 
fig_elev

ggsave(file = "fig_elev.svg", plot = fig_elev,  dpi = 300)


# Add in population codes if we want them?
pop_x <- dat %>%
  group_by(Pop) %>%
  summarise(
    x = mean(Elevation, na.rm = TRUE),
    .groups = "drop"
  )
ggplot() +
  geom_point(data = dat, aes(x = Elevation, y = Prop.Viable, color = Treatment), alpha = 0.6) +
  geom_line(data = preds, aes(x = Elevation, y = fit, color = Treatment), size = 1.2) +
  geom_ribbon(data = preds, aes(x = Elevation, ymin = lower, ymax = upper, fill = Treatment), alpha = 0.2) +
  scale_color_manual(values = c("C" = "#6c95ad", "H" = "#8a404e"))+
  scale_fill_manual(values = c("C" = "#6c95ad", "H" = "#8a404e"))+
  ggrepel::geom_text_repel(
    data = pop_x,
    aes(x = x, y = 1.02, label = Pop),
    inherit.aes = FALSE,
    direction = "x",
    nudge_y = 0,
    min.segment.length = 0,
    box.padding = 0.2,
    size = 3,
    color = "black"
  )+
  labs(
    x = "Elevation (m)",
    y = "Predicted pollen viability"
  ) +
  theme_bw()

## =======================
## Total pollen models
## ======================
# No significant elevation affect
m_total_elev <- glmmTMB(Total.Pollen ~ s_elev * Treatment + (1 | Treat_Block) + (1| Population), 
                   data = dat,
                   family = nbinom2())
summary(m_total_elev)
Anova(m_total_elev, type = "III")

# Significant treatment and some pop x treatment effects
m_total <- glmmTMB(Total.Pollen ~ Population * Treatment + (1 | Treat_Block), 
                   data = dat,
                   family = nbinom2())
summary(m_total)
Anova(m_total, type = "III")

# Total pollen
emm_total <- emmeans(m_total, ~ Treatment | Population, type = "link")
heat_total <- contrast(emm_total, method = "revpairwise")  # H - C on log scale = log(H/C)

# Prop pollen heat ratio using populations instead of s_elev

m_beta_pop <- glmmTMB(
  cbind(Viable.Pollen, Inviable.Pollen) ~ Population * Treatment + (1 | Treat_Block),
  data = dat,
  family = betabinomial()
)

emm_prop <- emmeans(m_beta_pop, ~ Treatment | Population, type = "link")
heat_prop <- contrast(emm_prop, method = "revpairwise")  # H - C on log scale = log(H/C)
ht <- as.data.frame(summary(heat_prop, infer =c(T,T)))


# Plot rate ratio with confidence intervals
# Dataframe from emmeans and Holm adjustment
ht <- as.data.frame(summary(heat_total, infer = c(TRUE, TRUE)))
ht$p_holm <- p.adjust(ht$p.value, method = "holm")

heat_ratio <- ht %>%
  transmute(
    Population,
    rr = exp(estimate),
    rr_low  = exp(asymp.LCL),
    rr_high = exp(asymp.UCL),
    p.value = p.value,
    p_holm = p_holm,
    sig = p_holm < 0.05
  )

# Plotted on log scale
hr <- ggplot(heat_ratio, aes(x = Population, y = rr, ymin = rr_low, ymax = rr_high, colour = sig)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(linewidth = 0.8) +
  scale_y_log10() +
  coord_flip() +
  scale_colour_manual(
    values = c(`TRUE` = "black", `FALSE` = "grey75"),
    breaks = c(TRUE, FALSE),
    labels = c("p < 0.05", "ns"),
    name = NULL
  ) +
  theme_bw() +
  labs(x = "Population", y = "Total Pollen Count Rate ratio (Heat / Control)")
hr 
ggsave(filename = "heat_ratio-TotalPollen.svg", hr, height = 5.5, width = 7, units = "in")

################################################
### M. guttatus and tilingii experiments 
###########################################
# ------------------------------------------------------------
# Maps of the two experiments
# ------------------------------------------------------------

us_map <- maps::map_data("state")

p1_sites <- pollen1 %>%
  group_by(Population) %>%
  summarise(
    Latitude = mean(Latitude, na.rm = TRUE),
    Longitude = mean(Longitude, na.rm = TRUE),
    Species_raw = first(as.character(Updated_species)),
    in_p1 = TRUE,
    .groups = "drop"
  ) %>%
  mutate(
    Species_clean = standardize_species(Species_raw),
    Species_plot = species_plot_label(Species_clean)
  )

p2_sites <- pollen2 %>%
  group_by(Population) %>%
  summarise(
    Latitude = mean(Latitude, na.rm = TRUE),
    Longitude = mean(Longitude, na.rm = TRUE),
    Species_raw = first(as.character(Species)),
    in_p2 = TRUE,
    .groups = "drop"
  ) %>%
  mutate(
    Species_clean = standardize_species(Species_raw),
    Species_plot = species_plot_label(Species_clean)
  )

pops_map <- full_join(
  p1_sites,
  p2_sites,
  by = "Population",
  suffix = c("_p1", "_p2")
) %>%
  mutate(
    Latitude = coalesce(Latitude_p1, Latitude_p2),
    Longitude = coalesce(Longitude_p1, Longitude_p2),
    Species_plot = coalesce(Species_plot_p1, Species_plot_p2),
    in_p1 = coalesce(in_p1, FALSE),
    in_p2 = coalesce(in_p2, FALSE),
    dataset_shape = case_when(
      in_p1 & in_p2 ~ "Experiment 1 & 2",
      in_p1 & !in_p2 ~ "Experiment 1 Only",
      !in_p1 & in_p2 ~ "Experiment 2 Only",
      TRUE ~ NA_character_
    ),
    region = case_when(
      Longitude < -113 &
        Latitude < 43 ~ "California",

      Longitude > -110 &
        Longitude < -102 &
        Latitude > 36 &
        Latitude < 41.5 ~ "Colorado",

      Longitude < -117 &
        Latitude > 45.5 ~ "Washington",

      TRUE ~ "Other"
    ),
    Species_plot = factor(
      Species_plot,
      levels = species_order
    )
  )

ca_pops <- pops_map %>%
  filter(region == "California")

co_pops <- pops_map %>%
  filter(region == "Colorado")

wa_pops <- pops_map %>%
  filter(region == "Washington")

# ------------------------------------------------------------
# Shared map scales
# ------------------------------------------------------------

experiment_shapes <- c(
  "Experiment 1 & 2" = 19,
  "Experiment 1 Only" = 1,
  "Experiment 2 Only" = 2
)

experiment_labels <- c(
  "Experiment 1 & 2" = "Experiments 1 + 2",
  "Experiment 1 Only" = "Experiment 1 only",
  "Experiment 2 Only" = "Experiment 2 only"
)

# ------------------------------------------------------------
# Colorado map
# ------------------------------------------------------------

co_xlim <- c(-107.35, -106.55)
co_ylim <- c(38.75, 39.65)

co_map <- ggplot() +
  geom_polygon(
    data = us_map,
    aes(x = long, y = lat, group = group),
    fill = "white",
    color = "grey80"
  ) +
  geom_path(
    data = us_map,
    aes(x = long, y = lat, group = group),
    color = "grey70",
    linewidth = 0.4
  ) +
  geom_point(
    data = co_pops,
    aes(
      x = Longitude,
      y = Latitude,
      color = Species_plot,
      shape = dataset_shape
    ),
    size = 3.5,
    stroke = 1
  ) +
  geom_text_repel(
    data = co_pops,
    aes(
      x = Longitude,
      y = Latitude,
      label = Population
    ),
    size = 2.6,
    box.padding = 0.7,
    point.padding = 0.35,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey40"
  ) +
  scale_shape_manual(
    values = experiment_shapes,
    labels = experiment_labels,
    drop = FALSE
  ) +
  scale_color_manual(
    values = species_colors,
    drop = FALSE
  ) +
  coord_fixed(
    xlim = co_xlim,
    ylim = co_ylim,
    expand = FALSE
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.8
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    )
  )

# ------------------------------------------------------------
# California map
# ------------------------------------------------------------

if (nrow(ca_pops) > 0) {
  ca_xlim <- c(
    min(ca_pops$Longitude, na.rm = TRUE) - 0.35,
    max(ca_pops$Longitude, na.rm = TRUE) + 0.35
  )

  ca_ylim <- c(
    min(ca_pops$Latitude, na.rm = TRUE) - 0.35,
    max(ca_pops$Latitude, na.rm = TRUE) + 0.35
  )
} else {
  ca_xlim <- c(-125, -114)
  ca_ylim <- c(32, 43)
}

ca_map <- ggplot() +
  geom_polygon(
    data = us_map,
    aes(x = long, y = lat, group = group),
    fill = "white",
    color = "grey80"
  ) +
  geom_path(
    data = us_map,
    aes(x = long, y = lat, group = group),
    color = "grey70",
    linewidth = 0.4
  ) +
  geom_point(
    data = ca_pops,
    aes(
      x = Longitude,
      y = Latitude,
      color = Species_plot,
      shape = dataset_shape
    ),
    size = 3.5,
    stroke = 1
  ) +
  geom_text_repel(
    data = ca_pops,
    aes(
      x = Longitude,
      y = Latitude,
      label = Population
    ),
    size = 2.6,
    box.padding = 0.2,
    point.padding = 0.1,
    max.overlaps = Inf,
    segment.color = NA
  ) +
  scale_shape_manual(
    values = experiment_shapes,
    labels = experiment_labels,
    drop = FALSE
  ) +
  scale_color_manual(
    values = species_colors,
    drop = FALSE
  ) +
  coord_fixed(
    xlim = ca_xlim,
    ylim = ca_ylim,
    expand = FALSE
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.8
    )
  )

# ------------------------------------------------------------
# Western overview map
# ------------------------------------------------------------

wa_pops <- wa_pops %>%
  mutate(
    label_x = if_else(
      Population == "BAG",
      Longitude + 0.5,
      Longitude
    ),
    label_y = if_else(
      Population == "BAG",
      Latitude + 0.15,
      Latitude
    )
  )

west_map <- ggplot() +
  geom_polygon(
    data = us_map,
    aes(x = long, y = lat, group = group),
    fill = "white",
    color = "grey80"
  ) +
  geom_path(
    data = us_map,
    aes(x = long, y = lat, group = group),
    color = "grey70",
    linewidth = 0.4
  ) +
  geom_point(
    data = pops_map,
    aes(
      x = Longitude,
      y = Latitude,
      color = Species_plot,
      shape = dataset_shape
    ),
    size = 3
  ) +
  geom_text(
    data = wa_pops,
    aes(
      x = label_x,
      y = label_y,
      label = Population
    ),
    size = 3
  ) +
  annotate(
    "rect",
    xmin = ca_xlim[1],
    xmax = ca_xlim[2],
    ymin = ca_ylim[1],
    ymax = ca_ylim[2],
    fill = NA,
    color = "black",
    linewidth = 0.5
  ) +
  annotate(
    "rect",
    xmin = co_xlim[1],
    xmax = co_xlim[2],
    ymin = co_ylim[1],
    ymax = co_ylim[2],
    fill = NA,
    color = "black",
    linewidth = 0.5
  ) +
  scale_shape_manual(
    values = experiment_shapes,
    labels = experiment_labels,
    drop = FALSE
  ) +
  scale_color_manual(
    values = species_colors,
    drop = FALSE
  ) +
  coord_fixed(
    xlim = c(-125.5, -104.5),
    ylim = c(31.5, 49.8),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.8
    )
  )

# ------------------------------------------------------------
# Legend
# ------------------------------------------------------------

legend_plot <- ggplot(
  pops_map,
  aes(
    x = 0,
    y = 0,
    color = Species_plot,
    shape = dataset_shape
  )
) +
  geom_point(size = 3, alpha = 0) +
  scale_color_manual(
    values = species_colors,
    name = "Species",
    breaks = species_order,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = experiment_shapes,
    name = "Experiment",
    labels = experiment_labels,
    drop = FALSE
  ) +
  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(alpha = 1, shape = 16)
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(alpha = 1, color = "black")
    )
  ) +
  coord_cartesian(
    xlim = c(-1, 1),
    ylim = c(-1, 1),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.spacing.x = unit(0.2, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )

legend_box <- cowplot::ggdraw() +
  cowplot::draw_plot(
    legend_plot,
    x = 0,
    y = 0,
    width = 1,
    height = 1
  )

top_row <- cowplot::plot_grid(
  west_map,
  ca_map,
  co_map,
  ncol = 3,
  rel_widths = c(1, 1, 1),
  labels = c("A.", "B.", "C."),
  label_x = 0.01,
  label_y = 1,
  hjust = 0,
  vjust = 1,
  label_fontface = "bold",
  label_size = 12
)

final_map <- cowplot::plot_grid(
  top_row,
  legend_box,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

save_plot(
  final_map,
  "fig_maps.svg",
  width = 9,
  height = 4.2
)

# Only in Pollen 2
pollen2 <- pollen2 %>% mutate(inviable_pollen = total_pollen - viable_pollen, 
                            prop_viable = viable_pollen / total_pollen)
# Only in Pollen 1
pollen1 <- pollen1 %>% mutate(inviable_pollen = total_pollen - num_viable,
                            viable_pollen = total_pollen - inviable_pollen,
                            prop_viable = num_viable/ total_pollen)

str(pollen1)
summary(pollen1)

str(pollen2)


# Pollen 2 only
guttatus <- pollen2 %>% filter(Group != "tilingii") %>% select(Population) %>% unique()


######################
## Data visualization
####################

# 36 C experiment = Pollen 1
treatment_labels <- c("0" = "Control", "1" = "Heat")

# Reorder by elevation for plotting
#dat <- pollen_env %>%
dat <- pollen1 %>%
 mutate(Pop = reorder(Population, Elevation, FUN = mean),
         Elev_num = as.numeric(as.character(Elevation)),   
         Elev_round = round(Elev_num),
         Elev_label = paste0(Elev_round),
         Population_ord = reorder(Elev_label, Elev_num, FUN = mean),
         elev_group = factor(if_else(Updated_species == "gutattus", "widespread", "restricted")),
         Treatment = factor(Treatment, levels = c("0","1"), labels = c("Control","Heat")))



ggplot(dat, aes(x = Pop, y = prop_viable, color = Treatment, fill = Treatment)) +
  # geom_violin(trim = FALSE, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
  #geom_jitter(aes(size = 0.5, alpha = 0.5)) +
  # facet_wrap(~ Treatment, labeller = labeller(Treatment = treatment_labels)) +
  labs(
    x = "Population",
    y = "Proportion Viable"
  ) + scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Control" = "#6c95ad", "Heat" = "#8a404e"),
                     labels = c("Control" = "Control", "Heat" = "Heat"))+
  scale_fill_manual(values = c("Control" = "#6c95ad", "Heat" = "#8a404e"),
                    labels = c("Control" = "Control", "Heat" = "Heat"))+
  theme_bw()

# Group by updated species
ggplot(dat, aes(x = Pop, y = prop_viable, color = Treatment, fill = Treatment)) +
  # geom_violin(trim = FALSE, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
  #geom_jitter(aes(size = 0.5, alpha = 0.5)) +
   facet_wrap(~ Updated_species, scales = "free_x") +
  labs(
    x = "Population",
    y = "Proportion Viable"
  ) + scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Control" = "#6c95ad", "Heat" = "#8a404e"),
                     labels = c("Control" = "Control", "Heat" = "Heat"))+
  scale_fill_manual(values = c("Control" = "#6c95ad", "Heat" = "#8a404e"),
                    labels = c("Control" = "Control", "Heat" = "Heat"))+
  theme_bw()

# Group by widespread vs restricted

ggplot(dat, aes(x = Pop, y = prop_viable, color = Treatment, fill = Treatment)) +
  # geom_violin(trim = FALSE, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
  #geom_jitter(aes(size = 0.5, alpha = 0.5)) +
  facet_wrap(~Species, scales = "free_x",
labeller = labeller(
  Species = c(
    gutattus = " M. guttatus",
    tilingii = "M. tilingii complex"
  )))+
  labs(
    x = "Population",
    y = "Proportion Viable"
  ) + scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Control" = "#6c95ad", "Heat" = "#8a404e"),
                     labels = c("Control" = "Control", "Heat" = "Heat"))+
  scale_fill_manual(values = c("Control" = "#6c95ad", "Heat" = "#8a404e"),
                    labels = c("Control" = "Control", "Heat" = "Heat"))+
  theme_bw()

# Indiviual species for presentation
gut <- dat %>% filter(Species == "gutattus")
til <- dat %>% filter(!Species == "gutattus")

til <- til %>% filter(!Pop =="CWC")

gut_plot <- ggplot(gut, aes(x = Pop, y = prop_viable, color = Treatment, fill = Treatment)) +
  # geom_violin(trim = FALSE, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
  geom_jitter(size = 1.5, alpha = 0.5) +
  labs(
    x = "Population",
    y = "Proportion Viable"
  ) + scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Control" = "#6c95ad", "Heat" = "#8a404e"),
                     labels = c("Control" = "Control", "Heat" = "Heat"))+
  scale_fill_manual(values = c("Control" = "#6c95ad", "Heat" = "#8a404e"),
                    labels = c("Control" = "Control", "Heat" = "Heat"))+
  theme_bw()

ggsave(filename = "pollen1_guttatus_pop_Viable.svg", plot = gut_plot, dpi = 300)

til_plot <- ggplot(til, aes(x = Pop, y = prop_viable, color = Treatment, fill = Treatment)) +
  # geom_violin(trim = FALSE, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5) +
  geom_jitter(size = 1.5, alpha = 0.5) +
  labs(
    x = "Population",
    y = "Proportion Viable"
  ) + scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Control" = "#6c95ad", "Heat" = "#8a404e"),
                     labels = c("Control" = "Control", "Heat" = "Heat"))+
  scale_fill_manual(values = c("Control" = "#6c95ad", "Heat" = "#8a404e"),
                    labels = c("Control" = "Control", "Heat" = "Heat"))+
  theme_bw()

ggsave(filename = "pollen1_tilingii_pop_Viable.svg", plot = til_plot, dpi = 300)


#################
## Models
################

# 36C Treatment pollen 1
# Elevational grouping (restricted vs widespread species) is poorly named
# Significant interaction between restricted vs widespread species and treatment
m_beta_p1 <- glmmTMB(
  cbind(viable_pollen, inviable_pollen) ~ elev_group * Treatment  + (1 + Treatment| Population),
  data = dat,
  family =betabinomial()
)
summary(m_beta_p1)
Anova(m_beta_p1, type = "III")
simulationOutput <- simulateResiduals(fittedModel = m_beta_p1, plot = T, re.form = NULL)
emmip(m_beta_p1, elev_group ~ Treatment, type = "response", CIs = TRUE) + theme_bw()

# Extract means from model
emmeans(m_beta_p1, ~ elev_group * Treatment, type = "response")

# Test random effects
m_beta_p1_no_pop <- glmmTMB(
  cbind(viable_pollen, inviable_pollen) ~ elev_group * Treatment + (1 + Treatment | Population),
  data = dat,
  family = betabinomial()
)

m_beta_p1_no_slope <- glmmTMB(
  cbind(viable_pollen, inviable_pollen) ~ elev_group * Treatment + (1 | Population),
  data = dat,
  family = betabinomial()
)

anova(m_beta_p1_no_pop, m_beta_p1_no_slope, m_beta_p1)

# Test population instead of elev group
m_beta_pop <- glmmTMB(
  cbind(viable_pollen, inviable_pollen) ~ Population * Treatment,
  data = gut,
  family = betabinomial(),
  control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4))
)
summary(m_beta_pop)

Anova(m_beta_pop, type = "III")

# Emmeans plot by species
emm_p1 <- emmeans(m_beta_p1, ~ elev_group * Treatment, type = "response")
emm_p1_df <- as.data.frame(emm_p1)

emm_p1_df

# Guttatus first
dat$elev_group <- factor(dat$elev_group, levels = c("widespread", "restricted"))

pollen_species <- ggplot(dat, aes(x = elev_group, y = prop_viable, fill = Treatment, color = Treatment)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.5,position = position_dodge(width = 0.7)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7),
    size = 1.4, alpha = 0.6
  ) +
  geom_point(
    data = emm_p1_df,
    aes(x = elev_group, y = prob, group = Treatment),
    inherit.aes = FALSE,
    position = position_dodge(width = 0.7),
    shape = 21, size = 3, fill = "white", colour = "black"
  ) +
  geom_errorbar(
    data = emm_p1_df,
    aes(x = elev_group, ymin = asymp.LCL, ymax = asymp.UCL, group = Treatment),
    inherit.aes = FALSE,
    position = position_dodge(width = 0.7),
    width = 0.1
  ) +
  scale_color_manual(values = c("Control" = "#6c95ad", "Heat" = "#8a404e")) +
  scale_fill_manual(values = c("Control" = "#6c95ad", "Heat" = "#8a404e"))+
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  labs(
    x = "Species",
    y = "Proportion viable pollen",
    fill = "Treatment",
    color = "Treatment"
  )+
  scale_x_discrete(labels = c(
    restricted = "M. tilingii complex",
    widespread = "M. guttatus"
  ))

ggsave(filename = "pollen_species.svg", plot = pollen_species, dpi = 300 )


# Population contrasts graph
emm_beta <- emmeans(m_beta_pop, ~ Treatment | Population, type = "link")
beta_contrasts <- contrast(emm_beta, method = "revpairwise")
beta_contrasts

pop_contrasts_df <- as.data.frame(summary(beta_contrasts, infer = c(TRUE, TRUE)))
pop_contrasts_df$p_holm <- p.adjust(pop_contrasts_df$p.value, method = "holm")
pop_contrasts_df

pop_or <- pop_contrasts_df %>%
  transmute(
    Population,
    or = exp(estimate),
    or_low = exp(asymp.LCL),
    or_high = exp(asymp.UCL),
    p.value,
    p_holm,
    sig = p.value < 0.05
  )

pop_or_plot <- pop_or %>%
  filter(is.finite(or), is.finite(or_low), is.finite(or_high),
         or > 0, or_low > 0, or_high > 0)

guttatus_p1 <- pollen1 %>% filter(Species == "gutattus") %>% select(Population) %>% unique() %>% pull()
pop_or_plot <- pop_or_plot %>% mutate(grouping = if_else(Population %in% guttatus_p1, "guttatus", "tilingii"))


# Pollen 1 population odds ratio plot
ggplot(pop_or, aes(x = Population, y = or, ymin = or_low, ymax = or_high, colour = sig)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(linewidth = 0.8) +
  scale_y_log10(limits = c(1e-5, 1.5)) +
  coord_flip() +
  scale_colour_manual(
    values = c(`TRUE` = "black", `FALSE` = "grey75"),
    breaks = c(TRUE, FALSE),
    labels = c("p < 0.05", "ns"),
    name = NULL
  ) +
  theme_bw() +
  labs(
    x = "Population",
    y = "Odds ratio for viable pollen (Heat / Control)"
  )




####################################
### Pollen 2 dataset
##################################
pollen2 <- pollen2 %>%
  mutate(Population = reorder(Population, Elevation, FUN = mean),
         Treatment = factor(Treatment, levels = c("C","H"), labels = c("Control","Heat")))


#Species grouping @ 30C
# Grouping all restricted tilingii complex spp. into one group
pollen2 <- pollen2 %>%
  mutate(type = ifelse(Species == "guttatus","guttatus", "tilingii"))
pollen2$type <- factor(pollen2$type, levels = c("guttatus","tilingii"))

pol2_pop_plot <- ggplot(pollen2, aes(x = Population, y = prop_viable, color = Treatment, fill = Treatment)) +
  #geom_violin(trim = FALSE, alpha = 0.3) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
  #geom_sina() +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7),
    size = 1.5,
    alpha = 0.7
  )+
  facet_wrap(~ type,scales = "free_x" ) +
  labs(
    title = "30C Heat",
    x = "Population",
    y = "Proportion Viable"
  ) + scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("C" = "#6c95ad", "H" = "#8a404e"))+
  scale_fill_manual(values = c("C" = "#6c95ad", "H" = "#8a404e"))+
  theme_bw()
pol2_pop_plot

ggsave(filename= "pollen2_pop_viable.svg", plot = pol2_pop_plot, dpi = 300)

#################
## Models Pollen2
################
# type = guttatus vs tilingii
# Group is elevation category
m_beta_group <- glmmTMB(cbind(viable_pollen, inviable_pollen) ~ Group * Treatment + (1 + Treatment | Population),
                     dispformula = ~ Treatment,
                     data = pollen2, family = betabinomial())

summary(m_beta_group)
Anova(m_beta_group, type = "III")

# Extract means
emmeans(m_beta_group, ~ Group * Treatment, type = "response")

simulationOutput2 <- simulateResiduals(fittedModel = m_beta_group, plot = TRUE, re.form = NULL)
testDispersion(simulationOutput2)
testUniformity(simulationOutput2)
testZeroInflation(simulationOutput2)

m_beta_species <- glmmTMB(cbind(viable_pollen, inviable_pollen) ~ type * Treatment + (1 + Treatment | Population),
                        dispformula = ~ Treatment,
                        data = pollen2, family = betabinomial())

summary(m_beta_species)
Anova(m_beta_species, type = "III")

simulationOutput2 <- simulateResiduals(fittedModel = m_beta_species, plot = TRUE, re.form = NULL)
testDispersion(simulationOutput2)
testUniformity(simulationOutput2)
testZeroInflation(simulationOutput2)

# Test random effects
m_beta_p2_no_pop <- glmmTMB(
  cbind(viable_pollen, inviable_pollen) ~ Group * Treatment + (1 + Treatment | Population),
  data = pollen2,
  family = betabinomial()
)

m_beta_p2_no_slope <- glmmTMB(
  cbind(viable_pollen, inviable_pollen) ~ Group * Treatment + (1 | Population),
  data = pollen2,
  family = betabinomial()
)

anova(m_beta_p2_no_pop, m_beta_p2_no_slope, m_beta_p2)

# Significant effect of heat but not by high elevation/low elevation/tilingii. 

# Plot the data and the emmeans
# Model-based estimates
emm_group <- emmeans(m_beta_group, ~ Group * Treatment, type = "response")
emm_group_df <- as.data.frame(emm_group)


plot_elev_group <- ggplot(pollen2, aes(x = Group, y = prop_viable, fill = Treatment, color = Treatment)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.5,
               position = position_dodge(width = 0.7)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7),
    size = 1.4, alpha = 0.6
  ) +
  geom_point(
    data = emm_group_df,
    aes(x = Group, y = prob, group = Treatment),
    inherit.aes = FALSE,
    position = position_dodge(width = 0.7),
    shape = 21, size = 3, fill = "white", colour = "black"
  ) +
  geom_errorbar(
    data = emm_group_df,
    aes(x = Group, ymin = asymp.LCL, ymax = asymp.UCL, group = Treatment),
    inherit.aes = FALSE,
    position = position_dodge(width = 0.7),
    width = 0.1
  ) +
  scale_color_manual(values = c("C" = "#6c95ad", "H" = "#8a404e")) +
  scale_fill_manual(values = c("C" = "#6c95ad", "H" = "#8a404e")) +
  theme_bw() +
  labs(
    x = "Group",
    y = "Proportion viable pollen",
    fill = "Treatment",
    color = "Treatment"
  ) +
  scale_x_discrete(labels = c(
    tilingii = "M. tilingii\ncomplex",
    high = "High-elevation\nM. guttatus",
    low = "Low-elevation\nM. guttatus"
  ))

ggsave(filename = "Pollen_elev_group.svg", plot_elev_group, dpi = 300)

# Guttatus only high vs low
g <- pollen2 %>% filter(Species == "guttatus")
m_gut <- glmmTMB(cbind(viable_pollen, inviable_pollen) ~ Group * Treatment + (1 + Treatment | Population),
                 dispformula = ~ Treatment,
                data = g, family = betabinomial())

summary(m_gut)
Anova(m_gut, type = "III")

# No significant difference between high and low guttatus

simulationOutput2 <- simulateResiduals(fittedModel = m_beta_p2, plot = TRUE, re.form = NULL)
testDispersion(simulationOutput2)
testUniformity(simulationOutput2)
testZeroInflation(simulationOutput2)

###################################
#### Flowers, Pollen, Heat
##################################

ggplot(pollen2, aes(x = Group, y = flowers_10days, fill = Treatment)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.6,
    outlier.shape = NA,
    color = "grey20"
  ) +
  geom_jitter(
    color = "lightgray",
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
    alpha = 0.35,
    size = 1
  ) +
  scale_x_discrete(labels = c(
    "guttatus" = "M. guttatus",
    "tilingii" = "M. tilingii complex"
  )) +
  theme_bw() +
  annotate("text", x = 1, y = 65, label = "ns", size = 5) +
  annotate("text", x = 2, y = 65, label = "**", size = 5) +
  labs(x = "Species", y = "Flower Number", fill = "Treatment", color = "Treatment") +
  scale_fill_manual(values = c("C" = "#6c95ad", "H" = "#8a404e")) +
  scale_color_manual(values = c("C" = "#6c95ad", "H" = "#8a404e"))


t.test(flowers_10days ~ Treatment, data = pollen2)
t.test(prop_viable ~ Treatment, data = pollen2)

# Flower # after 10 days
# type = species
m_flowers <- glmmTMB(
  flowers_10days ~ Treatment * Species + (1 | Population),
  data = pollen2,
  family = nbinom2()
)
summary(m_flowers)
Anova(m_flowers, type = "III")


# Group is high, low, tilingii
m_flowers_elev <- glmmTMB(
  flowers_10days ~ Treatment * Group + (1 | Population),
  data = pollen2,
  family = nbinom2()
)
summary(m_flowers_elev)
Anova(m_flowers_elev, type = "III")

AIC(m_flowers, m_flowers_elev)
# Difference in flower number is driven by species and exacerbated by heat


######################################
# Flower # to match species pollen 
##################################

# Emmeans plot by species
emm_p2_flwr <- emmeans(m_flowers_elev, ~ Group * Treatment, type = "response")
emm_p2_flwer_df <- as.data.frame(emm_p2_flwr)
emm_p2_flwer_df


flwr_species <- ggplot(pollen2, aes(x = Group, y = flowers_10days, fill = Treatment, color = Treatment)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.5,position = position_dodge(width = 0.7)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7),
    size = 1.4, alpha = 0.6
  ) +
  geom_point(
    data = emm_p2_flwer_df,
    aes(x = Group, y = response, group = Treatment),
    inherit.aes = FALSE,
    position = position_dodge(width = 0.7),
    shape = 21, size = 3, fill = "white", colour = "black"
  ) +
  geom_errorbar(
    data = emm_p2_flwer_df,
    aes(x = Group, ymin = asymp.LCL, ymax = asymp.UCL, group = Treatment),
    inherit.aes = FALSE,
    position = position_dodge(width = 0.7),
    width = 0.1
  ) +
  scale_color_manual(values = c("C" = "#6c95ad", "H" = "#8a404e")) +
  scale_fill_manual(values = c("C" = "#6c95ad", "H" = "#8a404e"))+
  theme_bw() +
  labs(
    x = "Group",
    y = "Flowers produced in 10 days",
    fill = "Treatment",
    color = "Treatment"
  )+
  scale_x_discrete(labels = c(
    tilingii = "M. tilingii\ncomplex",
    high = "High-elevation\nM. guttatus",
    low = "Low-elevation\nM. guttatus"
  ))

ggsave(filename = "pollen_species.svg", plot = pollen_species, dpi = 300 )

combined_p2 <- plot_elev_group + flwr_species +
  plot_layout(ncol = 2, guides = "collect")
combined_p2
ggsave(filename = "pollen2_plots.svg", plot = combined_p2, dpi = 300 )



#############################
### Table of pop codes 
############################

# A = Widespread comparison
pop1 <- pollen1 %>%
  transmute(
    Population,
    Species = Updated_species,
    Latitude,
    Longitude,
    Elevation,
    Experiment = "A"
  ) %>%
  distinct()

# B = elevation experiment for Guttatus
pop2 <- pollen2 %>%
  transmute(
    Population,
    Species,
    Latitude,
    Longitude,
    Elevation,
    Experiment = "B"
  ) %>%
  distinct()

# Sierra Nevada, Guttatus
pop3 <- p %>%
  transmute(
    Population = as.character(Pop),
    Species = "guttatus",
    Latitude,
    Longitude,
    Elevation,
    Experiment = "C"
  ) %>%
  distinct()

all_pops <- bind_rows(pop1, pop2, pop3) %>%
  group_by(Population, Species, Latitude, Longitude, Elevation) %>%
  summarise(
    Experiments = paste(sort(unique(Experiment)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(Species, Population)

all_pops


all_pops <- bind_rows(pop1, pop2, pop3) %>%
  mutate(
    Species = dplyr::recode(Species, "gutattus" = "guttatus")
  ) %>%
  arrange(Population, Experiment) %>%
  group_by(Population) %>%
  summarise(
    Species = first(Species),
    Latitude = first(Latitude),
    Longitude = first(Longitude),
    Elevation = first(Elevation),
    Experiments = paste(sort(unique(Experiment)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(Species, Population)
all_pops

# First tilingii species
all_pops <- all_pops %>%
  mutate(
    Species = factor(
      Species,
      levels = c("caespitosa", "minor", "tilingii", "guttatus")
    )
  ) %>%
  arrange(Species, Population)
all_pops %>% print(n=40)

############# 
# Format for output
#####################
all_pops %>%
  mutate(
    Species = case_when(
      Species == "guttatus" ~ "*M. guttatus*",
      Species == "tilingii" ~ "*M. tilingii*",
      Species == "minor" ~ "*M. minor*",
      Species == "caespitosa" ~ "*M. caespitosa*"
    )
  ) %>%
  gt() %>%
  fmt_markdown(columns = Species)

library(gt)

all_pops %>%
  gt() %>%
  cols_label(
    Population = "Population",
    Species = "Species",
    Latitude = "Latitude",
    Longitude = "Longitude",
    Elevation = "Elevation (m)",
    Experiments = "Experiment(s)"
  ) %>%
  fmt_number(
    columns = c(Latitude, Longitude),
    decimals = 2
  ) %>%
  fmt_number(
    columns = Elevation,
    decimals = 0
  ) %>%
  tab_header(
    title = "Populations used in growth chamber experiments"
  ) %>%
  tab_source_note(
    source_note = md(
      "**A** = Species comparison; **B** = Elevation-of-origin; **C** = Eastern Sierra Nevada"
    )
  )

all_pops


library(dplyr)
library(tidyr)
library(gt)

pop_table <- all_pops %>%
  mutate(
    A = ifelse(grepl("A", Experiments), "X", ""),
    B = ifelse(grepl("B", Experiments), "X", ""),
    C = ifelse(grepl("C", Experiments), "X", ""),
    Species = case_when(
      Species == "guttatus" ~ "*M. guttatus*",
      Species == "tilingii" ~ "*M. tilingii*",
      Species == "minor" ~ "*M. minor*",
      Species == "caespitosa" ~ "*M. caespitosa*"
    )
  ) %>%
  select(
    Population,
    Species,
    Latitude,
    Longitude,
    Elevation,
    A, B, C
  )

gt_table <- pop_table %>%
  gt() %>%
  fmt_markdown(columns = Species) %>%
  fmt_number(
    columns = c(Latitude, Longitude),
    decimals = 2
  ) %>%
  fmt_number(
    columns = Elevation,
    decimals = 0
  ) %>%
  cols_label(
    Population = "Population",
    Species = "Species",
    Latitude = "Latitude",
    Longitude = "Longitude",
    Elevation = "Elevation (m)",
    A = "A",
    B = "B",
    C = "C"
  ) %>%
  tab_header(
    title = md("**Populations used in growth chamber experiments**")
  ) %>%
  tab_source_note(
    source_note = md(
      "**A** = Species comparison experiment; **B** = Elevation-of-origin experiment; **C** = Eastern Sierra Nevada experiment."
    )
  ) %>%
  cols_align(
    align = "center",
    columns = c(A, B, C)
  )

library(writexl)
write_xlsx(pop_table, "Population_Table.xlsx")

gtsave(gt_table, "Population_Table.docx")
