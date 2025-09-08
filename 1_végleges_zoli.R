library(eurostat)
library(giscoR)
library(sf)
library(spdep)
library(dplyr)
library(tidyr)
library(stringr)
library(stringi)
library(janitor)
library(fuzzyjoin)
library(ggplot2)
library(tmap)
library(readr)
library(readxl)
library(mice)
library(writexl)
library(purrr)
library(rlang)
library(plm)
library(lmtest)
library(sandwich)
library(spatialreg)
library(splm)
library(ranger)
library(car)
library(factoextra)
library(mclust)
library(broom)
library(tibble)
library(pROC)
library(patchwork)
library(MatchIt)
library(cobalt)
library(urca)
library(vars)
library(conflicted)
library(grid)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")

#kép előkészítók

THEME_WORD <- theme_minimal(base_size = 11) +
  theme(
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.box         = "horizontal",
    legend.title       = element_text(size = 9),
    legend.text        = element_text(size = 8),
    legend.key.width   = unit(1.0, "lines"),
    legend.key.height  = unit(0.6, "lines"),
    plot.title         = element_text(face = "bold", hjust = 0.5),
    plot.subtitle      = element_text(color = "grey30"),
    plot.caption       = element_text(size = 9, color = "grey40"),
    plot.margin        = margin(t = 10, r = 10, b = 40, l = 10)
  )

theme_set(THEME_WORD)

bottom_guides <- guides(
  fill     = guide_legend(nrow = 1, byrow = TRUE),
  colour   = guide_legend(nrow = 1, byrow = TRUE),
  color    = guide_legend(nrow = 1, byrow = TRUE),
  linetype = guide_legend(nrow = 1, byrow = TRUE),
  shape    = guide_legend(nrow = 1, byrow = TRUE),
  size     = guide_legend(nrow = 1, byrow = TRUE)
)

with_bottom_legend <- function(p) {
  p +
    bottom_guides +
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      plot.margin      = margin(t = 10, r = 10, b = 40, l = 10)
    )
}

ggsave_word <- function(filename, plot = ggplot2::last_plot(), path = getwd(),
                        width = 6.5, height = 5.0, dpi = 300, ...) {

  ggplot2::ggsave(
    filename = filename,
    plot     = with_bottom_legend(plot),
    path     = path,
    units    = "in",
    width    = width,
    height   = height,
    dpi      = dpi,
    device   = "png",
    bg       = "white",
    ...
  )
}

save_base_word <- function(filename, path = getwd(),
                           width = 6.5, height = 5.0, dpi = 300, expr) {
  full <- file.path(path, filename)
  png(full, width = width, height = height, units = "in", res = dpi, bg = "white")
  on.exit(dev.off(), add = TRUE)
  force(expr)
  invisible(full)
}

#ADATELŐKÉSZÍTÉS
clean_key <- function(x) {
  x %>%
    str_replace_all("\\s+", " ") %>%
    str_replace_all("\\s*/\\s*", " / ") %>%
    str_remove_all("\\s*\\(statistical region 2016\\)\\s*") %>%
    str_trim() %>%
    str_to_lower() %>%
    stringi::stri_trans_general("Any-Latin; Latin-ASCII")
}
options(timeout = 600)
cache_dir <- file.path(Sys.getenv("LOCALAPPDATA"), "giscoR-cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
bad <- file.path(tempdir(), "giscoR", "NUTS_RG_20M_2024_4326_LEVL_2.geojson")
if (file.exists(bad)) unlink(bad, force = TRUE)

nuts2 <- get_eurostat_geospatial(
  output_class = "sf",
  resolution   = "20",
  nuts_level   = 2,
  year         = "2024",
  cache        = TRUE,
  update_cache = TRUE,     
  cache_dir    = cache_dir
) %>% 
  dplyr::filter(!CNTR_CODE %in% c("UK","GB"))

nuts_lookup <- nuts2 |>
  st_drop_geometry() |>
  transmute(NUTS_ID, CNTR_CODE,
            NAME_LATN = clean_key(NAME_LATN),
            NUTS_NAME = clean_key(NUTS_NAME)) |>
  pivot_longer(c(NAME_LATN, NUTS_NAME), names_to = "name_type",
               values_to = "name_clean") |>
  filter(!is.na(name_clean)) |>
  distinct(NUTS_ID, CNTR_CODE, name_clean)

#geoadat
read_indicator_long <- function(path) {
  df <- readxl::read_excel(path) %>%
    janitor::clean_names() %>%
    rename(region_label = time) %>%
    filter(region_label != "GEO (Labels)") %>%
    mutate(region_clean = clean_key(region_label))
  
  
  year_cols <- names(df)[stringr::str_detect(names(df), "^x?\\d{4}$")]
  
  df %>%
    
    mutate(across(all_of(year_cols), as.character)) %>%
    tidyr::pivot_longer(
      cols = all_of(year_cols),
      names_to = "year",
      values_to = "value",
      values_transform = list(value = as.character)  
    ) %>%
    mutate(
      year  = suppressWarnings(as.integer(stringr::str_remove(year, "^x"))),
      value = readr::parse_number(value, na = c(":", "", "NA"))
    ) %>%
    filter(!is.na(year))
}

#Namematchh
attach_nuts_keys <- function(long_tbl, max_dist = 0.08) {
  exact <- long_tbl %>%
    inner_join(nuts_lookup, by = c("region_clean" = "name_clean"))
  
  left <- long_tbl %>%
    anti_join(nuts_lookup, by = c("region_clean" = "name_clean")) %>%
    distinct(region_label, region_clean)
  
  if (nrow(left)) {
    fuzzy_hits <- fuzzyjoin::stringdist_inner_join(
      left, nuts_lookup,
      by = c("region_clean" = "name_clean"),
      method = "jw", max_dist = max_dist, distance_col = "dist"
    ) %>%
      group_by(region_label) %>%
      slice_min(dist, with_ties = FALSE) %>%
      ungroup() %>%
      select(region_clean, NUTS_ID, CNTR_CODE)
    
    bind_rows(
      exact,
      long_tbl %>% inner_join(fuzzy_hits, by = "region_clean")
    )
  } else {
    exact
  }
}

process_indicator <- function(path, indicator_name) {
  dat <- read_indicator_long(path) %>% attach_nuts_keys()
  
  dat %>%
    arrange(NUTS_ID, year) %>%
    group_by(NUTS_ID, year) %>%
    summarise(!!indicator_name := dplyr::first(na.omit(value), default = NA_real_),
              .groups = "drop")
}

#fájlok
files <- c(
  early_leavers   = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/early_leavers_total.xlsx",
  econ_act_rate   = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/econ_act_rates_total_total.xlsx",
  emp_rate        = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/emp_rate_total.xlsx",
  gdp_curr        = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/gdp_at_curr_mark_prices.xlsx",
  gerd            = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/GERD_all.xlsx",
  high_tech_emp   = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/high_tech_employment.xlsx",
  hh_income       = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/income_of_hh.xlsx",
  labour_prod_yoy = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/labour_prod_yearly_change.xlsx",
  life_exp_birth  = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/life_exp_at_birth_total.xlsx",
  migration       = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/migration.xlsx",
  urbanisation    = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/no_of_hh_degree_of_urb_total.xlsx",
  part_rate_educ  = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/part_rate_in_educ_total.xlsx",
  tertiary        = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/tertiary_total.xlsx",
  youth_unemp     = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/youth_unemp_total.xlsx",
  arop            = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/arop.xlsx"
)

# NAGY PANEL ÉPÍTÉS
indicator_tables <- lapply(names(files), function(nm) {
  message("Processing: ", nm)
  process_indicator(files[[nm]], nm)
}) %>% setNames(names(files))

wide_panel <- Reduce(function(x, y) full_join(x, y, by = c("NUTS_ID","year")),
                     indicator_tables)

unmatched_per_file <- lapply(names(files), function(nm) {
  raw <- read_indicator_long(files[[nm]])
  left <- raw %>%
    anti_join(nuts_lookup, by = c("region_clean" = "name_clean")) %>%
    distinct(region_label) %>%
    arrange(region_label)
  tibble(indicator = nm, n_unmatched = nrow(left), sample = paste(head(left$region_label, 5), collapse = "; "))
}) %>% bind_rows()

print(unmatched_per_file)

nuts2_with_panel <- nuts2 %>%
  left_join(wide_panel, by = "NUTS_ID")

final_full_dataset <- nuts2_with_panel[,c(3:5,12:28)]
final2023 <- final_full_dataset %>%
  filter(year == 2023) %>%
  st_transform(3035) %>%
  arrange(NUTS_ID)

#underlying háló megnézése, nem kapcsoltak kidobása
nb2023 <- poly2nb(final2023, queen = TRUE)


isolates <- which(card(nb2023) == 0)

isolated_ids <- final2023$NUTS_ID[isolates]
isolated_ids

final_full_dataset_clean <- final_full_dataset %>%
  filter(!NUTS_ID %in% isolated_ids)

nrow(final_full_dataset)-nrow(final_full_dataset_clean)
final_full_dataset <- final_full_dataset_clean
final_full_dataset$early_leavers_good <- 1/final_full_dataset$early_leavers
remove(nb2023, isolates, isolated_ids, final_full_dataset_clean, nuts_lookup, nuts2_with_panel, wide_panel, nuts2, unmatched_per_file, final2023)
#######################################################################
# MICE
final_no_geom <- final_full_dataset %>% st_drop_geometry()

id_vars <- c("NUTS_ID","NAME_LATN", "CNTR_CODE", "year")
num_vars <- setdiff(names(final_no_geom), id_vars)
final_no_geom[num_vars] <- lapply(final_no_geom[num_vars], as.numeric)
num_vars

pm <- quickpred(final_no_geom, exclude = id_vars)

set.seed(123)
imp <- mice(
  final_no_geom,
  m = 5,
  method = "pmm",
  maxit = 10,
  printFlag = TRUE,
  predictorMatrix = pm
)

imputed_once <- complete(imp, 1)

geom_key <- final_full_dataset %>%
  select(NUTS_ID, geometry) %>%
  distinct()                 

final <- imputed_once %>%
  left_join(geom_key, by = "NUTS_ID") %>%
  st_as_sf()                 

final_tbl <- final %>% st_drop_geometry()
colSums(is.na(final))
write_xlsx(final_tbl, "final_dataset.xlsx")

##### plot
tertiary_education_plot_2023 <- ggplot(final) +
  geom_sf(aes(fill = tertiary)) +
  scale_fill_viridis_c(option = "plasma", name = "Ráta (%)") +
  labs(title = "Felsőoktatási végzettséggel rendelkezők aránya (2023)") +
  coord_sf(xlim = c(-13, 45), ylim = c(32, 75), expand = FALSE) +
  bottom_guides 

ggsave_word(
  "tert_educ.png", plot = tertiary_education_plot_2023,
  path = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/figures"
)


###
# SIGMA ÉS BETA CONVERGENCE
###

#SIGMA
sigma_convergence <- final %>%
  st_drop_geometry() %>%
  group_by(year) %>%
  summarise(
    sd_early_leavers = sd(early_leavers, na.rm = TRUE),
    cv_early_leavers = sd(early_leavers, na.rm = TRUE) / mean(early_leavers, na.rm = TRUE),
    sd_part_rate     = sd(part_rate_educ, na.rm = TRUE),
    cv_part_rate     = sd(part_rate_educ, na.rm = TRUE) / mean(part_rate_educ, na.rm = TRUE),
    sd_tertiary      = sd(tertiary, na.rm = TRUE),
    cv_tertiary      = sd(tertiary, na.rm = TRUE) / mean(tertiary, na.rm = TRUE)
  )
sigma_convergence_long <- sigma_convergence %>%
  select(year, starts_with("cv")) %>%
  pivot_longer(-year, names_to = "indicator", values_to = "cv")

sigma_convergence_plot <- ggplot(sigma_convergence_long,
                                 aes(x = year, y = cv, color = indicator)) +
  geom_line() +
  labs(title = "σ-konvergencia (NUTS 2-es régiók közti szóródás)",
       y = "Variancia Koefficiense", x = "Év") +
  bottom_guides

ggsave_word(
  "sigma_conv.png", plot = sigma_convergence_plot,
  path = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/figures"
)

sigma_convergence_plot
#beta

beta_early <- final %>%
  st_drop_geometry() %>%
  filter(year %in% c(2013, 2023)) %>%
  select(NUTS_ID, year, early_leavers) %>%
  pivot_wider(names_from = year, values_from = early_leavers,
              names_prefix = "y") %>%
  mutate(
    growth = (log(y2023) - log(y2013)) / (2023 - 2013),
    log_initial = log(y2013)
  ) %>%
  filter(!is.na(growth), !is.na(log_initial))

beta_early_model <- lm(growth ~ log_initial, data = beta_early)
summary(beta_early_model)


beta_tertiary <- final %>%
  st_drop_geometry() %>%
  filter(year %in% c(2013, 2023)) %>%
  select(NUTS_ID, year, tertiary) %>%
  pivot_wider(names_from = year, values_from = tertiary,
              names_prefix = "y") %>%
  mutate(
    growth = (log(y2023) - log(y2013)) / (2023 - 2013),
    log_initial = log(y2013)
  ) %>%
  filter(!is.na(growth), !is.na(log_initial))

beta_tertiary_model <- lm(growth ~ log_initial, data = beta_tertiary)
summary(beta_tertiary_model)

#plot
beta_convergence_plot <- ggplot(beta_tertiary,
                                aes(x = log_initial, y = growth)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Beta konvergencia a felsőoktatási végzettséget szerzők arányában (2013–2023)",
    x = "A kezdeti arány logaritmusa (2013)",
    y = "Évesített növekedési ráta (logaritmizálva)"
  ) +
  bottom_guides

ggsave_word(
  "beta_conv.png", plot = beta_convergence_plot,
  path = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/figures"
)

beta_participation <- final %>%
  st_drop_geometry() %>%
  filter(year %in% c(2013, 2023)) %>%
  select(NUTS_ID, year, part_rate_educ) %>%
  pivot_wider(names_from = year, values_from = part_rate_educ,
              names_prefix = "y") %>%
  mutate(
    growth = (log(y2023) - log(y2013)) / (2023 - 2013),
    log_initial = log(y2013)
  ) %>%
  filter(!is.na(growth), !is.na(log_initial))

beta_participation_model <- lm(growth ~ log_initial, data = beta_participation)
summary(beta_participation_model)

###
# Spatial Autocorrelacion : Moran's I, Geary's C, Local Moran I
###
final$early_leavers <- final$early_leavers_good
colnames(final)
final <- final[,-20]
#2013
colSums(is.na(final))
final2013 <- final %>%
  filter(year == 2013) %>%
  st_transform(3035) %>%
  arrange(NUTS_ID)

nb2013 <- poly2nb(final2013, queen = TRUE)
lw2013 <- nb2listw(nb2013, style = "W", zero.policy = TRUE)

g2013 <- nb2lines(nb2013, coords = st_coordinates(st_centroid(final2013)), as_sf = TRUE)
st_crs(g2013) <- st_crs(final2013)

ggplot() +
  geom_sf(data = final2013, fill = "white", color = "grey60", linewidth = 0.3) +
  geom_sf(data = g2013, color = "black", linewidth = 0.4, alpha = 0.8) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "white"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
ggsave("network.png", plot = last_plot(), path = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/figures")

vals_early2013 <- final2013$early_leavers
vals_part2013  <- final2013$part_rate_educ
vals_tert2013  <- final2013$tertiary

#Moran’s I ÉS Geary’s C
moran.test(vals_early2013, lw2013, zero.policy = TRUE)
geary.test(vals_early2013, lw2013, zero.policy = TRUE)

moran.test(vals_part2013, lw2013, zero.policy = TRUE)
geary.test(vals_part2013, lw2013, zero.policy = TRUE)

moran.test(vals_tert2013, lw2013, zero.policy = TRUE)
geary.test(vals_tert2013, lw2013, zero.policy = TRUE)

#2018
final2018 <- final %>%
  filter(year == 2018) %>%
  st_transform(3035) %>%
  arrange(NUTS_ID)

nb2018 <- poly2nb(final2018, queen = TRUE)
lw2018 <- nb2listw(nb2018, style = "W", zero.policy = TRUE)

vals_early2018 <- final2018$early_leavers
vals_part2018  <- final2018$part_rate_educ
vals_tert2018  <- final2018$tertiary

moran.test(vals_early2018, lw2018, zero.policy = TRUE)
geary.test(vals_early2018, lw2018, zero.policy = TRUE)

moran.test(vals_part2018, lw2018, zero.policy = TRUE)
geary.test(vals_part2018, lw2018, zero.policy = TRUE)

moran.test(vals_tert2018, lw2018, zero.policy = TRUE)
geary.test(vals_tert2018, lw2018, zero.policy = TRUE)

#2023
final2023 <- final %>%
  filter(year == 2023) %>%
  st_transform(3035) %>%
  arrange(NUTS_ID)

nb2023 <- poly2nb(final2023, queen = TRUE)
lw2023 <- nb2listw(nb2023, style = "W", zero.policy = TRUE)

vals_early2023 <- final2023$early_leavers
vals_part2023  <- final2023$part_rate_educ
vals_tert2023  <- final2023$tertiary

moran.test(vals_early2023, lw2023, zero.policy = TRUE)
geary.test(vals_early2023, lw2023, zero.policy = TRUE)

moran.test(vals_part2023, lw2023, zero.policy = TRUE)
geary.test(vals_part2023, lw2023, zero.policy = TRUE)

moran.test(vals_tert2023, lw2023, zero.policy = TRUE)
geary.test(vals_tert2023, lw2023, zero.policy = TRUE)

# LOCMORAN
vals_early <- final2023$early_leavers
vals_early[is.na(vals_early)] <- mean(vals_early, na.rm = TRUE)

lmi_early <- localmoran(vals_early, lw2023, zero.policy = TRUE)
final2023$LMI_early_leavers <- lmi_early[, "Ii"]

vals_part <- final2023$part_rate_educ
vals_part[is.na(vals_part)] <- mean(vals_part, na.rm = TRUE)

lmi_part <- localmoran(vals_part, lw2023, zero.policy = TRUE)
final2023$LMI_part_rate_educ <- lmi_part[, "Ii"]

#tertiary
vals <- final2023$tertiary
vals[is.na(vals)] <- mean(vals, na.rm = TRUE)

z <- scale(vals)[,1]
lag_z <- lag.listw(lw2023, z, zero.policy = TRUE)


lmi <- localmoran(z, lw2023, zero.policy = TRUE)
pvals <- pvals <- p.adjust(lmi[,"Pr(z != E(Ii))"], method = "BH")
alpha <- 0.15

cluster <- case_when(
  pvals >= alpha                     ~ "Nem szignifikáns (10%-os szignifikancia szint)",
  z > 0  & lag_z > 0                 ~ "Hotspot (Magas-Magas)",
  z < 0  & lag_z < 0                 ~ "Coldspot (Alacsony-Alacsony)",
  z > 0  & lag_z < 0                 ~ "Magas-Alacsony (magas, alacsony környezetben)",
  z < 0  & lag_z > 0                 ~ "Alacsony-magas (alacsony, magas környezetben)",
  TRUE                               ~ "Nem szignifikáns (10%-os szignifikancia szint)"
)

final2023 <- final2023 %>%
  mutate(
    LMI_tertiary = lmi[,"Ii"],
    LMI_padj     = pvals,
    z_tertiary   = z,
    z_lag        = lag_z,
    lisa_cluster = factor(cluster,
                          levels = c("Hotspot (Magas-Magas)",
                                     "Coldspot (Alacsony-Alacsony)",
                                     "Magas-Alacsony (magas, alacsony környezetben)",
                                     "Alacsony-magas (alacsony, magas környezetben)",
                                     "Nem szignifikáns (10%-os szignifikancia szint)")
  )
)

p_tert_lisa <- ggplot(final2023) +
  geom_sf(aes(fill = lisa_cluster), color = "grey60", linewidth = 0.1) +
  scale_fill_manual(
    name = "LISA besorolás",
    values = c(
      "Hotspot (Magas-Magas)" = "#b2182b",
      "Coldspot (Alacsony-Alacsony)" = "#2166ac",
      "Magas-Alacsony (magas, alacsony környezetben)" = "#ef8a62",
      "Alacsony-magas (alacsony, magas környezetben)" = "green",
      "Nem szignifikáns (10%-os szignifikancia szint)" = "#f0f0f0"
    ),
    drop = FALSE
  ) +
  coord_sf(crs = 4326, xlim = c(-13, 45), ylim = c(32, 75), expand = FALSE, datum = NA) +
  labs(title = "LISA klaszterek — Felsőoktatási végzettséget szerzők aránya, 2023") +
  bottom_guides

ggsave_word(
  "tert_lisa_clusters.png", plot = p_tert_lisa,
  path = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/figures"
)


###
# Lineáris Panelmodellek
###

final_df <- final %>% st_drop_geometry()
final_df$NUTS_ID <- as.factor(final_df$NUTS_ID)

num_vars <- final_df %>%
  select(-NUTS_ID, -year) %>%
  select(where(is.numeric)) %>%
  names()

const_cols <- names(which(sapply(final_df[num_vars], function(x) sd(x, na.rm = TRUE) == 0)))
if (length(const_cols)) {
  message("Dropping constant columns: ", paste(const_cols, collapse = ", "))
  num_vars <- setdiff(num_vars, const_cols)
}

rhs_gdp <- setdiff(num_vars, "gdp_curr")
fmla_gdp <- as.formula(paste("gdp_curr ~", paste(rhs_gdp, collapse = " + ")))

rhs_emp <- setdiff(num_vars, "emp_rate")
fmla_emp <- as.formula(paste("emp_rate ~", paste(rhs_emp, collapse = " + ")))


pool_gdp <- plm(fmla_gdp, data = final_df, index = c("NUTS_ID","year"), model = "pooling")
pool_emp <- plm(fmla_emp, data = final_df, index = c("NUTS_ID","year"), model = "pooling")

fe_gdp <- plm(fmla_gdp, data = final_df, index = c("NUTS_ID","year"), model = "within", effect = "individual")
fe_emp <- plm(fmla_emp, data = final_df, index = c("NUTS_ID","year"), model = "within", effect = "individual")

twfe_gdp <- plm(fmla_gdp, data = final_df, index = c("NUTS_ID","year"), model = "within", effect = "twoways")
twfe_emp <- plm(fmla_emp, data = final_df, index = c("NUTS_ID","year"), model = "within", effect = "twoways")

re_gdp <- plm(fmla_gdp, data = final_df, index = c("NUTS_ID","year"), model = "random")
re_emp <- plm(fmla_emp, data = final_df, index = c("NUTS_ID","year"), model = "random")


vcov_pool_gdp <- vcovHC(pool_gdp, type = "HC1", cluster = "group")
vcov_pool_emp <- vcovHC(pool_emp, type = "HC1", cluster = "group")
vcov_fe_gdp   <- vcovHC(fe_gdp,   type = "HC1", cluster = "group")
vcov_fe_emp   <- vcovHC(fe_emp,   type = "HC1", cluster = "group")
vcov_re_gdp   <- vcovHC(re_gdp,   type = "HC1", cluster = "group")
vcov_re_emp   <- vcovHC(re_emp,   type = "HC1", cluster = "group")

vcov_twfe_gdp <- vcovSCC(twfe_gdp, type = "HC1") 
vcov_twfe_emp <- vcovSCC(twfe_emp, type = "HC1")

cat("\n=== Pooled OLS (GDP) ===\n"); print(coeftest(pool_gdp, vcov. = vcov_pool_gdp))
cat("\n=== Pooled OLS (EMP) ===\n"); print(coeftest(pool_emp, vcov. = vcov_pool_emp))

cat("\n=== FE (individual) GDP ===\n"); print(coeftest(fe_gdp, vcov. = vcov_fe_gdp))
cat("\n=== FE (individual) EMP ===\n"); print(coeftest(fe_emp, vcov. = vcov_fe_emp))

cat("\n=== Kétirányú FE GDP (DK SEs) ===\n"); print(coeftest(twfe_gdp, vcov. = vcov_twfe_gdp))
cat("\n=== Kétirányú FE EMP (DK SEs) ===\n"); print(coeftest(twfe_emp, vcov. = vcov_twfe_emp))

cat("\n=== Random Effects GDP ===\n"); print(coeftest(re_gdp, vcov. = vcov_re_gdp))
cat("\n=== Random Effects EMP ===\n"); print(coeftest(re_emp, vcov. = vcov_re_emp))


cat("\n--- F-test: pooled vs FE (GDP) ---\n"); print(pFtest(fe_gdp, pool_gdp))
cat("\n--- F-test: pooled vs FE (EMP) ---\n"); print(pFtest(fe_emp, pool_emp))

cat("\n--- LM test (BP): pooled vs RE (GDP) ---\n"); print(plmtest(pool_gdp, type = "bp"))
cat("\n--- LM test (BP): pooled vs RE (EMP) ---\n"); print(plmtest(pool_emp, type = "bp"))

cat("\n--- Hausman: FE vs RE (GDP) ---\n"); print(phtest(fe_gdp, re_gdp))
cat("\n--- Hausman: FE vs RE (EMP) ---\n"); print(phtest(fe_emp, re_emp))


cat("\n--- Wooldridge tes AR(1) (GDP FE) ---\n"); print(pdwtest(fe_gdp))
cat("\n--- Wooldridge tes AR(1) (EMP FE) ---\n"); print(pdwtest(fe_emp))
cat("\n--- Pesaran CD test (GDP FE) ---\n"); print(pcdtest(fe_gdp, test = "cd"))
cat("\n--- Pesaran CD test (EMP FE) ---\n"); print(pcdtest(fe_emp, test = "cd"))

# további teszteks

car::vif(lm(update(fmla_gdp, . ~ .), data = final_df)) %>% sort(decreasing = TRUE) -> vif_gdp
car::vif(lm(update(fmla_emp, . ~ .), data = final_df)) %>% sort(decreasing = TRUE) -> vif_emp

cat("\nBP (GDP FE):\n"); print(bptest(fe_gdp, studentize = TRUE))
cat("\nBP (EMP FE):\n"); print(bptest(fe_emp, studentize = TRUE))

cat("\nWooldridge AR(1) GDP FE:\n"); print(pdwtest(fe_gdp))
cat("\nWooldridge AR(1) EMP FE:\n"); print(pdwtest(fe_emp))
cat("\nPesaran CD GDP FE:\n"); print(pcdtest(fe_gdp, test = "cd"))
cat("\nPesaran CD EMP FE:\n"); print(pcdtest(fe_emp, test = "cd"))


### 
# Spatial panel modelek
###
years_full <- 2013:2023

panel_df <- final %>% st_drop_geometry() %>% arrange(NUTS_ID, year)

all_num <- panel_df %>%
  select(-NUTS_ID, -year) %>%
  select(where(is.numeric)) %>%
  names()

all_num <- setdiff(
  all_num,
  names(which(sapply(panel_df[all_num], function(x) sd(x, na.rm = TRUE) == 0)))
)

rhs_gdp  <- setdiff(all_num, "gdp_curr")
fmla_gdp <- as.formula(paste("gdp_curr ~", paste(rhs_gdp, collapse = " + ")))

data_gdp_bal <- panel_df %>%
  select(NUTS_ID, year, gdp_curr, all_of(rhs_gdp)) %>%
  filter(year %in% years_full) %>%
  group_by(NUTS_ID) %>%
  filter(n() == length(years_full)) %>%  
  ungroup() %>%
  tidyr::drop_na() %>%
  arrange(NUTS_ID, year)

rhs_emp  <- setdiff(all_num, "emp_rate")
fmla_emp <- as.formula(paste("emp_rate ~", paste(rhs_emp, collapse = " + ")))

data_emp_bal <- panel_df %>%
  select(NUTS_ID, year, emp_rate, all_of(rhs_emp)) %>%
  filter(year %in% years_full) %>%
  group_by(NUTS_ID) %>%
  filter(n() == length(years_full)) %>%
  ungroup() %>%
  tidyr::drop_na() %>%
  arrange(NUTS_ID, year)

ids_gdp <- unique(data_gdp_bal$NUTS_ID)
ids_emp <- unique(data_emp_bal$NUTS_ID)

data_gdp_bal$NUTS_ID <- factor(data_gdp_bal$NUTS_ID, levels = ids_gdp)
data_emp_bal$NUTS_ID <- factor(data_emp_bal$NUTS_ID, levels = ids_emp)

regions_sf <- final %>%
  filter(year == 2023) %>%
  select(NUTS_ID, geometry) %>%
  distinct()

reg_gdp <- regions_sf %>%
  filter(NUTS_ID %in% ids_gdp) %>%
  mutate(NUTS_ID = factor(NUTS_ID, levels = ids_gdp)) %>%
  arrange(NUTS_ID) %>%
  st_transform(3035)

reg_emp <- regions_sf %>%
  filter(NUTS_ID %in% ids_emp) %>%
  mutate(NUTS_ID = factor(NUTS_ID, levels = ids_emp)) %>%
  arrange(NUTS_ID) %>%
  st_transform(3035)


nb_gdp <- poly2nb(reg_gdp, queen = TRUE)

if (any(card(nb_gdp) == 0)) {
  coords <- st_coordinates(st_centroid(reg_gdp))
  nb_knn <- knearneigh(coords, k = 1) |> knn2nb()
  nb_gdp <- include.self(nb_gdp)  
  empty <- which(card(nb_gdp) == 0)
  nb_gdp[empty] <- nb_knn[empty]
}
lw_gdp <- nb2listw(nb_gdp, style = "W", zero.policy = TRUE)

#súly
nb_emp <- poly2nb(reg_emp, queen = TRUE)
if (any(card(nb_emp) == 0)) {
  coords <- st_coordinates(st_centroid(reg_emp))
  nb_knn <- knearneigh(coords, k = 1) |> knn2nb()
  nb_emp <- include.self(nb_emp)
  empty <- which(card(nb_emp) == 0)
  nb_emp[empty] <- nb_knn[empty]
}
lw_emp <- nb2listw(nb_emp, style = "W", zero.policy = TRUE)

stopifnot(length(lw_gdp$neighbours) == length(ids_gdp))
stopifnot(length(lw_emp$neighbours) == length(ids_emp))

# SAR
sar_panel_gdp <- spml(fmla_gdp, data = data_gdp_bal,
                      index = c("NUTS_ID","year"),
                      listw = lw_gdp,
                      model = "within", effect = "individual",
                      lag = TRUE)

sar_panel_emp <- spml(fmla_emp, data = data_emp_bal,
                      index = c("NUTS_ID","year"),
                      listw = lw_emp,
                      model = "within", effect = "individual",
                      lag = TRUE)

# SEM
sem_panel_gdp <- spml(fmla_gdp, data = data_gdp_bal,
                      index = c("NUTS_ID","year"),
                      listw = lw_gdp,
                      model = "within", effect = "individual",
                      spatial.error = "b")

sem_panel_emp <- spml(fmla_emp, data = data_emp_bal,
                      index = c("NUTS_ID","year"),
                      listw = lw_emp,
                      model = "within", effect = "individual",
                      spatial.error = "b")

#SDM
sdm_panel_gdp <- spml(fmla_gdp, data = data_gdp_bal,
                      index = c("NUTS_ID","year"),
                      listw = lw_gdp,
                      model = "within", effect = "individual",
                      lag = TRUE, Durbin = TRUE)

sdm_panel_emp <- spml(fmla_emp, data = data_emp_bal,
                      index = c("NUTS_ID","year"),
                      listw = lw_emp,
                      model = "within", effect = "individual",
                      lag = TRUE, Durbin = TRUE)

cat("\n=== PANEL SAR (gdp_curr) ===\n");  print(summary(sar_panel_gdp))
cat("\n=== PANEL SEM (gdp_curr) ===\n");  print(summary(sem_panel_gdp))
cat("\n=== PANEL SDM (gdp_curr) ===\n");  print(summary(sdm_panel_gdp))
cat("\n=== PANEL SAR (emp_rate) ===\n");  print(summary(sar_panel_emp))
cat("\n=== PANEL SEM (emp_rate) ===\n");  print(summary(sem_panel_emp))
cat("\n=== PANEL SDM (emp_rate) ===\n");  print(summary(sdm_panel_emp))


res_gdp <- residuals(sar_panel_gdp) 
res_gdp_2023 <- data_gdp_bal %>%
  mutate(.resid = as.numeric(res_gdp)) %>%
  filter(year == 2023) %>%
  arrange(NUTS_ID)  
cat("\nRezidumok Moran's I — SAR panel GDP (2023):\n")
print(moran.test(res_gdp_2023$.resid, lw_gdp, zero.policy = TRUE))

res_sem_g <- residuals(sem_panel_gdp)
res_sem_g_2023 <- data_gdp_bal %>%
  mutate(.resid = as.numeric(res_sem_g)) %>%
  filter(year == 2023) %>%
  arrange(NUTS_ID)
cat("\nRezidumokl Moran's I — SEM panel GDP (2023):\n")
print(moran.test(res_sem_g_2023$.resid, lw_gdp, zero.policy = TRUE))

res_emp <- residuals(sar_panel_emp)
res_emp_2023 <- data_emp_bal %>%
  mutate(.resid = as.numeric(res_emp)) %>%
  filter(year == 2023) %>%
  arrange(NUTS_ID)
cat("\nRezidumok Moran's I — SAR panel EMP (2023):\n")
print(moran.test(res_emp_2023$.resid, lw_emp, zero.policy = TRUE))


T_gdp <- length(unique(data_gdp_bal$year)) 
T_emp <- length(unique(data_emp_bal$year))

imp_sdm_panel_gdp <- impacts(sdm_panel_gdp, listw = lw_gdp, time = T_gdp, R = 500)
cat("\n=== Hatás — SDM PANEL (gdp_curr) ===\n")
print(summary(imp_sdm_panel_gdp, zstats = TRUE))

imp_sar_panel_gdp <- impacts(sar_panel_gdp, listw = lw_gdp, time = T_gdp, R = 500)
cat("\n=== Hatás — SAR PANEL (gdp_curr) ===\n")
print(summary(imp_sar_panel_gdp, zstats = TRUE))

imp_sdm_panel_emp <- impacts(sdm_panel_emp, listw = lw_emp, time = T_emp, R = 500)
cat("\n=== Hatáss — SDM PANEL (emp_rate) ===\n")
print(summary(imp_sdm_panel_emp, zstats = TRUE))

imp_sar_panel_emp <- impacts(sar_panel_emp, listw = lw_emp, time = T_emp, R = 500)
cat("\n=== Hatás — SAR PANEL (emp_rate) ===\n")
print(summary(imp_sar_panel_emp, zstats = TRUE))


####
# Random Forests, GDP és Employment 
###
set.seed(123)

stopifnot(exists("final"))

rf_df <- tryCatch({
  if (inherits(final, "sf")) sf::st_drop_geometry(final) else final
}, error = function(e) final)

if ("NUTS_ID" %in% names(rf_df)) rf_df$NUTS_ID <- as.factor(rf_df$NUTS_ID)

num_vars <- rf_df %>% select(where(is.numeric)) %>% names()
const_cols <- names(which(sapply(rf_df[num_vars], function(x) sd(x, na.rm = TRUE) == 0)))
if (length(const_cols)) {
  message("Dropping constant numeric columns: ", paste(const_cols, collapse = ", "))
  rf_df <- rf_df %>% select(-all_of(const_cols))
  num_vars <- setdiff(num_vars, const_cols)
}

stopifnot(all(c("gdp_curr", "emp_rate") %in% names(rf_df)))
stopifnot("year" %in% names(rf_df))


rmse  <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
r2cor <- function(y, yhat) cor(y, yhat, use = "complete.obs")^2

print_rf_report <- function(rf, name = "") {
  cat("\n=== Random Forest:", name, "===\n")
  cat("Fák      :", rf$num.trees, "\n")
  cat("OOB MSE     :", rf$prediction.error, "\n")
  cat("OOB RMSE    :", sqrt(rf$prediction.error), "\n")
  if (!is.null(rf$r.squared)) cat("OOB R^2     :", rf$r.squared, "\n")
  if (!is.null(rf$mtry)) cat("mtry        :", rf$mtry, "\n")
  if (!is.null(rf$min.node.size)) cat("min.node    :", rf$min.node.size, "\n")
  if (!is.null(rf$sample.fraction)) cat("sample.frac :", rf$sample.fraction, "\n")
  if (!is.null(rf$splitrule)) cat("splitrule   :", rf$splitrule, "\n")
  if (!is.null(rf$variable.importance)) {
    cat("\nLegfontosabb változók (permutatációs):\n")
    print(sort(rf$variable.importance, decreasing = TRUE)[1:min(15, length(rf$variable.importance))])
  }
}

make_formula <- function(target, rhs_vars) {
  as.formula(paste(target, "~", paste(rhs_vars, collapse = " + ")))
}


num_vars_all <- rf_df %>% select(where(is.numeric)) %>% names()

# GDP
rhs_gdp_full <- setdiff(num_vars_all, "gdp_curr")
fmla_gdp_full <- make_formula("gdp_curr", rhs_gdp_full)
rf_gdp_full <- ranger(
  formula    = fmla_gdp_full,
  data       = rf_df,
  importance = "impurity",
  num.trees  = 500,
  mtry       = floor(sqrt(length(rhs_gdp_full))),
  seed       = 123
)
cat("\n=== Alapmodell — GDP ===\n")
print_rf_report(rf_gdp_full, "GDP (full, impurity imp)")

# EMP
rhs_emp_full <- setdiff(num_vars_all, "emp_rate")
fmla_emp_full <- make_formula("emp_rate", rhs_emp_full)
rf_emp_full <- ranger(
  formula    = fmla_emp_full,
  data       = rf_df,
  importance = "impurity",
  num.trees  = 500,
  mtry       = floor(sqrt(length(rhs_emp_full))),
  seed       = 123
)
cat("\n=== Alapmodell — EMP ===\n")
print_rf_report(rf_emp_full, "Employment (full, imp)")


train <- rf_df %>% filter(year <= 2019)
test  <- rf_df %>% filter(year >= 2020)

stopifnot(nrow(train) > 0, nrow(test) > 0)

# GDP time-split
rhs_gdp <- setdiff(names(train %>% select(where(is.numeric))), "gdp_curr")
fmla_gdp <- make_formula("gdp_curr", rhs_gdp)

rf_gdp_time <- ranger(
  formula    = fmla_gdp,
  data       = train,
  num.trees  = 1000,
  importance = "permutation",
  seed       = 123
)
pred_gdp <- predict(rf_gdp_time, data = test)$predictions
cat("\n--- Tt split (train<=2019, test>=2020) ---\n")
print_rf_report(rf_gdp_time, "GDP (train)")
cat("\nGDP Test RMSE:", rmse(test$gdp_curr, pred_gdp))
cat("\nGDP Test R^2 :", r2cor(test$gdp_curr, pred_gdp), "\n")

# EMP time-split
rhs_emp <- setdiff(names(train %>% select(where(is.numeric))), "emp_rate")
fmla_emp <- make_formula("emp_rate", rhs_emp)

rf_emp_time <- ranger(
  formula    = fmla_emp,
  data       = train,
  num.trees  = 1000,
  importance = "permutation",
  seed       = 123
)
pred_emp <- predict(rf_emp_time, data = test)$predictions
cat("\n--- Tt split (train<=2019, test>=2020) ---\n")
print_rf_report(rf_emp_time, "Employment (train)")
cat("\nEMP Test RMSE:", rmse(test$emp_rate, pred_emp))
cat("\nEMP Test R^2 :", r2cor(test$emp_rate, pred_emp), "\n")



tune_rf <- function(formula, data_train, num_trees = 800,
                    mtry_grid, min_node_grid = c(5, 10, 20, 50),
                    sample_frac_grid = c(0.632, 0.8, 1.0),
                    seed = 123) {
  grid <- expand.grid(mtry = unique(pmax(1, mtry_grid)),
                      min.node.size = min_node_grid,
                      sample.fraction = sample_frac_grid)
  
  fit_one <- function(mtry, min.node.size, sample.fraction) {
    ranger(
      formula          = formula,
      data             = data_train,
      num.trees        = num_trees,
      mtry             = mtry,
      min.node.size    = min.node.size,
      sample.fraction  = sample.fraction,
      replace          = TRUE,
      importance       = "permutation",
      seed             = seed
    )
  }
  
  res <- lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    mdl <- fit_one(g$mtry, g$min.node.size, g$sample.fraction)
    cbind(grid[i, ], oob_mse = mdl$prediction.error, oob_rmse = sqrt(mdl$prediction.error),
          oob_r2 = mdl$r.squared, model = I(list(mdl)))
  })
  res_df <- do.call(rbind, res)
  res_df <- res_df[order(res_df$oob_rmse), ]
  rownames(res_df) <- NULL
  res_df
}

p_gdp <- length(rhs_gdp)
mtry_grid_gdp <- unique(c(floor(sqrt(p_gdp)), floor(p_gdp/4), floor(p_gdp/3), floor(p_gdp/2)))
mtry_grid_gdp <- mtry_grid_gdp[mtry_grid_gdp >= 1]

tune_gdp <- tune_rf(fmla_gdp, train,
                    num_trees = 800,
                    mtry_grid = mtry_grid_gdp,
                    min_node_grid = c(5, 10, 20, 50),
                    sample_frac_grid = c(0.632, 0.8, 1.0),
                    seed = 123)

cat("\n=== GDP javítás (top 5 var OOB RMSE szerint) ===\n")
print(head(tune_gdp[, c("mtry","min.node.size","sample.fraction","oob_rmse","oob_r2")], 5))

best_gdp <- tune_gdp[1, ]
rf_gdp_best <- best_gdp$model[[1]]
print_rf_report(rf_gdp_best, "GDP (best OOB)")

pred_gdp_best <- predict(rf_gdp_best, data = test)$predictions
cat("\nGDP (BEST) Test RMSE:", rmse(test$gdp_curr, pred_gdp_best))
cat("\nGDP (BEST) Test R^2 :", r2cor(test$gdp_curr, pred_gdp_best), "\n")

p_emp <- length(rhs_emp)
mtry_grid_emp <- unique(c(floor(sqrt(p_emp)), floor(p_emp/4), floor(p_emp/3), floor(p_emp/2)))
mtry_grid_emp <- mtry_grid_emp[mtry_grid_emp >= 1]

tune_emp <- tune_rf(fmla_emp, train,
                    num_trees = 800,
                    mtry_grid = mtry_grid_emp,
                    min_node_grid = c(5, 10, 20, 50),
                    sample_frac_grid = c(0.632, 0.8, 1.0),
                    seed = 123)

cat("\n=== EMP javítás (top 5 var OOB RMSE szerint) ===\n")
print(head(tune_emp[, c("mtry","min.node.size","sample.fraction","oob_rmse","oob_r2")], 5))

best_emp <- tune_emp[1, ]
rf_emp_best <- best_emp$model[[1]]
print_rf_report(rf_emp_best, "Employment (best OOB)")

pred_emp_best <- predict(rf_emp_best, data = test)$predictions
cat("\nEMP (BEST) Test RMSE:", rmse(test$emp_rate, pred_emp_best))
cat("\nEMP (BEST) Test R^2 :", r2cor(test$emp_rate, pred_emp_best), "\n")



###
# Unsupervised klaszter
###
educ2023 <- final %>%
  filter(year == 2023) %>%
  select(NUTS_ID, geometry, early_leavers, part_rate_educ, tertiary)

educ_scaled <- educ2023 %>%
  st_drop_geometry() %>%
  select(early_leavers, part_rate_educ, tertiary) %>%
  scale()

#Elbow és silhouette 

p_elbow <- factoextra::fviz_nbclust(
  educ_scaled, kmeans, method = "wss",
  k.max = 10       
) +
  labs(
    title   = "Klaszterszám kiválasztás — Könyök (WSS)",
    subtitle= "Standardizált oktatási indikátorok",
    x       = "Klaszterszám (K)",
    y       = "Within-cluster SSE (WSS)",
    caption = "Forrás: saját számítás, k-means"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "grey30"),
    plot.caption    = element_text(size = 9, color = "grey40"),
    panel.grid.minor= element_blank()
  )

p_sil <- factoextra::fviz_nbclust(
  educ_scaled, kmeans, method = "silhouette",
  k.max = 10
) +
  labs(
    title   = "Klaszterszám kiválasztás — Sziluett",
    subtitle= "Standardizált oktatási indikátorok",
    x       = "Klaszterszám (K)",
    y       = "Átlagos sziluett-érték",
    caption = "Forrás: saját számítás, k-means"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "grey30"),
    plot.caption    = element_text(size = 9, color = "grey40"),
    panel.grid.minor= element_blank()
  )

# Individual panels
p_elbow <- p_elbow + bottom_guides
p_sil   <- p_sil   + bottom_guides

ggsave_word(
  "elbow_wss.png", plot = p_elbow,
  path = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/figures"
)
ggsave_word(
  "silhouette.png", plot = p_sil,
  path = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/figures"
)

# Combined (collect legends to one bottom bar)
combined <- (p_elbow + p_sil + patchwork::plot_layout(guides = "collect")) &
  theme(legend.position = "bottom")

ggsave_word(
  "cluster_selection_combined.png", plot = combined, height = 3.8,
  path = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/figures"
)

#ARI
set.seed(123)
k <- 2
labels_list <- replicate(50, kmeans(educ_scaled, centers = k, nstart = 10)$cluster, simplify = FALSE)
ref <- labels_list[[1]]
ari_vals <- sapply(labels_list[-1], function(z) mclust::adjustedRandIndex(ref, z))
cat(sprintf("\nKlaszterek stabilitása (k=%d): átlag ARI = %.3f\n", k, mean(ari_vals)))

educ_mat <- educ2023 %>%
  st_drop_geometry() %>%
  select(early_leavers, part_rate_educ, tertiary) %>%
  scale()   


set.seed(123)
k <- 3
km_res <- kmeans(educ_mat, centers = k, nstart = 25)
summary(km_res)

educ2023$cluster <- factor(km_res$cluster)

cluster_plot <- ggplot(educ2023) +
  geom_sf(aes(fill = cluster), color = "grey60", linewidth = 0.1) +
  scale_fill_brewer(name = "Klaszter", palette = "Set1") +
  coord_sf(crs = 4326, xlim = c(-13, 45), ylim = c(32, 75), expand = FALSE, datum = NA) +
  labs(title = paste0("Oktatási indikátorok klaszterezése (", k, " csoport) — 2023")) +
  bottom_guides

ggsave_word(
  "cluster_map_good_2023.png", plot = cluster_plot,
  path = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/figures"
)



### SAVING MINDENT ###

dir.create("out/tables", recursive = TRUE, showWarnings = FALSE)

tidy_coeftest <- function(ct, model_name) {

  df <- tryCatch(as.data.frame(ct, stringsAsFactors = FALSE),
                 error = function(e) as.data.frame(as.matrix(ct)))
  if (!nrow(df)) {
    return(tibble::tibble(term = character(), estimate = numeric(),
                          std_error = numeric(), statistic = numeric(),
                          p_value = numeric(), model = character()))
  }
  
  rn <- rownames(df)
  if (is.null(rn)) rn <- paste0("param_", seq_len(nrow(df)))

  cn <- tolower(colnames(df))

  pick <- function(pat) {
    hit <- which(grepl(pat, cn))
    if (length(hit)) as.numeric(df[[hit[1]]]) else rep(NA_real_, nrow(df))
  }
  
  est  <- pick("^(estimate|coef|value)$")
  se   <- pick("^(std(\\.?\\s*error)?|se)$")
  stat <- pick("(t\\s*value|z\\s*value|stat)")
  pval <- pick("^(pr\\(>\\|[tz]\\|\\)|p(>|-)?value|^p$)")
  
  tibble::tibble(
    term      = rn,
    estimate  = est,
    std_error = se,
    statistic = stat,
    p_value   = pval,
    model     = model_name
  )
}


tidy_lm_coef <- function(lm_obj, model_name) {
  sm <- summary(lm_obj)$coefficients
  as.data.frame(sm) %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    rename(estimate = Estimate, std_error = `Std. Error`,
           statistic = `t value`, p_value = `Pr(>|t|)`) %>%
    mutate(model = model_name, .before = 1)
}

tidy_test <- function(name, stat, p, df = NA, param = NA) {
  tibble(test = name, statistic = as.numeric(stat), p_value = as.numeric(p),
         df = as.character(df), parameter = as.character(param))
}

safe_tidy_spatial <- function(summ, model_name) {

  out <-
    tryCatch({
      if (!is.null(summ$Coef)) as.data.frame(summ$Coef)
      else if (!is.null(summ$coefficients)) as.data.frame(summ$coefficients)
      else if (!is.null(summ$coef)) as.data.frame(summ$coef)
      else stop("No coefficient table found")
    }, error = function(e) NULL)
  
  if (is.null(out)) return(tibble())
  out %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%

    rename_with(~sub("Std\\.Error|Std\\. Error", "std_error", .x)) %>%
    rename_with(~sub("Pr\\(>|z\\)|Pr\\(>|t\\)|p\\.value", "p_value", .x)) %>%
    rename_with(~sub("^Estimate$", "estimate", .x)) %>%
    mutate(model = model_name, .before = 1)
}

###
# DATA MICE
###
sheet_data_assembly <- list(
  sigma_convergence  = sigma_convergence,
  sigma_cv_long      = sigma_convergence_long
)

###
# Beta konv
###
sheet_beta <- bind_rows(
  tidy_lm_coef(beta_early_model,       "beta_early (2013→2023)"),
  tidy_lm_coef(beta_tertiary_model,    "beta_tertiary (2013→2023)"),
  tidy_lm_coef(beta_participation_model,"beta_participation (2013→2023)")
)

###
# spatial autocorr
###
compute_mg <- function(dat_sf, varname, listw, year) {
  v <- dat_sf[[varname]]
  tibble(
    year = year, variable = varname,
    moran_i = unname(moran.test(v, listw, zero.policy = TRUE)$estimate[["Moran I statistic"]]),
    moran_p = moran.test(v, listw, zero.policy = TRUE)$p.value,
    geary_c = unname(geary.test(v, listw, zero.policy = TRUE)$estimate[["Geary C statistic"]]),
    geary_p = geary.test(v, listw, zero.policy = TRUE)$p.value
  )
}

vars_sp <- c("early_leavers","part_rate_educ","tertiary")

sheet_spatial_global <- bind_rows(
  map_df(vars_sp, ~compute_mg(final2013, .x, lw2013, 2013)),
  map_df(vars_sp, ~compute_mg(final2018, .x, lw2018, 2018)),
  map_df(vars_sp, ~compute_mg(final2023, .x, lw2023, 2023))
)

lisa2023 <- final2023 %>%
  st_drop_geometry() %>%
  select(NUTS_ID, NAME_LATN,
         LMI_early_leavers, LMI_part_rate_educ, LMI_tertiary)

###
# linmodell, fe, re és tesztek
###
ct_pool_gdp <- coeftest(pool_gdp, vcov. = vcov_pool_gdp)
ct_pool_emp <- coeftest(pool_emp, vcov. = vcov_pool_emp)
ct_fe_gdp   <- coeftest(fe_gdp,   vcov. = vcov_fe_gdp)
ct_fe_emp   <- coeftest(fe_emp,   vcov. = vcov_fe_emp)
ct_twfe_gdp <- coeftest(twfe_gdp, vcov. = vcov_twfe_gdp)
ct_twfe_emp <- coeftest(twfe_emp, vcov. = vcov_twfe_emp)
ct_re_gdp   <- coeftest(re_gdp,   vcov. = vcov_re_gdp)
ct_re_emp   <- coeftest(re_emp,   vcov. = vcov_re_emp)

sheet_panel_linear <- bind_rows(
  tidy_coeftest(ct_pool_gdp, "pooled_OLS_gdp_curr"),
  tidy_coeftest(ct_pool_emp, "pooled_OLS_emp_rate"),
  tidy_coeftest(ct_fe_gdp,   "FE_gdp_curr"),
  tidy_coeftest(ct_fe_emp,   "FE_emp_rate"),
  tidy_coeftest(ct_twfe_gdp, "TWFE_gdp_curr"),
  tidy_coeftest(ct_twfe_emp, "TWFE_emp_rate"),
  tidy_coeftest(ct_re_gdp,   "RE_gdp_curr"),
  tidy_coeftest(ct_re_emp,   "RE_emp_rate")
)

sheet_panel_tests <- bind_rows(
  tidy_test("F-test pooled vs FE (GDP)",  pFtest(fe_gdp, pool_gdp)$statistic,  pFtest(fe_gdp, pool_gdp)$p.value),
  tidy_test("F-test pooled vs FE (EMP)",  pFtest(fe_emp, pool_emp)$statistic,  pFtest(fe_emp, pool_emp)$p.value),
  tidy_test("LM BP pooled vs RE (GDP)",   plmtest(pool_gdp, type = "bp")$statistic, plmtest(pool_gdp, type = "bp")$p.value),
  tidy_test("LM BP pooled vs RE (EMP)",   plmtest(pool_emp, type = "bp")$statistic, plmtest(pool_emp, type = "bp")$p.value),
  tidy_test("Hausman FE vs RE (GDP)",     phtest(fe_gdp, re_gdp)$statistic,    phtest(fe_gdp, re_gdp)$p.value),
  tidy_test("Hausman FE vs RE (EMP)",     phtest(fe_emp, re_emp)$statistic,    phtest(fe_emp, re_emp)$p.value),
  tidy_test("Wooldridge AR(1) (GDP FE)",  pdwtest(fe_gdp)$statistic,           pdwtest(fe_gdp)$p.value),
  tidy_test("Wooldridge AR(1) (EMP FE)",  pdwtest(fe_emp)$statistic,           pdwtest(fe_emp)$p.value),
  tidy_test("Pesaran CD (GDP FE)",        pcdtest(fe_gdp, test = "cd")$statistic, pcdtest(fe_gdp, test = "cd")$p.value),
  tidy_test("Pesaran CD (EMP FE)",        pcdtest(fe_emp, test = "cd")$statistic, pcdtest(fe_emp, test = "cd")$p.value),
  tidy_test("BP heterosk (GDP FE)",       bptest(fe_gdp, studentize = TRUE)$statistic, bptest(fe_gdp, studentize = TRUE)$p.value),
  tidy_test("BP heterosk (EMP FE)",       bptest(fe_emp, studentize = TRUE)$statistic, bptest(fe_emp, studentize = TRUE)$p.value)
)


###
# SAR/SEM/SDM
###
sheet_spatial_panel <- bind_rows(
  safe_tidy_spatial(summary(sar_panel_gdp), "Panel_SAR_gdp_curr"),
  safe_tidy_spatial(summary(sem_panel_gdp), "Panel_SEM_gdp_curr"),
  safe_tidy_spatial(summary(sdm_panel_gdp), "Panel_SDM_gdp_curr"),
  safe_tidy_spatial(summary(sar_panel_emp), "Panel_SAR_emp_rate"),
  safe_tidy_spatial(summary(sem_panel_emp), "Panel_SEM_emp_rate"),
  safe_tidy_spatial(summary(sdm_panel_emp), "Panel_SDM_emp_rate")
)

###
# RF
###

rf_importance_df <- function(rf, model_label) {
  vi <- rf$variable.importance
  if (is.null(vi)) return(tibble(model = model_label, variable = character(), importance = numeric()))
  tibble(
    model      = model_label,
    variable   = names(vi),
    importance = as.numeric(vi)
  ) %>% arrange(desc(importance))
}

pred_df <- function(test_df, y_true, y_pred, id_cols = c("NUTS_ID","year"), label = "") {
  tibble(!!!test_df[id_cols], y_true = y_true, y_pred = as.numeric(y_pred)) %>%
    mutate(model = label, .before = 1)
}

rf_imp_gdp_full <- rf_importance_df(rf_gdp_full, "RF_gdp_curr_full_impurity")
rf_imp_emp_full <- rf_importance_df(rf_emp_full, "RF_emp_rate_full_impurity")

rf_imp_gdp_time <- rf_importance_df(rf_gdp_time, "RF_gdp_curr_time_permutation")
rf_imp_emp_time <- rf_importance_df(rf_emp_time, "RF_emp_rate_time_permutation")

rf_imp_gdp_best <- rf_importance_df(rf_gdp_best, "RF_gdp_curr_tuned_permutation")
rf_imp_emp_best <- rf_importance_df(rf_emp_best, "RF_emp_rate_tuned_permutation")

rf_importances <- bind_rows(
  rf_imp_gdp_full, rf_imp_emp_full,
  rf_imp_gdp_time, rf_imp_emp_time,
  rf_imp_gdp_best, rf_imp_emp_best
)

preds_gdp_time <- pred_df(test, test$gdp_curr, pred_gdp, label = "RF_gdp_curr_time")
preds_emp_time <- pred_df(test, test$emp_rate, pred_emp, label = "RF_emp_rate_time")

preds_gdp_best <- pred_df(test, test$gdp_curr, pred_gdp_best, label = "RF_gdp_curr_tuned")
preds_emp_best <- pred_df(test, test$emp_rate, pred_emp_best, label = "RF_emp_rate_tuned")

rf_predictions_test <- bind_rows(
  preds_gdp_time, preds_emp_time,
  preds_gdp_best, preds_emp_best
)

rf_perf <- bind_rows(
  tibble(
    model = c("RF_time_gdp_curr","RF_time_emp_rate"),
    RMSE  = c(sqrt(mean((pred_gdp - test$gdp_curr)^2, na.rm = TRUE)),
              sqrt(mean((pred_emp - test$emp_rate)^2,  na.rm = TRUE))),
    R2    = c(cor(pred_gdp, test$gdp_curr, use = "complete.obs")^2,
              cor(pred_emp, test$emp_rate,  use = "complete.obs")^2)
  ),
  tibble(
    model = c("RF_tuned_gdp_curr","RF_tuned_emp_rate"),
    RMSE  = c(sqrt(mean((pred_gdp_best - test$gdp_curr)^2, na.rm = TRUE)),
              sqrt(mean((pred_emp_best - test$emp_rate)^2,  na.rm = TRUE))),
    R2    = c(cor(pred_gdp_best, test$gdp_curr, use = "complete.obs")^2,
              cor(pred_emp_best, test$emp_rate,  use = "complete.obs")^2)
  )
)

tuning_gdp_full <- tune_gdp %>%
  mutate(model = "RF_gdp_curr", .before = 1) %>%
  select(model, mtry, min.node.size, sample.fraction, oob_mse, oob_rmse, oob_r2)

tuning_emp_full <- tune_emp %>%
  mutate(model = "RF_emp_rate", .before = 1) %>%
  select(model, mtry, min.node.size, sample.fraction, oob_mse, oob_rmse, oob_r2)

tuning_gdp_top  <- head(tuning_gdp_full[order(tuning_gdp_full$oob_rmse), ], 20)
tuning_emp_top  <- head(tuning_emp_full[order(tuning_emp_full$oob_rmse), ], 20)

dir.create("C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/models", recursive = TRUE, showWarnings = FALSE)
saveRDS(rf_gdp_full, file = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/models/rf_gdp_full.rds")
saveRDS(rf_emp_full, file = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/models/rf_emp_full.rds")
saveRDS(rf_gdp_time, file = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/models/rf_gdp_time.rds")
saveRDS(rf_emp_time, file = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/models/rf_emp_time.rds")
saveRDS(rf_gdp_best, file = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/models/rf_gdp_tuned.rds")
saveRDS(rf_emp_best, file = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/models/rf_emp_tuned.rds")

###
# kmeans klaszterek
###
cluster_labels <- educ2023 %>%
  st_drop_geometry() %>%
  select(NUTS_ID) %>%
  mutate(cluster = as.integer(educ2023$cluster))

cluster_centers <- as_tibble(km_res$centers, .name_repair = "minimal") %>%
  mutate(cluster = 1:nrow(.), .before = 1)

cluster_summary <- tibble(
  k = nrow(cluster_centers),
  tot_withinss = sum(km_res$withinss),
  betweenss    = km_res$betweenss,
  totss        = km_res$totss
)

cluster_stability <- tibble(
  k = k,
  mean_ARI = mean(ari_vals)
)

###
# EXCEL
###
sheets <- list(
  "Sigma_convergence"        = sigma_convergence,
  "Beta_convergence"         = sheet_beta,
  "Spatial_global"           = sheet_spatial_global,
  "LISA_2023_values"         = lisa2023,
  "Panel_linear_coefs"       = sheet_panel_linear,
  "Panel_spec_tests"         = sheet_panel_tests,
  "Spatial_panel_coefs"      = sheet_spatial_panel,
  "Cluster_labels_2023"      = cluster_labels,
  "Cluster_centers_2023"     = cluster_centers,
  "Cluster_summary"          = cluster_summary,
  "Cluster_stability"        = cluster_stability,
  "RF_importance"          = rf_importances,
  "RF_perf_summary"        = rf_perf,
  "RF_predictions_test"    = rf_predictions_test,
  "RF_tuning_GDP_top20"    = tuning_gdp_top,
  "RF_tuning_EMP_top20"    = tuning_emp_top,
  "RF_tuning_GDP_full"     = tuning_gdp_full,
  "RF_tuning_EMP_full"     = tuning_emp_full
)

write_xlsx(sheets, path = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/results_master.xlsx")
message("Minden mentve")

options(timeout = 600)


###
# PSM nemekre
###

clean_key <- function(x) {
  x %>%
    str_replace_all("\\s+", " ") %>%
    str_replace_all("\\s*/\\s*", " / ") %>%
    str_remove_all("\\s*\\(statistical region 2016\\)\\s*") %>%
    str_trim() %>%
    str_to_lower() %>%
    stringi::stri_trans_general("Any-Latin; Latin-ASCII")
}

z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

cache_dir <- file.path(Sys.getenv("LOCALAPPDATA"), "giscoR-cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
bad <- file.path(tempdir(), "giscoR", "NUTS_RG_20M_2024_4326_LEVL_2.geojson")
if (file.exists(bad)) unlink(bad, force = TRUE)

nuts2 <- get_eurostat_geospatial(
  output_class = "sf",
  resolution   = "20",
  nuts_level   = 2,
  year         = "2024",
  cache        = TRUE,
  update_cache = TRUE,
  cache_dir    = cache_dir
) %>% 
  dplyr::filter(!CNTR_CODE %in% c("UK","GB"))

nuts_lookup <- nuts2 |>
  st_drop_geometry() |>
  transmute(NUTS_ID, CNTR_CODE,
            NAME_LATN = clean_key(NAME_LATN),
            NUTS_NAME = clean_key(NUTS_NAME)) |>
  pivot_longer(c(NAME_LATN, NUTS_NAME), names_to = "name_type",
               values_to = "name_clean") |>
  filter(!is.na(name_clean)) |>
  distinct(NUTS_ID, CNTR_CODE, name_clean)

read_indicator_long <- function(path) {
  df <- readxl::read_excel(path) %>%
    janitor::clean_names() %>%
    rename(region_label = time) %>%
    filter(region_label != "GEO (Labels)") %>%
    mutate(region_clean = clean_key(region_label))
  
  year_cols <- names(df)[stringr::str_detect(names(df), "^x?\\d{4}$")]
  
  df %>%
    mutate(across(all_of(year_cols), as.character)) %>%
    tidyr::pivot_longer(
      cols = all_of(year_cols),
      names_to = "year",
      values_to = "value",
      values_transform = list(value = as.character)
    ) %>%
    mutate(
      year  = suppressWarnings(as.integer(stringr::str_remove(year, "^x"))),
      value = readr::parse_number(value, na = c(":", "", "NA"))
    ) %>%
    filter(!is.na(year))
}

attach_nuts_keys <- function(long_tbl, max_dist = 0.08) {
  exact <- long_tbl %>%
    inner_join(nuts_lookup, by = c("region_clean" = "name_clean"))
  
  left <- long_tbl %>%
    anti_join(nuts_lookup, by = c("region_clean" = "name_clean")) %>%
    distinct(region_label, region_clean)
  
  if (nrow(left)) {
    fuzzy_hits <- fuzzyjoin::stringdist_inner_join(
      left, nuts_lookup,
      by = c("region_clean" = "name_clean"),
      method = "jw", max_dist = max_dist, distance_col = "dist"
    ) %>%
      group_by(region_label) %>%
      slice_min(dist, with_ties = FALSE) %>%
      ungroup() %>%
      select(region_clean, NUTS_ID, CNTR_CODE)
    
    bind_rows(
      exact,
      long_tbl %>% inner_join(fuzzy_hits, by = "region_clean")
    )
  } else {
    exact
  }
}

process_indicator <- function(path, indicator_name) {
  dat <- read_indicator_long(path) %>% attach_nuts_keys()
  dat %>%
    arrange(NUTS_ID, year) %>%
    group_by(NUTS_ID, year) %>%
    summarise(!!indicator_name := dplyr::first(na.omit(value), default = NA_real_),
              .groups = "drop")
}

build_sex_panel <- function(files_named_vec) {

  indicator_tables <- lapply(names(files_named_vec), function(nm) {
    message("Processing: ", nm)
    process_indicator(files_named_vec[[nm]], nm)
  }) %>% setNames(names(files_named_vec))
  
  wide_panel <- Reduce(function(x, y) full_join(x, y, by = c("NUTS_ID","year")),
                       indicator_tables)

  nuts2_with_panel <- nuts2 %>% left_join(wide_panel, by = "NUTS_ID")

  final_full_dataset <- nuts2_with_panel %>%
    dplyr::select(
      NUTS_ID, CNTR_CODE, year,

      tidyselect::any_of(c(
        "early_leavers","early_leavers_males","early_leavers_females",
        "part_rate_in_educ","part_rate_in_educ_males","part_rate_in_educ_females",
        "part_rate_educ_males","part_rate_educ_females",
        "tertiary","tertiary_males","tertiary_females",
        "econ_act_rate","econ_act_rate_males","econ_act_rate_females",
        "econ_act_rates","econ_act_rates_males_total","econ_act_rates_females_total",
        "emp_rate","emp_rate_males","emp_rate_females",
        "employment_rate","employment_rate_males","employment_rate_females",
        "life_exp_at_birth","life_exp_at_birth_males","life_exp_at_birth_females",
        "life_exp_birth","life_exp_birth_males","life_exp_birth_females"
      ))
    )

  final2023 <- final_full_dataset %>%
    dplyr::filter(year == 2023) %>%
    sf::st_transform(3035) %>%
    dplyr::arrange(NUTS_ID)
  
  nb2023 <- spdep::poly2nb(final2023, queen = TRUE)
  isolates <- which(spdep::card(nb2023) == 0)
  isolated_ids <- final2023$NUTS_ID[isolates]

  final_full_dataset %>%
    dplyr::filter(!NUTS_ID %in% isolated_ids)
}

files_male <- c(
  early_leavers_males    = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/males/early_leavers_males.xlsx",
  part_rate_educ_males   = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/males/part_rate_in_educ_males.xlsx",
  tertiary_males         = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/males/tertiary_males.xlsx",
  econ_act_rate_males    = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/males/econ_act_rates_males_total.xlsx",
  emp_rate_males         = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/males/emp_rate_males.xlsx",
  life_exp_birth_males   = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/males/life_exp_at_birth_males.xlsx"
)

files_female <- c(
  early_leavers_females  = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/females/early_leavers_females.xlsx",
  part_rate_educ_females = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/females/part_rate_in_educ_females.xlsx",
  tertiary_females       = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/females/tertiary_females.xlsx",
  econ_act_rate_females  = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/females/econ_act_rates_females_total.xlsx",
  emp_rate_females       = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/females/emp_rate_females.xlsx",
  life_exp_birth_females = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/females/life_exp_at_birth_females.xlsx"
)

final_full_dataset_males   <- build_sex_panel(files_male)
final_full_dataset_females <- build_sex_panel(files_female)

rename_to_common <- function(df) {
  df %>%
    dplyr::rename(
      early_leavers = tidyselect::any_of(c(
        "early_leavers","early_leavers_males","early_leavers_females"
      )),
      part_rate_in_educ = tidyselect::any_of(c(
        "part_rate_in_educ","part_rate_educ_males","part_rate_in_educ_males",
        "part_rate_educ_females","part_rate_in_educ_females"
      )),
      tertiary_attain = tidyselect::any_of(c(
        "tertiary","tertiary_males","tertiary_females"
      )),
      econ_act_rate = tidyselect::any_of(c(
        "econ_act_rate","econ_act_rate_males","econ_act_rate_females",
        "econ_act_rates","econ_act_rates_males_total","econ_act_rates_females_total"
      )),
      emp_rate = tidyselect::any_of(c(
        "emp_rate","emp_rate_males","emp_rate_females",
        "employment_rate","employment_rate_males","employment_rate_females",
        "emp_rates","employment_rates"
      )),
      life_exp_at_birth = tidyselect::any_of(c(
        "life_exp_at_birth","life_exp_at_birth_males","life_exp_at_birth_females",
        "life_exp_birth","life_exp_birth_males","life_exp_birth_females"
      ))
    )
}

men23 <- final_full_dataset_males %>%
  rename_to_common() %>%
  filter(year == 2023) %>%
  mutate(gender = "male")

women23 <- final_full_dataset_females %>%
  rename_to_common() %>%
  filter(year == 2023) %>%
  mutate(gender = "female")

common_ids <- intersect(men23$NUTS_ID, women23$NUTS_ID)
men23   <- filter(men23,   NUTS_ID %in% common_ids)
women23 <- filter(women23, NUTS_ID %in% common_ids)

needed_covars <- c("econ_act_rate","emp_rate","life_exp_at_birth")
for (nm in needed_covars) if (!nm %in% names(men23))   men23[[nm]]   <- NA_real_
for (nm in needed_covars) if (!nm %in% names(women23)) women23[[nm]] <- NA_real_


edu23 <- bind_rows(men23, women23) %>%
  dplyr::select(tidyselect::any_of(c(
    "NUTS_ID","CNTR_CODE","gender",
    "early_leavers","part_rate_in_educ","tertiary_attain",
    "econ_act_rate","emp_rate","life_exp_at_birth"
  )))

edu23 <- edu23 %>%
  mutate(
    early_leavers     = as.numeric(early_leavers),
    part_rate_in_educ = as.numeric(part_rate_in_educ),
    tertiary_attain   = as.numeric(tertiary_attain)
  ) %>%
  tidyr::drop_na(early_leavers, part_rate_in_educ, tertiary_attain) %>%
  mutate(

    z_early_leavers    =  z(early_leavers),
    z_participation    =  z(part_rate_in_educ),
    z_tertiary         =  z(tertiary_attain),
    z_econ_act_rate    =  z(econ_act_rate),
    z_emp_rate         =  z(emp_rate),
    z_life_exp_at_birth=  z(life_exp_at_birth),
    edu_profile_eqw    = rowMeans(cbind(z_early_leavers, z_participation, z_tertiary), na.rm = TRUE)
  )

psm_df_full <- edu23 %>%
  mutate(treat_female = as.integer(gender == "female")) %>%
  sf::st_drop_geometry()

covars <- c("econ_act_rate","emp_rate","life_exp_at_birth","CNTR_CODE")
psm_df_model <- psm_df_full %>%
  dplyr::select(treat_female, edu_profile_eqw, dplyr::all_of(covars)) %>%
  mutate(CNTR_CODE = as.factor(CNTR_CODE)) %>%
  filter(complete.cases(.[, c("treat_female", covars)]))

form_ps <- treat_female ~ econ_act_rate + emp_rate + life_exp_at_birth

logit <- glm(form_ps, data = psm_df_model, family = binomial(), na.action = na.exclude)
psm_df_model$psscore <- predict(logit, newdata = psm_df_model, type = "response")
summary(logit)

matching <- matchit(
  form_ps,
  data     = psm_df_model,
  method   = "nearest",
  distance = "logit",      
  exact    = ~ CNTR_CODE, 
  replace  = TRUE,
  ratio    = 1
)

summary(matching)
matched <- match.data(matching)

stopifnot("edu_profile_eqw" %in% names(matched))

treatment <- matched$edu_profile_eqw[matched$treat_female == 1]
control   <- matched$edu_profile_eqw[matched$treat_female == 0]

t.test(treatment, control, conf.level = 0.95)

m1 <- lm(edu_profile_eqw ~ econ_act_rate + emp_rate + life_exp_at_birth + CNTR_CODE,
         data = psm_df_model)
betas <- summary(m1)$coefficients
betas_df <- data.frame(
  Term = rownames(betas),
  Estimate = betas[, "Estimate"]
) %>% filter(Term != "(Intercept)")

coef_plot <- ggplot(betas_df, aes(x = reorder(Term, Estimate), y = Estimate)) +
  geom_col() +
  coord_flip() +
  labs(title = "Regressziós bétaértékek (Gender blokk)",
       x = "Változók", y = "Becslés") +
  bottom_guides

ggsave_word(
  "gender_coeffs.png", plot = coef_plot,
  path = "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/figures"
)



###
# VAR magyar
###

magyar <- read_excel("C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/magyar/magyar_okt_idősor.xlsx")

# Csak: Tanév (1), Középisk. összesen (6), Érettségi összesen (12)
magyar_var <- magyar[, c(1,6,12)] %>%
  rename(
    total_students = `Középiskolában tanuló összesen`,
    graduates      = `Érettségi vizsgát tett tanuló összesen`
  ) %>%
  arrange(Tanév)

# Éves ts (szintek)
ts_magyar <- ts(magyar_var[, c("total_students","graduates")],
                start = min(magyar_var$Tanév), frequency = 1)

adf_students  <- ur.df(ts_magyar[,"total_students"], type = "drift", selectlags = "AIC")
adf_graduates <- ur.df(ts_magyar[,"graduates"],      type = "drift", selectlags = "AIC")
cat("\n--- ADF total_students ---\n");  print(summary(adf_students))
cat("\n--- ADF graduates ---\n");       print(summary(adf_graduates))

#Johansen kointegrációs teszt
jo <- ca.jo(ts_magyar, type = "trace", K = 2, ecdet = "const", spec = "transitory")
cat("\n--- Johansen trace test ---\n"); print(summary(jo))

crit_cols <- colnames(jo@cval)
five_col  <- if ("5pct" %in% crit_cols) "5pct" else if ("5%" %in% crit_cols) "5%" else crit_cols[grep("5", crit_cols)[1]]
rank_5pct <- sum(jo@teststat > jo@cval[, five_col], na.rm = TRUE)
cat("\n>>> Becsült kointegrációs rang (5%): r =", rank_5pct, "\n")

diff_magyar <- diff(ts_magyar)                     
diff_df <- as.data.frame(diff_magyar)              
colnames(diff_df) <- c("y1","y2")                 

lagselect <- VARselect(diff_df, lag.max = 10, type = "const")
p_aic <- suppressWarnings(as.integer(lagselect$selection[["AIC(n)"]]))
if (is.na(p_aic) || p_aic < 1) p_aic <- 2L
cat(">>> AIC szerinti késleltetés (p):", p_aic, "\n")

#van-e kointegráció?
if (rank_5pct >= 1) {
  cat("\n### Kointegrációt találtunk → SVEC modell (szinteken) ###\n")

  vars_lvls <- colnames(ts_magyar)
 
  LR <- matrix(NA_real_, 2, 2, dimnames = list(vars_lvls, vars_lvls))
  SR <- matrix(NA_real_, 2, 2, dimnames = list(vars_lvls, vars_lvls))
  diag(SR) <- 1

  LR[vars_lvls[2], vars_lvls[1]] <- 0

  svec_model <- SVEC(x = jo, LR = LR, SR = SR, r = 1, lrtest = FALSE)
  print(summary(svec_model))

  plot(irf(svec_model, n.ahead = 10, boot = TRUE, runs = 1000, ci = 0.95))
  print(fevd(svec_model, n.ahead = 10))
  
} else {
  cat("\n### NINCS kointegráció → VAR (Δ) + SVAR (AB) + opcionális BQ ###\n")

  var_model <- VAR(diff_df, p = p_aic, type = "const")
  print(summary(var_model))

  cat("\nRoots |abs|:\n"); print(Mod(roots(var_model)))

  cat("\n--- Diagnosztika ---\n")
  print(serial.test(var_model, lags.pt = 12, type = "PT.asymptotic"))
  print(arch.test(var_model, lags.multi = 5))
  print(normality.test(var_model))

  nm <- colnames(var_model$y)    
  stopifnot(length(nm) == 2)

  var_model$y <- var_model$y[, nm, drop = FALSE]
  
  Amat <- matrix(c(1,    NA,
                   0,     1), nrow = 2, byrow = TRUE)
  Bmat <- diag(2)

  svar_ab <- try(SVAR(var_model, estmethod = "direct",  Amat = Amat, Bmat = Bmat), silent = TRUE)
  if (inherits(svar_ab, "try-error")) {

    svar_ab <- SVAR(var_model, estmethod = "scoring", Amat = Amat, Bmat = Bmat)
  }
}
summary(svar_ab)


plot(irf(svar_ab, n.ahead = 10, boot = TRUE, runs = 1000, ci = 0.95))
fevd_svar <- fevd(svar_ab, n.ahead = 10); fevd_svar

###
#maradek mentese
###


figures_dir <- "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/figures"
excel_path  <- "C:/Users/molna/Documents/KUTATÁSOK/EU_oktatás/ADATOK/model_eval/results_master.xlsx"
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

tidy_test <- function(name, stat, p, df = NA, param = NA) {
  tibble::tibble(test = name, statistic = as.numeric(stat), p_value = as.numeric(p),
                 df = as.character(df), parameter = as.character(param))
}

tidy_coefs_list <- function(cf_list) {

  out <- lapply(names(cf_list), function(eq) {
    tibble::tibble(
      equation = eq,
      term     = names(cf_list[[eq]]),
      estimate = as.numeric(cf_list[[eq]])
    )
  })
  dplyr::bind_rows(out)
}

###
# PSM SAVES
###

psm_logit_coefs  <- tryCatch(broom::tidy(logit), error = function(e) tibble::tibble())
psm_logit_glance <- tryCatch(broom::glance(logit), error = function(e) tibble::tibble())

psm_matchit_summary <- tryCatch({
  s <- summary(matching)
  tibble::tibble(
    n_treated  = s$nn[["Treated"]],
    n_control  = s$nn[["Control"]],
    method     = s$call$method %||% NA_character_,
    ratio      = s$call$ratio %||% NA_real_
  )
}, error = function(e) tibble::tibble())

psm_balance <- tryCatch({
  bt <- cobalt::bal.tab(matching, un = TRUE)
  as_tibble(bt$Balance, rownames = "covariate")
}, error = function(e) tibble::tibble())

psm_ttest <- tryCatch({
  matched <- match.data(matching)
  t.res <- t.test(matched$edu_profile_eqw[matched$treat_female == 1],
                  matched$edu_profile_eqw[matched$treat_female == 0])
  tidy_test("PSM matched outcome t-test",
            t.res$statistic, t.res$p.value,
            df = t.res$parameter, param = "diff of means")
}, error = function(e) tibble::tibble())

try({
  lp <- cobalt::love.plot(matching, stats = "mean.diffs", thresholds = 0.1, abs = TRUE,
                          var.order = "unadjusted", stars = "raw", line = TRUE)
  ggsave_word("psm_loveplot.png", plot = lp, path = figures_dir)
}, silent = TRUE)


###
#var svar mentesek
###

adf_tbl <- tryCatch({
  tibble::tibble(
    series    = c("total_students","graduates"),
    teststat  = c(adf_students@teststat[1], adf_graduates@teststat[1]),
    cval_5pct = c(adf_students@cval[1,"5pct"], adf_graduates@cval[1,"5pct"]),
    reject_5  = teststat < cval_5pct
  )
}, error = function(e) tibble::tibble())

johansen_tbl <- tryCatch({
  crit_cols <- colnames(jo@cval)
  five_col  <- if ("5pct" %in% crit_cols) "5pct" else if ("5%" %in% crit_cols) "5%" else crit_cols[grep("5", crit_cols)[1]]
  tibble::tibble(
    statistic = as.numeric(jo@teststat),
    cval_5pct = as.numeric(jo@cval[, five_col]),
    reject_5  = statistic > cval_5pct
  )
}, error = function(e) tibble::tibble())

var_coefs_tbl <- tryCatch({
  cf <- coef(var_model)              
  tidy_coefs_list(cf) %>% mutate(model = "VAR(diff)", .before = 1)
}, error = function(e) tibble::tibble())

var_diag_tbl <- tryCatch({
  s <- serial.test(var_model, lags.pt = 12, type = "PT.asymptotic")
  a <- arch.test(var_model, lags.multi = 5)
  n <- normality.test(var_model)
  dplyr::bind_rows(
    tidy_test("VAR serial (Portmanteau)", s$statistic, s$p.value),
    tidy_test("VAR ARCH LM (multi)",      a$statistic, a$p.value),
    tidy_test("VAR normality (JB)",       n$jb.stat[,"JB"], n$jb.stat[,"JB p-value"])
  )
}, error = function(e) tibble::tibble())

var_roots_tbl <- tryCatch({
  r <- Mod(roots(var_model))
  tibble::tibble(root_modulus = as.numeric(r),
                 stable = root_modulus < 1)
}, error = function(e) tibble::tibble())

var_irf_tbl <- tryCatch({
  ir <- irf(var_model, n.ahead = 10, boot = TRUE, runs = 1000, ci = 0.95)
  save_base_word("var_irf.png", path = figures_dir, expr = plot(ir))

  comp <- function(arr, name) {

    tibble::tibble()
  }
  tibble::tibble(note = "IRFs saved to PNG (var_irf.png)")
}, error = function(e) tibble::tibble())

var_fevd_tbl <- tryCatch({
  fv <- fevd(var_model, n.ahead = 10)

  out <- lapply(names(fv), function(eq) {
    as.data.frame(fv[[eq]]) |>
      tibble::rownames_to_column("horizon") |>
      tidyr::pivot_longer(-horizon, names_to = "shock", values_to = "share") |>
      mutate(equation = eq, .before = 1)
  })
  dplyr::bind_rows(out)
}, error = function(e) tibble::tibble())

svar_coefs_tbl <- tryCatch({

  Ahat <- svar_ab$A; Bhat <- svar_ab$B
  rbind(
    tibble::as_tibble(as.data.frame(Ahat)) |>
      mutate(row = rownames(Ahat), .before = 1) |>
      tidyr::pivot_longer(-row, names_to = "col", values_to = "estimate") |>
      mutate(matrix = "A"),
    tibble::as_tibble(as.data.frame(Bhat)) |>
      mutate(row = rownames(Bhat), .before = 1) |>
      tidyr::pivot_longer(-row, names_to = "col", values_to = "estimate") |>
      mutate(matrix = "B")
  ) |>
    mutate(model = "SVAR(AB)", .before = 1)
}, error = function(e) tibble::tibble())

svar_irf_fevd_note <- tryCatch({
  ir_svar <- irf(svar_ab, n.ahead = 10, boot = TRUE, runs = 1000, ci = 0.95)
  save_base_word("svar_irf.png", path = figures_dir, expr = plot(ir_svar))
  fv_svar <- fevd(svar_ab, n.ahead = 10)

  save_base_word("svar_fevd.png", path = figures_dir, expr = plot(fv_svar))
  tibble::tibble(note = "SVAR IRF/FEVD plots saved (svar_irf.png, svar_fevd.png)")
}, error = function(e) tibble::tibble())

###
# excelbe
###

if (!exists("sheets") || !is.list(sheets)) sheets <- list()  

sheets$"PSM_logit_coefs"    <- psm_logit_coefs
sheets$"PSM_logit_glance"   <- psm_logit_glance
sheets$"PSM_matchit_summary"<- psm_matchit_summary
sheets$"PSM_balance_table"  <- psm_balance
sheets$"PSM_matched_ttest"  <- psm_ttest

sheets$"VAR_ADF"            <- adf_tbl
sheets$"Johansen_trace"     <- johansen_tbl
sheets$"VAR_coefficients"   <- var_coefs_tbl
sheets$"VAR_diagnostics"    <- var_diag_tbl
sheets$"VAR_roots"          <- var_roots_tbl
sheets$"VAR_FEVD"           <- var_fevd_tbl
sheets$"SVAR_AB_coefs"      <- svar_coefs_tbl
sheets$"SVAR_plots_note"    <- svar_irf_fevd_note

writexl::write_xlsx(sheets, path = excel_path)
message("PSM + VAR/SVAR mentve: táblák Excelben, ábrák a 'figures' mappában.")

