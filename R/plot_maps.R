
filter_boot_mid_end_century <- function(nben_ann) {
  return(nben_ann[year %in% c(2050,2100) &
                    c5model != "CMIP5" &
                    dmg_spec %chin% "BHM SR BOOT"])
}


# Map plots
# Bootstrap run
# Temp cluster
# CMIP5 models
# quantity = {nbenefit_value, nbenefit_share}
map_nben <- function(year0, temp, quantity = "nbenefit_value", show_legend = TRUE, output_file) {
  
  loadd(sr15c_info)
  
  dd <- readd(nben_boot_plot)
  
  # Select scenarios and add info
  scen_ref <- unique(sr15c_info$scenario_ref)
  
  nben <- dd[!scenario %in% scen_ref &
               !is.na(nbenefit) &
               region %in% r5_regions &
               year %in% c(2050,2100) &
               c5model == 1 &
               dmg_spec %in% "BHM SR BOOT" ]
  
  nben <- merge(nben, sr15c_info, by = c("model","scenario"))
  
  nben <- nben[temp_2100 <= 2]
  nben[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nben[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nben[, temp_2100 := NULL]
  
  
  # Median over scenario
  nbensum <- nben[,.(nbenefit_value = quantile(nbenefit,0.5),
                     nbenefit_share = sum(ifelse(nbenefit > 0, 1, 0)) / .N * 100),
                     by = "year,region,temp_clus"]

  # Compare clusters
  dstat <- dcast(nbensum, year + region ~ temp_clus, value.var = quantity)

  r5_iso3 <- fread('data/r5_iso3_mapping.csv')
  
  # nbenefits
  if (quantity == "nbenefit_value") {
    lim0 <- min(nbensum$nbenefit_value)
    lim1 <- max(nbensum$nbenefit_value)
  }
  
  # nbenefits share
  if (quantity == "nbenefit_share") {
    lim0 <- 0
    lim1 <- 100
  }
  
  world <- ne_countries(scale = "small", returnclass = "sf")
  world <- subset(world,!adm0_a3 %in% c("ATA","FJI")) 
  
  world <- merge(world,r5_iso3, by.x = "adm0_a3", by.y = "iso3")
  world <- merge(world,dstat[year == year0, .(region,nben = get(temp))], by = "region")
  
  if (quantity == "nbenefit_value") {
    world$nben <- world$nben * 1e2 #T$
  }
  
  # Remove small islands
  target_crs <- '+proj=eqearth +wktext'
  world0 <- st_transform(world, crs = target_crs)
  
  #world1 <- ms_filter_islands(world0, min_area = 1.5e10)
  
  p <- ggplot(data = world0) +
    geom_sf(aes(fill = nben)) +
    coord_sf(datum = target_crs,
             expand = FALSE,
             clip = "off") +
    theme_void()
  
  if (quantity == "nbenefit_value") {
    p <- p + scale_fill_gradient2(limits = c(lim0,lim1) * 1e2, 
                                  name = "[%]",
                                  guide = ifelse(show_legend,"colourbar","none"))
  }
  
  if (quantity == "nbenefit_share") {
    p <- p + scale_fill_viridis_c(option = "viridis",
                                  name = "[%]",
                                  guide = ifelse(show_legend,"colourbar","none"),
                                  limits = c(lim0,lim1)) 
  }
  
  if (show_legend) {
    p <- p + theme(legend.position = "bottom",
                   legend.key.width = unit(2, "cm"),
                   legend.text = element_text(size = 16),
                   legend.title = element_text(size = 18),
                   legend.direction = "horizontal")
  }
  
  ggsave(output_file, plot = p, width = 6, height = 4, dpi = 120)
  
  data.table(file = output_file)


}


