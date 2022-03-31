# Plot damages

# Plot GDP loss on map in 2100
plot_dmg_2100 <- function(output_file, aggreg = "iso3") {
  
  loadd(sr15c_info)
  loadd(ssp_gdpcap)
  r5_iso3 <- fread('data/r5_iso3_mapping.csv')
  
  dmg_spec_val <- data.table(dmg_spec = c(bhm_spec,djo_spec,"HS2017 NCAT","HS2017 TOT","DICE2016","KW2020"),
                             dmg_spec_label = c(bhm_spec_label,djo_spec_label,"HS NCAT","HS TOT","DICE","KW LVL"),
                             dmg_spec_star = c(bhm_spec_star,djo_spec_star,1,1,1,1))

  sspec <- c('BHM SR',
             'BHM LR ORIG',
             'BHM LR RP',
             'PRETIS2018',
             'DJO 0L ORIG',
             'DJO 0L',
             'DJO 5L ORIG',
             'DJO 5L')
  
  dmg_map <- list()
  
  world <- ne_countries(scale = "small", returnclass = "sf")
  world <- subset(world,!adm0_a3 %in% c("ATA","FJI")) 
  
  world <- merge(world,r5_iso3, by.x = "adm0_a3", by.y = "iso3")
  # Remove small islands
  target_crs <- '+proj=eqearth +wktext'
  world0 <- st_transform(world, crs = target_crs)

  dmg_plot = list()
  
  gdp_gwt <- readd("gdp_gwt", character_only = T)
  
  for (f in sspec) {
    
    
    dd <- gdp_gwt[c5model == 1 & year %in% dmg_diag_years &
              scenario %in% dmg_diag_scenario]
    dd <- dd[,.(gdp_cc = mean(gdp_cc)), by = "scenario,region,year"]
    
    dmg_spec_label <- dmg_spec_val[dmg_spec == gdp_gwt$dmg_spec[1], dmg_spec_label]
    
    map_ssp <- data.table(scenario = unique(dd$scenario),
                    ssp = sapply(unique(dd$scenario),scen_ssp, USE.NAMES = F))
    
    ss <- merge(ssp_gdpcap[year %in% dmg_diag_years], map_ssp, by = "ssp")
    
    
    #R5 regions
    ss0 <- merge(ss,r5_iso3, by = "iso3")
    ss0 <- ss0[,.(gdp_nocc = sum(gdpcap_nocc * pop * 1e-3)), # -> Billions
        by = c("region","scenario","year")]
    #World
    ss1 <- ss[,.(region = "World", gdp_nocc = sum(gdpcap_nocc * pop * 1e-3)), # -> Billions
              by = c("scenario","year")]
    #ISO3
    ss2 <- ss[,.(gdp_nocc = sum(gdpcap_nocc * pop * 1e-3)), # -> Billions
              by = c("iso3","scenario","year")]
    setnames(ss2,"iso3","region")
    
    ss <- rbindlist(list(ss0,ss1,ss2), use.names = TRUE)
    
    dd <- merge(dd,ss, by = c("region","scenario","year"))
    
    dd[, gdp_loss := ((gdp_cc / gdp_nocc) - 1) * 100]
    dd[,dmg_spec := gdp_gwt$dmg_spec[1]]

    dmg_map <- c(dmg_map, list(dd))
    
    dd[gdp_loss > 200, gdp_loss := 200]
    if (aggreg == "iso3") {
      world1 <- merge(world0,dd[scenario == dmg_diag_scenario[1] & year == dmg_diag_years[1], 
                                .(adm0_a3=region,gdp_loss)], 
                      by = "adm0_a3") 
    }
    if (aggreg == "r5") {
      world1 <- merge(world0,dd[scenario == dmg_diag_scenario[1] & year == dmg_diag_years[1], 
                                .(region,gdp_loss)], 
                      by = "region")
      
    }

    p1 <- ggplot(data = world1) +
      geom_sf(aes(fill = gdp_loss)) +
      coord_sf(datum = target_crs,
               expand = FALSE,
               clip = "off") +
      theme_void() + 
        scale_fill_gradient2(limits = c(-100,200), breaks = c(-100,-10,-1,1,10,100), 
                                  name = "[%]",
                             trans = modulus_trans(0.25)) +
      labs(title = paste(dmg_spec_label, dmg_diag_scenario[1]),
           subtitle = paste0("World: ", dd[scenario == dmg_diag_scenario[1] & region == "World", round(gdp_loss,2)], " %"))
    
    
    if (aggreg == "iso3") {
      world2 <- merge(world0,dd[scenario == dmg_diag_scenario[2] & year == dmg_diag_years[1], 
                                .(adm0_a3=region,gdp_loss)], 
                      by = "adm0_a3") 
    }
    if (aggreg == "r5") {
      world2 <- merge(world0,dd[scenario == dmg_diag_scenario[2] & year == dmg_diag_years[1], 
                                .(region,gdp_loss)], 
                      by = "region")
      
    }
    p2 <- ggplot(data = world2) +
      geom_sf(aes(fill = gdp_loss)) +
      coord_sf(datum = target_crs,
               expand = FALSE,
               clip = "off") +
      theme_void() + 
      scale_fill_gradient2(limits = c(-100,200), breaks = c(-100,-10,-1,1,10,100), 
                           name = "[%]",
                           trans = modulus_trans(0.25)) +
      labs(title = paste(dmg_spec_label, dmg_diag_scenario[2]),
           subtitle = paste0("World: ", dd[scenario == dmg_diag_scenario[2] & region == "World", round(gdp_loss,2)], " %"))
    
    dmg_plot <- c(dmg_plot, list(p1), list(p2))
  
  }
  
  dmg_map <- rbindlist(dmg_map)
  
  pp1 <- dmg_plot[[1]] + dmg_plot[[2]] + dmg_plot[[3]] + dmg_plot[[4]] +
    dmg_plot[[5]] + dmg_plot[[6]] + dmg_plot[[7]] + dmg_plot[[8]] + 
    plot_layout(ncol = 2) +
    plot_annotation(
      title = "Economic impacts from climate change in 2100 (GDP loss [%]) - part I")
  
  pp2 <- dmg_plot[[9]] + dmg_plot[[10]] + dmg_plot[[11]] + dmg_plot[[12]] +
    dmg_plot[[13]] + dmg_plot[[14]] + dmg_plot[[15]] + dmg_plot[[16]] + 
    plot_layout(ncol = 2) +
    plot_annotation(
      title = "Economic impacts from climate change in 2100 (GDP loss [%]) - part II")
  
  
  ggsave(output_file, plot = pp1, width = 10, height = 12, dpi = 120)
  ggsave(str_replace(output_file,'_dmg_','_dmg2_'), plot = pp2, width = 10, height = 12, dpi = 120)
  
  return(dmg_map)
  
}
