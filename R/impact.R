# Generic impact functions

read_dmg_param <- function(ssp_gdp,wdi_stat) {
  return(c(read_param_bhm(),
           read_param_djo(),
           read_param_levels(),
           param_rp(ssp_gdp,wdi_stat)))
}

# Project country-level damage function
# Aggregate at region-level
project_cc_boot_impact <- function(dmg_spec, climate_data = "climate_med") {
    
  loadd(dmg_param)
  loadd(ssp_gdpcap) 
  climate <- readd(climate_data, character_only = TRUE)
  loadd(clim_hist)
  
  # Mappings
  r5_iso3 <- fread('data/r5_iso3_mapping.csv')
  
  gdp_cc <- NULL
  
  cmodels00 <- unique(climate$c5model)
  sw_c <- unique(climate$sw)
  
  for (i in cmodels00) {
    
    for (j in sw_c) {
      
      if (dmg_spec %in% c(bhm_spec_boot)) {
        dd <- project_cc_impact_bhm(ssp_gdpcap,climate[c5model == i & sw == j],
                                    clim_hist,dmg_param,dmg_spec)
      }
  
      # Regional aggregation
      map_ssp <- data.table(scenario = unique(dd$scenario),
                            ssp = sapply(unique(dd$scenario),scen_ssp, USE.NAMES = F))
      dd <- merge(dd,map_ssp, by = "scenario")
      dd <- merge(dd,ssp_gdpcap[,.(ssp,iso3,year,pop)],
                      by = c("ssp","iso3","year"))
      dd <- merge(dd,r5_iso3, by = "iso3")
      
      # World
      dd0 <- dd[,.(gdp_cc = sum(gdpcap_cc * pop * 1e-3)), # -> Billions (1e6 * 1e-9)
                by = c("dmg_spec","runid","c5model","model","scenario","year")]
      dd0[,region := "World"]
      
      # R5 regions
      dd1 <- dd[,.(gdp_cc = sum(gdpcap_cc * pop * 1e-3)), # -> Billions
                by = c("dmg_spec","runid","c5model","model","scenario","region","year")]
      
      dd <- rbindlist(list(dd1,dd0),use.names = TRUE)
      
      dd[, sw := j]
      
      gdp_cc <- rbind(gdp_cc,dd)
      
      climate <- climate[!(c5model == i & sw == j)]
      
      #Clean memory
      rm(dd0,dd1)
      gc()
    
    }

  }

  return(gdp_cc)
}

# Project country-level damage function
# Aggregate at region-level
project_cc_gwt_impact <- function(dmg_spec, climate_data = "climate_med") {
  
  loadd(ssp_gdpcap,clim_hist)
  loadd(dmg_param)
  climate = readd(climate_data, character_only = TRUE)
  
  # Mappings
  r5_iso3 <- fread('data/r5_iso3_mapping.csv')
  
  gdp_cc <- NULL
  
  cmodels <- unique(climate$c5model)
  sw_c <- unique(climate$sw)
  
  for (i in cmodels) {
    
    for (j in sw_c) {
    
    if (dmg_spec %in% c(bhm_spec)) {
      dd <- project_cc_impact_bhm(ssp_gdpcap,climate[c5model == i & sw == j],
                                  clim_hist,dmg_param,dmg_spec)
    }
    
    if (dmg_spec %in% djo_spec) {
      dd <- project_cc_impact_djo(ssp_gdpcap,climate[c5model == i & sw == j],
                                  clim_hist,dmg_param,dmg_spec)
    }
    
    # Regional aggregation
    map_ssp <- data.table(scenario = unique(dd$scenario),
                          ssp = sapply(unique(dd$scenario),scen_ssp, USE.NAMES = F))
    dd <- merge(dd,map_ssp, by = "scenario")
    dd <- merge(dd,ssp_gdpcap[,.(ssp,iso3,year,pop)],
                by = c("ssp","iso3","year"))
    
    # ISO3 for c5model = CMIP5 = 1 and a selection of dmg_spec and scenario
    if (i == 1) {
      dd0 <- dd[scenario %in% dmg_diag_scenario & year %in% dmg_diag_years,
                .(gdp_cc = sum(gdpcap_cc * pop * 1e-3)), # -> Billions (1e6 * 1e-9)
                by = c("dmg_spec","runid","c5model","model","scenario","iso3","year")]
      setnames(dd0,"iso3","region")
    } else {
      dd0 <- NULL
    }
    
    dd <- merge(dd,r5_iso3, by = "iso3")
    
    # World
    dd1 <- dd[,.(gdp_cc = sum(gdpcap_cc * pop * 1e-3)), # -> Billions (1e6 * 1e-9)
              by = c("dmg_spec","runid","c5model","model","scenario","year")]
    dd1[,region := "World"]
    
    # R5 regions
    dd2 <- dd[,.(gdp_cc = sum(gdpcap_cc * pop * 1e-3)), # -> Billions
              by = c("dmg_spec","runid","c5model","model","scenario","region","year")]
    
    dd <- rbindlist(list(dd1,dd2,dd0),use.names = TRUE)
    
    dd[, sw := j]
    
    gdp_cc <- rbind(gdp_cc,dd)
    
    climate <- climate[!(c5model == i & sw == j)]
    
    #Clean memory
    rm(dd0,dd1,dd2)
    gc()
    
    }
    
  }
  
  return(gdp_cc)
}

# Rich-Poor Threshold
param_rp <- function(ssp_gdp,wdi_stat) {
  
  # Get iso3 list from SSP
  valid_iso3 <- unique(ssp_gdp$iso3)
  
  # Compute median from WDI statistics (USD2018)
  # The separating value y∗ is the median per capita GDP in 2010, 
  # i.e. at the end of the historical period in BHM2015.
  rp_threshold <- mean(wdi_stat[iso3 %in% valid_iso3 & year %in% 2010, 
                                median(gdp / pop), 
                                by = "year"]$V1)
  
  return(list(rp = rp_threshold))
  
}