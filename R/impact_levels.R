# DICE-2016 
# HS-2017
# KW-2020
# damage functions with level quadratic function

read_param_levels <- function() {
  
  p_list <- list()
  # DICE2016 
  p_list = c(p_list, list(DICE2016 = c(alpha1 = 0, alpha2 = -0.00236)))
  # HS-2017 (preferred specification * 1.25 to include missing impacts = non-catastrophic) = -0.007438
  p_list = c(p_list, list(`HS2017 NCAT` = c(alpha1 = 0, alpha2 = -0.59504 * 1.25 * 1e-2)))
  # HS-2017 (Only market) = -0.001218 (16%)
  p_list = c(p_list, list(`HS2017 MKT` = c(alpha1 = 0, alpha2 = -(0.59504 * 1.25 - 0.622) * 1e-2)))
  # HS-2017 (preferred specification * 1.25 to include missing impacts = total damages) = -0.010038
  p_list = c(p_list, list(`HS2017 TOT` = c(alpha1 = 0, alpha2 = -(0.59504 * 1.25 + 0.260) * 1e-2)))
  # HS-2017 (preferred specification * 1.25 to include missing impacts = total damages + productivity) = -0.01145
  p_list = c(p_list, list(`HS2017 TOTP` = c(alpha1 = 0, alpha2 = -(0.59504 * 1.25 + 0.260 + 0.113 * 1.25) * 1e-2)))
  # KW-2020 (Panel-based (population weighted) global damages)
  p_list = c(p_list, list(KW2020 = c(alpha1 = -0.0357, alpha2 = -0.0018 / 2)))
  # TAKAKURA2019 (SSP-based fitting) rows = ssp{1..5} cols = a x + b x + c
  p_list = c(p_list, list(TAKAKURA2019 = matrix(c( 0.09169, 0.20563, -0.12881,
                                                   0.07625, 0.21465, -0.11746,
                                                  -0.32617, 0.33814,  0.19575,
                                                   0.05442, 0.24189, -0.13341,
                                                   0.14312, 0.17533, -0.15604),
                                               nrow = 5,  byrow = T)))
  
  return(p_list)
  
}

project_cc_impact_levels <- function() {
  
  loadd(ssp_gdpcap,
        temp_magicc,
        temp_magicc_p05,
        temp_magicc_p95,
        gmt_hist,
        dmg_param)
  
  # split scenarios according to SSP
  scen <- data.table(scenario = unique(temp_magicc$scenario))
  scen[, ssp := sapply(scenario,scen_ssp)]
  
  # SSP GWP
  ssp_gwp <- ssp_gdpcap[,.(gdp_ref = sum(gdpcap_nocc * pop) * 1e-3), by = c("ssp","year")]

  # Year should be numbers
  temp_magicc[, year := as.numeric(year)]
  temp_magicc_p05[, year := as.numeric(year)]
  temp_magicc_p95[, year := as.numeric(year)]
  
  # Unbiasied global mean temperature
  hist_magicc <- temp_magicc[year >= 2005 & year <= 2019, .(gmt = mean(value)), by = year]
  hist_cru <- gmt_hist[year >= 2005 & year <= 2019, gmt]
  
  diffgmt <- function(x) {
    return(abs(sum(hist_magicc$gmt - x - hist_cru)))
  }
  search <- optimize(diffgmt, c(0,1)) # minimizing distance between MAGICC and CRU
  d2 <- search$minimum # approx = 0.078
  
  # Year should be numbers
  temp_magicc[, value := value - d2]
  temp_magicc_p05[, value := value - d2]
  temp_magicc_p95[, value := value - d2]
  
  # Restrict to keep_years
  temp_magicc <- temp_magicc[year %in% keep_years]
  temp_magicc_p05 <- temp_magicc_p05[year %in% keep_years]
  temp_magicc_p95 <- temp_magicc_p95[year %in% keep_years]
  
  
  #  Compute GWP with CC for all damages functions
  proj_gdp <- NULL
  for (dmgf in c("DICE2016","HS2017 NCAT","HS2017 TOT","HS2017 TOTP","HS2017 MKT","KW2020")) {
    
    #median
    dd <- merge(temp_magicc,scen,by = "scenario")
    dd <- merge(dd,ssp_gwp,by = c("ssp","year"))
    dd[, gdp_cc := gdp_ref * (1 + (dmg_param[[dmgf]][1] * value + 
                                   dmg_param[[dmgf]][2] * value * value))]
    dd[, dmg_spec := dmgf]
    dd[, c5model := 1]
    dd[, sw := 1]
    proj_gdp <- rbindlist(list(proj_gdp,
                               dd[, .(dmg_spec, c5model, sw, model, scenario, year, gdp_cc)]))

    #p05
    dd <- merge(temp_magicc_p05,scen,by = "scenario")
    dd <- merge(dd,ssp_gwp,by = c("ssp","year"))
    dd[, gdp_cc := gdp_ref * (1 + (dmg_param[[dmgf]][1] * value + 
                                     dmg_param[[dmgf]][2] * value * value))]
    dd[, dmg_spec := dmgf]
    dd[, c5model := 1]
    dd[, sw := 2]
    proj_gdp <- rbindlist(list(proj_gdp,
                               dd[, .(dmg_spec, c5model, sw, model, scenario, year, gdp_cc)]))
    
    #p95
    dd <- merge(temp_magicc_p95,scen,by = "scenario")
    dd <- merge(dd,ssp_gwp,by = c("ssp","year"))
    dd[, gdp_cc := gdp_ref * (1 + (dmg_param[[dmgf]][1] * value + 
                                     dmg_param[[dmgf]][2] * value * value))]
    dd[, dmg_spec := dmgf]
    dd[, c5model := 1]
    dd[, sw := 3]
    proj_gdp <- rbindlist(list(proj_gdp,
                               dd[, .(dmg_spec, c5model, sw, model, scenario, year, gdp_cc)]))
    
  }
  
  for (dmgf in c("TAKAKURA2019")) {
    
    #median
    dd <- merge(temp_magicc,scen,by = "scenario")
    dd <- merge(dd,ssp_gwp,by = c("ssp","year"))
    dd[, ssp_num := as.numeric(str_extract(ssp,"\\d"))]
    dd[, gdp_cc := gdp_ref * (1 - (dmg_param[[dmgf]][ssp_num,1] * value + 
                                     dmg_param[[dmgf]][ssp_num,2] * value * value +
                                     dmg_param[[dmgf]][ssp_num,3]) / 100)]
    dd[, dmg_spec := dmgf]
    dd[, c5model := 1]
    dd[, sw := 1]
    proj_gdp <- rbindlist(list(proj_gdp,
                               dd[, .(dmg_spec, c5model, sw, model, scenario, year, gdp_cc)]))
    
    #p05
    dd <- merge(temp_magicc_p05,scen,by = "scenario")
    dd <- merge(dd,ssp_gwp,by = c("ssp","year"))
    dd[, ssp_num := as.numeric(str_extract(ssp,"\\d"))]
    dd[, gdp_cc := gdp_ref * (1 - (dmg_param[[dmgf]][ssp_num,1] * value + 
                                     dmg_param[[dmgf]][ssp_num,2] * value * value +
                                     dmg_param[[dmgf]][ssp_num,3]) / 100)]
    dd[, dmg_spec := dmgf]
    dd[, c5model := 1]
    dd[, sw := 2]
    proj_gdp <- rbindlist(list(proj_gdp,
                               dd[, .(dmg_spec, c5model, sw, model, scenario, year, gdp_cc)]))
    
    #p95
    dd <- merge(temp_magicc_p95,scen,by = "scenario")
    dd <- merge(dd,ssp_gwp,by = c("ssp","year"))
    dd[, ssp_num := as.numeric(str_extract(ssp,"\\d"))]
    dd[, gdp_cc := gdp_ref * (1 - (dmg_param[[dmgf]][ssp_num,1] * value + 
                                     dmg_param[[dmgf]][ssp_num,2] * value * value +
                                     dmg_param[[dmgf]][ssp_num,3]) / 100)]
    dd[, dmg_spec := dmgf]
    dd[, c5model := 1]
    dd[, sw := 3]
    proj_gdp <- rbindlist(list(proj_gdp,
                               dd[, .(dmg_spec, c5model, sw, model, scenario, year, gdp_cc)]))
    
  }
  
  proj_gdp[, runid := 0]
  proj_gdp[, region := "World"]
  
  return(proj_gdp)
  
}
