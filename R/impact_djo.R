# DJO damage function

# Estimates for DJO based on (DJO, 2012)
read_param_djo <- function() {
  
  # param: c(rich coeff, poor coeff)
  # Non-significant coefficients are set to zero, refer to Table
  
  p_list <- list()
  # 0-lag specification (coeff from Table 3, col 1)
  p_list = c(p_list, list(`DJO 0L` = c(0 / 100, -1.394 / 100)))
  p_list = c(p_list, list(`DJO 0L ORIG` = c(0.261 / 100, -1.394 / 100)))
  # 5-lag specification (coeff from Table 3, col 3)
  p_list = c(p_list, list(`DJO 5L` = c(0 / 100, -1.235 / 100)))
  p_list = c(p_list, list(`DJO 5L ORIG` = c(-0.180 / 100, -1.235 / 100)))
  # 10-lag specification (coeff from Table 3, col 4)
  p_list = c(p_list, list(`DJO 10L` = c(0 / 100, -1.171 / 100)))
  p_list = c(p_list, list(`DJO 10L ORIG` = c(-0.152 / 100, -1.171 / 100)))
  
  return(p_list)
}

warming_effect_djo <- function(temp, temp_baseline, gdpcap_tm1, param, y_star){
  if (gdpcap_tm1 > y_star) {
    return(param[1] * (temp - temp_baseline))
  } else {
    return(param[2] * (temp - temp_baseline))
  }
}

project_djo <- function(tas, gdpr, gdpc_2010, tas_base, param, y_star){
  stopifnot(length(tas) == 19) # 2010-2100 [5-year]
  stopifnot(length(gdpr) == 91) # 2010-2100

  .gdpcap <- rep(gdpc_2010,91)
  idx_tas <- ceiling((1:91) / 5)
  idx_tas1 <- pmin(19,ceiling((1:91) / 5) + 1)
  for (i in 2:91) {
    tas_i <- tas[idx_tas[i]] + ((i - 1) %% 5) * (tas[idx_tas1[i]] - tas[idx_tas[i]]) / 5
    .delta <- warming_effect_djo(tas_i - tas[1] + tas_base, tas_base, 
                                    .gdpcap[i - 1], param, y_star)
    .gdpcap[i] <- .gdpcap[i-1] * (1 + gdpr[i] + .delta)
  }
  
  return(list(year = keep_years, 
              gdpcap_cc = .gdpcap[keep_years_idx]))
}


project_cc_impact_djo <- function(ssp_gdpcap,climate,clim_hist,dmg_param,spec) {
  
  bhm_dta <- dmg_param[['BHM DATASET']]
  bhm_baseline <- bhm_dta[year >= 2000, .(tas = mean(temp)), 
                          by = "iso3"]
  
  # split scenarios according to SSP
  scen <- data.table(scenario = unique(climate$scenario))
  
  scen[, ssp := sapply(scenario,scen_ssp)]

  # check iso3 list
  all_iso3 <- intersect(unique(climate$iso3),
                        unique(ssp_gdpcap$iso3))
  all_ssp <- unique(scen$ssp)
  all_comb <- expand.grid(all_iso3,all_ssp,stringsAsFactors = F)
  
  # add missing country in bhm baseline
  miss <- clim_hist[iso3 %in% all_iso3[!all_iso3 %in% unique(bhm_dta$iso3)] &
                      year >= 2000,
                    .(tas = mean(tas)),
                    by = "iso3"]
  bhm_baseline <- rbind(bhm_baseline,miss)

  proj_gdpcap <- NULL
  for (i in 1:nrow(all_comb)) {
                             
     i_iso3 <- all_comb[i,1]    
     i_ssp <- all_comb[i,2]
     
     .gdpr <- ssp_gdpcap[ssp == i_ssp & iso3 == i_iso3, gdpr]
     .bclim <- bhm_baseline[iso3 == i_iso3, tas]
     .gdpc_2010 <- ssp_gdpcap[ssp == i_ssp & iso3 == i_iso3 & year == 2010, gdpcap_nocc]
     
     xx <- climate[scenario %chin% scen[ssp == i_ssp, scenario] & iso3 == i_iso3, 
             project_djo(tas, .gdpr, .gdpc_2010, 
                         .bclim, 
                         dmg_param[[spec]], 
                         dmg_param[['rp']]),
               by = c("c5model,model,scenario,iso3")]
     
     proj_gdpcap <- rbindlist(list(proj_gdpcap,xx))
                           
  }
  
  proj_gdpcap[, dmg_spec := spec]
  proj_gdpcap[, runid := 1]
  
  return(proj_gdpcap)
  
}