# Compute NPV of net benefits
compute_net_benef_npv <- function(gdp_boot, discount_rates, extrap = "cc") {
    
  gdp_mit <- readd('gdp_mit')
  sr15c_runs <- readd('sr15c_runs')
  
  nben_ann <- compute_net_benefit_ann(gdp_mit, gdp_boot, sr15c_runs, 
                                      disc_years_keep,
                                      extrap)
  
  nben_ann <- nben_ann[!is.na(nbenefit)]
  
  dd <- lapply(discount_rates, 
               function(x) nben_ann[,.(dr = x, 
                                    npv_nbenefit = npvalue(.SD$year, .SD$nbenefit, x, 2020:2300),
                                    npv_gdp_cc_ref = npvalue(.SD$year, .SD$gdp_cc_ref, x, 2020:2300)),
                                    by = c("dmg_spec","runid","c5model","sw","model","scenario","region")])
  dd <- rbindlist(dd)
  
  dd[, nbenefit := npv_nbenefit / npv_gdp_cc_ref]
  
  return(dd)
  
}

# Compute net benefits for the period 2020-2100
compute_net_benef_years <- function(gdp_boot) {
    
  gdp_mit <- readd('gdp_mit')
  sr15c_runs <- readd('sr15c_runs')
  
  nben_ann <- compute_net_benefit_ann(gdp_mit, gdp_boot, sr15c_runs,
                                      c(seq(2020,2100,by = 10)))
  nben_ann[, nbenefit := nbenefit / gdp_cc_ref]

  return(nben_ann)
  
}

# Compute annual net benefits
compute_net_benefit_ann <- function(gdp_mit, gdp_cc, sr15c_runs, years, extrap = "cc") {

  # Add reference scenarios
  scen_ref <- unique(sr15c_runs$scenario_ref)
  gdp_cc_ref <- gdp_cc[scenario %in% scen_ref, 
                       .(dmg_spec,runid,c5model,sw,model,
                         scenario_ref = scenario,
                         region,year,
                         gdp_cc_ref = gdp_cc)]
  
  benef <- merge(gdp_cc,sr15c_runs,by = c("model","scenario"))
  benef <- merge(gdp_cc_ref, benef, 
                 by = c("dmg_spec","runid","c5model","sw","scenario_ref",
                        "model","region","year"))
  
  # Annual avoided damages 
  benef <- benef[year %in% years, 
                 .(admg = gdp_cc - gdp_cc_ref, gdp_cc_ref),
                 by = c("dmg_spec","runid","c5model","sw",
                        "model","scenario","region","year")]
  
  # Extrapolation
  if (max(years) > 2100) {
    extra_years <- years[years %in% 2101:2300]
    # constant extrapolation
    if (str_sub(extrap,1,1) == "c") {
      extra2100 <- benef[year == 2100]
      extra2100 <- extra2100[, .(year = extra_years, 
                                 admg = admg, 
                                 gdp_cc_ref = gdp_cc_ref), 
                             by = "dmg_spec,runid,c5model,sw,model,scenario,region"]
    }
    # linear trend extrapolation
    if (str_sub(extrap,1,1) == "t") {
      extra2100 <- benef[year %in% c(2090,2100)]
      extra2100 <- extra2100[, .(year = extra_years, 
                                 admg = pmax(0, admg[2] + (admg[2] - admg[1]) / 10 * (extra_years - 2100)),
                                 gdp_cc_ref = gdp_cc_ref[2] + (gdp_cc_ref[2] - gdp_cc_ref[1]) / 10 * (extra_years - 2100)), 
                             by = "dmg_spec,runid,c5model,sw,model,scenario,region"]
    }
    benef <- rbindlist(list(benef,extra2100))
  }
  
  setkey(benef,dmg_spec,runid,c5model,model,scenario,region,year)
  
  # Annual mitigation cost [% of GDP REF] 
  cost <- merge(gdp_mit,sr15c_runs,by = c("model","scenario"))
  cost <- cost[gdp_mit > gdp_ref, gdp_mit := gdp_ref]
  cost <- cost[year %in% years,
                  .(cmit = sum(gdp_ref - gdp_mit) / sum(gdp_ref)),
                  by = c("scenario","model","region","year")]
  cost[, year := as.numeric(year)]
  
  # Negative values, as very small, are set to zero
  cost[cmit < 0, cmit := 0]

  # Baseline
  cost <- merge(cost,sr15c_runs,by = c("model","scenario"))
  gdp_cc_ref <- gdp_cc[scenario %chin% unique(cost$scenario_ref)]
  setnames(gdp_cc_ref, "scenario", "scenario_ref")
  cost <- merge(gdp_cc_ref, cost, 
                by = c("scenario_ref","model","region","year"), 
                allow.cartesian = TRUE)
  cost[, cmit := cmit * gdp_cc]
  cost[, c("scenario_ref","gdp_cc","rcp_clus","temp_2100") := NULL]

  # Extrapolation
  if (max(years) > 2100) {
    extra_years <- years[years %in% 2101:2300]
    # constant extrapolation
    if (str_sub(extrap,2,2) == "c") {
      extra2100 <- cost[year == 2100]
      extra2100 <- extra2100[, .(year = extra_years, cmit = cmit), 
                             by = "model,region,dmg_spec,runid,c5model,sw,scenario"]
    }
    # linear trend extrapolation
    if (str_sub(extrap,2,2) == "t") {
      extra2100 <- cost[year %in% c(2090,2100)]
      extra2100 <- extra2100[, .(year = extra_years, 
                                 cmit = pmax(0, cmit[2] + (cmit[2] - cmit[1]) / 10 * (extra_years - 2100))), 
                             by = "model,region,dmg_spec,runid,c5model,sw,scenario"]
    }
    # decrease to 0
    if (str_sub(extrap,2,2) == "0") {
      extra2100 <- cost[year == 2100]
      extra2100 <- extra2100[, .(year = extra_years, 
                                 cmit = pmax(0, cmit + (0 - cmit) / (2200 - 2100) * (extra_years - 2100))), 
                             by = "model,region,dmg_spec,runid,c5model,sw,scenario"]
    }
    cost <- rbindlist(list(cost,extra2100),use.names=TRUE)
  }
    
  setkey(cost,dmg_spec,runid,c5model,model,scenario,region,year)
  
  net_benef <- merge(benef,cost,
                     by = c("dmg_spec","runid","c5model","sw",
                            "model","scenario","region","year"), 
                     all = TRUE)
  
  net_benef[, nbenefit := admg - cmit]

  net_benef <- net_benef[!is.na(nbenefit)]
  
  return(net_benef)
}

# Compute delta CEBGE
compute_diff_cebge <- function(gdp_mit, gdp_cc, years, prtps, extrap = "tc", eta = 1.001) {
  
  loadd(sr15c_runs)
  loadd(ssp_pop)
  
  # Add reference scenarios
  scen_ref <- unique(sr15c_runs$scenario_ref)
  gdp_cc_ref <- gdp_cc[scenario %in% scen_ref & region == "World", 
                       .(dmg_spec,runid,c5model,sw,model,
                         scenario_ref = scenario,
                         region,year,
                         gdp_cc_ref = gdp_cc)]
  map_ssp <- data.table(scenario = unique(gdp_cc$scenario),
                        ssp = sapply(unique(gdp_cc$scenario),scen_ssp, USE.NAMES = F))
  wssp <- ssp_pop[, .(pop = sum(value)) , by = "ssp,year"]
  
  welf <- merge(gdp_cc,sr15c_runs,by = c("model","scenario"))
  welf <- merge(gdp_cc_ref, welf, 
                 by = c("dmg_spec","runid","c5model","sw","scenario_ref",
                        "model","region","year"))
  welf <- merge(map_ssp, welf, by = "scenario")
  welf <- merge(wssp, welf, by = c("ssp","year"))
  
  # Add mitigation costs
  
  # utility
  welf[, u_cc := (gdp_cc / pop * 1e3)^(1 - eta) / (1 - eta)]
  welf[, u_cc_ref := (gdp_cc_ref / pop * 1e3)^(1 - eta) / (1 - eta)]
  
  # remove useless columns
  welf[, c("scenario_ref","gdp_cc_ref","gdp_cc","rcp_clus","temp_2100","pop") := NULL]
  
  # Extrapolation
  if (max(years) > 2100) {
    extra_years <- years[years %in% 2101:2300]
    # constant extrapolation
    if (str_sub(extrap,1,1) == "c") {
      extra2100 <- welf[year == 2100]
      extra2100 <- extra2100[, .(year = extra_years, u_cc = u_cc, u_cc_ref = u_cc_ref), by = "dmg_spec,runid,c5model,sw,model,scenario,region"]
    }
    # linear trend extrapolation
    if (str_sub(extrap,1,1) == "t") {
      extra2100 <- welf[year %in% c(2095,2100)]
      extra2100 <- extra2100[, .(year = extra_years, 
                                 u_cc = u_cc[2] + (u_cc[2] - u_cc[1]) / 5 * (extra_years - 2100),
                                 u_cc_ref = u_cc_ref[2] + (u_cc_ref[2] - u_cc_ref[1]) / 5 * (extra_years - 2100)), 
                             by = "dmg_spec,runid,c5model,sw,model,ssp,scenario,region"]
    }
    welf <- rbindlist(list(welf,extra2100),use.names=TRUE)
  }
  
  # compute welfare
  welf <- merge(wssp, welf, by = c("ssp","year"))
  welf[, lu_cc := pop * 1e6 * u_cc]
  welf[, lu_cc_ref := pop * 1e6 * u_cc_ref]
  
  welf <- lapply(prtps, 
               function(x) welf[,.(prtp = x,
                                       w_cc = npvalue(.SD$year, .SD$lu_cc, x, 2020:2300),
                                       w_cc_ref = npvalue(.SD$year, .SD$lu_cc_ref, x, 2020:2300)),
                                    by = c("dmg_spec","runid","c5model","sw","model","scenario","region")])
  welf <- rbindlist(welf)
  
  # Compute delta CEBGE
  welf[, dcebge := (w_cc / w_cc_ref)^(1/(1 - eta)) - 1]
  
  return(welf)
}
