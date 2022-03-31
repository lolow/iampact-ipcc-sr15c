
fig_extrap_admg <- function(output_file){

  gdp_cc <- readd(gdp_gwt_all_BHM.SR)
loadd(gdp_mit)
loadd(sr15c_runs)
years <- 2020:2300

gdp_cc <- gdp_cc[region == "World" & c5model == 1]

all_benef <- NULL
for ( extrap in c("cc","tc")) {
  
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
                 .(admg = gdp_cc - gdp_cc_ref),
                 by = c("dmg_spec","runid","c5model","sw",
                        "model","scenario","region","year")]
  
  # Extrapolation
  if (max(years) > 2100) {
    extra_years <- years[years %in% 2101:2300]
    # constant extrapolation
    if (str_sub(extrap,1,1) == "c") {
      extra2100 <- benef[year == 2100]
      extra2100 <- extra2100[, .(year = extra_years, admg = admg), by = "dmg_spec,runid,c5model,sw,model,scenario,region"]
    }
    # linear trend extrapolation
    if (str_sub(extrap,1,1) == "t") {
      extra2100 <- benef[year %in% c(2095,2100)]
      extra2100 <- extra2100[, .(year = extra_years, 
                                 admg = pmax(0, admg[2] + (admg[2] - admg[1]) / 5 * (extra_years - 2100))), 
                             by = "dmg_spec,runid,c5model,sw,model,scenario,region"]
    }
    benef <- rbindlist(list(benef,extra2100))
  }
  
  if (extrap=="cc") {benef[, extrap := "Constant"]}
  if (extrap=="tc") {benef[, extrap := "Trend"]}
  
  all_benef <- rbindlist(list(all_benef, benef))

}

p <- ggplot() +
  geom_line(aes(x = year, y = admg * 1e-3, group = paste(model, scenario)), 
            data = all_benef[year <= 2100]) +
  geom_line(aes(x = year, y = admg * 1e-3, color = extrap, group = paste(model, scenario,extrap)), 
            data = all_benef[year >= 2100], alpha = 0.5) +
  scale_color_d3(name = "Extrapolation\nScenario", guide = "none") +
  #guides(color = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(~ extrap) +
  theme_bw() +
  labs(x = "", y = "[T$]")

ggsave(output_file, plot = p, width = 8, height = 3)

}

fig_extrap_cmit <- function(output_file){
  
  gdp_cc <- readd(gdp_gwt_all_BHM.SR)
  loadd(gdp_mit)
  loadd(sr15c_runs)
  years <- 2020:2300
  
  gdp_cc <- gdp_cc[region == "World" & c5model == 1]
  
  all_cmit <- NULL
  for ( extrap in c("cc","ct","c0")) {
    
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
    
    if (extrap=="cc") {cost[, extrap := "Constant"]}
    if (extrap=="ct") {cost[, extrap := "Trend"]}
    if (extrap=="c0") {cost[, extrap := "Decreasing"]}
    
    all_cmit <- rbindlist(list(all_cmit, cost))
    
  }
  
  p <- ggplot() +
    geom_line(aes(x = year, y = cmit * 1e-3, group = paste(model, scenario)), 
              data = all_cmit[year <= 2100]) +
    geom_line(aes(x = year, y = cmit * 1e-3, color = extrap, group = paste(model, scenario,extrap)), 
              data = all_cmit[year >= 2100], alpha = 0.5) +
    scale_color_d3(name = "Extrapolation\nScenario", guide = "none") +
    #guides(color = guide_legend(override.aes = list(alpha = 1))) +
    facet_wrap(~ extrap) +
    theme_bw() +
    labs(x = "", y = "[T$]")
  
  ggsave(output_file, plot = p, width = 8, height = 3)
  
}
