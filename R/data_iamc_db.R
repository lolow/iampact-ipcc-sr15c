# Read IAMC-format database
# SR1.5C IAMC database

get_sr15c_outliers <- function() {
  dd <- fread("scenario,model,why
                  ADVANCE_INDC,WITCH-GLOBIOM 4.2,Very high npc cmit in R5ASIA 85% + R5OECD90+EU 55%
                  ADVANCE_2030_1.5C-2100,POLES ADVANCE,Very High npv cmit in R5REF 86%
                  ADVANCE_2020_1.5C-2100,POLES ADVANCE,Very High npv cmit in R5REF 49%
                  ADVANCE_2030_Price1.5C,POLES ADVANCE,Very High npv cmit in R5REF 49%")
  setkey(dd, model, scenario)
  return(dd)
}

sr15c_query <- function(sr15c_csv, query_var, sr15c_mo_no, sr15c_scen_pat, sr15c_outliers) {
  sr15c_header <- tolower(unlist(fread(sr15c_csv, nrows = 1), use.names = F))
  sr15c_header_ids <- sr15c_header[1:5]
  sr15c_coltype <- c(rep("character", 5), rep("numeric", length(sr15c_header) - 5))

  fread(str_glue("grep -e ',{query_var},' {sr15c_csv}"),
    col.names = sr15c_header,
    colClasses = sr15c_coltype
  ) %>%
    long_table(ids = sr15c_header_ids) %>%
    setkey(model, scenario, region, variable, unit, year) %>%
    .[!is.na(value)] %>%
    .[!model %chin% sr15c_mo_no] %>%
    .[str_detect(scenario, sr15c_scen_pat)] %>%
    .[!sr15c_outliers[, .(model, scenario)]] %>%
    .[, .(model, scenario, region, year, value)] %>%
    setkey(model, scenario, region)
}

get_scen_info <- function(trf_magicc, temp_magicc, gmt_hist) {
  find_baseline <- function(x) {
    if (str_detect(x, "SSP\\d-(\\d\\d|Baseline)")) {
      return(paste0(str_extract(x, "SSP\\d"), "-Baseline"))
    }
    if (str_detect(x, "SFCM_")) {
      return(str_replace(x, "1p5Degree|2Degree", "Baseline"))
    }
    if (str_detect(x, "^EN_")) {
      return("EN_NPi2100")
    }
    if (str_detect(x, "CD-LINKS_")) {
      return("CD-LINKS_NoPolicy")
    }
    if (str_detect(x, "ADVANCE_")) {
      return("ADVANCE_NoPolicy")
    }
    if (str_detect(x, "EMF33_")) {
      return("EMF33_Baseline")
    }
    if (str_detect(x, "TERL_")) {
      return(str_replace(x, "15D|2D", "Baseline"))
    }
    if (str_detect(x, "Ratchet")) {
      return("Reference")
    }
    if (str_detect(x, "GEA_Mix")) {
      return("GEA_Mix_base")
    }
    if (str_detect(x, "GEA_Eff")) {
      return("GEA_Eff_base")
    }
    if (str_detect(x, "CEMICS-")) {
      return("CEMICS-Ref")
    }
    if (str_detect(x, "SMP_\\w+_Def")) {
      return("SMP_REF_Def")
    }
    if (str_detect(x, "SMP_*")) {
      return("SMP_REF_Sust")
    }
    return(NA)
  }

  hist_magicc <- temp_magicc[year >= 2005 & year <= 2019, .(gmt = mean(value)), by = year]
  hist_cru <- gmt_hist[year >= 2005 & year <= 2019, gmt]
  
  diffgmt <- function(x) {
    return(abs(sum(hist_magicc$gmt - x - hist_cru)))
  }
  search <- optimize(diffgmt, c(0,1)) # minimizing distance between MAGICC and CRU
  d2 <- search$minimum # approx = 0.078
  
  trf_2100 <- trf_magicc[year == 2100]
  temp_2100 <- temp_magicc[year == 2100]
  temp_2100[, value := value - d2] # unbiased temperature

  dd <- trf_magicc[, .(model, scenario)] %>%
    unique() %>%
    .[, scenario_ref := sapply(scenario, find_baseline)]

  rcp_forc <- c(1.9, 2.6, 3.4, 4.5, 6.0, 8.5)
  rcp_name <- c("19", "26", "34", "45", "60", "85")

  clust <- trf_2100[, rcp_clus := sapply(value, function(x) rcp_name[which.min(abs(x - rcp_forc))])][, .(model, scenario, rcp_clus)]

  dd <- merge(dd, clust, by = c("model", "scenario"))

  dd <- merge(dd, temp_2100[, .(model, scenario, temp_2100 = value)], by = c("model", "scenario"))

  return(dd)
}

select_scen_runs <- function(gdp_mit, sr15c_info) {

  # Keep only RCP19,26,34 and Reference scenario
  dd <- sr15c_info[rcp_clus %in% c("19", "26", "34") | scenario == scenario_ref]

  # Keep only models providing mitigation costs
  miti_models <- unique(gdp_mit$model)
  dd <- dd[model %in% miti_models]
  dd <- dd[!scenario %in% c("GEA_Mix_1p5C_AdvNCO2_PartialDelay2020","ADVANCE_2030_Price1.5C")]

  return(dd)
}

# Consolidate GDP mitigation from SR15C database
consolidate_gdp_mit <- function(gdp_mit_mer, gdp_mit_ppp, gdp_mit_mac, gdp_mit_nrg, sr15c_info, sr15c_scen) {

  # From ADVANCE and CD-LINKS take GDP|MER
  gdp1 <- gdp_mit_mer[str_detect(scenario, sr15c_scen)]

  # Compute GDP loss
  # basepol <- sr15c_runs[scenario %chin% gdp1$scenario]
  gdp1_ref <- gdp1[
    scenario %chin% unique(sr15c_info$scenario_ref),
    .(model, scenario_ref = scenario, region, year, value_ref = value)
  ]
  gdp1 <- gdp1[!scenario %chin% unique(sr15c_info$scenario_ref)]
  gdp1 <- merge(gdp1, sr15c_info, by = c("model", "scenario"))
  gdp1 <- merge(gdp1, gdp1_ref, by = c("model", "scenario_ref", "region", "year"))
  gdp1[, gdp_loss := value - value_ref]

  # Detect model which do not provide GDP loss
  mo_nrg <- gdp1[gdp_loss == 0 & !str_detect(scenario, "NoPolicy") & year == 2050, unique(model)]
  gdp1 <- gdp1[!model %chin% mo_nrg]

  # From SSPDB take GDP|PPP (all computation will have to be done relatively then)
  gdp2 <- gdp_mit_ppp[str_detect(scenario, "^SSP\\d*-")]
  gdp2_ref <- gdp2[
    scenario %chin% unique(sr15c_info$scenario_ref),
    .(model, scenario_ref = scenario, region, year, value_ref = value)
  ]
  gdp2 <- gdp2[!scenario %chin% unique(sr15c_info$scenario_ref)]
  gdp2 <- merge(gdp2, sr15c_info, by = c("model", "scenario"))
  gdp2 <- merge(gdp2, gdp2_ref, by = c("model", "scenario_ref", "region", "year"))
  gdp2[, gdp_loss := value - value_ref]

  # Detect model which do not provide GDP loss
  mo_nrg <- gdp2[gdp_loss == 0 & !str_detect(scenario, "NoPolicy") & year == 2050, unique(model)]
  gdp2 <- gdp2[!model %chin% mo_nrg]

  # Harmonize tables
  gdp1 <- gdp1[, .(model, scenario, region, year, gdp_ref = value_ref, gdp_mit = value)]
  # gdp1a <- gdp1a[, .(model, scenario, region, year, gdp_ref, gdp_mit = gdp)]
  gdp2 <- gdp2[, .(model, scenario, region, year, gdp_ref = value_ref, gdp_mit = value)]

  return(rbindlist(list(gdp1, gdp2)))
}
