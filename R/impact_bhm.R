# Burke et al. damage function

# Load coefficient for estimates and bootstrap
read_param_bhm <- function() {
  
  mb_list <- list()
  
  # Original BHM2015 dataset
  dta <- setDT(read.dta13('data/BDD2018/data/input/GrowthClimateDataset_Stata13.dta'))
  dta <- dta[, .(iso3 = iso, year, pop = Pop, gdp = TotGDP * 1.12538, #(USD2010->USD2018) 
                growth = growthWDI,temp = UDel_temp_popweight)]
  dta <- dta[!iso3 %in% c('COD', 'ROU')]
  dta <- dta[iso3 == 'ZAR', iso3 := 'COD']
  dta <- dta[iso3 == 'ROM', iso3 := 'ROU']
  dta <- na.omit(dta)
  mb_list = c(mb_list, list(`BHM DATASET` = dta))
  
  # Non-significant coefficients are set to zero for main estimates
  # (ref. replication code)
  # Pooled bootstrap_noLag 
  pb <- read_fst('data/BDD2018/data/output/bhm_sr_param.fst', as.data.table = T)
  mb <- as.matrix(pb[,.(b1,b2,optimal,normalize)])
  mb_list = c(mb_list, list(`BHM SR` = mb))
  # Rich/Poor noLag specification
  pb <- fread("data/BurkeHsiangMiguel2015_Replication/data/output/bootstrap/bootstrap_richpoor.csv")
  mb <- as.matrix(pb[,.(temp,temppoor,temp2,temp2poor)])
  mb_list = c(mb_list, list(`BHM SR RP` = mb))
  # Pooled bootstrap_5Lag
  pb <- read_fst('data/BDD2018/data/output/bhm_lr_param.fst', as.data.table = T)
  mb <- as.matrix(pb[,.(b1,b2,optimal,normalize)])
  mb[1,2] <- 0 # Quadratic term not significant
  mb_list = c(mb_list, list(`BHM LR` = mb))
  pb <- read_fst('data/BDD2018/data/output/bhm_lr_param.fst', as.data.table = T)
  mb <- as.matrix(pb[,.(b1,b2,optimal,normalize)])
  mb_list = c(mb_list, list(`BHM LR ORIG` = mb))
  # Rich/Poor 5-lag specification 
  pb <- fread("data/BurkeHsiangMiguel2015_Replication/data/output/bootstrap/bootstrap_richpoor_5lag.csv")
  mb <- as.matrix(pb[,.(tlin,tlinpoor,tsq,tsqpoor)])
  mb[1,4] <- 0 # Quadratic term not significant
  mb[1,c(1,3)] <- 0 # Rich coefficients not significant
  mb_list = c(mb_list, list(`BHM LR RP` = mb))
  pb <- fread("data/BurkeHsiangMiguel2015_Replication/data/output/bootstrap/bootstrap_richpoor_5lag.csv")
  mb <- as.matrix(pb[,.(tlin,tlinpoor,tsq,tsqpoor)])
  mb_list = c(mb_list, list(`BHM LR RP ORIG` = mb))
  # KW2020 (panel regression Tab4 (2) - reproducing BHM SR, but at higher resolution)
  mb_list = c(mb_list, list(`KW2020 BHM` = matrix(c(0.00947,-7.09e-04),nrow = 1)))
  # BT2019 (Tab S3 panel A Column (1-2), Baseline specification)
  # BT2019 STATE - linear not significant
  mb_list = c(mb_list, list(`BT2019 STATE` = matrix(c(0.0033392,-3.013e-04),nrow = 1)))
  mb_list = c(mb_list, list(`BT2019 CNTRY` = matrix(c(0.0108371,-5.752e-04),nrow = 1)))
  # Pretis2018 (Tab 1 Specification M2)
  # Baseline is 2006-2014 in paper
  mb_list = c(mb_list, list(`PRETIS2018` = matrix(c(0.0115,-4e-04),nrow = 1)))
  # Newell2018 (Tab 4 Specification DJO*+T2)
  mb_list = c(mb_list, list(`NEWELL2018 DJO T2` = matrix(c(0.013219,-3.85e-04),nrow = 1)))
  # HS2019 (Tab 4 Column 1)
  mb_list = c(mb_list, list(`HS2019` = matrix(c(0.00928,-4.535e-04),nrow = 1)))
  # ACEVODO2020 (Tab 1 Column 5)
  mb_list = c(mb_list, list(`ACEVODO2020` = matrix(c(1.347 / 100,-0.051 / 100),nrow = 1)))
  
  return(mb_list)
}

g_rich <- function(par, temp) { return(par[1] * temp + par[3] * temp^2) }
g_poor <- function(par, temp) { return(par[2] * temp + par[4] * temp^2) }

g_pool <- function(par, temp) { return(par[1] * temp + par[2] * temp^2) }

# rid = 1 estimates
# rid = 2:1000 bootstrap coefficients
warming_effect <- function(temp, temp_baseline, param){
    return(g_pool(param,temp) - g_pool(param,temp_baseline))
}

project_bhm <- function(tas, gdpr, gdpc_2010, tas_base, param){
  stopifnot(length(tas) == 19) # 2010-2100 [5-year]
  stopifnot(length(gdpr) == 91) # 2010-2100

  .gdpcap <- rep(gdpc_2010,91)
  idx_tas <- ceiling((1:91) / 5)
  idx_tas1 <- pmin(19,ceiling((1:91) / 5) + 1)
  for (i in 2:91) {
    tas_i <- tas[idx_tas[i]] + ((i - 1) %% 5) * (tas[idx_tas1[i]] - tas[idx_tas[i]]) / 5
    .delta <- warming_effect(tas_i - tas[1] + tas_base, tas_base, param)
    .gdpcap[i] <- .gdpcap[i - 1] * (1 + gdpr[i] + .delta)
  }

  return(list(year = keep_years, 
              gdpcap_cc = .gdpcap[keep_years_idx]))
}

warming_effect_rp <- function(temp, temp_baseline, gdpcap_tm1, param, y_star){
  if (gdpcap_tm1 > y_star) {
    return(g_rich(param,temp) - g_rich(param,temp_baseline))
  } else {
    return(g_poor(param,temp) - g_poor(param,temp_baseline))
  }
}

project_bhm_richpoor <- function(tas, gdpr, gdpc_2010, tas_base, param, y_star){
  stopifnot(length(tas) == 19) # 2010-2100 [5-year]
  stopifnot(length(gdpr) == 91) # 2010-2100

  .gdpcap <- rep(gdpc_2010,91)
  idx_tas <- ceiling((1:91) / 5)
  idx_tas1 <- pmin(19,ceiling((1:91) / 5) + 1)
  for (i in 2:21) {
    tas_i <- tas[idx_tas[i]] + ((i - 1) %% 5) * (tas[idx_tas1[i]] - tas[idx_tas[i]]) / 5
    .delta <- warming_effect_rp(tas_i - tas[1] + tas_base, tas_base, .gdpcap[i-1], param, y_star)
    .gdpcap[i] <- .gdpcap[i-1] * (1 + gdpr[i] + .delta)
  }
  
  return(list(year = keep_years, 
              gdpcap_cc = .gdpcap[keep_years_idx]))
}

project_cc_impact_bhm <- function(ssp_gdpcap,climate,clim_hist,dmg_param,spec) {
  
  runid <- 1
  richpoor <- FALSE
  if (str_detect(spec,"RP")) {
    richpoor <- TRUE
  }
  if (str_detect(spec,"BOOT")) {
    runid <- eval(parse(text = str_sub(spec,str_locate(spec,"BOOT")[1,"end"] + 1))) + 1
    par <- dmg_param[[str_trim(str_replace(str_extract(spec,"(\\s|\\w)*BOOT"),"BOOT",""))]]
    spec <- str_trim(str_extract(spec,"(\\s|\\w)*BOOT"))
  } else {
    runid <- 1
    par <- dmg_param[[spec]]
  }
  bhm_dta <- dmg_param[['BHM DATASET']]
  bhm_baseline <- bhm_dta[year >= 2000, .(tas = mean(temp)), 
                          by = "iso3"]

  # split scenarios according to SSP
  scen <- data.table(scenario = unique(climate$scenario))
  
  scen[, ssp := sapply(scenario,scen_ssp)]

  # check iso3 list
  all_iso3 <- intersect(unique(climate$iso3), unique(ssp_gdpcap$iso3))
  all_ssp <- unique(scen[,ssp])
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

    if (richpoor) {
      xx <- climate[scenario %in% scen[ssp == i_ssp, scenario] & iso3 == i_iso3, 
              project_bhm_richpoor(tas, .gdpr, .gdpc_2010, .bclim, par[runid,], dmg_param[['rp']]),
              by = c("c5model,model,scenario,iso3")]
    } else {
      xx <- climate[scenario %in% scen[ssp == i_ssp, scenario] & iso3 == i_iso3,
                project_bhm(tas, .gdpr, .gdpc_2010, .bclim, par[runid,]),
                by = c("c5model,model,scenario,iso3")]
    }
    
    proj_gdpcap <- rbindlist(list(proj_gdpcap,xx))
    
  }

  proj_gdpcap[, dmg_spec := spec]
  proj_gdpcap[, runid := runid]

  return(proj_gdpcap)
  
}

# check countries

if (F) {
  loadd(climate_med)
  climate <- climate_med[model == "WITCH-GLOBIOM 4.4" & str_detect(scenario,"CD-LINKS") & c5model == 1]

  loadd(ssp_gdpcap,clim_hist,dmg_param)
  spec <- "BHM SR"
  
  res <- project_cc_impact_bhm(ssp_gdpcap,climate,clim_hist,dmg_param,spec) 
  
  ssps <- ssp_gdpcap[,.(ssp,iso3,year=as.numeric(year),pop)]
  map_ssp <- data.table(scenario = unique(res$scenario),
                        ssp = sapply(unique(res$scenario),scen_ssp, USE.NAMES = F))
  dd <- merge(res,map_ssp, by = "scenario")
  dd <- merge(dd,ssps,by = c("ssp","iso3","year"))
  dd[, gdp_cc := gdpcap_cc * pop * 1e-3]
  
  # ALL GDP_CAP
  ggplot(res[scenario == "CD-LINKS_NPi2020_1000" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=gdpcap_cc,group=iso3), color = "blue", alpha = 0.1)
  ggplot(dd[scenario == "CD-LINKS_NPi2020_1000" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=gdp_cc,group=iso3), color = "blue", alpha = 0.1)
  ggplot(res[scenario == "CD-LINKS_NoPolicy" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=gdpcap_cc,group=iso3), color = "blue", alpha = 0.1)
  ggplot(dd[scenario == "CD-LINKS_NoPolicy" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=gdp_cc,group=iso3), color = "blue", alpha = 0.1)
  
  # ALL GDP_R
  ggplot(res[scenario == "CD-LINKS_NPi2020_1000" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=deltax,group=iso3), color = "blue", alpha = 0.1)
  ggplot(res[scenario == "CD-LINKS_NPi2020_1600" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=deltax,group=iso3), color = "blue", alpha = 0.1)
  ggplot(res[scenario == "CD-LINKS_NoPolicy" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=deltax,group=iso3), color = "blue", alpha = 0.1)
  
  
  # GDP_CAP
  library(ggplot2)
  cnt <- "FRA"
  ggplot(res[iso3 == cnt & scenario == "CD-LINKS_NPi2020_1000" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=gdpcap_cc,color = "POL")) + 
    geom_line(aes(x=as.numeric(year),y=gdpcap_cc,color = "REF_CC"), 
              data = res[iso3 == cnt & scenario == "CD-LINKS_NoPolicy" & year <= 2300]) + 
    geom_line(aes(x=year,y=V1,color = "REF"),data=ssp_gdpcap[iso3 == cnt & ssp == "ssp2"  & year <= 2300, sum(gdpcap_nocc), by = c("year")])
  cnt <- "FRA"
  ggplot(dd[iso3 == cnt & scenario == "CD-LINKS_NPi2020_1000" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=gdp_cc,color = "POL")) + 
    geom_line(aes(x=as.numeric(year),y=gdp_cc,color = "REF_CC"), 
              data = dd[iso3 == cnt & scenario == "CD-LINKS_NoPolicy" & year <= 2300]) + 
    geom_line(aes(x=year,y=V1,color = "REF"),data=ssp_gdpcap[iso3 == cnt & ssp == "ssp2"  & year <= 2300, sum(gdpcap_nocc* pop * 1e-3), by = c("year")])
  
  cnt <- "USA"
  ggplot(res[iso3 == cnt & scenario == "CD-LINKS_NPi2020_1000" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=gdpcap_cc,color = "POL")) + 
    geom_line(aes(x=as.numeric(year),y=gdpcap_cc,color = "REF_CC"), 
              data = res[iso3 == cnt & scenario == "CD-LINKS_NoPolicy" & year <= 2300]) + 
    geom_line(aes(x=year,y=V1,color = "REF"),data=ssp_gdpcap[iso3 == cnt & ssp == "ssp2"  & year <= 2300, sum(gdpcap_nocc), by = c("year")])
  ggplot(dd[iso3 == cnt & scenario == "CD-LINKS_NPi2020_1000" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=gdp_cc,color = "POL")) + 
    geom_line(aes(x=as.numeric(year),y=gdp_cc,color = "REF_CC"), 
              data = dd[iso3 == cnt & scenario == "CD-LINKS_NoPolicy" & year <= 2300]) + 
    geom_line(aes(x=year,y=V1,color = "REF"),data=ssp_gdpcap[iso3 == cnt & ssp == "ssp2"  & year <= 2300, sum(gdpcap_nocc* pop * 1e-3), by = c("year")])
 
  cnt <- "ISL"
  ggplot(res[iso3 == cnt & scenario == "CD-LINKS_NPi2020_1000" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=gdpcap_cc,color = "POL")) + 
    geom_line(aes(x=as.numeric(year),y=gdpcap_cc,color = "REF_CC"), 
              data = res[iso3 == cnt & scenario == "CD-LINKS_NoPolicy" & year <= 2300]) + 
    geom_line(aes(x=year,y=V1,color = "REF"),data=ssp_gdpcap[iso3 == cnt & ssp == "ssp2"  & year <= 2300, sum(gdpcap_nocc), by = c("year")]) 
  ggplot(dd[iso3 == cnt & scenario == "CD-LINKS_NPi2020_1000" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=gdp_cc,color = "POL")) + 
    geom_line(aes(x=as.numeric(year),y=gdp_cc,color = "REF_CC"), 
              data = dd[iso3 == cnt & scenario == "CD-LINKS_NoPolicy" & year <= 2300]) + 
    geom_line(aes(x=year,y=V1,color = "REF"),data=ssp_gdpcap[iso3 == cnt & ssp == "ssp2"  & year <= 2300, sum(gdpcap_nocc* pop * 1e-3), by = c("year")])
  
  cnt <- "RUS"
  ggplot(res[iso3 == cnt & scenario == "CD-LINKS_NPi2020_1000" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=gdpcap_cc), color = "blue") + 
    geom_line(aes(x=as.numeric(year),y=gdpcap_cc), 
              data = res[iso3 == cnt & scenario == "CD-LINKS_NoPolicy" & year <= 2300], color = "red") + 
    geom_line(aes(x=year,y=V1),data=ssp_gdpcap[iso3 == cnt & ssp == "ssp2"  & year <= 2300, sum(gdpcap_nocc), by = c("year")],color="green")
  ggplot(dd[iso3 == cnt & scenario == "CD-LINKS_NPi2020_1000" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=gdp_cc,color = "POL")) + 
    geom_line(aes(x=as.numeric(year),y=gdp_cc,color = "REF_CC"), 
              data = dd[iso3 == cnt & scenario == "CD-LINKS_NoPolicy" & year <= 2300]) + 
    geom_line(aes(x=year,y=V1,color = "REF"),data=ssp_gdpcap[iso3 == cnt & ssp == "ssp2"  & year <= 2300, sum(gdpcap_nocc* pop * 1e-3), by = c("year")])
  
  cnt <- "FIN"
  ggplot(res[iso3 == cnt & scenario == "CD-LINKS_NPi2020_1000" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=gdpcap_cc), color = "blue") + 
    geom_line(aes(x=as.numeric(year),y=gdpcap_cc), 
              data = res[iso3 == cnt & scenario == "CD-LINKS_NoPolicy" & year <= 2300], color = "red") + 
    geom_line(aes(x=year,y=V1),data=ssp_gdpcap[iso3 == cnt & ssp == "ssp2"  & year <= 2300, sum(gdpcap_nocc), by = c("year")],color="green")
  ggplot(dd[iso3 == cnt & scenario == "CD-LINKS_NPi2020_1000" & year <= 2100]) + 
    geom_line(aes(x=as.numeric(year),y=gdp_cc,color = "POL")) + 
    geom_line(aes(x=as.numeric(year),y=gdp_cc,color = "REF_CC"), 
              data = dd[iso3 == cnt & scenario == "CD-LINKS_NoPolicy" & year <= 2100]) + 
    geom_line(aes(x=year,y=V1,color = "REF"),data=ssp_gdpcap[iso3 == cnt & ssp == "ssp2"  & year <= 2100, sum(gdpcap_nocc* pop * 1e-3), by = c("year")])
  
  
  #WORLD
  wdd <- dd[!iso3 %in% c("MNG","FIN","RUS","ISL","CAN","NOR","SWE"), .(gdp_cc = sum(gdp_cc)), by = "runid,dmg_spec,c5model,model,scenario,year"]
  ggplot(wdd[scenario == "CD-LINKS_NPi2020_1000" & year <= 2300]) + 
    geom_line(aes(x=as.numeric(year),y=gdp_cc,color = "POL")) + 
    geom_line(aes(x=as.numeric(year),y=gdp_cc,color = "REF_CC"), 
              data = wdd[scenario == "CD-LINKS_NoPolicy" & year <= 2300]) + 
    geom_line(aes(x=year,y=V1,color = "REF"),data=ssp_gdpcap[ssp == "ssp2"  & year <= 2300, sum(gdpcap_nocc* pop * 1e-3), by = c("year")])
  
  
  dd_share <- dd[,gdp_cc_share := gdp_cc / sum(gdp_cc), 
              by = "runid,dmg_spec,c5model,model,year"]
  
  
  # table
  rest <- dcast(res, runid + dmg_spec + c5model + model + iso3 + year ~ scenario, value.var = "gdpcap_cc", fun.aggregate = sum)
  
  res_share <- res[, gdpcap_cc_share := gdpcap_cc / sum(gdpcap_cc)  , by = "runid,dmg_spec,c5model,model,year"]
  rest_share <- dcast(res_share, runid + dmg_spec + c5model + model + iso3 + year ~ scenario, 
                      value.var = "gdpcap_cc_share", fun.aggregate = sum)
  
}