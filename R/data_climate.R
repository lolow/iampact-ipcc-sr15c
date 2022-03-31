# CLIMATE

downscale_climate <- function(clim_dwnscl, 
                              temp_magicc, 
                              cmodel,
                              sr15c_runs) {
  
  estim <- clim_dwnscl[mo == cmodel][r.squared > 0.5]
  
  tt <- temp_magicc[year >= 2010 & year <= 2100]
  tt <- tt[scenario %in% c(sr15c_runs$scenario,sr15c_runs$scenario_ref)]
  
  # rcp_temp 2005 0.7548798 historical
  # gmt_2005 = 0.8079997
  gmt_2005 <- mean(temp_magicc[year == 2005,value])
  
  tas <- lapply(tt$value, function(gmt) estim$temp * (gmt - gmt_2005) + estim$inter) %>%
    unlist(use.names = FALSE)
  
  # Add indexes
  clim = data.table(tas = tas)
  clim[, c5model := which(cmodels == cmodel)]
  clim[, iso3 := rep(estim$reg, times = nrow(tt))]
  clim[, model := rep(tt$model, each = nrow(estim))]
  clim[, scenario := rep(tt$scenario, each = nrow(estim))]
  clim[, year := rep(tt$year, each = nrow(estim))]
  clim[, sw := 1]
  
  # Keep less years to save memory
  clim <- clim[year %in% seq(2010,2100,5)]
  
  setkey(clim, c5model, model, scenario, iso3, year, sw)

  return(clim)
  
}


read_hist_climate <- function(clim_hist_csv, climate) {
  
  # load historical temp 1980-2017
  dd <- fread(clim_hist_csv)
  
  # select iso3 and columns
  dd <- dd[iso3 %in% unique(climate$iso3), .(year,iso3,tas = tas_pop_wgt_mean_hist)]

  return(dd)
    
}

read_hist_gmt <- function() {
  
  hadcrut4gl <- fread("https://crudata.uea.ac.uk/cru/data/temperature/HadCRUT4-gl.dat", fill = TRUE)
  
  # keep year and annual value
  cru_ann <- hadcrut4gl[!is.na(V14), c(1, 14)]
  setnames(cru_ann, c("V1", "V14"), c("year", "gmt"))
  
  # Keep only years only 2019
  cru_ann <- cru_ann[year < 2020]
  
  # average 1850-1880
  ref_preind <- cru_ann[year >= 1850 & year <= 1880, mean(gmt)]
  
  # compute gmt from preind
  cru_ann[, gmt := gmt - ref_preind]
  
  return(cru_ann)
}
