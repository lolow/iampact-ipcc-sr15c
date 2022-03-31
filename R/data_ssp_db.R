# SSP original assumption database

# Extract variable from SSP database from OECD projection
ssp_query <- function(ssp_csv,src,variable) {
  ssp_header <- tolower(unlist(fread(ssp_csv, nrows = 1), use.names = F))
  ssp_header_ids <- ssp_header[1:5]
  ssp_coltype <- c(rep("character",5),rep("numeric",length(ssp_header) - 5))
  
  dd <- fread(str_glue("grep -e '\"{variable}\"' {ssp_csv}"), 
        col.names = ssp_header,
        colClasses = ssp_coltype) %>%
    long_table(ids = ssp_header_ids) %>%
    .[!is.na(value)] %>%
    .[model == src] %>%
    .[, year := as.numeric(year)] %>%
    .[, ssp := tolower(str_extract(scenario,'SSP\\d'))] %>%
    setnames('region','iso3') %>%
    .[, .(ssp,iso3,year,value)] %>%
    setkey(year)
  
  return(dd)
  
}

exp_interp <- function(x,y,x0 = 2018, y0 = NULL) {
  xout = 2010:2100
  # interpolate linearly yearly gdp
  yout <- approx(x,y,xout)$y
  # original rate
  rate <- c()
  for (i in 1:(length(x) - 1)) {
    rate[i] <- (y[i + 1] / y[i])^(1 / (x[i + 1] - x[i])) - 1
  }
  rate_out <- approx(x,c(0,rate),xout)$y
  yout2 <- yout
  for (i in 2:(length(xout))) {
    yout2[i] <- yout2[i - 1] * (1 + rate_out[i])
  }
  # Ensure 2018 value is correct, keeping same growth rate
  yout2 <- yout2 * y0 / yout2[which(xout==x0)]
  return(yout2)
}

# Compute GDP per cap
compute_gdp_cap <- function(ssp_gdp,ssp_pop,wdi_stat) {
  
  # Collect SSP GDP PPP (Original in billion US$2005/yr)
  # and GDP MER historical value in 2018 from WDI  
  .gdp <- ssp_gdp[year >= 2005,.(ssp,iso3,year,gdp = value)]
  .gdp <- merge(.gdp,
                wdi_stat[year == 2018, .(iso3,gdp_2018 = gdp * 1e-9)],
                by = 'iso3', all.x = TRUE)
  .gdp[iso3 == "IRN", gdp_2018 := wdi_stat[year == 2017 & iso3 == "IRN",gdp] * 1e-9] # Missing data, estimated from 2017 WDI
  miss <- .gdp[year %in% c(2015,2020) & is.na(gdp_2018) & !iso3 %in% c("IRN")]
  miss <- miss[,.(gdp = mean(gdp)), by = "iso3,year"]
  miss <- miss[,.(gdp_2018_estimate = approx(year,gdp,2018)$y), by = "iso3"]
  .gdp <- merge(.gdp, miss, by = 'iso3', all.x = TRUE)
  .gdp[is.na(gdp_2018), gdp_2018 := gdp_2018_estimate]
  .gdp[, gdp_2018_estimate := NULL]
  # Country with GDP == 0 [TWN and SOM] in 2010
  .gdp[year == 2005 & iso3 == "TWN", 
       gdp := 2 * .gdp[year == 2010 & iso3 == "TWN", mean(gdp)] - 
         .gdp[year == 2015 & iso3 == "TWN", mean(gdp)]]
  .gdp[year == 2005 & iso3 == "SOM", 
       gdp := 2 * .gdp[year == 2010 & iso3 == "SOM", mean(gdp)] - 
         .gdp[year == 2015 & iso3 == "SOM", mean(gdp)]]
  # WDI includes TWN in China, but not SSP so remove it
  .gdp[iso3 == "CHN", 
       gdp_2018 := gdp_2018 - miss[iso3 == "TWN"]$gdp_2018_estimate]

  # Recalibrate to 2018 historical GDP MER, keeping based on SSP growth
  .gdp <- .gdp[,.(year = 2010:2100, 
                  gdp = exp_interp(year,gdp,2018,gdp_2018[1])), 
               by = "ssp,iso3"]

  # average SSP until 2018
  .gdp[year <= 2018, gdp := mean(gdp), by = c("iso3,year")]

  # Convert SSP POP (Original in millions)
  .pop <- ssp_pop[year >= 2005,.(ssp,iso3,year,pop = value)]
  .pop <- merge(.pop,
                wdi_stat[year == 2018, .(iso3,pop_2018 = pop * 1e-6)],
                by = 'iso3', all.x = TRUE)
  .pop[iso3 == "IRN", pop_2018 := 81.79936] # Missing data, estimated from 2016-2017 WDI
  miss <- .pop[year %in% c(2015,2020) & is.na(pop_2018) & !iso3 %in% c("IRN")]
  miss <- miss[,.(pop = mean(pop)), by = "iso3,year"]
  miss <- miss[,.(pop_2018_estimate = approx(year,pop,2018)$y), by = "iso3"]
  .pop <- merge(.pop, miss, by = 'iso3', all.x = TRUE)
  .pop[is.na(pop_2018), pop_2018 := pop_2018_estimate]
  .pop[, pop_2018_estimate := NULL]
  # Country with POP == 0 [TWN]
  .pop[year == 2005 & iso3 == "TWN", 
       pop := 2 * .pop[year == 2010 & iso3 == "TWN", mean(pop)] - 
         .pop[year == 2015 & iso3 == "TWN", mean(pop)]]
  # WDI includes TWN in China, but not SSP so remove it
  .pop[iso3 == "CHN", 
       pop_2018 := pop_2018 - miss[iso3 == "TWN"]$pop_2018_estimate]
  
  # Recalibrate to 2018 historical value
  .pop <- .pop[,.(year = 2010:2100, 
                  pop = exp_interp(year,pop,2018,pop_2018[1])), 
               by = "ssp,iso3"]
  
  # average SSP until 2018
  .pop[year <= 2018, pop := mean(pop), by = c("iso3,year")]

  dd <- merge(.gdp,.pop,by = c("ssp","iso3","year"))
  dd[, gdpcap_nocc := gdp / pop * 1e3]
  dd[, gdpr := (gdpcap_nocc / shift(gdpcap_nocc, 1)) - 1, by = "ssp,iso3"]
  dd <- dd[year >= 2010, .(ssp,iso3,year,gdpcap_nocc,gdpr,pop)]
  setkey(dd, ssp, iso3, year)

  return(dd)
}