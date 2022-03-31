#load data from World Development Indicators
library(WDI)
library(data.table)
library(fst)

# Check units with
#WDIsearch("GDP, PPP")
#WDIsearch("Population, total")

wdi_varlist <- fread('wdi_identifier,description,varname
SP.POP.TOTL,Population,pop
NY.GDP.MKTP.PP.KD,"GDP, PPP (constant ???? international $)",gdp_ppp
NY.GDP.MKTP.KD,"GDP (constant 2010 US$)",gdp
') # changing year!! as now 2017
#NY.GDP.MKTP.IN,"GDP Deflator",gdp_deflator# cannot be downloaded

allwdi <- data.table(WDI(country = "all", 
                         start = 1960, 
                         end = 2018,
                         indicator = wdi_varlist$wdi_identifier, 
                         extra = TRUE))
setnames(allwdi, "iso3c", "iso3")
setnames(allwdi, wdi_varlist$wdi_identifier, wdi_varlist$varname)

wdi_stat <- allwdi[iso3 != "<NA>" & !is.na(pop) & !is.na(gdp), .(year, iso3, gdp, gdp_ppp, pop)]

# Compute GDP in USD2019
deflat <- fread('data/gdp_deflator_usa.csv')
deflat_2005 <- deflat[year == 2005, value]
deflat_2010 <- deflat[year == 2011, value]
deflat_2011 <- deflat[year == 2011, value]
deflat_2017 <- deflat[year == 2017, value]
deflat_2018 <- deflat[year == 2018, value]
deflat_2015 <- deflat[year == 2015, value]

wdi_stat[, gdp := gdp / deflat_2010 * deflat_2018]
wdi_stat[, gdp_ppp := gdp_ppp / deflat_2017 * deflat_2018]

write_fst(wdi_stat,'data/wdi_stat.fst')
