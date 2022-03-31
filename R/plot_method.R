


plot_method_temp <- function(outfile) {
  
  temp <- readd('temp_magicc')
  gmt_hist <- readd('gmt_hist')
  sr15c_runs <- readd('sr15c_runs')
  
  temp[, year := as.numeric(year)]
  
  hist_magicc <- temp[year >= 2005 & year <= 2019, .(gmt = mean(value)), by = year]
  hist_cru <- gmt_hist[year >= 2005 & year <= 2019, gmt]
  
  diffgmt <- function(x) {
    return(abs(sum(hist_magicc$gmt - x - hist_cru)))
  }
  search <- optimize(diffgmt, c(0,1)) # minimizing distance between MAGICC and CRU
  d2 <- search$minimum # approx = 0.078
  
  # Unbiased temperature from MAGICC
  temp <- temp[,value := value - d2]
  
  temp <- temp[year >= 2010]
  temp <- merge(temp,sr15c_runs, by = c("model","scenario"))
  
  temp <- temp[temp_2100 <= 2 | (scenario_ref == scenario)]
  temp[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  temp[, temp_clus := ifelse(scenario_ref == scenario, 3, temp_clus)]
  temp[, temp_clus := factor(temp_clus, levels = c(3,2,1), labels = c("Baseline","2\u00B0C","1.5\u00B0C"))]
  temp[, temp_2100 := NULL]
  
  end_temp <- temp[year == 2100, .(ymin = min(value), ymax = max(value)), by = "temp_clus"]
  
  end_temp_rib <- rbind(end_temp[,.(year = 2000, temp_clus, ymin, ymax)],
                        end_temp[,.(year = 2120, temp_clus, ymin, ymax)])
  
  p <- ggplot(temp) +
    geom_ribbon(aes(x = year, 
                       ymin = ymin, 
                       ymax = ymax, 
                       fill = temp_clus),
                   size = 1,
                   alpha = 0.2,
                   data = end_temp_rib) +
    geom_line(aes(x = year, 
                  y = value,
                  color = temp_clus,
                  group = paste(model,scenario)),
              size = 0.5,
              alpha = 0.5) +
    geom_line(aes(x = year,
                  y = gmt),
              data = gmt_hist) +
    geom_text(aes(x = ifelse(temp_clus=="Baseline",2065,2101), y = (ymin + ymax) / 2, 
                  label = temp_clus),
              hjust = 0,
              size = 3,
              data = end_temp) +
    labs(x = "", y = "[\u00B0C]") +
    scale_color_npg(guide = "none") +
    scale_fill_npg(guide = "none") +
    scale_y_continuous(breaks = c(0,1,2,3,4,5,6)) +
    coord_cartesian(xlim = c(2010, 2108), ylim = c(0.66,4.5)) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 10)) +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "gray"))
  p
  
  ggsave(outfile, width = 5, height = 3.5)
  
}



plot_method_temp_map <- function(outfile) {
  
  dd <- readd('climate_med_092603c2')
  
  dd <- dd[c5model == 1 & model == "WITCH-GLOBIOM 4.4" & scenario == "CD-LINKS_NPi2020_1000" & year == 2100]
  
  world <- ne_countries(scale = "small", returnclass = "sf")
  world <- subset(world,!adm0_a3 %in% c("ATA","FJI")) 
  
  world <- merge(world,dd, by.x = "adm0_a3", by.y = "iso3")

  # Remove small islands
  target_crs <- '+proj=eqearth +wktext'
  world0 <- st_transform(world, crs = target_crs)
  
  p <- ggplot(data = world0) +
    geom_sf(aes(fill = tas)) +
    coord_sf(datum = target_crs,
             expand = FALSE,
             clip = "off") +
    theme_void() + 
    scale_fill_viridis_c(option = "magma", name = "[\u00B0C]") + 
    theme(legend.position = c(0.7,0),
          legend.direction = "horizontal",
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 12))
  p
  
  ggsave(outfile, width = 5, height = 3.5)
  
}


plot_method_lvl_dmg <- function(outfile) {
  
  loadd('dmg_param')
  
  # compute damage for warming between 0.5 to 2.5
  ddmg <- list()
  tt <- seq(0.5,6.2,by = 0.2)
  for (dmgf in c("DICE2016","HS2017 NCAT","HS2017 TOT","KW2020")) {
  dmg <- dmg_param[[dmgf]][1] * tt + dmg_param[[dmgf]][2] * tt * tt
    ddmg <- c(ddmg,list(data.table(dmgf = dmgf, id = 1, gmt = tt, dmg = dmg)))
  }
  for (dmgf in c("TAKAKURA2019")) {
    for (ssp_num in 1:5) {
      dmg <- - (dmg_param[[dmgf]][ssp_num,1] * tt +
        dmg_param[[dmgf]][ssp_num,2] * tt * tt +
        dmg_param[[dmgf]][ssp_num,3]) / 100
      ddmg <- c(ddmg,list(data.table(dmgf = dmgf, id = ssp_num, gmt = tt, dmg = dmg)))
    }
  }
  
  ddmg <- rbindlist(ddmg)
  
  
  p <- ggplot(data = ddmg) +
    geom_hline(yintercept = 0, color = "#888888") +
    geom_line(aes(x = gmt, y = dmg * 100, color = dmgf, linetype = factor(id), group = paste(id,dmgf)), size = 1) +
    scale_color_npg(name = "") +
    theme_bw() +
    theme(legend.position = c(0.2,0.4)) +
    theme(legend.text = element_text(size = 8)) +
    theme(legend.background = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "gray")) +
    scale_x_continuous(breaks = c(0,1,2,3,4,5,6)) +
    scale_linetype(guide = "none") +
    labs(x = "[\u00B0C]", y = "GDP loss [%]")
  p
  
  ggsave(outfile, width = 5, height = 3.5)
  
}

plot_method_gwt_dmg <- function(outfile) {
  
  loadd('dmg_param')
  
  # compute damage for warming between 0.5 to 2.5
  ddmg <- list()
  tt <- seq(-4, 30, 0.2)
  for (dmgf in c("BHM SR")) {
    for (i in 2:999) {
      dmg <- dmg_param[[dmgf]][i,1] + dmg_param[[dmgf]][i,2] * 2 * tt
      ddmg <- c(ddmg,list(data.table(dmgf = dmgf, id = i, vrt = "all", analysis = "1.", temp = tt, dmg = dmg, opt = dmg_param[[dmgf]][i,3], size = 0.6)))
    }
  }
  for (dmgf in c(bhm_spec[bhm_spec_star==1])) {
      if (!str_detect(dmgf,"RP")) {
  
        b1 <- dmg_param[[dmgf]][1,1]
        b2 <- dmg_param[[dmgf]][1,2]
        
          opt = -b1/(2*b2)
          dmg <- b1 + b2 * 2 * tt
          ddmg <- c(ddmg,list(data.table(dmgf = dmgf, id = 1, 
                                         vrt = "all", analysis = "2.", 
                                         temp = tt, dmg = dmg, 
                                         opt = opt, 
                                         size = 0.8)))
          

      } else {

        b1 <- dmg_param[[dmgf]][1,1]
        b2 <- dmg_param[[dmgf]][1,2]
        b3 <- dmg_param[[dmgf]][1,3]
        b4 <- dmg_param[[dmgf]][1,4]
                
        vrt <- "rich"
          opt = -b1/(2*b3)
          dmg <- b1 * tt + b3 * tt * tt - (b3*opt^2+b1*opt)
          dmg <- b1  + b3 * 2 * tt
          ddmg <- c(ddmg,list(data.table(dmgf = dmgf, id = 1, 
                                         vrt = vrt, analysis = "2.", 
                                         temp = tt, dmg = dmg, 
                                         opt = opt, 
                                         size = 0.8)))

        vrt <- "poor"
        opt = -b2/(2*b4)
        dmg <- b2 + b4 * 2 * tt
        ddmg <- c(ddmg,list(data.table(dmgf = dmgf, id = 1, 
                                       vrt = vrt, analysis = "2.", 
                                       temp = tt, dmg = dmg, 
                                       opt = opt, 
                                       size = 0.8)))
    }
  }
  for (dmgf in c(djo_spec[djo_spec_star==1])) {

      
      b1 <- dmg_param[[dmgf]][1]
      b2 <- dmg_param[[dmgf]][2]

      vrt <- "rich"

      if( b1 != 0) {
        opt = NA
        dmg <- b1 
        ddmg <- c(ddmg,list(data.table(dmgf = dmgf, id = 1, 
                                       vrt = vrt, analysis = "2.", 
                                       temp = tt, dmg = dmg, 
                                       opt = opt, 
                                       size = 0.8)))
      }

      vrt <- "poor"
      if( b2 != 0) {
        opt = NA
        dmg <- b2 
        ddmg <- c(ddmg,list(data.table(dmgf = dmgf, id = 1, 
                                       vrt = vrt, analysis = "2.", 
                                       temp = tt, dmg = dmg, 
                                       opt = opt, 
                                       size = 0.8)))
      }

  }
  
  
  ddmg <- rbindlist(ddmg)
  ddmg[opt < -4, opt := NA]
  ddmg[opt > 30, opt := NA]
  
  p <- ggplot(data = ddmg) +
    geom_line(aes(x = temp, y = dmg * 100, 
                  color = as.numeric(opt), 
                  #color = vrt,
                  group = paste(dmgf,id,analysis,vrt),
                  size = size,
                  linetype = vrt,
                  alpha = ifelse(analysis == 1, 0.1, 0.5))) +
    geom_hline(yintercept = 0, color = "#888888") +
    scale_color_viridis_c(name = "opt.") +
    scale_alpha(guide = "none") +
    scale_linetype(name = "") +
    scale_size_identity(guide = "none") +
    facet_wrap(~ analysis) +
    theme_bw() +
    theme(legend.position = c(0.85,0.85)) +
    theme(legend.text = element_text(size = 8)) +
    theme(legend.background = element_blank()) +
    theme(legend.box = "horizontal") +
    theme(strip.background = element_blank()) +
    theme(strip.text = element_text(hjust = 0)) +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "gray")) +
    labs(x = "[\u00B0C]", y = "Growth effect [%/\u00B0C]")
  p
  
  ggsave(outfile, width = 5, height = 3.5)
  
}

plot_method_admg <- function(outfile) {
  
  dd <- readd('ssp_gdpcap')
  dd[, year := as.numeric(year)]
  dd <- dd[ssp == "ssp2" & iso3 == "FRA" & year >= 2020 & year <= 2100]
  dd[, coeff := seq(0,1, length.out = 81)]
  dd <- dd[, gdpcap_pol_cc := gdpcap_nocc * (1 - 0.1 * coeff)]
  dd <- dd[, gdpcap_cc := gdpcap_nocc * (1 - 0.33 * coeff)]
  
  llab <- fread("year,value,label,color,fontface
                2097,126218,+ policy climate,b,bold
                2095,81199,+ baseline climate,a,bold")
  #2057,114687,SSP Baseline,a,plain
  
  p <- ggplot(data = dd) +
    geom_ribbon(aes(x = year, ymin = gdpcap_pol_cc, ymax = gdpcap_cc, fill = 'Avoided damages'), alpha = 0.7) +
    geom_line(aes(x = year, y = gdpcap_nocc), linetype = 2) +
    geom_line(aes(x = year, y = gdpcap_cc, color = "a"), size = 1) +
    geom_line(aes(x = year, y = gdpcap_pol_cc, color = "b"), size = 1) +
    geom_label(aes(x = year, y = value, label = label, color = color, fontface = fontface),
                    hjust = 0,
                    data = llab) +
    annotate('label', x = 2062, y = 110000, label = 'SSP Baseline') +
    scale_x_continuous(expand = c(0,0), limits = c(2020, 2130), breaks = c(2025, 2050, 2075, 2100)) +
    scale_color_npg(guide = "none") +
    scale_fill_manual(values = "gray", name = "") +
    labs(x = "", y = "GDP [T$]") +
    theme_bw() +
    theme(legend.position = c(0.8,0.2)) +
    theme(legend.background = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "gray")) +
    theme(axis.text.y = element_blank())
  p
  
  ggsave(outfile, width = 5, height = 3.5)
  
}

plot_method_cmit <- function(outfile) {
  
  gdp_mit <- readd('gdp_mit')
  sr15c_runs <- readd('sr15c_runs')
  
  dd <- copy(gdp_mit)
  dd[, year := as.numeric(year)]
  
  cmit <- dd[year >= 2010]
  cmit <- merge(cmit,sr15c_runs, by = c("model","scenario"))
  cmit[, cmit_share := gdp_mit / gdp_ref - 1]
  cmit <- cmit[temp_2100 <= 2]
  cmit[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  cmit[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5\u00B0C", "2\u00B0C"))]
  cmit[, temp_2100 := NULL]
  cmit[cmit_share > 0, cmit_share := 0]
  
  cmit_stat <-   cmit[region == "World", .(cmin = quantile(cmit_share,0.05),
                                           cmax = quantile(cmit_share,0.95),
                                           c83 = quantile(cmit_share,2.5 / 3),
                                           c16 = quantile(cmit_share,0.5 / 3),
                                           cmed = quantile(cmit_share,0.5),
                                           N = .N),
                      by = c("temp_clus","year")]

  llab <- fread("year,value,label,color,fontface
                2075,3.3,Policy + baseline climate,b,bold
                2075,0.5,No policy + baseline climate,a,bold")
  #2057,114687,SSP Baseline,a,plain
  
  cmed_smooth <- loess(cmed ~ year, data=cmit_stat[temp_clus == "2\u00B0C"])
  
  dd <- data.table(year = 2020:2100, cmed = predict(cmed_smooth, newdata = data.frame(year=2020:2100)))
                   
  p <- ggplot(data = dd) +
    geom_ribbon(aes(x = year, ymin = 0, ymax = - cmed * 100, fill = 'Mitigation costs'), alpha = 0.7) +
    geom_line(aes(x = year, y = 0, color = "a"), linetype = 1, size = 1) +
    geom_line(aes(x = year, y = - cmed * 100, color = "b"), size = 1) +
    #geom_line(aes(x = year, y = gdpcap_pol_cc, color = "b"), size = 1) +
    geom_label(aes(x = year, y = value, label = label, color = color, fontface = fontface),
               hjust = 0,
               data = llab) +
    #annotate('label', x = 2062, y = 110000, label = 'SSP Baseline') +
    scale_x_continuous(expand = c(0,0), limits = c(2020, 2130), breaks = c(2025, 2050, 2075, 2100)) +
    scale_y_continuous(limits = c(0,5)) +
    scale_color_npg(guide = "none") +
    scale_fill_manual(values = "gray", name = "") +
    labs(x = "", y = "Costs [T$]") +
    theme_bw() +
    theme(legend.position = c(0.7,0.9)) +
    theme(legend.background = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "gray")) +
    theme(axis.text.y = element_blank())
  p
  
  ggsave(outfile, width = 5, height = 3.5)
  
}