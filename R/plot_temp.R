plot_temp <- function(outfile) {
  
  temp <- readd('temp_magicc')
  gmt_hist <- readd('gmt_hist')
  sr15c_runs <- readd('sr15c_runs')
  
  scen_ref <- sr15c_runs$scenario_ref
  
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
  
  #temp <- temp[temp_2100 <= 2]
  temp[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  temp[, temp_clus := ifelse(temp_2100 > 2, 0, temp_clus)]
  temp <- temp[temp_clus >= 1 | scenario %in% scen_ref ]
  temp[, temp_clus := factor(temp_clus, levels = 0:2, labels = c("Reference","1.5\u00B0C", "2\u00B0C"))]
  temp[, temp_2100 := NULL]

  end_temp <- temp[year == 2100, .(ymin = min(value), ymax = max(value)), by = "temp_clus"]
  
  p <- ggplot(temp) +
    geom_line(aes(x = year, 
                  y = value,
                  color = temp_clus,
                  group = paste(model,scenario)),
              size = 0.5,
              alpha = 0.5) +
    geom_line(aes(x = year,
                  y = gmt),
              data = gmt_hist) +
    geom_linerange(aes(x = 2101, 
                       ymin = ymin, 
                       ymax = ymax, 
                       color = temp_clus),
                   size = 1,
                   data = end_temp) +
    geom_text(aes(x = 2102, y = (ymin + ymax) / 2, 
                  label = temp_clus),
              hjust = 0,
              data = end_temp) +
    labs(x = "", y = "[\u00B0C]") +
    scale_color_brewer(palette = "Set2", guide = "none") +
    scale_y_continuous(breaks = seq(0,7,0.5)) +
    coord_cartesian(xlim = c(2010, 2108), ylim = c(0.66,5.25)) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 10)) +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "gray"))
  p
  
  ggsave(outfile, width = 8, height = 5)
  
}

plot_temp_quantile <- function(outfile) {
  
  temp <- readd('temp_magicc')
  temp_p05 <- readd('temp_magicc_p05')
  temp_p95 <- readd('temp_magicc_p95')
  
  temp[, qscen := "50"]
  temp_p05[, qscen := "05"]
  temp_p95[, qscen := "95"]
  
  temp <- rbindlist(list(temp,temp_p05,temp_p95))
  
  gmt_hist <- readd('gmt_hist')
  sr15c_runs <- readd('sr15c_runs')
  
  scen_ref <- sr15c_runs$scenario_ref
  
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
  
  #temp <- temp[temp_2100 <= 2]
  temp[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  temp[, temp_clus := ifelse(temp_2100 > 2, 0, temp_clus)]
  temp <- temp[temp_clus >= 1 | scenario %in% scen_ref ]
  temp[, temp_clus := factor(temp_clus, levels = c(0,2,1), labels = c("Reference","2\u00B0C","1.5\u00B0C"))]
  temp[, temp_2100 := NULL]

  temp_ribbon <- temp[, .(temp_min = min(value),
                          temp_max = max(value),
                          temp_med = quantile(value,0.5)),
                      by = "year,temp_clus,qscen"]
  
  temp_ribbon[qscen == "50", qscen_label := "median"]
  temp_ribbon[qscen == "05", qscen_label := "5th"]
  temp_ribbon[qscen == "95", qscen_label := "95th"]
  
  temp_ribbon[, qscen_label := factor(qscen_label, 
                                      levels = c("5th",
                                                 "median",
                                                 "95th"))]

  p <- ggplot(temp_ribbon) +
    geom_ribbon(aes(x = year, 
                    ymin = temp_min,
                    ymax = temp_max,
                    fill = qscen_label),
                alpha = 0.8) +
    geom_line(aes(x = year,
                  y = gmt),
              data = gmt_hist) +
    facet_wrap(~ temp_clus) +
    labs(x = "", y = "[\u00B0C]") +
    scale_fill_brewer(palette = "Set1",
                       name = "Temperature\npercentile") +
    scale_y_continuous(breaks = seq(0,10,0.5)) +
    coord_cartesian(xlim = c(2010, 2100), ylim = c(0.66,9)) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 10)) +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "gray")) +
    theme(legend.position = c(.85,.7))
  p
  
  ggsave(outfile, width = 8, height = 5)
  
}

  