plot_cmit <- function(outfile) {
  
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
  cmit_stat[, N := max(N), by = "temp_clus"]

  p <- ggplot(cmit[region == "World"]) +
    geom_line(aes(x = year, 
                  y = - cmit_share * 100,
                  group = paste(model,scenario)), 
              color = 'gray',
              size = 0.5,
              alpha = 0.5) +
    geom_errorbar(aes(x = year,
                      ymin = - cmin * 100,
                      ymax = - cmax * 100),
                  size = 1,
                  data = cmit_stat[year %in% c(2050,2100)],
                  width = 3) +
    geom_linerange(aes(x = year,
                       ymin = - c16 * 100,
                       ymax = - c83 * 100),
                   size = 2,
                   data = cmit_stat[year %in% c(2050,2100)]) +
    geom_point(aes(x = year,
                   y = - cmed * 100),
               data = cmit_stat[year %in% c(2050,2100)],
               size = 4) +
    geom_point(aes(x = year,
                   y = - cmed * 100),
               data = cmit_stat[year %in% c(2050,2100)],
               color = "white",
               size = 3) +
    geom_text(aes(x = 2010, 
                  y = max(cmit[region == "World"]$cmit_share*-100),
                  label = paste("N =",N)),
              data = cmit_stat[year %in% c(2010)],
              hjust = 0.25,
              vjust = 0) +
    #scale_x_continuous(expand = c(0.01,0.01)) +
    #scale_y_continuous(expand = c(0.01,0.01)) +
    facet_wrap(~ temp_clus) +
    labs(x = "", y = "[% GDP]") +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 10))
  p
  
  ggsave(outfile, width = 8, height = 5)
  
}


  