
# Select only decade year
filter_boot_sample <- function(nben_ann) {
  return(nben_ann[year %in% seq(2020,2100,by = 10) &
                    !is.na(nbenefit) &
                    region == "World" &
                    dmg_spec %chin% "BHM SR BOOT"])
}

plot_nben_boot_rcp <- function(dr0,dmg_spec0,output_file) {
  
  nben_npv_boot <- readd('nben_npv_boot')
  sr15c_runs <- readd('sr15c_runs')
  
  disc_rate = as.numeric(dr0) / 100
  
  nbft <- copy(nben_npv_boot[dr == disc_rate &
                               region == "World" &
                               dmg_spec == dmg_spec0 &
                               c5model == 1])
  nbft <- merge(nbft, sr15c_runs, by = c("model","scenario"))
  nbft <- nbft[,.(model,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft <- nbft[temp_2100 <= 2]
  nbft[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft[, temp_2100 := NULL]
  
  setkey(nbft,dmg_spec,model,scenario,temp_clus)

  # stats
  
  #Virtually certain: 99 to 100 per cent probability
  #Extremely likely: Over 95 per cent. According to the IPCC, it is now “extremely likely” that human activity is to blame for climate change
  #Very likely: Above 90 per cent
  #Likely: Above 66 per cent
  
  nbftstat <- nbft[, .(cmin = quantile(nbenefit, 0.5 / 3 ), # likely range - IPCC [66%]
                       cmax = quantile(nbenefit, 2.5 / 3), 
                       cmed = quantile(nbenefit, 0.5),      # median
                       c05  = quantile(nbenefit, 0.05),     # very likely range - IPCC [90%]
                       c95  = quantile(nbenefit, 0.95),
                       c01  = quantile(nbenefit, 0.01),     # 
                       c99  = quantile(nbenefit, 0.99),
                       mmax  = max(nbenefit),     # 
                       mmin  = min(nbenefit),
                       mmm  = mean(nbenefit)),              # mean
                   by = c("dmg_spec","temp_clus")]

  p <- ggplot(nbftstat, aes(x = temp_clus)) +
    geom_hline(yintercept = 0, color = "black") +
    geom_violin(mapping = aes(x = temp_clus, 
                              y = nbenefit * 1e2,
                              fill = dmg_spec,
                              group = paste(temp_clus,dmg_spec)), 
                data = nbft, 
                trim = F,
                position = position_dodge(width = 0.33),
                fill = 'gray', alpha = 0.5, color = NA) +
    geom_errorbar(aes(ymin = c05 * 1e2, 
                       ymax = c95 * 1e2,
                       group = paste(dmg_spec)),
                   size = 1,
                  width = 0.33,
                   position = position_dodge(width = 0.33)) +
    geom_linerange(aes(ymin = cmin * 1e2, 
                       ymax = cmax * 1e2,
                       group = paste(dmg_spec)),
                   size = 5,
                   position = position_dodge(width = 0.33)) +
    geom_point(aes(y = cmed * 1e2,
                   group = paste(dmg_spec)),
               position = position_dodge(width = 0.33),
               size = 5) +
    geom_point(aes(y = cmed * 1e2,
                   group = paste(dmg_spec)), 
               color = "white",
               position = position_dodge(width = 0.33),
               size = 4) +
    scale_color_npg(name = "") +
    scale_y_continuous(labels = scales::percent_format(scale = 1), expand = c(0,0)) + 
    coord_cartesian(ylim = c(nbftstat[,min(c01) * 1e2 * 1.1],
                             nbftstat[,max(c99) * 1e2 * 1.1])) +
    labs(
      y = "Policy benefits [%]", 
      x = "",
      title = str_glue("{round(dr0)}% dr.")
    ) +
    theme_light() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text = element_text(size = 14),
          axis.line.x = element_blank())
  
  ggsave(output_file, plot = p, width = 6, height = 6)
  
  output_file
  
}

plot_nben_boot_rcp_nocluster <- function(dr0,dmg_spec0,output_file, main = FALSE) {
  
  nben_npv_boot <- readd('nben_npv_boot')
  nben_npv_boot_p05 <- readd('nben_npv_boot_p05')
  nben_npv_boot_p95 <- readd('nben_npv_boot_p95')
  
  sr15c_runs <- readd('sr15c_runs')
  
  disc_rate = as.numeric(dr0) / 100
  
  nbft <- copy(nben_npv_boot[dr == disc_rate &
                               region == "World" &
                               dmg_spec == dmg_spec0 &
                               c5model == 1])
  nbft_p05 <- copy(nben_npv_boot_p05[dr == disc_rate &
                               region == "World" &
                               dmg_spec == dmg_spec0 &
                               c5model == 1])
  nbft_p95 <- copy(nben_npv_boot_p95[dr == disc_rate &
                                       region == "World" &
                                       dmg_spec == dmg_spec0 &
                                       c5model == 1])
  nbft[, qscen := "50"]
  nbft_p05[, qscen := "05"]
  nbft_p95[, qscen := "95"]
  nbft <- rbindlist(list(nbft,nbft_p05,nbft_p95))
  nbft <- merge(nbft, sr15c_runs, by = c("model","scenario"))
  nbft <- nbft[,.(model,scenario,dmg_spec,qscen,nbenefit,temp_2100)]
  
  nbft <- nbft[temp_2100 <= 2]
  nbft[, temp_2100 := NULL]
  
  nbft[qscen == "05", qscen_label := "5th"]
  nbft[qscen == "50", qscen_label := "median"]
  nbft[qscen == "95", qscen_label := "95th"]
  nbft[, qscen_label := factor(qscen_label, levels = c("5th",
                                                           "median",
                                                           "95th"))]
  
  setkey(nbft,qscen_label,dmg_spec,qscen,model,scenario)
  
  # stats
  
  #Virtually certain: 99 to 100 per cent probability
  #Extremely likely: Over 95 per cent. According to the IPCC, it is now “extremely likely” that human activity is to blame for climate change
  #Very likely: Above 90 per cent
  #Likely: Above 66 per cent
  
  nbftstat <- nbft[, .(cmin = quantile(nbenefit, 0.5 / 3 ), # likely range - IPCC [66%]
                       cmax = quantile(nbenefit, 2.5 / 3), 
                       cmed = quantile(nbenefit, 0.5),      # median
                       c05  = quantile(nbenefit, 0.05),     # very likely range - IPCC [90%]
                       c95  = quantile(nbenefit, 0.95),
                       c01  = quantile(nbenefit, 0.01),     # 
                       c99  = quantile(nbenefit, 0.99),
                       mmax  = max(nbenefit),     # 
                       mmin  = min(nbenefit),
                       mmm  = mean(nbenefit)),              # mean
                   by = c("dmg_spec","qscen")]

  nbftstat[qscen == "05", qscen_label := "5th"]
  nbftstat[qscen == "50", qscen_label := "median"]
  nbftstat[qscen == "95", qscen_label := "95th"]
  nbftstat[, qscen_label := factor(qscen_label, levels = c("5th",
                                                           "median",
                                                           "95th"))]
  
  if (dr0 == 1) {
    breaks <- seq(-60000,20000, by = 5000)
  }  
  if (dr0 == 2) {
    breaks <- seq(-16000,5000, by = 2000)
  }
  if (dr0 == 3) {
    breaks <- seq(-5000,1500, by = 500)
  }
  if (dr0 == 4) {
    breaks <- seq(-600,400, by = 200)
  }
  
  d05 = density(nbft[qscen=="05"]$nbenefit * 1e2)
  d50 = density(nbft[qscen=="50"]$nbenefit * 1e2)
  d95 = density(nbft[qscen=="95"]$nbenefit * 1e2)
  scl05 <- max(d05[['y']]) #* length(d05[['y']])
  scl50 <- max(d50[['y']]) #* length(d50[['y']])
  scl95 <- max(d95[['y']]) #* length(d95[['y']])
  scl <- max(scl05,scl50,scl95)
  show_legend <- TRUE
  p <- ggplot(nbftstat) +
    geom_density(mapping = aes(x = nbenefit * 1e2, 
                               fill = qscen_label, 
                               color = qscen_label),
                 data = nbft, 
                 inherit.aes = FALSE, 
                 alpha = 0.2) +
    geom_jitter(mapping = aes(x = nbenefit * 1e2, 
                              y = -0.133 * scl - 0.66 * 0.5 / 3 * scl - as.numeric(qscen) / 100 * scl,
                              color = qscen_label),
                data = nbft,
                size = 0.1,
                alpha = 0.03, 
                height = 0.66 * 0.5 / 3 * scl) +
    geom_vline(xintercept = 0, color = "black") +
    geom_segment(aes(y = -0.10 * scl - as.numeric(qscen) / 100 * scl, 
                     yend = -0.40 * scl - as.numeric(qscen) / 100 * scl,
                     x = c05 * 1e2, xend = c05 * 1e2, color = qscen_label),
                 size = 0.8) +
    geom_segment(aes(y = -0.10 * scl - as.numeric(qscen) / 100 * scl, 
                     yend = -0.40 * scl - as.numeric(qscen) / 100 * scl,
                     x = c95 * 1e2, xend = c95 * 1e2, color = qscen_label),
                 size = 0.8) +
    geom_segment(aes(y = -0.25 * scl - as.numeric(qscen) / 100 * scl, 
                     yend = -0.25 * scl - as.numeric(qscen) / 100 * scl,
                     x = c05 * 1e2, xend = cmin * 1e2, color = qscen_label),
                 size = 0.8) +
    geom_segment(aes(y = -0.25 * scl - as.numeric(qscen) / 100 * scl,
                     yend = -0.25 * scl - as.numeric(qscen) / 100 * scl,
                     x = cmax * 1e2, xend = c95 * 1e2, color = qscen_label),
                 size = 0.8) +
    geom_rect(aes(ymin = -0.075 * scl - as.numeric(qscen) / 100 * scl, 
                  ymax = -0.425 * scl - as.numeric(qscen) / 100 * scl,
                  xmin = cmin * 1e2, xmax = cmax * 1e2, 
                  color = qscen_label),
                 fill = "white", alpha = 0.5, size = 0.8) +
    geom_segment(aes(y = -0.075 * scl - as.numeric(qscen) / 100 * scl, 
                     yend = -0.425 * scl - as.numeric(qscen) / 100 * scl,
                  x = cmed * 1e2, xend = cmed * 1e2, 
                  color = qscen_label),
              size = 0.8) +
    scale_x_continuous(labels = scales::percent_format(scale = 1), expand = c(0,0)) +
    scale_color_brewer(palette = "Set1",
                       guide = "none",
                       name = "Temperature\npercentile") +
    scale_fill_brewer(palette = "Set1",
                       guide = ifelse(show_legend,"legend","none"),
                       name = "Temperature\npercentile") +
    guides(fill = guide_legend(override.aes = list(alpha = 1, color = NA) ) ) +
    coord_cartesian(xlim = c(min(nbftstat$c01),max(nbftstat$c99)) * 1e2) +
    labs(x = str_glue("Policy benefits [%]"), y = "") +
    theme_light() +
    theme(axis.text.y = element_blank()) +
    theme(axis.ticks.y = element_blank())
  
    ggsave(output_file, plot = p, width = 8, height = 3, dpi = 200)
  
  output_file
  
}

plot_nben_boot_rcp_nocluster_main <- function(dr0,dmg_spec0,output_file, main = FALSE) {
  
  nben_npv_boot <- readd('nben_npv_boot')
  nben_npv_boot_p05 <- readd('nben_npv_boot_p05')
  nben_npv_boot_p95 <- readd('nben_npv_boot_p95')
  
  sr15c_runs <- readd('sr15c_runs')
  
  disc_rate = as.numeric(dr0) / 100
  
  nbft <- copy(nben_npv_boot[dr == disc_rate &
                               region == "World" &
                               dmg_spec == dmg_spec0 &
                               c5model == 1])
  nbft_p05 <- copy(nben_npv_boot_p05[dr == disc_rate &
                                       region == "World" &
                                       dmg_spec == dmg_spec0 &
                                       c5model == 1])
  nbft_p95 <- copy(nben_npv_boot_p95[dr == disc_rate &
                                       region == "World" &
                                       dmg_spec == dmg_spec0 &
                                       c5model == 1])
  nbft[, qscen := "50"]
  nbft_p05[, qscen := "05"]
  nbft_p95[, qscen := "95"]
  nbft <- rbindlist(list(nbft,nbft_p05,nbft_p95))
  nbft <- merge(nbft, sr15c_runs, by = c("model","scenario"))
  nbft <- nbft[,.(model,scenario,dmg_spec,qscen,nbenefit,temp_2100)]
  
  nbft <- nbft[temp_2100 <= 2]
  nbft[, temp_2100 := NULL]
  
  nbft[qscen == "05", qscen_label := "5th"]
  nbft[qscen == "50", qscen_label := "median"]
  nbft[qscen == "95", qscen_label := "95th"]
  nbft[, qscen_label := factor(qscen_label, levels = c("5th",
                                                       "median",
                                                       "95th"))]
  
  setkey(nbft,qscen_label,dmg_spec,qscen,model,scenario)
  
  # stats
  
  #Virtually certain: 99 to 100 per cent probability
  #Extremely likely: Over 95 per cent. According to the IPCC, it is now “extremely likely” that human activity is to blame for climate change
  #Very likely: Above 90 per cent
  #Likely: Above 66 per cent
  
  nbftstat <- nbft[, .(cmin = quantile(nbenefit, 0.5 / 3 ), # likely range - IPCC [66%]
                       cmax = quantile(nbenefit, 2.5 / 3), 
                       cmed = quantile(nbenefit, 0.5),      # median
                       c05  = quantile(nbenefit, 0.05),     # very likely range - IPCC [90%]
                       c95  = quantile(nbenefit, 0.95),
                       c01  = quantile(nbenefit, 0.01),     # 
                       c99  = quantile(nbenefit, 0.99),
                       mmax  = max(nbenefit),     # 
                       mmin  = min(nbenefit),
                       mmm  = mean(nbenefit)),              # mean
                   by = c("dmg_spec","qscen")]
  
  nbftstat[qscen == "05", qscen_label := "5th"]
  nbftstat[qscen == "50", qscen_label := "median"]
  nbftstat[qscen == "95", qscen_label := "95th"]
  nbftstat[, qscen_label := factor(qscen_label, levels = c("5th",
                                                           "median",
                                                           "95th"))]
  
  if (dr0 == 1) {
    breaks <- seq(-60000,20000, by = 5000)
  }  
  if (dr0 == 2) {
    breaks <- seq(-16000,5000, by = 2000)
  }
  if (dr0 == 3) {
    breaks <- seq(-5000,1500, by = 500)
  }
  if (dr0 == 4) {
    breaks <- seq(-600,400, by = 200)
  }
  
  d05 = density(nbft[qscen=="05"]$nbenefit * 1e2)
  d50 = density(nbft[qscen=="50"]$nbenefit * 1e2)
  d95 = density(nbft[qscen=="95"]$nbenefit * 1e2)
  scl05 <- max(d05[['y']]) #* length(d05[['y']])
  scl50 <- max(d50[['y']]) #* length(d50[['y']])
  scl95 <- max(d95[['y']]) #* length(d95[['y']])
  scl <- max(scl05,scl50,scl95)
  show_legend <- TRUE
  p <- ggplot(nbftstat) +
    geom_density(mapping = aes(x = nbenefit * 1e2, 
                               fill = qscen_label, 
                               color = qscen_label),
                 data = nbft, 
                 inherit.aes = FALSE, 
                 alpha = 0.2) +
    geom_jitter(mapping = aes(x = nbenefit * 1e2, 
                              y = -0.133 * scl - 0.66 * 0.5 / 3 * scl - as.numeric(qscen) / 100 * scl,
                              color = qscen_label),
                data = nbft,
                size = 0.1,
                alpha = 0.03, 
                height = 0.66 * 0.5 / 3 * scl) +
    geom_vline(xintercept = 0, color = "black") +
    geom_segment(aes(y = -0.10 * scl - as.numeric(qscen) / 100 * scl, 
                     yend = -0.40 * scl - as.numeric(qscen) / 100 * scl,
                     x = c05 * 1e2, xend = c05 * 1e2, color = qscen_label),
                 size = 0.8) +
    geom_segment(aes(y = -0.10 * scl - as.numeric(qscen) / 100 * scl, 
                     yend = -0.40 * scl - as.numeric(qscen) / 100 * scl,
                     x = c95 * 1e2, xend = c95 * 1e2, color = qscen_label),
                 size = 0.8) +
    geom_segment(aes(y = -0.25 * scl - as.numeric(qscen) / 100 * scl, 
                     yend = -0.25 * scl - as.numeric(qscen) / 100 * scl,
                     x = c05 * 1e2, xend = cmin * 1e2, color = qscen_label),
                 size = 0.8) +
    geom_segment(aes(y = -0.25 * scl - as.numeric(qscen) / 100 * scl,
                     yend = -0.25 * scl - as.numeric(qscen) / 100 * scl,
                     x = cmax * 1e2, xend = c95 * 1e2, color = qscen_label),
                 size = 0.8) +
    geom_rect(aes(ymin = -0.075 * scl - as.numeric(qscen) / 100 * scl, 
                  ymax = -0.425 * scl - as.numeric(qscen) / 100 * scl,
                  xmin = cmin * 1e2, xmax = cmax * 1e2, 
                  color = qscen_label),
              fill = "white", alpha = 0.5, size = 0.8) +
    geom_segment(aes(y = -0.075 * scl - as.numeric(qscen) / 100 * scl, 
                     yend = -0.425 * scl - as.numeric(qscen) / 100 * scl,
                     x = cmed * 1e2, xend = cmed * 1e2, 
                     color = qscen_label),
                 size = 0.8) +
    scale_x_continuous(labels = scales::percent_format(scale = 1), expand = c(0,0)) +
    scale_color_brewer(palette = "Set1",
                       guide = "none",
                       name = "Temperature\npercentile") +
    scale_fill_brewer(palette = "Set1",
                      guide = ifelse(show_legend,"legend","none"),
                      name = "Temperature\npercentile") +
    guides(fill = guide_legend(override.aes = list(alpha = 1, color = NA) ) ) +
    coord_cartesian(xlim = c(min(nbftstat$c01),max(nbftstat$c99)) * 1e2) +
    labs(x = str_glue("Net present policy benefits [%]"), y = "") +
    theme_light() +
    theme(axis.text.y = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(legend.position = c(0.15,0.82))

    ggsave(output_file, plot = p, width = 5.5, height = 4.5, dpi = 200)

  output_file
  
}

plot_nben_boot_temp <- function(dr0,dmg_spec0,output_file) {
  
  nben_npv_boot <- readd('nben_npv_boot')
  sr15c_runs <- readd('sr15c_runs')
  
  disc_rate = as.numeric(dr0) / 100
  
  nbft <- copy(nben_npv_boot[dr == disc_rate &
                               region == "World" &
                               dmg_spec == dmg_spec0 &
                               c5model == 1])
  nbft <- merge(nbft, sr15c_runs, by = c("model","scenario"))
  nbft <- nbft[,.(model,scenario,dmg_spec,nbenefit,temp_2100)]

  # stats
  nbftstat <- nbft[, .(cmin = quantile(nbenefit, 0.5 / 3 ), # likely range - IPCC
                       cmax = quantile(nbenefit, 2.5 / 3), 
                       cmed = quantile(nbenefit, 0.5),      # median
                       c01  = quantile(nbenefit, 0.01),     # 
                       c99  = quantile(nbenefit, 0.99),
                       mmm  = mean(nbenefit)),              # mean
                   by = c("model","scenario")]
  nbftstat <- merge(nbftstat, sr15c_runs, by = c("model","scenario"))
  
  p <- ggplot(nbftstat, aes(x = temp_2100)) +
            geom_hline(yintercept = 0, color = "black") +
            geom_linerange(aes(ymin = cmin * 1e2, ymax = cmax * 1e2,
                               color = cmed * 1e2,
                               group = paste(model,scenario)),
                           size = 1) +
    geom_point(aes(y = cmed * 1e2,
                   group = paste(model,scenario)), 
               size = 2) +
    geom_point(aes(y = cmed * 1e2,
                   group = paste(model,scenario)), 
               color = "white",
               size = 1) +
    geom_smooth(aes(y = cmed * 1e2), 
                color = muted('green'),
                data = nbftstat, 
                method = 'glm',
                se = TRUE) +
    scale_y_continuous(labels = scales::percent_format(scale = 1), expand = c(0,0)) +
    #scale_y_continuous(labels = comma, expand = c(0,0)) +
    scale_x_continuous(breaks = c(1.25,1.5,1.75,2,2.25), limits = c(1,2.5)) +
    #scale_color_gradient2(name = "Net benefits\nmedian [T$]") +
    scale_color_gradient2(name = "Policy benefits\nmedian [%]") +
    labs(
      #y = "Net Present Value [Trillion USD2018]", 
      y = "Policy benefits [%]", 
      x = "Global Mean Temperature in 2100, median estimate [C]",
      #title = "Policy benefits",
      subtitle = str_glue("{round(dr0)}% discount rate")
    ) +
    theme_minimal() +
    theme(
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          #legend.position = "right",
          legend.text = element_text(size = 8),
          axis.text = element_text(size = 11),
          axis.line.x = element_blank())
  
  ggsave(output_file, plot = p, width = 10, height = 5)
    
  output_file

}

plot_nben_boot_ann <- function(output_file, data_file) {
  
  sr15c_runs <- readd(sr15c_runs)
  nben <- readd(data_file, character_only = TRUE)
  
  # Select scenarios and add info
  scen_ref <- unique(sr15c_runs$scenario_ref)
  
  nben <- nben[!scenario %in% scen_ref &
               !is.na(nbenefit) &
               region == "World" &
               scenario != "GEA_Mix_1p5C_AdvNCO2_PartialDelay2020" & # wrong value in 2020
               c5model == 1,
               .(dmg_spec,c5model,model,scenario,year,nbenefit,runid)]
  
  nben <- merge(nben, sr15c_runs, by = c("model","scenario"))
  nben[, c("scenario_ref","dmg_spec","rcp_clus") := NULL]
  
  nben <- nben[temp_2100 <= 2]
  nben[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nben[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nben[, temp_2100 := NULL]

  nben_stat <- nben[, .(cmed = quantile(nbenefit,0.5),
                        cmean = mean(nbenefit),
                        gt0 = length(nbenefit[nbenefit>0])/length(nbenefit)), 
                      by = "temp_clus,year"]
  
  allcomb <- unique(nben[,paste(c5model,runid)])
  
  # select randomly 1000 
  idx <- sample(1:length(allcomb), min(length(allcomb),1000))

  nben_smp <- nben[paste(c5model,runid) %in% allcomb[idx]]
  
  interp_spline <- function(SD, step = 2) {
    
    return (list(year = seq(2020,2100,by = step),
                nbenefit = spline(x = SD$year,
                                  y = SD$nbenefit,
                                  xout = seq(2020,2100,by = step))[['y']]
    ))
    
  }
  interp_spline2 <- function(SD, step = 2) {
    
    return (list(year = seq(2020,2100,by = step),
                 gt0 = spline(x = SD$year,
                                   y = SD$gt0,
                                   xout = seq(2020,2100,by = step))[['y']]
    ))
    
  }
  
  nben_smp2 <- nben_smp[, interp_spline(.SD), 
                        by = "model,scenario,c5model,runid,temp_clus"]

  nben_stat2 <- nben_stat[,.(temp_clus,year,nbenefit = cmed)]
  nben_stat2 <- nben_stat2[,interp_spline(.SD,1), by = "temp_clus"]
  nben_stat2 <- nben_stat2[, .(year,nbenefit,nbenefitp1 = shift(nbenefit, -1L)), by = "temp_clus"]
  nben_stat2 <- nben_stat2[year != 2100]
  
  
  nben_stat3 <- nben_stat[,.(temp_clus,year,gt0)]
  nben_stat3 <- nben_stat3[,interp_spline2(.SD,1), by = "temp_clus"]
  
  nben_stat00 <- merge(nben_stat2,nben_stat3, by = c("temp_clus","year"))

  
  half_line <- 11 / 2
  
  p <- ggplot(data = nben_smp2) +
    geom_line(aes(x = year, y = nbenefit * 1e2,
                  group = paste(model,c5model,scenario,runid)),
              alpha = 0.01, size = 0.5) +
    geom_hline(yintercept = 0, size = 1) +
    geom_segment(aes(x = year, xend = year + 1,
                   y = nbenefit * 1e2, yend = nbenefitp1 * 1e2),
               size = 1.6,
               color = "white",
               lineend = "round",
               linejoin = "bevel",
               alpha = 1,
               data = nben_stat00) +
    geom_segment(aes(x = year, xend = year + 1,
                     y = nbenefit * 1e2, yend = nbenefitp1 * 1e2,
                     color = gt0 * 100),
                 size = 1.5,
                 lineend = "round",
                 linejoin = "bevel",
                 alpha = 1,
                 data = nben_stat00) +
    geom_tile(aes(x = year, y = -50, fill = gt0 * 100),
              height = 12,
              colour = "grey",
              data = nben_stat[year > 2020]) +
    geom_text(aes(x = year, y = -50,
                  label = paste(round(gt0 * 100))),
              colour = "darkgray",
              size = 3,
              data = nben_stat[year > 2020 & gt0 <= 0.25]) +
    geom_text(aes(x = year, y = -50,
                label = paste(round(gt0 * 100))),
            colour = "black",
            size = 3,
            data = nben_stat[year > 2020 & gt0 > 0.25]) +
    scale_color_viridis_c(name = "% scenarios > 0", 
                       guide = "colorbar",
                       limits = c(NA,100)) +
    scale_fill_viridis_c(name = "[%] > 0", 
                         guide = "none") +
    scale_y_continuous(labels = scales::percent_format(scale = 1), expand = c(0,0)) +
    #scale_y_continuous(breaks = seq(-500,500,by = 100), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(2020,2100,by = 10),
                       labels = c("2020","","2040","","2060","","2080","","2100")) +
    coord_cartesian(ylim = c(-57,100), 
                    xlim = c(2015,2105)) +
    facet_wrap( ~ temp_clus) +
    #labs(y = "Net benefits [T$/yr]", x = "") +
    labs(y = "Policy benefits [%]", x = "") +
    theme_light() +
    guides(color = guide_colourbar(title.position="bottom", title.hjust = 0.5)) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "right",
      legend.box.margin = margin(215),
      axis.text = element_text(size = 11),
      axis.line.x = element_blank(),
      strip.text = element_text(size = 12, color = "black"),
      strip.background = element_blank())

  ggsave(output_file, plot = p, width = 8, height = 5, dpi = 120)
  
  data.table(file = output_file)
  
}

plot_nben_boot_ann_nocluster <- function(output_file, data_file) {
  
  sr15c_runs <- readd(sr15c_runs)
  nben <- readd(data_file, character_only = TRUE)
  
  # Select scenarios and add info
  scen_ref <- unique(sr15c_runs$scenario_ref)
  
  nben <- nben[!scenario %in% scen_ref &
                 !is.na(nbenefit) &
                 region == "World" &
                 scenario != "GEA_Mix_1p5C_AdvNCO2_PartialDelay2020" & # wrong value in 2020
                 c5model == 1,
               .(dmg_spec,c5model,model,scenario,year,nbenefit,runid)]
  
  nben <- merge(nben, sr15c_runs, by = c("model","scenario"))
  nben[, c("scenario_ref","dmg_spec","rcp_clus") := NULL]
  
  nben <- nben[temp_2100 <= 2]
  nben[, temp_2100 := NULL]
  
  nben_stat <- nben[, .(cmed = quantile(nbenefit,0.5),
                        cmean = mean(nbenefit),
                        gt0 = length(nbenefit[nbenefit>0])/length(nbenefit)), 
                    by = "year"]
  
  allcomb <- unique(nben[,paste(c5model,runid)])
  
  # select randomly 1000 
  idx <- sample(1:length(allcomb), min(length(allcomb),1000))
  
  nben_smp <- nben[paste(c5model,runid) %in% allcomb[idx]]
  
  interp_spline <- function(SD, step = 2) {
    
    return (list(year = seq(2020,2100,by = step),
                 nbenefit = spline(x = SD$year,
                                   y = SD$nbenefit,
                                   xout = seq(2020,2100,by = step))[['y']]
    ))
    
  }
  interp_spline2 <- function(SD, step = 2) {
    
    return (list(year = seq(2020,2100,by = step),
                 gt0 = spline(x = SD$year,
                              y = SD$gt0,
                              xout = seq(2020,2100,by = step))[['y']]
    ))
    
  }
  
  nben_smp2 <- nben_smp[, interp_spline(.SD), 
                        by = "model,scenario,c5model,runid"]
  
  nben_stat2 <- nben_stat[,.(year,nbenefit = cmed)]
  nben_stat2 <- nben_stat2[,interp_spline(.SD,1)]
  nben_stat2 <- nben_stat2[, .(year,nbenefit,nbenefitp1 = shift(nbenefit, -1L))]
  nben_stat2 <- nben_stat2[year != 2100]
  
  
  nben_stat3 <- nben_stat[,.(year,gt0)]
  nben_stat3 <- nben_stat3[,interp_spline2(.SD,1)]
  
  nben_stat00 <- merge(nben_stat2,nben_stat3, by = c("year"))
  
  
  half_line <- 11 / 2
  
  p <- ggplot(data = nben_smp2) +
    geom_line(aes(x = year, y = nbenefit * 1e2,
                  group = paste(model,c5model,scenario,runid)),
              alpha = 0.01, size = 0.5) +
    geom_hline(yintercept = 0, size = 1) +
    geom_segment(aes(x = year, xend = year + 1,
                     y = nbenefit * 1e2, yend = nbenefitp1 * 1e2),
                 size = 1.6,
                 color = "white",
                 lineend = "round",
                 linejoin = "bevel",
                 alpha = 1,
                 data = nben_stat00) +
    geom_segment(aes(x = year, xend = year + 1,
                     y = nbenefit * 1e2, yend = nbenefitp1 * 1e2,
                     color = gt0 * 100),
                 size = 1.5,
                 lineend = "round",
                 linejoin = "bevel",
                 alpha = 1,
                 data = nben_stat00) +
    geom_tile(aes(x = year, y = -50, fill = gt0 * 100),
              height = 12,
              colour = "grey",
              data = nben_stat[year > 2020]) +
    geom_text(aes(x = year, y = -50, 
                  label = paste(round(gt0 * 100),"%")),
              colour = "darkgray",
              size = 3,
              data = nben_stat[year > 2020 & gt0 <= 0.25]) +
    geom_text(aes(x = year, y = -50, 
                  label = paste(round(gt0 * 100),"%")),
              colour = "black",
              size = 3,
              data = nben_stat[year > 2020 & gt0 > 0.25]) +
    scale_color_viridis_c(name = "% scenarios > 0", 
                          guide = "colorbar",
                          limits = c(NA,100)) +
    scale_fill_viridis_c(name = "[%] > 0", 
                         guide = "none") +
    scale_y_continuous(labels = scales::percent_format(scale = 1), expand = c(0,0)) +
#    scale_y_continuous(breaks = seq(-500,500,by = 100), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(2020,2100,by = 10),
                       labels = c("2020","2030","2040","2050","2060",
                                  "2070","2080","2090","2100"),
                       expand = c(0,0)) +
    coord_cartesian(ylim = c(-57,100), 
                    xlim = c(2015, 2105)) +
#    labs(y = "Net benefits [T$/yr]", x = "") +
    labs(y = "Policy benefits [%]", x = "") +
    theme_light() +
    guides(color = guide_colourbar(title.position="bottom", title.hjust = 0.5)) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = c(0,1),
      legend.justification = c(-0.1,1.1),
      axis.text = element_text(size = 11),
      axis.line.x = element_blank(),
      strip.text = element_text(size = 12, color = "black"),
      strip.background = element_blank())
  
  ggsave(output_file, plot = p, width = 6, height = 5, dpi = 200)
  
  data.table(file = output_file)
  
}

plot_nben_rcp_v2_x <- function(dr0,output_file,show_legend = TRUE) {
  
  nben_npv_gwt <- readd('nben_npv_gwt')
  nben_npv_gwt_p05 <- readd('nben_npv_gwt_p05')
  nben_npv_gwt_p95 <- readd('nben_npv_gwt_p95')
  nben_npv_lvl <- readd('nben_npv_lvl')
  sr15c_runs <- readd('sr15c_runs')
  
  disc_rate = as.numeric(dr0) / 100
  
  nbft1 <- rbindlist(list(nben_npv_gwt[dr == disc_rate & region == "World" & c5model == 1]))
  
  nbft1 <- merge(nbft1, sr15c_runs, by = c("model","scenario"))
  nbft1 <- nbft1[,.(model,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft1 <- nbft1[temp_2100 <= 2]
  nbft1[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft1[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft1[, temp_2100 := NULL]
  nbft1[, dmg_type := "gwt"]
  nbft1[str_detect(dmg_spec, "LVL"), dmg_type := "lvl"]
  
  setkey(nbft1,dmg_spec,model,scenario,temp_clus)
  
  nbftstat1 <- rbindlist(list(
    nbft1[, .(value = quantile(nbenefit,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1[, .(value = quantile(nbenefit,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1[, .(value = quantile(nbenefit,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat1[, qtemp := "50"]
  
  nbft1_p05 <- rbindlist(list(nben_npv_gwt_p05[dr == disc_rate & region == "World" & c5model == 1]))
  
  nbft1_p05 <- merge(nbft1_p05, sr15c_runs, by = c("model","scenario"))
  nbft1_p05 <- nbft1_p05[,.(model,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft1_p05 <- nbft1_p05[temp_2100 <= 2]
  nbft1_p05[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft1_p05[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft1_p05[, temp_2100 := NULL]
  nbft1_p05[, dmg_type := "gwt"]
  nbft1_p05[str_detect(dmg_spec, "LVL"), dmg_type := "lvl"]
  
  setkey(nbft1_p05,dmg_spec,model,scenario,temp_clus)
  
  nbftstat1_p05 <- rbindlist(list(
    nbft1_p05[, .(value = quantile(nbenefit,0.5), qscen = "50"),
              by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1_p05[, .(value = quantile(nbenefit,0.05), qscen = "05"),
              by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1_p05[, .(value = quantile(nbenefit,0.95), qscen = "95"),
              by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat1_p05[, qtemp := "05"]
  
  nbft1_p95 <- rbindlist(list(nben_npv_gwt_p95[dr == disc_rate & region == "World" & c5model == 1]))
  
  nbft1_p95 <- merge(nbft1_p95, sr15c_runs, by = c("model","scenario"))
  nbft1_p95 <- nbft1_p95[,.(model,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft1_p95 <- nbft1_p95[temp_2100 <= 2]
  nbft1_p95[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft1_p95[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft1_p95[, temp_2100 := NULL]
  nbft1_p95[, dmg_type := "gwt"]
  nbft1_p95[str_detect(dmg_spec, "LVL"), dmg_type := "lvl"]
  
  setkey(nbft1_p95,dmg_spec,model,scenario,temp_clus)
  
  nbftstat1_p95 <- nbft1_p95[, .(cmax = quantile(nbenefit,0.5),
                                 cmax_05 = quantile(nbenefit,0.05),
                                 cmax_95 = quantile(nbenefit,0.95)),
                             by = c("temp_clus","dmg_spec","dmg_type")]
  nbftstat1_p95 <- rbindlist(list(
    nbft1_p95[, .(value = quantile(nbenefit,0.5), qscen = "50"),
              by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1_p95[, .(value = quantile(nbenefit,0.05), qscen = "05"),
              by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1_p95[, .(value = quantile(nbenefit,0.95), qscen = "95"),
              by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat1_p95[, qtemp := "95"]
  
  nbftstat1 <- rbindlist(list(nbftstat1, nbftstat1_p05, nbftstat1_p95))
  
  nbft2 <- copy(nben_npv_lvl[dr == disc_rate &
                               region == "World"])
  
  nbft2 <- merge(nbft2, sr15c_runs, by = c("model","scenario"))
  
  nbft2 <- nbft2[,.(model,c5model,sw,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft2 <- nbft2[temp_2100 <= 2]
  nbft2[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft2[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft2[, temp_2100 := NULL]
  nbft2[, dmg_type := "lvl"]
  
  setkey(nbft2,dmg_spec,c5model,sw,model,scenario,temp_clus)
  
  nbftstat2 <- rbindlist(list(
    nbft2[sw == 1, .(value = quantile(nbenefit,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 1, .(value = quantile(nbenefit,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 1, .(value = quantile(nbenefit,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat2[, qtemp := "50"]
  nbftstat2_p05 <- rbindlist(list(
    nbft2[sw == 2, .(value = quantile(nbenefit,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 2, .(value = quantile(nbenefit,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 2, .(value = quantile(nbenefit,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat2_p05[, qtemp := "05"]
  nbftstat2_p95 <- rbindlist(list(
    nbft2[sw == 3, .(value = quantile(nbenefit,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 3, .(value = quantile(nbenefit,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 3, .(value = quantile(nbenefit,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat2_p95[, qtemp := "95"]
  
  nbftstat2 <- rbindlist(list(nbftstat2, nbftstat2_p05, nbftstat2_p95))
  
  nbftstat <- rbind(nbftstat1,nbftstat2)
  
  
  # label dmg spec
  
  dmg_spec_val <- data.table(dmg_spec = c(bhm_spec,djo_spec,lvl_spec),
                             dmg_spec_label = c(bhm_spec_label,djo_spec_label,lvl_spec_label),
                             dmg_spec_star = c(bhm_spec_star,djo_spec_star,lvl_spec_star))
  
  nbftstatx <- merge(nbftstat,dmg_spec_val,by = c("dmg_spec"))
  
  nbftstatx <- dcast(nbftstatx, dmg_spec + temp_clus + dmg_type + qtemp + dmg_spec_label + dmg_spec_star ~ qscen)
  
  dmg_spec_order <- nbftstatx[qtemp == "50",mean(`50`),by = "dmg_spec_label"][order(V1), dmg_spec_label]
  
  nbftstatx[, dmg_spec_label := factor(dmg_spec_label, 
                                       levels = dmg_spec_order)]
  
  nbftstatx[, dmg_type := factor(dmg_type, 
                                 levels = c("gwt","lvl"),
                                 labels = c("growth-based","level-based"))]
  
  wdd <- 0.33
  
  nbftstatx[qtemp == "50", size_point := 4]
  nbftstatx[qtemp != "50", size_point := 3]
  
  nbftstatx[qtemp == "50", size_line := 1]
  nbftstatx[qtemp != "50", size_line := 0.8]
  
  
  nbftstatx[qtemp == "05", qscen_label := "5th"]
  nbftstatx[qtemp == "50", qscen_label := "median"]
  nbftstatx[qtemp == "95", qscen_label := "95th"]
  nbftstatx[, qscen_label := factor(qscen_label, levels = c("5th",
                                                            "median",
                                                            "95th"))]
  
  p <- ggplot(nbftstatx) +
    geom_hline(yintercept = 0, color = "black")  +
    geom_linerange(aes(x = dmg_spec_label,
                       ymin = `05` * 1e2,
                       ymax = `95` * 1e2,
                       color = temp_clus,
                       alpha = dmg_spec_star,
                       group = temp_clus,
                       size = size_line),
                   position = position_dodge(width = wdd)) +
    geom_point(aes(x = dmg_spec_label,
                   y = `50` * 1e2,
                   shape = dmg_type,
                   alpha = dmg_spec_star,
                   group = temp_clus,
                   size = size_point),
               position = position_dodge(width = wdd)) +
    geom_point(aes(x = dmg_spec_label,
                   y = `50` * 1e2,
                   color = temp_clus,
                   shape = dmg_type,
                   alpha = dmg_spec_star,
                   group = temp_clus,
                   size = size_point - 1),
               position = position_dodge(width = wdd)) +
    scale_size_identity() +
    scale_color_brewer(palette = "Set2",
                       guide = ifelse(show_legend,"legend","none"),
                       name = "Temperature\ncluster") +
    scale_alpha_continuous(range = c(0.25,1), guide = "none") +
    scale_shape_discrete(name = "Damage function", guide = ifelse(show_legend,"legend","none")) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    facet_wrap(~ qscen_label) +
    coord_flip() +
    labs(
      #y = "Net Present Value [Trillion USD2018]", x = ""
      y = "Policy benefits [%]", x = ""
    ) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "right",
      axis.text = element_text(size = 11),
      axis.line.x = element_blank(),
      strip.text = element_text(size = 12, color = "black"),
      strip.background = element_blank())
  
  ggsave(output_file, plot = p, width = 14, height = 7, dpi = 120)
  
  data.table(file = output_file)
  
}

plot_nben_rcp_v2_hs2017 <- function(dr0,output_file,show_legend = TRUE) {
  
  nben_npv_gwt <- readd('nben_npv_gwt')
  nben_npv_gwt_p05 <- readd('nben_npv_gwt_p05')
  nben_npv_gwt_p95 <- readd('nben_npv_gwt_p95')
  nben_npv_lvl <- readd('nben_npv_lvl')
  sr15c_runs <- readd('sr15c_runs')
  
  disc_rate = as.numeric(dr0) / 100
  
  nbft1 <- rbindlist(list(nben_npv_gwt[dr == disc_rate & region == "World" & c5model == 1]))
  
  nbft1 <- merge(nbft1, sr15c_runs, by = c("model","scenario"))
  nbft1 <- nbft1[,.(model,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft1 <- nbft1[temp_2100 <= 2]
  nbft1[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft1[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft1[, temp_2100 := NULL]
  nbft1[, dmg_type := "gwt"]
  nbft1[str_detect(dmg_spec, "LVL"), dmg_type := "lvl"]
  
  setkey(nbft1,dmg_spec,model,scenario,temp_clus)
  
  nbftstat1 <- rbindlist(list(
    nbft1[, .(value = quantile(nbenefit,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1[, .(value = quantile(nbenefit,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1[, .(value = quantile(nbenefit,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat1[, qtemp := "50"]
  
  nbft1_p05 <- rbindlist(list(nben_npv_gwt_p05[dr == disc_rate & region == "World" & c5model == 1]))
  
  nbft1_p05 <- merge(nbft1_p05, sr15c_runs, by = c("model","scenario"))
  nbft1_p05 <- nbft1_p05[,.(model,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft1_p05 <- nbft1_p05[temp_2100 <= 2]
  nbft1_p05[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft1_p05[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft1_p05[, temp_2100 := NULL]
  nbft1_p05[, dmg_type := "gwt"]
  nbft1_p05[str_detect(dmg_spec, "LVL"), dmg_type := "lvl"]
  
  setkey(nbft1_p05,dmg_spec,model,scenario,temp_clus)
  
  nbftstat1_p05 <- rbindlist(list(
    nbft1_p05[, .(value = quantile(nbenefit,0.5), qscen = "50"),
              by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1_p05[, .(value = quantile(nbenefit,0.05), qscen = "05"),
              by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1_p05[, .(value = quantile(nbenefit,0.95), qscen = "95"),
              by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat1_p05[, qtemp := "05"]
  
  nbft1_p95 <- rbindlist(list(nben_npv_gwt_p95[dr == disc_rate & region == "World" & c5model == 1]))
  
  nbft1_p95 <- merge(nbft1_p95, sr15c_runs, by = c("model","scenario"))
  nbft1_p95 <- nbft1_p95[,.(model,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft1_p95 <- nbft1_p95[temp_2100 <= 2]
  nbft1_p95[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft1_p95[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft1_p95[, temp_2100 := NULL]
  nbft1_p95[, dmg_type := "gwt"]
  nbft1_p95[str_detect(dmg_spec, "LVL"), dmg_type := "lvl"]
  
  setkey(nbft1_p95,dmg_spec,model,scenario,temp_clus)
  
  nbftstat1_p95 <- nbft1_p95[, .(cmax = quantile(nbenefit,0.5),
                                 cmax_05 = quantile(nbenefit,0.05),
                                 cmax_95 = quantile(nbenefit,0.95)),
                             by = c("temp_clus","dmg_spec","dmg_type")]
  nbftstat1_p95 <- rbindlist(list(
    nbft1_p95[, .(value = quantile(nbenefit,0.5), qscen = "50"),
              by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1_p95[, .(value = quantile(nbenefit,0.05), qscen = "05"),
              by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1_p95[, .(value = quantile(nbenefit,0.95), qscen = "95"),
              by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat1_p95[, qtemp := "95"]
  
  nbftstat1 <- rbindlist(list(nbftstat1, nbftstat1_p05, nbftstat1_p95))
  
  nbft2 <- copy(nben_npv_lvl[dr == disc_rate &
                               region == "World"])
  
  nbft2 <- merge(nbft2, sr15c_runs, by = c("model","scenario"))
  
  nbft2 <- nbft2[,.(model,c5model,sw,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft2 <- nbft2[temp_2100 <= 2]
  nbft2[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft2[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft2[, temp_2100 := NULL]
  nbft2[, dmg_type := "lvl"]
  
  setkey(nbft2,dmg_spec,c5model,sw,model,scenario,temp_clus)
  
  nbftstat2 <- rbindlist(list(
    nbft2[sw == 1, .(value = quantile(nbenefit,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 1, .(value = quantile(nbenefit,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 1, .(value = quantile(nbenefit,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat2[, qtemp := "50"]
  nbftstat2_p05 <- rbindlist(list(
    nbft2[sw == 2, .(value = quantile(nbenefit,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 2, .(value = quantile(nbenefit,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 2, .(value = quantile(nbenefit,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat2_p05[, qtemp := "05"]
  nbftstat2_p95 <- rbindlist(list(
    nbft2[sw == 3, .(value = quantile(nbenefit,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 3, .(value = quantile(nbenefit,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 3, .(value = quantile(nbenefit,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat2_p95[, qtemp := "95"]
  
  nbftstat2 <- rbindlist(list(nbftstat2, nbftstat2_p05, nbftstat2_p95))
  
  nbftstat <- rbind(nbftstat1,nbftstat2)
  
  
  # label dmg spec
  
  dmg_spec_val <- data.table(dmg_spec = c(bhm_spec,djo_spec,lvl_spec),
                             dmg_spec_label = c(bhm_spec_label,djo_spec_label,lvl_spec_label),
                             dmg_spec_star = c(bhm_spec_star,djo_spec_star,lvl_spec_star))
  
  nbftstatx <- merge(nbftstat,dmg_spec_val,by = c("dmg_spec"))
  
  nbftstatx <- dcast(nbftstatx, dmg_spec + temp_clus + dmg_type + qtemp + dmg_spec_label + dmg_spec_star ~ qscen)
  
  dmg_spec_order <- c("HS MKT","HS NCAT","HS TOT","HS TOT+P")
  
  nbftstatx[, dmg_spec_label := factor(dmg_spec_label, 
                                       levels = dmg_spec_order)]
  
  nbftstatx[, dmg_type := factor(dmg_type, 
                                 levels = c("gwt","lvl"),
                                 labels = c("growth-based","level-based"))]
  
  wdd <- 0.33
  
  nbftstatx[qtemp == "50", size_point := 4]
  nbftstatx[qtemp != "50", size_point := 3]
  
  nbftstatx[qtemp == "50", size_line := 1]
  nbftstatx[qtemp != "50", size_line := 0.8]
  
  
  nbftstatx[qtemp == "05", qscen_label := "5th"]
  nbftstatx[qtemp == "50", qscen_label := "median"]
  nbftstatx[qtemp == "95", qscen_label := "95th"]
  nbftstatx[, qscen_label := factor(qscen_label, levels = c("5th",
                                                            "median",
                                                            "95th"))]
  
  p <- ggplot(nbftstatx[dmg_spec_label %in% dmg_spec_order]) +
    geom_hline(yintercept = 0, color = "black")  +
    geom_linerange(aes(x = dmg_spec_label,
                       ymin = `05` * 1e2,
                       ymax = `95` * 1e2,
                       color = qscen_label,
                       alpha = dmg_spec_star,
                       group = qtemp,
                       size = size_line),
                   position = position_dodge(width = wdd)) +
    geom_point(aes(x = dmg_spec_label,
                   y = `50` * 1e2,
                   shape = dmg_type,
                   alpha = dmg_spec_star,
                   group = qtemp,
                   size = size_point),
               position = position_dodge(width = wdd)) +
    geom_point(aes(x = dmg_spec_label,
                   y = `50` * 1e2,
                   color = qscen_label,
                   shape = dmg_type,
                   alpha = dmg_spec_star,
                   group = qtemp,
                   size = size_point - 1),
               position = position_dodge(width = wdd)) +
    scale_size_identity() +
    scale_color_brewer(palette = "Set1",
                       name = "Temperature\npercentile") +
    scale_alpha_continuous(range = c(0.25,1), guide = "none") +
    scale_shape_discrete(name = "Damage function", guide = "none") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    facet_wrap(~ temp_clus) +
    coord_flip() +
    labs(
      #y = "Net Present Value [Trillion USD2018]", x = ""
      y = "Policy benefits [%]", x = ""
    ) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "right",
      axis.text = element_text(size = 11),
      axis.line.x = element_blank(),
      strip.text = element_text(size = 12, color = "black"),
      strip.background = element_blank())
  
  ggsave(output_file, plot = p, width = 10, height = 4, dpi = 120)
  
  data.table(file = output_file)
  
}

plot_nben_rcp_v2_med <- function(dr0,output_file,show_legend = TRUE) {
  
  nben_npv_gwt <- readd('nben_npv_gwt')
  nben_npv_gwt_p05 <- readd('nben_npv_gwt_p05')
  nben_npv_gwt_p95 <- readd('nben_npv_gwt_p95')
  nben_npv_lvl <- readd('nben_npv_lvl')
  sr15c_runs <- readd('sr15c_runs')
  
  disc_rate = as.numeric(dr0) / 100
  
  nbft1 <- rbindlist(list(nben_npv_gwt[dr == disc_rate & region == "World" & c5model == 1]))
  
  nbft1 <- merge(nbft1, sr15c_runs, by = c("model","scenario"))
  nbft1 <- nbft1[,.(model,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft1 <- nbft1[temp_2100 <= 2]
  nbft1[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft1[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft1[, temp_2100 := NULL]
  nbft1[, dmg_type := "gwt"]
  nbft1[str_detect(dmg_spec, "LVL"), dmg_type := "lvl"]
  
  setkey(nbft1,dmg_spec,model,scenario,temp_clus)
  
  nbftstat1 <- rbindlist(list(
    nbft1[, .(value = quantile(nbenefit,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1[, .(value = quantile(nbenefit,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1[, .(value = quantile(nbenefit,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat1[, qtemp := "50"]
  
  nbft1_p05 <- rbindlist(list(nben_npv_gwt_p05[dr == disc_rate & region == "World" & c5model == 1]))
  
  nbft1_p05 <- merge(nbft1_p05, sr15c_runs, by = c("model","scenario"))
  nbft1_p05 <- nbft1_p05[,.(model,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft1_p05 <- nbft1_p05[temp_2100 <= 2]
  nbft1_p05[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft1_p05[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft1_p05[, temp_2100 := NULL]
  nbft1_p05[, dmg_type := "gwt"]
  nbft1_p05[str_detect(dmg_spec, "LVL"), dmg_type := "lvl"]
  
  setkey(nbft1_p05,dmg_spec,model,scenario,temp_clus)
  
  nbftstat1_p05 <- rbindlist(list(
    nbft1_p05[, .(value = quantile(nbenefit,0.5), qscen = "50"),
              by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1_p05[, .(value = quantile(nbenefit,0.05), qscen = "05"),
              by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1_p05[, .(value = quantile(nbenefit,0.95), qscen = "95"),
              by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat1_p05[, qtemp := "05"]
  
  nbft1_p95 <- rbindlist(list(nben_npv_gwt_p95[dr == disc_rate & region == "World" & c5model == 1]))
  
  nbft1_p95 <- merge(nbft1_p95, sr15c_runs, by = c("model","scenario"))
  nbft1_p95 <- nbft1_p95[,.(model,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft1_p95 <- nbft1_p95[temp_2100 <= 2]
  nbft1_p95[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft1_p95[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft1_p95[, temp_2100 := NULL]
  nbft1_p95[, dmg_type := "gwt"]
  nbft1_p95[str_detect(dmg_spec, "LVL"), dmg_type := "lvl"]
  
  setkey(nbft1_p95,dmg_spec,model,scenario,temp_clus)
  
  nbftstat1_p95 <- nbft1_p95[, .(cmax = quantile(nbenefit,0.5),
                                 cmax_05 = quantile(nbenefit,0.05),
                                 cmax_95 = quantile(nbenefit,0.95)),
                             by = c("temp_clus","dmg_spec","dmg_type")]
  nbftstat1_p95 <- rbindlist(list(
    nbft1_p95[, .(value = quantile(nbenefit,0.5), qscen = "50"),
              by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1_p95[, .(value = quantile(nbenefit,0.05), qscen = "05"),
              by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1_p95[, .(value = quantile(nbenefit,0.95), qscen = "95"),
              by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat1_p95[, qtemp := "95"]
  
  nbftstat1 <- rbindlist(list(nbftstat1, nbftstat1_p05, nbftstat1_p95))
  
  nbft2 <- copy(nben_npv_lvl[dr == disc_rate &
                               region == "World"])
  
  nbft2 <- merge(nbft2, sr15c_runs, by = c("model","scenario"))
  
  nbft2 <- nbft2[,.(model,c5model,sw,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft2 <- nbft2[temp_2100 <= 2]
  nbft2[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft2[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft2[, temp_2100 := NULL]
  nbft2[, dmg_type := "lvl"]
  
  setkey(nbft2,dmg_spec,c5model,sw,model,scenario,temp_clus)
  
  nbftstat2 <- rbindlist(list(
    nbft2[sw == 1, .(value = quantile(nbenefit,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 1, .(value = quantile(nbenefit,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 1, .(value = quantile(nbenefit,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat2[, qtemp := "50"]
  nbftstat2_p05 <- rbindlist(list(
    nbft2[sw == 2, .(value = quantile(nbenefit,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 2, .(value = quantile(nbenefit,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 2, .(value = quantile(nbenefit,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat2_p05[, qtemp := "05"]
  nbftstat2_p95 <- rbindlist(list(
    nbft2[sw == 3, .(value = quantile(nbenefit,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 3, .(value = quantile(nbenefit,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 3, .(value = quantile(nbenefit,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat2_p95[, qtemp := "95"]
  
  nbftstat2 <- rbindlist(list(nbftstat2, nbftstat2_p05, nbftstat2_p95))
  
  nbftstat <- rbind(nbftstat1,nbftstat2)
  
  
  # label dmg spec
  
  dmg_spec_val <- data.table(dmg_spec = c(bhm_spec,djo_spec,lvl_spec),
                             dmg_spec_label = c(bhm_spec_label,djo_spec_label,lvl_spec_label),
                             dmg_spec_star = c(bhm_spec_star,djo_spec_star,lvl_spec_star))
  
  nbftstatx <- merge(nbftstat,dmg_spec_val,by = c("dmg_spec"))
  
  nbftstatx <- dcast(nbftstatx, dmg_spec + temp_clus + dmg_type + qtemp + dmg_spec_label + dmg_spec_star ~ qscen)
  
  dmg_spec_order <- nbftstatx[qtemp == "50",mean(`50`),by = "dmg_spec_label"][order(V1), dmg_spec_label]
  
  nbftstatx[, dmg_spec_label := factor(dmg_spec_label, 
                                       levels = dmg_spec_order)]
  
  nbftstatx[, dmg_type := factor(dmg_type, 
                                 levels = c("gwt","lvl"),
                                 labels = c("growth-based","level-based"))]
  
  wdd <- 0.33
  
  nbftstatx[qtemp == 50, size_point := 4]
  nbftstatx[qtemp != 50, size_point := 3]
  
  nbftstatx[qtemp == 50, size_line := 1]
  nbftstatx[qtemp != 50, size_line := 0.8]
  
  
  p <- ggplot(nbftstatx[qtemp==50]) +
    geom_hline(yintercept = 0, color = "black")  +
    geom_linerange(aes(x = dmg_spec_label,
                       ymin = `05` * 1e2,
                       ymax = `95` * 1e2,
                       alpha = dmg_spec_star,
                       group = qtemp,
                       size = size_line),
                   color = "black",
                   position = position_dodge(width = wdd)) +
    geom_point(aes(x = dmg_spec_label,
                   y = `50` * 1e2,
                   shape = dmg_type,
                   alpha = dmg_spec_star,
                   group = qtemp,
                   size = size_point),
               position = position_dodge(width = wdd)) +
    geom_point(aes(x = dmg_spec_label,
                   y = `50` * 1e2,
                   color = dmg_type,
                   shape = dmg_type,
                   alpha = dmg_spec_star,
                   group = qtemp,
                   size = size_point - 1),
               position = position_dodge(width = wdd)) +
    scale_size_identity() +
    scale_color_npg(guide = ifelse(show_legend,"legend","none"),
                       name = "Damage function") +
    scale_alpha_continuous(range = c(0.25,1), guide = "none") +
    scale_shape_discrete(name = "Damage function", guide = ifelse(show_legend,"legend","none")) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    guides(color = guide_legend(override.aes = list(size = 4) ) ) +
    facet_wrap(~ temp_clus) +
    coord_flip() +
    labs(
#      y = "Net Present Value [Trillion USD2018]", x = ""
       y = "Policy benefits [%]", x = ""
    ) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "right",
      axis.text = element_text(size = 11),
      axis.line.x = element_blank(),
      strip.text = element_text(size = 12, color = "black"),
      strip.background = element_blank())
  
  ggsave(output_file, plot = p, width = 10, height = 7, dpi = 120)
  
  data.table(file = output_file)
  
}

plot_nben_gwt_annual <- function(row_i,output_file) {
  
  nben_ann_gwt <- readd('nben_ann_gwt')
  nben_ann_gwt_p05 <- readd('nben_ann_gwt_p05')
  nben_ann_gwt_p95 <- readd('nben_ann_gwt_p95')
  nben_ann_gwt[, qscen := "50"]
  nben_ann_gwt_p05[, qscen := "05"]
  nben_ann_gwt_p95[, qscen := "95"]
  
  nben_ann_gwt <- rbindlist(list(nben_ann_gwt,nben_ann_gwt_p05,nben_ann_gwt_p95))
  
  sr15c_runs <- readd('sr15c_runs')
  
  scen_ref <- unique(sr15c_runs$scenario_ref)
  
  nbft <- nben_ann_gwt[!scenario %in% scen_ref &
                         !is.na(nbenefit) &
                         region == "World" & 
                         c5model != 1 &
                         dmg_spec %in% c(bhm_spec[bhm_spec_star==1],
                                         djo_spec[djo_spec_star==1]) 
                         ]
  nbft <- nbft[,.(model,scenario,dmg_spec,qscen,runid,c5model,year,nbenefit)]
  nbft <- merge(nbft, sr15c_runs, by = c("model","scenario"))

  nbft <- nbft[temp_2100 <= 2]
  nbft[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft[, temp_2100 := NULL]
  
  setkey(nbft,dmg_spec,c5model,model,scenario,temp_clus,qscen,year)
  
  # stats
  nbftstat <- nbft[, .(cmed = quantile(nbenefit,0.5),
                       cmean = mean(nbenefit),
                       cmin = min(nbenefit),
                       cmax = max(nbenefit),
                       c05 = quantile(nbenefit,0.05),
                       c95 = quantile(nbenefit,0.95),
                       c33 = quantile(nbenefit,0.33),
                       c66 = quantile(nbenefit,0.66),
                       gt0 = length(nbenefit[nbenefit>0])/length(nbenefit)),
                   by = c("qscen","temp_clus","dmg_spec","year")]
  
  dmg_spec_order <- nbftstat[,mean(cmed),by = "dmg_spec"][order(V1), dmg_spec]
  
  nbftstat[, dmg_spec := factor(dmg_spec, 
                                levels = dmg_spec_order)]

  sel_dmg_spec <- dmg_spec_order[((row_i - 1) * 4) + 1:4]
  
  nbftstat[qscen == "05", qscen_label := "5th"]
  nbftstat[qscen == "50", qscen_label := "median"]
  nbftstat[qscen == "95", qscen_label := "95th"]
  nbftstat[, qscen_label := factor(qscen_label, levels = c("5th",
                                                            "median",
                                                            "95th"))]
  if (row_i == 3) {
    show_legend <- "legend"  
  } else {
    show_legend <- "none"
  }
  
  p <- ggplot(data = nbftstat[dmg_spec %in% sel_dmg_spec]) +
    geom_hline(yintercept = 0, size = 0.5) +
    geom_ribbon(aes(x = year, ymin = c05 * 1e2, 
                    ymax = c95 * 1e2,
                    fill = qscen_label),
                alpha = 0.15, size = 0.5) +
    geom_line(aes(x = year, 
                  y = cmed * 1e2,
                  color = qscen_label),
              alpha = 0.6,
              size = 1) +
    #scale_fill_viridis_c(name = "% scenarios > 0", guide = "none") +
    scale_color_brewer(palette = "Set1",
                       guide = "none",
                       name = "Temperature\npercentile") +
    scale_fill_brewer(palette = "Set1",
                      guide = show_legend,
                      name = "Temperature\npercentile") +
    #scale_y_continuous(breaks = seq(-500,500,by = 100), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(-50,200,by = 50),
                       labels = scales::percent_format(scale = 1)) +
    scale_x_continuous(breaks = seq(2020,2100,by = 20), 
                       guide = guide_axis(check.overlap = T),
                       expand = c(0,0)) +
    coord_cartesian(ylim = c(-50,150), xlim = c(2020,2105)) +
    facet_grid(temp_clus ~ dmg_spec) +
#    labs(y = "Net benefits [T$/yr]", x = "") +
    labs(y = "Net benefits [%]", x = "") +
    theme_bw() +
    theme(panel.spacing.x = unit(4, "mm"),
          panel.spacing.y = unit(4, "mm"),
          strip.background = element_blank())
    
    if (row_i == 3) {
      p <- p +guides(fill = guide_legend(override.aes = list(alpha = 1, color = NA) ) )
    }

  ggsave(output_file, plot = p, width = 12, height = 6, dpi = 120)
  
  data.table(file = output_file)
  
}

plot_nben_lvl_annual <- function(row_i,output_file) {
  
  nben_lvl_dmg <- readd('nben_lvl_dmg')
  sr15c_runs <- readd('sr15c_runs')
  
  scen_ref <- unique(sr15c_runs$scenario_ref)
  
  nbft <- nben_lvl_dmg[!scenario %in% scen_ref &
                         !is.na(nbenefit) &
                         region == "World"]
  nbft <- nbft[,.(model,scenario,dmg_spec,c5model,sw,year,nbenefit)]
  nbft <- merge(nbft, sr15c_runs, by = c("model","scenario"))
  
  nbft <- nbft[temp_2100 <= 2]
  nbft[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft[, temp_2100 := NULL]

  nbft[sw == 1, qscen := "50"]
  nbft[sw == 2, qscen := "05"]
  nbft[sw == 3, qscen := "95"]
  
  setkey(nbft,dmg_spec,c5model,model,scenario,temp_clus,sw,year)
  
  
  nbftstat <- nbft[, .(c05 = quantile(nbenefit,0.05),
                       c95 = quantile(nbenefit,0.95),
                       cmed = quantile(nbenefit,0.5),
                       gt0 = length(nbenefit[nbenefit>0])/length(nbenefit)),
                   by = c("temp_clus","dmg_spec","qscen","year")]
  
  dmg_spec_order <- nbftstat[,mean(cmed),by = "dmg_spec"][order(V1), dmg_spec]
  
  nbftstat[, dmg_spec := factor(dmg_spec, 
                                levels = dmg_spec_order)]
  
  nbftstat[qscen == "05", qscen_label := "5th"]
  nbftstat[qscen == "50", qscen_label := "median"]
  nbftstat[qscen == "95", qscen_label := "95th"]
  nbftstat[, qscen_label := factor(qscen_label, levels = c("5th",
                                                           "median",
                                                           "95th"))]
  if (row_i == 2) {
    show_legend <- "legend"  
  } else {
    show_legend <- "none"
  }
  
  sel_dmg_spec <- dmg_spec_order[((row_i - 1) * 4) + 1:4]
  
  p <- ggplot(data = nbftstat[dmg_spec %in% sel_dmg_spec]) +
    geom_hline(yintercept = 0, size = 0.5) +
    geom_ribbon(aes(x = year, ymin = c05 * 1e2, 
                    ymax = c95 * 1e2,
                    fill = qscen_label),
                alpha = 0.15, size = 0.5) +
    geom_line(aes(x = year, 
                  y = cmed * 1e2,
                  color = qscen_label),
              alpha = 0.6,
              size = 1) +
    #scale_fill_viridis_c(name = "% scenarios > 0", guide = "none") +
    scale_color_brewer(palette = "Set1",
                       guide = "none",
                       name = "Temperature\npercentile") +
    scale_fill_brewer(palette = "Set1",
                      guide = show_legend,
                      name = "Temperature\npercentile") +
    scale_x_continuous(breaks = seq(2020,2100,by = 20), 
                       guide = guide_axis(check.overlap = T),
                       expand = c(0,0)) +
    scale_y_continuous(breaks = seq(-50,200,by = 50),
                       labels = scales::percent_format(scale = 1)) +
    coord_cartesian(ylim = c(-50,150), xlim = c(2020,2105)) +
    facet_grid(temp_clus ~ dmg_spec) +
#    labs(y = "Net benefits [T$/yr]", x = "") +
    labs(y = "Net benefits [%]", x = "") +
    theme_bw() +
    theme(panel.spacing.x = unit(4, "mm"),
          panel.spacing.y = unit(4, "mm"),
          strip.background = element_blank())
  
  if (row_i == 2) {
    p <- p +guides(fill = guide_legend(override.aes = list(alpha = 1, color = NA) ) )
  }
  
  
  ggsave(output_file, plot = p, width = 12, height = 6, dpi = 120)
  
  data.table(file = output_file)
  
}

plot_alt_npv <- function(dr0, outfile) {
  
  sr15c_runs <- readd('sr15c_runs')
  
  # Load alternatives
  qscen <- c("05","50","95")
  qscen_fix <- c("_p05","","_p95")
  alts <- c("cc","ct","c0","tc","tt","t0")
  dd <- NULL
  for (p in seq_along(qscen)) {
    for (a in alts) {
     xx <- readd(paste0('nben_npv_gwt',qscen_fix[p],"_alt_",a), character_only = TRUE)
     xx <- xx[region == "World"]
     xx[, qscen := qscen[p]]
     xx[, extrap := a]
     xx[, dmg_type := "gwt"]
     dd <- rbindlist(list(dd,xx))
    }
  }
  for (a in alts) {
    xx <- readd(paste0('nben_npv_lvl_alt_',a), character_only = TRUE)
    xx <- xx[region == "World"]
    xx[sw == 1, qscen := "50"]
    xx[sw == 2, qscen := "05"]
    xx[sw == 3, qscen := "95"]
    xx[, extrap := a]
    xx[, dmg_type := "lvl"]
    dd <- rbindlist(list(dd,xx))
  }

  
  disc_rate = as.numeric(dr0) / 100
  
  dd <- rbindlist(list(dd[dr == disc_rate & region == "World"]))
  
  dd <- merge(dd, sr15c_runs, by = c("model","scenario"))
  dd <- dd[,.(extrap,model,scenario,dmg_spec,dmg_type,qscen,nbenefit,temp_2100)]
  
  # Compute stats
  npv_stat <- dd[, .(cmin = quantile(nbenefit,0.05),
                     cmax = quantile(nbenefit,0.95),
                     c83 = quantile(nbenefit,2.5 / 3),
                     c16 = quantile(nbenefit,0.5 / 3),
                     cmed = quantile(nbenefit,0.5)),
                 by = c("dmg_spec","dmg_type","qscen","extrap")]
  
  npv_stat[, extrap_admg := str_sub(extrap,1,1)]
  npv_stat[, extrap_admg_label := ifelse(extrap_admg == "c", "constant", "trend")]
  npv_stat[, extrap_cmit := str_sub(extrap,2,2)]
  npv_stat[, extrap_cmit_label := ifelse(extrap_cmit == "c", "constant", "trend")]
  npv_stat[, extrap_cmit_label := ifelse(extrap_cmit == "0", "decreasing", extrap_cmit_label)]

  npv_stat[qscen == "05", qscen_label := "temp = 5th perc."]
  npv_stat[qscen == "50", qscen_label := "temp = median"]
  npv_stat[qscen == "95", qscen_label := "temp = 95th perc."]
  npv_stat[, qscen_label := factor(qscen_label, levels = c("temp = 5th perc.",
                                                           "temp = median",
                                                           "temp = 95th perc."))]
  
  
  # label dmg spec
  
  dmg_spec_val <- data.table(dmg_spec = c(bhm_spec,djo_spec,lvl_spec),
                             dmg_spec_label = c(bhm_spec_label,djo_spec_label,
                                               lvl_spec_label),
                             dmg_spec_star = c(bhm_spec_star,djo_spec_star,
                                               lvl_spec_star))
  
  npv_statx <- merge(npv_stat,dmg_spec_val,by = c("dmg_spec"))
  
  dmg_spec_order <- npv_statx[qscen == "50",mean(cmed),by = "dmg_spec_label"][order(V1), dmg_spec_label]
  
  npv_statx[, dmg_spec_label := factor(dmg_spec_label, 
                                       levels = dmg_spec_order)]
  
  
  p <- ggplot(npv_statx) +
    geom_hline(yintercept = 0, color = "black") +
    geom_point(aes(x = dmg_spec_label, 
                   y = cmed * 1e2, 
                   color = extrap_cmit_label,
                   shape = extrap_admg_label,
                   group = extrap_admg_label),
               position = position_dodge(width = 0.5),
               size = 3, alpha = 0.75) +
    geom_linerange(aes(x = dmg_spec_label, 
                   ymin = cmin * 1e2,
                   ymax = cmax * 1e2,
                   color = extrap_cmit_label,
                   group = extrap_admg_label),
               position = position_dodge(width = 0.5),
               size = 0.5, alpha = 0.75) +
    facet_wrap(~ qscen_label) +
    scale_color_d3(name = "Extrapolation\nMitigation cost") +
    scale_shape_discrete(name = "Extrapolation\nAvoided damage") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    coord_flip() +
    labs(
#      y = "Net Present Value [Trillion USD2018]", x = ""
      y = "Policy benefits [%]", x = ""
    ) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "right",
      axis.text = element_text(size = 11),
      axis.line.x = element_blank(),
      strip.text = element_text(size = 12, color = "black"),
      strip.background = element_blank())

  ggsave(outfile, plot = p, width = 12, height = 6, dpi = 120)
  
}

plot_welfare_all <- function(prtp0,output_file,show_legend = TRUE) {
  
  nben_npv_gwt <- readd('nben_welf_gwt')
  nben_npv_gwt_p05 <- readd('nben_welf_gwt_p05')
  nben_npv_gwt_p95 <- readd('nben_welf_gwt_p95')
  nben_npv_lvl <- readd('nben_welf_lvl')
  sr15c_runs <- readd('sr15c_runs')
  
  nbft1 <- rbindlist(list(nben_npv_gwt[prtp == prtp0 & region == "World" & c5model == 1]))
  
  nbft1 <- merge(nbft1, sr15c_runs, by = c("model","scenario"))
  nbft1 <- nbft1[,.(model,scenario,dmg_spec,dcebge,temp_2100)]
  
  nbft1 <- nbft1[temp_2100 <= 2]
  nbft1[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft1[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft1[, temp_2100 := NULL]
  nbft1[, dmg_type := "gwt"]
  nbft1[str_detect(dmg_spec, "LVL"), dmg_type := "lvl"]
  
  setkey(nbft1,dmg_spec,model,scenario,temp_clus)
  
  nbftstat1 <- rbindlist(list(
    nbft1[, .(value = quantile(dcebge,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1[, .(value = quantile(dcebge,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1[, .(value = quantile(dcebge,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat1[, qtemp := "50"]
  
  nbft1_p05 <- rbindlist(list(nben_npv_gwt_p05[prtp == prtp0 & region == "World" & c5model == 1]))
  
  nbft1_p05 <- merge(nbft1_p05, sr15c_runs, by = c("model","scenario"))
  nbft1_p05 <- nbft1_p05[,.(model,scenario,dmg_spec,dcebge,temp_2100)]
  
  nbft1_p05 <- nbft1_p05[temp_2100 <= 2]
  nbft1_p05[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft1_p05[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft1_p05[, temp_2100 := NULL]
  nbft1_p05[, dmg_type := "gwt"]
  nbft1_p05[str_detect(dmg_spec, "LVL"), dmg_type := "lvl"]
  
  setkey(nbft1_p05,dmg_spec,model,scenario,temp_clus)
  
  nbftstat1_p05 <- rbindlist(list(
    nbft1_p05[, .(value = quantile(dcebge,0.5), qscen = "50"),
              by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1_p05[, .(value = quantile(dcebge,0.05), qscen = "05"),
              by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1_p05[, .(value = quantile(dcebge,0.95), qscen = "95"),
              by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat1_p05[, qtemp := "05"]
  
  nbft1_p95 <- rbindlist(list(nben_npv_gwt_p95[prtp == prtp0 & region == "World" & c5model == 1]))
  
  nbft1_p95 <- merge(nbft1_p95, sr15c_runs, by = c("model","scenario"))
  nbft1_p95 <- nbft1_p95[,.(model,scenario,dmg_spec,dcebge,temp_2100)]
  
  nbft1_p95 <- nbft1_p95[temp_2100 <= 2]
  nbft1_p95[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft1_p95[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft1_p95[, temp_2100 := NULL]
  nbft1_p95[, dmg_type := "gwt"]
  nbft1_p95[str_detect(dmg_spec, "LVL"), dmg_type := "lvl"]
  
  setkey(nbft1_p95,dmg_spec,model,scenario,temp_clus)
  
  nbftstat1_p95 <- nbft1_p95[, .(cmax = quantile(dcebge,0.5),
                                 cmax_05 = quantile(dcebge,0.05),
                                 cmax_95 = quantile(dcebge,0.95)),
                             by = c("temp_clus","dmg_spec","dmg_type")]
  nbftstat1_p95 <- rbindlist(list(
    nbft1_p95[, .(value = quantile(dcebge,0.5), qscen = "50"),
              by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1_p95[, .(value = quantile(dcebge,0.05), qscen = "05"),
              by = c("temp_clus","dmg_spec","dmg_type")],
    nbft1_p95[, .(value = quantile(dcebge,0.95), qscen = "95"),
              by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat1_p95[, qtemp := "95"]
  
  nbftstat1 <- rbindlist(list(nbftstat1, nbftstat1_p05, nbftstat1_p95))
  
  nbft2 <- copy(nben_npv_lvl[prtp == prtp0 &
                               region == "World"])
  
  nbft2 <- merge(nbft2, sr15c_runs, by = c("model","scenario"))
  
  nbft2 <- nbft2[,.(model,c5model,sw,scenario,dmg_spec,dcebge,temp_2100)]
  
  nbft2 <- nbft2[temp_2100 <= 2]
  nbft2[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft2[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft2[, temp_2100 := NULL]
  nbft2[, dmg_type := "lvl"]
  
  setkey(nbft2,dmg_spec,c5model,sw,model,scenario,temp_clus)
  
  nbftstat2 <- rbindlist(list(
    nbft2[sw == 1, .(value = quantile(dcebge,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 1, .(value = quantile(dcebge,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 1, .(value = quantile(dcebge,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat2[, qtemp := "50"]
  nbftstat2_p05 <- rbindlist(list(
    nbft2[sw == 2, .(value = quantile(dcebge,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 2, .(value = quantile(dcebge,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 2, .(value = quantile(dcebge,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat2_p05[, qtemp := "05"]
  nbftstat2_p95 <- rbindlist(list(
    nbft2[sw == 3, .(value = quantile(dcebge,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 3, .(value = quantile(dcebge,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type")],
    nbft2[sw == 3, .(value = quantile(dcebge,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type")]))
  nbftstat2_p95[, qtemp := "95"]
  
  nbftstat2 <- rbindlist(list(nbftstat2, nbftstat2_p05, nbftstat2_p95))
  
  nbftstat <- rbind(nbftstat1,nbftstat2)
  
  
  # label dmg spec
  
  dmg_spec_val <- data.table(dmg_spec = c(bhm_spec,djo_spec,lvl_spec),
                             dmg_spec_label = c(bhm_spec_label,djo_spec_label,lvl_spec_label),
                             dmg_spec_star = c(bhm_spec_star,djo_spec_star,lvl_spec_star))
  
  nbftstatx <- merge(nbftstat,dmg_spec_val,by = c("dmg_spec"))
  
  nbftstatx <- dcast(nbftstatx, dmg_spec + temp_clus + dmg_type + qtemp + dmg_spec_label + dmg_spec_star ~ qscen)
  
  dmg_spec_order <- nbftstatx[qtemp == "50",mean(`50`),by = "dmg_spec_label"][order(V1), dmg_spec_label]
  
  nbftstatx[, dmg_spec_label := factor(dmg_spec_label, 
                                       levels = dmg_spec_order)]
  
  nbftstatx[, dmg_type := factor(dmg_type, 
                                 levels = c("gwt","lvl"),
                                 labels = c("growth-based","level-based"))]
  
  wdd <- 0.33
  
  nbftstatx[qtemp == "50", size_point := 4]
  nbftstatx[qtemp != "50", size_point := 3]
  
  nbftstatx[qtemp == "50", size_line := 1]
  nbftstatx[qtemp != "50", size_line := 0.8]
  
  
  nbftstatx[qtemp == "05", qscen_label := "5th"]
  nbftstatx[qtemp == "50", qscen_label := "median"]
  nbftstatx[qtemp == "95", qscen_label := "95th"]
  nbftstatx[, qscen_label := factor(qscen_label, levels = c("5th",
                                                            "median",
                                                            "95th"))]
  
  p <- ggplot(nbftstatx) +
    geom_hline(yintercept = 0, color = "black")  +
    geom_linerange(aes(x = dmg_spec_label,
                       ymin = `05` * 1e2,
                       ymax = `95` * 1e2,
                       color = qscen_label,
                       alpha = dmg_spec_star,
                       group = qtemp,
                       size = size_line),
                   position = position_dodge(width = wdd)) +
    geom_point(aes(x = dmg_spec_label,
                   y = `50` * 1e2,
                   shape = dmg_type,
                   alpha = dmg_spec_star,
                   group = qtemp,
                   size = size_point),
               position = position_dodge(width = wdd)) +
    geom_point(aes(x = dmg_spec_label,
                   y = `50` * 1e2,
                   color = qscen_label,
                   shape = dmg_type,
                   alpha = dmg_spec_star,
                   group = qtemp,
                   size = size_point - 1),
               position = position_dodge(width = wdd)) +
    scale_size_identity() +
    scale_color_brewer(palette = "Set1",
                       guide = ifelse(show_legend,"legend","none"),
                       name = "Temperature\npercentile") +
    scale_alpha_continuous(range = c(0.25,1), guide = "none") +
    scale_shape_discrete(name = "Damage function", guide = ifelse(show_legend,"legend","none")) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    facet_wrap(~ temp_clus) +
    coord_flip() +
    labs(
      y = "Delta CEGBE [%]", x = ""
    ) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "right",
      axis.text = element_text(size = 11),
      axis.line.x = element_blank(),
      strip.text = element_text(size = 12, color = "black"),
      strip.background = element_blank())
  
  ggsave(output_file, plot = p, width = 12, height = 7, dpi = 120)
  
  data.table(file = output_file)
  
}

plot_nben_rcp_v2_model <- function(dr0,output_file,show_legend = TRUE) {
  
  nben_npv_gwt <- readd('nben_npv_gwt')
  nben_npv_gwt_p05 <- readd('nben_npv_gwt_p05')
  nben_npv_gwt_p95 <- readd('nben_npv_gwt_p95')
  nben_npv_lvl <- readd('nben_npv_lvl')
  sr15c_runs <- readd('sr15c_runs')
  
  disc_rate = as.numeric(dr0) / 100
  
  nbft1 <- rbindlist(list(nben_npv_gwt[dr == disc_rate & region == "World" & c5model == 1]))
  
  nbft1 <- merge(nbft1, sr15c_runs, by = c("model","scenario"))
  nbft1 <- nbft1[,.(model,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft1 <- nbft1[temp_2100 <= 2]
  nbft1[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft1[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft1[, temp_2100 := NULL]
  nbft1[, dmg_type := "gwt"]
  nbft1[str_detect(dmg_spec, "LVL"), dmg_type := "lvl"]
  
  nbft1[, modtype := ifelse(str_detect(model,"MESSAGE"), "BU", "TD")]
  
  setkey(nbft1,dmg_spec,model,scenario,temp_clus,modtype)
  
  nbftstat1 <- rbindlist(list(
    nbft1[, .(value = quantile(nbenefit,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type","modtype")],
    nbft1[, .(value = quantile(nbenefit,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type","modtype")],
    nbft1[, .(value = quantile(nbenefit,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type","modtype")]))
  nbftstat1[, qtemp := "50"]
  
  nbft1_p05 <- rbindlist(list(nben_npv_gwt_p05[dr == disc_rate & region == "World" & c5model == 1]))
  
  nbft1_p05 <- merge(nbft1_p05, sr15c_runs, by = c("model","scenario"))
  nbft1_p05 <- nbft1_p05[,.(model,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft1_p05 <- nbft1_p05[temp_2100 <= 2]
  nbft1_p05[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft1_p05[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft1_p05[, temp_2100 := NULL]
  nbft1_p05[, dmg_type := "gwt"]
  nbft1_p05[str_detect(dmg_spec, "LVL"), dmg_type := "lvl"]
  
  
  nbft1_p05[, modtype := ifelse(str_detect(model,"MESSAGE"), "BU", "TD")]
  
  setkey(nbft1_p05,dmg_spec,model,scenario,temp_clus,modtype)
  
  nbftstat1_p05 <- rbindlist(list(
    nbft1_p05[, .(value = quantile(nbenefit,0.5), qscen = "50"),
              by = c("temp_clus","dmg_spec","dmg_type","modtype")],
    nbft1_p05[, .(value = quantile(nbenefit,0.05), qscen = "05"),
              by = c("temp_clus","dmg_spec","dmg_type","modtype")],
    nbft1_p05[, .(value = quantile(nbenefit,0.95), qscen = "95"),
              by = c("temp_clus","dmg_spec","dmg_type","modtype")]))
  nbftstat1_p05[, qtemp := "05"]
  
  nbft1_p95 <- rbindlist(list(nben_npv_gwt_p95[dr == disc_rate & region == "World" & c5model == 1]))
  
  nbft1_p95 <- merge(nbft1_p95, sr15c_runs, by = c("model","scenario"))
  nbft1_p95 <- nbft1_p95[,.(model,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft1_p95 <- nbft1_p95[temp_2100 <= 2]
  nbft1_p95[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft1_p95[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft1_p95[, temp_2100 := NULL]
  nbft1_p95[, dmg_type := "gwt"]
  nbft1_p95[str_detect(dmg_spec, "LVL"), dmg_type := "lvl"]
  
  nbft1_p95[, modtype := ifelse(str_detect(model,"MESSAGE"), "BU", "TD")]
  
  setkey(nbft1_p95,dmg_spec,model,scenario,temp_clus,modtype)
  
  nbftstat1_p95 <- nbft1_p95[, .(cmax = quantile(nbenefit,0.5),
                                 cmax_05 = quantile(nbenefit,0.05),
                                 cmax_95 = quantile(nbenefit,0.95)),
                             by = c("temp_clus","dmg_spec","dmg_type")]
  nbftstat1_p95 <- rbindlist(list(
    nbft1_p95[, .(value = quantile(nbenefit,0.5), qscen = "50"),
              by = c("temp_clus","dmg_spec","dmg_type","modtype")],
    nbft1_p95[, .(value = quantile(nbenefit,0.05), qscen = "05"),
              by = c("temp_clus","dmg_spec","dmg_type","modtype")],
    nbft1_p95[, .(value = quantile(nbenefit,0.95), qscen = "95"),
              by = c("temp_clus","dmg_spec","dmg_type","modtype")]))
  nbftstat1_p95[, qtemp := "95"]
  
  nbftstat1 <- rbindlist(list(nbftstat1, nbftstat1_p05, nbftstat1_p95))
  
  nbft2 <- copy(nben_npv_lvl[dr == disc_rate &
                               region == "World"])
  
  nbft2 <- merge(nbft2, sr15c_runs, by = c("model","scenario"))
  
  nbft2 <- nbft2[,.(model,c5model,sw,scenario,dmg_spec,nbenefit,temp_2100)]
  
  nbft2 <- nbft2[temp_2100 <= 2]
  nbft2[, temp_clus := ifelse(temp_2100 <= 1.5, 1, 2)]
  nbft2[, temp_clus := factor(temp_clus, levels = 1:2, labels = c("1.5C", "2C"))]
  nbft2[, temp_2100 := NULL]
  nbft2[, dmg_type := "lvl"]
  
  nbft2[, modtype := ifelse(str_detect(model,"MESSAGE"), "BU", "TD")]
  
  setkey(nbft2,dmg_spec,model,scenario,temp_clus,modtype)
  
  nbftstat2 <- rbindlist(list(
    nbft2[sw == 1, .(value = quantile(nbenefit,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type","modtype")],
    nbft2[sw == 1, .(value = quantile(nbenefit,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type","modtype")],
    nbft2[sw == 1, .(value = quantile(nbenefit,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type","modtype")]))
  nbftstat2[, qtemp := "50"]
  nbftstat2_p05 <- rbindlist(list(
    nbft2[sw == 2, .(value = quantile(nbenefit,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type","modtype")],
    nbft2[sw == 2, .(value = quantile(nbenefit,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type","modtype")],
    nbft2[sw == 2, .(value = quantile(nbenefit,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type","modtype")]))
  nbftstat2_p05[, qtemp := "05"]
  nbftstat2_p95 <- rbindlist(list(
    nbft2[sw == 3, .(value = quantile(nbenefit,0.5), qscen = "50"),
          by = c("temp_clus","dmg_spec","dmg_type","modtype")],
    nbft2[sw == 3, .(value = quantile(nbenefit,0.05), qscen = "05"),
          by = c("temp_clus","dmg_spec","dmg_type","modtype")],
    nbft2[sw == 3, .(value = quantile(nbenefit,0.95), qscen = "95"),
          by = c("temp_clus","dmg_spec","dmg_type","modtype")]))
  nbftstat2_p95[, qtemp := "95"]
  
  nbftstat2 <- rbindlist(list(nbftstat2, nbftstat2_p05, nbftstat2_p95))
  
  nbftstat <- rbind(nbftstat1,nbftstat2)
  
  
  # label dmg spec
  
  dmg_spec_val <- data.table(dmg_spec = c(bhm_spec,djo_spec,lvl_spec),
                             dmg_spec_label = c(bhm_spec_label,djo_spec_label,lvl_spec_label),
                             dmg_spec_star = c(bhm_spec_star,djo_spec_star,lvl_spec_star))
  
  nbftstatx <- merge(nbftstat,dmg_spec_val,by = c("dmg_spec"))
  
  nbftstatx <- dcast(nbftstatx, modtype + dmg_spec + temp_clus + dmg_type + qtemp + dmg_spec_label + dmg_spec_star ~ qscen)
  
  dmg_spec_order <- nbftstatx[qtemp == "50",mean(`50`),by = "dmg_spec_label"][order(V1), dmg_spec_label]
  
  nbftstatx[, dmg_spec_label := factor(dmg_spec_label, 
                                       levels = dmg_spec_order)]
  
  nbftstatx[, dmg_type := factor(dmg_type, 
                                 levels = c("gwt","lvl"),
                                 labels = c("growth-based","level-based"))]
  
  wdd <- 0.33
  
  nbftstatx[qtemp == "50", size_point := 4]
  nbftstatx[qtemp != "50", size_point := 3]
  
  nbftstatx[qtemp == "50", size_line := 1]
  nbftstatx[qtemp != "50", size_line := 0.8]
  
  
  nbftstatx[qtemp == "05", qscen_label := "5th"]
  nbftstatx[qtemp == "50", qscen_label := "median"]
  nbftstatx[qtemp == "95", qscen_label := "95th"]
  nbftstatx[, qscen_label := factor(qscen_label, levels = c("5th",
                                                            "median",
                                                            "95th"))]
  
  p <- ggplot(nbftstatx[dmg_spec_star ==1 & qtemp == "50"]) +
    geom_hline(yintercept = 0, color = "black")  +
    geom_linerange(aes(x = dmg_spec_label,
                       ymin = `05` * 1e2,
                       ymax = `95` * 1e2,
                       color = modtype,
                       alpha = dmg_spec_star,
                       group = modtype,
                       size = size_line),
                   position = position_dodge(width = wdd)) +
    geom_point(aes(x = dmg_spec_label,
                   y = `50` * 1e2,
                   shape = modtype,
                   alpha = dmg_spec_star,
                   group = modtype,
                   size = size_point),
               position = position_dodge(width = wdd)) +
    geom_point(aes(x = dmg_spec_label,
                   y = `50` * 1e2,
                   color = modtype,
                   shape = modtype,
                   alpha = dmg_spec_star,
                   group = modtype,
                   size = size_point - 1),
               position = position_dodge(width = wdd)) +
    scale_size_identity() +
    scale_color_brewer(palette = "Set1",
                       guide = ifelse(show_legend,"legend","none"),
                       name = "Temperature\npercentile") +
    scale_alpha_continuous(range = c(0.25,1), guide = "none") +
    scale_shape_discrete(name = "Damage function", guide = ifelse(show_legend,"legend","none")) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    facet_grid( ~ temp_clus) +
    coord_flip() +
    labs(
      #y = "Net Present Value [Trillion USD2018]", x = ""
      y = "Policy benefits [%]", x = ""
    ) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "right",
      axis.text = element_text(size = 11),
      axis.line.x = element_blank(),
      strip.text = element_text(size = 12, color = "black"),
      strip.background = element_blank())
  
  ggsave(output_file, plot = p, width = 10, height = 7, dpi = 120)
  
  data.table(file = output_file)
  
}
