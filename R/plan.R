# 1 + 20 cmodels
cmodels <- c(
  "CMIP5",
  "ACCESS1-0", "BNU-ESM", "CCSM4", "CMCC-CMS", "GFDL-CM3",
  "GFDL-ESM2G", "GFDL-ESM2M", "GISS-E2-H", "GISS-E2-H-CC", "GISS-E2-R",
  "GISS-E2-R-CC", "HadGEM2-CC", "HadGEM2-ES", "IPSL-CM5A-LR", "IPSL-CM5A-MR",
  "IPSL-CM5B-LR", "MPI-ESM-LR", "MPI-ESM-MR", "NorESM1-M", "inmcm4"
)

bhm_spec <- c(
  "BHM SR", "BHM SR RP",
  "BHM LR", "BHM LR RP",
  "BHM LR ORIG", "BHM LR RP ORIG",
  "KW2020 BHM", "BT2019 STATE", "BT2019 CNTRY",
  "PRETIS2018",
  "HS2019","ACEVODO2020"
)
bhm_spec_label <- c(
  "BHM SR", "BHM SR RP",
  "BHM LR L", "BHM LR P",
  "BHM LR", "BHM LR RP",
  "KW BHM", "BT STATE", "BT CNTRY",
  "PRETIS",
  "HENSELER","ACEVODO"
)
bhm_spec_star <- c(
  1, 1,
  1, 1,
  0, 0,
  0, 0, 1,
  1,
  1, 1
)

set.seed(42)

bootrunid <- 1:1000

bhm_spec_boot <- c(
  paste("BHM SR BOOT", bootrunid)
)

djo_spec <- c("DJO 0L", "DJO 5L", "DJO 10L",
              "DJO 0L ORIG", "DJO 5L ORIG", "DJO 10L ORIG")
djo_spec_label <- c("DJO 0L P", "DJO 5L P", "DJO 10L P",
                   "DJO 0L RP", "DJO 5L RP", "DJO 10L RP")
djo_spec_star <- c(1, 1, 1,
                   0, 0, 0)

lvl_spec <- c("HS2017 NCAT","HS2017 TOT","DICE2016","KW2020","TAKAKURA2019","HS2017 TOTP","HS2017 MKT")
lvl_spec_label <- c("HS NCAT","HS TOT","DICE","KW LVL","TAKAKURA","HS TOT+P","HS MKT")
lvl_spec_star <- c(1,1,1,1,1,1,1)

dmg_spec_boot <- bhm_spec_boot
dmg_spec_gwth <- c(bhm_spec,djo_spec)

dmg_spec_all <- c(bhm_spec, djo_spec, dmg_spec_boot)
dmg_spec_plot <- c(bhm_spec, djo_spec)

discount_rates <- c(0.01, 0.02, 0.03)
drplot <- c(1, 2, 3)

keep_years <- c(seq(2020,2100,by = 5))
keep_years_idx <- which(2010:2100 %in% keep_years)

disc_years_keep <- c(seq(2020,2100,by = 5),seq(2120,2300,by = 20))

r5_regions <- c("R5OECD90+EU","R5ASIA","R5LAM","R5MAF","R5REF")

prtps <- c(0.00, 0.01, 0.02, 0.03)
prtpplot <- c(0, 1, 2)

dmg_diag_scenario <- c("SSP3-Baseline","SSP2-26")
dmg_diag_years <- c(2100)

zeus <- !Sys.info()[["nodename"]] %in% c("spleen")

extrap_nben <- c("cc","c0","ct","tc","t0","tt") # "cc" is default

plan <- drake_plan(
  
  # WDI Statistics
  wdi_stat = target(read_fst(file_in("data/wdi_stat.fst"),
    as.data.table = TRUE
  )),

  # ORIGINAL SSP ASSUMPTIONS
  # Source: https://tntcat.iiasa.ac.at/SspDb/dsd
  ssp_csv = target("data/SspDb_country_data_2013-06-12.csv",
    format = "file"
  ),
  
  ## OECD Population [millions]
  ssp_pop = ssp_query(ssp_csv, "OECD Env-Growth", "Population"),
  ## OECD GDP PPP [billion US$2005/yr]
  ssp_gdp = ssp_query(ssp_csv, "OECD Env-Growth", "GDP|PPP"),
  ## Compute gdp per capita for SSP scenario [USD2018/yr]
  ssp_gdpcap = compute_gdp_cap(ssp_gdp, ssp_pop, wdi_stat),

  # IAMC SCENARIOS: SR15C DB
  # Source: https://data.ene.iiasa.ac.at/iamc-1.5c-explorer/#/downloads
  sr15c_csv = target("data/iamc15_scenario_data_all_regions_r2_0.csv",
    format = "file"
  ),
  sr15c_mo_no = target(c(
    "GENeSYS-MOD 1.0", # <2100
    "IEA Energy Technology Perspective Model 2017", # <2100
    "Shell World Energy Model 2018", # <2100
    "IEA World Energy Model 2017" # no baseline
  ), format = "rds"),
  sr15c_scen = target("^ADVANCE_|^CD-LINKS_|^SSP\\d-|^SFCM_|^EMF33_|^TERL_|^GEA_|^CEMICS",
    format = "rds"
  ),
  sr15c_outliers = get_sr15c_outliers(),
  ## GDP|MER["billion US$2010/yr"]
  gdp_mit_mer = sr15c_query(
    sr15c_csv,
    "GDP|MER",
    sr15c_mo_no, sr15c_scen, sr15c_outliers
  ),
  ## GDP|PPP["billion US$2010/yr"]
  gdp_mit_ppp = sr15c_query(
    sr15c_csv,
    "GDP|PPP",
    sr15c_mo_no, sr15c_scen, sr15c_outliers
  ),
  ## Area under MAC Curve ["billion US$2010/yr"]
  gdp_mit_mac = sr15c_query(
    sr15c_csv,
    "Policy Cost|Area under MAC Curve",
    sr15c_mo_no, sr15c_scen, sr15c_outliers
  ),
  ## Area under MAC Curve ["billion US$2010/yr"]
  gdp_mit_nrg = sr15c_query(
    sr15c_csv,
    "Policy Cost|Additional Total Energy System Cost",
    sr15c_mo_no, sr15c_scen, sr15c_outliers
  ),
  ## Global Mean Temperature [C] median
  temp_magicc = sr15c_query(
    sr15c_csv,
    "AR5 climate diagnostics|Temperature|Global Mean|MAGICC6|MED",
    sr15c_mo_no, sr15c_scen, sr15c_outliers
  ),
  ## Global Mean Temperature [C] percentile 5%
  temp_magicc_p05 = sr15c_query(
    sr15c_csv,
    "AR5 climate diagnostics|Temperature|Global Mean|MAGICC6|P5",
    sr15c_mo_no, sr15c_scen, sr15c_outliers
  ),
  ## Global Mean Temperature [C] percentile 95%
  temp_magicc_p95 = sr15c_query(
    sr15c_csv,
    "AR5 climate diagnostics|Temperature|Global Mean|MAGICC6|P95",
    sr15c_mo_no, sr15c_scen, sr15c_outliers
  ),
  ## Total Radiative Forcing [W.m-2]
  trf_magicc = sr15c_query(
    sr15c_csv,
    "AR5 climate diagnostics|Forcing|MAGICC6|MED",
    sr15c_mo_no, sr15c_scen, sr15c_outliers
  ),
  # Scenario clusters and reference scenario
  sr15c_info = get_scen_info(trf_magicc, temp_magicc, gmt_hist),
  # Mitigation costs
  ## consolidate GDP REF and GDP MIT ["billion US$2010/yr"]
  gdp_mit = consolidate_gdp_mit(
    gdp_mit_mer,
    gdp_mit_ppp,
    gdp_mit_mac,
    gdp_mit_nrg,
    sr15c_info,
    sr15c_scen
  ),
  # Scenario selection
  sr15c_runs = select_scen_runs(gdp_mit, sr15c_info),
  
  # CLIMATE
  ## climate downscaling coefficients
  clim_dwnscl = fread(file_in("data/fit_estimates.csv")),
  ## Country-level temperature and precipitation
  climate_med = target(downscale_climate(clim_dwnscl, 
                                         temp_magicc, 
                                         cmodels,
                                         sr15c_runs),
     dynamic = map(cmodels)
  ),
  climate_p05 = target(downscale_climate(clim_dwnscl, 
                                         temp_magicc_p05, 
                                         cmodels,
                                         sr15c_runs),
                       dynamic = map(cmodels)
  ),
  climate_p95 = target(downscale_climate(clim_dwnscl, 
                                         temp_magicc_p95, 
                                         cmodels,
                                         sr15c_runs),
                       dynamic = map(cmodels)
  ),
  ## Climate historical timeseries extended until 2019
  clim_hist = read_hist_climate(file_in("data/hist_temp.csv"),
                                climate_med),
  ## CRU GMT until 2019 wrt 1850-1880
  gmt_hist = read_hist_gmt(),
  # ECONOMIC IMPACTS OF CC
  dmg_param = target(read_dmg_param(ssp_gdp, wdi_stat),
     format = "rds"
  ),

)

# BOOTSTRAP damage function 
plan_boot <- drake_plan(
  
  ## GROWTH - BOOTSTRAP SPECIFICATION
  ## Project GDP per capita with climate and aggregate GDP
  gdp_boot = target(project_cc_boot_impact(
    dmg_spec
  ),
    transform = cross(dmg_spec = !!dmg_spec_boot)
  ),
  nben_boot_dmg = target(compute_net_benef_years(gdp_boot),
    transform = map(gdp_boot)
  ),
  nben_npv_boot_dmg = target(compute_net_benef_npv(gdp_boot,discount_rates),
    transform = map(gdp_boot)
  ),
  nben_npv_boot = target(rbind(nben_npv_boot_dmg),
    transform = combine(nben_npv_boot_dmg)
  ),
  
  ## DATA for PLOTS
  nben_ann_boot_smp = target(filter_boot_sample(nben_boot_dmg),
                               transform = map(nben_boot_dmg)),
  data_boot_smp = target(rbind(nben_ann_boot_smp),
                           transform = combine(nben_ann_boot_smp)),
  nben_boot_dmg_plot = target(filter_boot_mid_end_century(nben_boot_dmg),
                                           transform = map(nben_boot_dmg)),
  nben_boot_plot = target(rbind(nben_boot_dmg_plot),
                             transform = combine(nben_boot_dmg_plot)),
)

## LVL
plan_lvl <- drake_plan(
  # LEVELS
  gdp_levels = target(project_cc_impact_levels()),
  nben_lvl_dmg = target(compute_net_benef_years(gdp_levels)),
  nben_npv_lvl = target(compute_net_benef_npv(
    gdp_levels,
    discount_rates)
  ),
  nben_npv_lvl_alt = target(compute_net_benef_npv(
    gdp_levels,
    discount_rates,
    extrap),
    transform = cross(extrap = !!extrap_nben)),
  nben_welf_lvl = target(compute_diff_cebge(gdp_mit,gdp_levels,disc_years_keep,prtps))
)

## GWT
plan_gwt <- drake_plan(
  # GROWTH - ESTIMATES
  # Project GDP per capita with climate and aggregate GDP
  gdp_gwt_all = target(project_cc_gwt_impact(dmg_spec),
     transform = cross(dmg_spec = !!dmg_spec_gwth)
  ),
  gdp_gwt = target(rbind(gdp_gwt_all),
     transform = combine(gdp_gwt_all)
  ),
  nben_ann_gwt = target(compute_net_benef_years(gdp_gwt)),
  nben_npv_gwt = target(compute_net_benef_npv(gdp_gwt,discount_rates)),
  nben_npv_gwt_alt = target(compute_net_benef_npv(gdp_gwt,discount_rates,extrap),
                            transform = cross(extrap = !!extrap_nben)),
  nben_welf_gwt = target(compute_diff_cebge(gdp_mit,gdp_gwt,disc_years_keep,prtps)),
  

  #p05
  gdp_gwt_all_p05 = target(project_cc_gwt_impact(dmg_spec,
                                                 climate_data = "climate_p05"),
                           transform = cross(dmg_spec = !!dmg_spec_gwth)
  ),
  gdp_gwt_p05 = target(rbind(gdp_gwt_all_p05),
                       transform = combine(gdp_gwt_all_p05)
  ),
  nben_ann_gwt_p05 = target(compute_net_benef_years(gdp_gwt_p05)),
  nben_npv_gwt_p05 = target(compute_net_benef_npv(gdp_gwt_p05,discount_rates)),
  nben_npv_gwt_p05_alt = target(compute_net_benef_npv(gdp_gwt_p05,discount_rates,extrap),
                            transform = cross(extrap = !!extrap_nben)),
  nben_welf_gwt_p05 = target(compute_diff_cebge(gdp_mit,gdp_gwt_p05,disc_years_keep,prtps)),
  
  #p95
  gdp_gwt_all_p95 = target(project_cc_gwt_impact(dmg_spec,
                                                 climate_data = "climate_p95"),
                           transform = cross(dmg_spec = !!dmg_spec_gwth)
  ),
  gdp_gwt_p95 = target(rbind(gdp_gwt_all_p95),
                       transform = combine(gdp_gwt_all_p95)
  ),
  nben_ann_gwt_p95 = target(compute_net_benef_years(gdp_gwt_p95)),
  nben_npv_gwt_p95 = target(compute_net_benef_npv(gdp_gwt_p95,discount_rates)),
  nben_npv_gwt_p95_alt = target(compute_net_benef_npv(gdp_gwt_p95,discount_rates,extrap),
                                transform = cross(extrap = !!extrap_nben)),
  nben_welf_gwt_p95 = target(compute_diff_cebge(gdp_mit,gdp_gwt_p05,disc_years_keep,prtps))
)

## BOOT - P05 
plan_boot_p05 <- drake_plan(
  ## GROWTH - BOOTSTRAP SPECIFICATION
  ## Project GDP per capita with climate and aggregate GDP
  gdp_boot_p05 = target(project_cc_boot_impact(dmg_spec, 
                                               climate_data = "climate_p05"),
        transform = cross(dmg_spec = !!dmg_spec_boot)),
  nben_boot_dmg_p05 = target(compute_net_benef_years(gdp_boot_p05),
        transform = map(gdp_boot_p05)),
  nben_npv_boot_dmg_p05 = target(compute_net_benef_npv(gdp_boot_p05,discount_rates),
        transform = map(gdp_boot_p05)),
  nben_npv_boot_p05 = target(rbind(nben_npv_boot_dmg_p05),
        transform = combine(nben_npv_boot_dmg_p05)),
  nben_ann_boot_smp_p05 = target(filter_boot_sample(nben_boot_dmg_p05),
        transform = map(nben_boot_dmg_p05)),
  data_boot_smp_p05 = target(rbind(nben_ann_boot_smp_p05),
        transform = combine(nben_ann_boot_smp_p05)),
  nben_boot_dmg_plot_p05 = target(filter_boot_mid_end_century(nben_boot_dmg_p05),
                              transform = map(nben_boot_dmg_p05)),
  nben_boot_plot_p05 = target(rbind(nben_boot_dmg_plot_p05),
                          transform = combine(nben_boot_dmg_plot_p05)),
)

## BOOT - P95
plan_boot_p95 <- drake_plan(
  ## GROWTH - BOOTSTRAP SPECIFICATION
  ## Project GDP per capita with climate and aggregate GDP
  gdp_boot_p95 = target(project_cc_boot_impact(dmg_spec, 
                                               climate_data = "climate_p95"),
                transform = cross(dmg_spec = !!dmg_spec_boot)),
  nben_boot_dmg_p95 = target(compute_net_benef_years(gdp_boot_p95),
                transform = map(gdp_boot_p95)),
  nben_npv_boot_dmg_p95 = target(compute_net_benef_npv(gdp_boot_p95,discount_rates),
                transform = map(gdp_boot_p95)),
  nben_npv_boot_p95 = target(rbind(nben_npv_boot_dmg_p95),
                transform = combine(nben_npv_boot_dmg_p95)),
  nben_ann_boot_smp_p95 = target(filter_boot_sample(nben_boot_dmg_p95),
                transform = map(nben_boot_dmg_p95)),
  data_boot_smp_p95 = target(rbind(nben_ann_boot_smp_p95),
                transform = combine(nben_ann_boot_smp_p95)),
  nben_boot_dmg_plot_p95 = target(filter_boot_mid_end_century(nben_boot_dmg_p95),
                                  transform = map(nben_boot_dmg_p95)),
  nben_boot_plot_p95 = target(rbind(nben_boot_dmg_plot_p95),
                              transform = combine(nben_boot_dmg_plot_p95)),
)

plan2 <- drake_plan(
  
  # Manuscript
  
  ## Fig1
  fig_method_1 = plot_method_temp(file_out(!!file.path('artifacts','plots',"fig_method_temp.pdf"))),
  fig_method_2 = plot_method_temp_map(file_out(!!file.path('artifacts','plots',"fig_method_temp_map.pdf"))),
  fig_method_3 = plot_method_lvl_dmg(file_out(!!file.path('artifacts','plots',"fig_method_lvl_dmg.pdf"))),
  fig_method_4 = plot_method_gwt_dmg(file_out(!!file.path('artifacts','plots',"fig_method_gwt_dmg.pdf"))),
  fig_method_5 = plot_method_admg(file_out(!!file.path('artifacts','plots',"fig_method_admg.pdf"))),
  fig_method_6 = plot_method_cmit(file_out(!!file.path('artifacts','plots',"fig_method_cmit.pdf"))),
  
  ## Fig2a + Supplementary Figure 5
  fig_dist_bhm_cluster = target(plot_nben_boot_rcp(dr, "BHM SR BOOT", 
    output_file = file_out(!!file.path("artifacts", "plots", paste0(.id_chr, ".pdf")))),
    transform = cross(dr = !!drplot)
  ),
  
  fig_dist_bhm_all = target(plot_nben_boot_rcp_nocluster(dr, "BHM SR BOOT", 
    output_file = file_out(!!file.path("artifacts", "plots", paste0(.id_chr, ".png")))),
    transform = cross(dr = !!drplot)
  ),
  
  fig_dist_bhm_all_main = target(plot_nben_boot_rcp_nocluster_main(2, "BHM SR BOOT", 
    output_file = file_out(!!file.path("artifacts", "plots", "fig_dist_bhm_all_main_2.png")))),
  
  ## Fig2b
  plot3 = target(plot_nben_boot_ann(
    output_file = file_out(!!file.path('artifacts','plots',"plot_nben_boot_ann_bhm_sr.png")), 
    data_file = "data_boot_smp")
  ),
  plot3_0 = target(plot_nben_boot_ann_nocluster(
    output_file = file_out(!!file.path('artifacts','plots',"plot_nben_boot_ann_bhm_sr_0.png")), 
    data_file = "data_boot_smp")
  ),
  
  
  ## Fig3
  map_nb1 = map_nben(year0 = 2050, temp = "1.5C", quantity = "nbenefit_share", 
                     show_legend = TRUE, 
                     output_file = file_out(!!file.path('artifacts','plots',"map_nben_shr_2050_15C.pdf"))),
  map_nb2 = map_nben(year0 = 2100, temp = "1.5C", quantity = "nbenefit_share", 
                     show_legend = TRUE, 
                     output_file = file_out(!!file.path('artifacts','plots',"map_nben_shr_2100_15C.pdf"))),
  map_nb3 = map_nben(year0 = 2050, temp = "2C", quantity = "nbenefit_share", 
                     show_legend = FALSE, 
                     output_file = file_out(!!file.path('artifacts','plots',"map_nben_shr_2050_2C.pdf"))),
  map_nb4 = map_nben(year0 = 2100, temp = "2C", quantity = "nbenefit_share", 
                     show_legend = FALSE, 
                     output_file = file_out(!!file.path('artifacts','plots',"map_nben_shr_2100_2C.pdf"))),
  map_nb5 = map_nben(year0 = 2100, temp = "1.5C", quantity = "nbenefit_value", 
                     show_legend = TRUE, 
                     output_file = file_out(!!file.path('artifacts','plots',"map_nben_val_2100_15C.pdf"))),
  map_nb6 = map_nben(year0 = 2100, temp = "2C", quantity = "nbenefit_value", 
                     show_legend = FALSE, 
                     output_file = file_out(!!file.path('artifacts','plots',"map_nben_val_2100_2C.pdf"))),
  
  ## Fig4
  all_plot3_v3 = plot_nben_rcp_v2_med(2,
    output_file = file_out(!!file.path('artifacts','plots',"plot_all_dr2.pdf")),
    show_legend = T
  ),
  
  # Supplementary Material
  
  ## Supplementary Figure 2

  fig_si_cmit = plot_cmit(file_out(!!file.path('artifacts','plots',"fig_si_cmit.pdf"))),
  
  ## Supplementary Figure 3
  
  fig_si_temp = plot_temp(file_out(!!file.path('artifacts','plots',"fig_si_temp.pdf"))),

  ## Supplementary Figure 4
  
  fig_si_temp_q = plot_temp_quantile(file_out(!!file.path('artifacts','plots',"fig_si_temp_quantile.pdf"))),
  
  # Supplementary Figure 5
  
  plot2 = target(plot_nben_boot_temp(dr, "BHM SR BOOT",
    output_file = file_out(!!file.path("artifacts", "plots", paste0(.id_chr, ".pdf")))),
    transform = cross(dr = !!drplot)
  ),
  
  ## Supplementary Figure 6

  plot3_05 = target(plot_nben_boot_ann(
    output_file = file_out(!!file.path('artifacts','plots',"plot_nben_boot_ann_bhm_sr_p05.png")), 
    data_file = "data_boot_smp_p05")
  ),
  plot3_95 = target(plot_nben_boot_ann(
    output_file = file_out(!!file.path('artifacts','plots',"plot_nben_boot_ann_bhm_sr_p95.png")), 
    data_file = "data_boot_smp_p95")
  ),
  
  ## Supplementary Figure 8
  
  all_plot2 = plot_nben_rcp_v2(1, 
   output_file = file_out(!!file.path('artifacts','plots',"plot_all_dr1.pdf")),
   show_legend = T
  ),
  all_plot4 = plot_nben_rcp_v2(3, 
   output_file = file_out(!!file.path('artifacts','plots',"plot_all_dr3.pdf")),
  show_legend = T
  ),
  
  ## Supplementary Figure 9

  dmg1 = plot_dmg_2100(output_file = file_out(!!file.path('artifacts','plots',"plot_dmg_2100.pdf"))),
  
  ## Supplementary Figure 10

  plot8a = target(plot_nben_gwt_annual(1,
                                       output_file = file_out(!!file.path('artifacts','plots',"plot_nben_gwt_ann1.pdf")))),
  plot8b = target(plot_nben_gwt_annual(2,
                                       output_file = file_out(!!file.path('artifacts','plots',"plot_nben_gwt_ann2.pdf")))),
  plot8c = target(plot_nben_gwt_annual(3,
                                       output_file = file_out(!!file.path('artifacts','plots',"plot_nben_gwt_ann3.pdf")))),
  plot9a = target(plot_nben_lvl_annual(1,
                                       output_file = file_out(!!file.path('artifacts','plots',"plot_nben_lvl_ann1.pdf")))),
  plot9b = target(plot_nben_lvl_annual(2,
                                       output_file = file_out(!!file.path('artifacts','plots',"plot_nben_lvl_ann2.pdf")))),
  
  ## Supplementary Methods
  
  plotsix1 = target(fig_extrap_admg(output_file = file_out(!!file.path('artifacts','plots',"plot_extrap_admg.pdf")))),
  plotsix2 = target(fig_extrap_cmit(output_file = file_out(!!file.path('artifacts','plots',"plot_extrap_cmit.pdf")))),
  
  plotext1 = target(plot_alt_npv(1,outfile = file_out(!!file.path('artifacts','plots',"plot_extrap_1.pdf")))),
  plotext2 = target(plot_alt_npv(2,outfile = file_out(!!file.path('artifacts','plots',"plot_extrap_2.pdf")))),
  plotext3 = target(plot_alt_npv(3,outfile = file_out(!!file.path('artifacts','plots',"plot_extrap_3.pdf")))),

  ploths1 = target(plot_nben_rcp_v2_hs2017(1,file_out(!!file.path('artifacts','plots',"plot_hs2017_1.pdf")))),
  ploths2 = target(plot_nben_rcp_v2_hs2017(2,file_out(!!file.path('artifacts','plots',"plot_hs2017_2.pdf")))),
  ploths3 = target(plot_nben_rcp_v2_hs2017(3,file_out(!!file.path('artifacts','plots',"plot_hs2017_3.pdf")))),
  
  plotmod = target(plot_nben_rcp_v2_model(2,file_out(!!file.path('artifacts','plots',"plot_modtype.pdf")))),
  
  ## Welfare analysis
  
  plotwelf0 = plot_welfare_all(0,file_out(!!file.path('artifacts','plots',"plot_dcegbe_0.pdf"))),
  plotwelf1 = plot_welfare_all(0.01,file_out(!!file.path('artifacts','plots',"plot_dcegbe_1.pdf"))),
  plotwelf2 = plot_welfare_all(0.02,file_out(!!file.path('artifacts','plots',"plot_dcegbe_2.pdf"))),
  plotwelf3 = plot_welfare_all(0.03,file_out(!!file.path('artifacts','plots',"plot_dcegbe_3.pdf"))),
  
  ## Alternative axis for a better comparison of temperature cluster
  all_plot3_vx = plot_nben_rcp_v2_x(2, 
                                  output_file = file_out(!!file.path('artifacts','plots',"plot_all_dr2x_all.pdf")),
                                  show_legend = T
  )
  
)
