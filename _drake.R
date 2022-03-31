# Basic drake file
# Needs to be adjusted according to the plan and the resources

source("R/packages.R")  # Load packages.
source("R/functions.R") # Load functions.
source("R/plan.R")      # Drake plan.

future::plan(future::multicore)

config <- drake_config(plan, 
                       #plan_boot,
                       #plan_lvl, 
                       #plan_gwt,
                       #plan_boot_p05,
                       #plan_boot_p95,
                       parallelism = "future",
                       format = 'fst_dt',
                       jobs = 4)
