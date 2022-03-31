# IAMPACT-SR1.5C - Net benefits of the IPCC SR1.5C Scenarios

This repository is released under the Apache License 2.0; see the
[LICENSE](LICENSE) for details.

Copyright 2022 Laurent Drouet.

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

## Usage

The code is written in R and relies on two packages for reproducibility:
`drake` for the workflow and `renv` for the dependency management. To
run the analysis, you need to install first the `renv` packages, as it
will be suggested if you open R from the main directory. Then execute
once the following command to install the missing dependencies:

``` r
renv::restore()
```

To launch the data analysis, run

    drake::r_make()

The default plan (`plan`) will be executed. As the data analysis is
heavily computing demanding, some additional plans need to run
sequentially.

For this, modify the `_drake.R` file and comment/uncomment the first
argument of the `drake_config` function so that only one plan name is
specified. The plans to run are:

-   `plan` (default)
-   `plan_boot`
-   `plan_lvl`
-   `plan_gwt`
-   `plan_boot_p05`
-   `plan_boot_p95`

Relaunch the command `r_make()` for each plan.

Moreover, the drake configuration should be adjusted to you computing
system. For example, using LSF job queues with 100 nodes of 36 cores:

``` r
options(clustermq.scheduler = "lsf", 
        clustermq.template = "lsf_clustermq.tmpl")

make(plan_boot,
     parallelism = "clustermq",
     memory_strategy = "autoclean",
     garbage_collection = TRUE,
     caching = "worker",
     jobs = 100,
     format = 'fst_dt',
     template = list(job_name = "iboot", 
                     log_file = "iamp_boot.log"
                     ))
```

Overall, it took 24+ hours to run the whole data on our HPC.

To build the figures, run

``` r
source("R/packages.R")
source("R/packages_plots.R")
source("R/functions_.R")
source("R/plan.R")

config <- drake_config(plan2)
    
make(plan2)
```

For any questions, please ask laurent.drouet [at] eiee.org.

## Funding acknowledgement

<img src="./data/Flag_of_Europe.svg.png" alt="EU logo" width="80" height="54" align="left"/>
This project has received funding from the European Union’s Horizon 2020
research and innovation programme under grant agreements Nos 821471
(ENGAGE) and 821124 (NAVIGATE).
