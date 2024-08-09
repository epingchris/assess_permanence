This script is the implementation of the method described in Rau et al. (2024).
It uses carbon loss values in the project area and in the counterfactual scenario of a REDD+ project to simulate yearly values of:
1. Additionality (Mg CO2e)
2. Anticipated release (Mg CO2e)
3. Equivalent permanence (EP, between 0 and 1)
4. Reversal risks of (between 0 and 1)

The carbon loss values can be either taken from observed time series, randomly drawn from a distribution fitted to observed values, or randomly drawn from exponential distributions for a theoretical project.

## Overview

1. SimulatePermanence() performs the simulation; its output is provided to the "outlist" argument of the SaveStandard() function
2. SaveStandard() performs basic summary and saves it to an RDS object; its output is provided to the "summary_out" argument of the PlotStandard() function
3. PlotStandard() generates and saves basic time series plots

## Explanation of the SimulatePermanence() function

The core function that performs the simulation is SimulatePermanence(), which  takes the following arguments:
1. type: character, "theo" for theoretical projects, "real" for real-life projects

2. mean_drawdown: numerical vector, mean drawdown rate(s) for theoretical project(s)
Its inverse will be used as the lambda parameter of an exponential distribution
When aggregate_type is NULL and the length of this argument > 1, it initialises a custom theoretical aggregated project

3. sites: character vector, name(s) for the real-life project(s)
The function will look for a data frame of carbon flux time series at the path "/project_input_data/sites.csv"
When aggregate_type is NULL and the length of this argument > 1, it initialises a custom real-life aggregated project

4. aggregate_type: character, specifying default settings for theoretical/real-life aggregated projects
This argument overrides both mean_drawdown and sites when not NULL:
"A": theoretical, mean_drawdown = c(1.1, 1.1, 1.1, 5)
"B": theoretical, mean_drawdown = c(1.1, 1.1, 5, 5)
"C": theoretical, mean_drawdown = c(1.1, 5, 5, 5)
"three": real-life, sites = c("Gola_country", "CIF_Alto_Mayo", "VCS_1396")
"four": real-life, sites = c("Gola_country", "CIF_Alto_Mayo", "VCS_1396", "VCS_934")

5. verbose: boolean, whether to print basic output at each timestep (default: FALSE)
6. runtime: boolean, whether to print (default: TRUE)
7. omega: numeric, threshold of acceptable reversal risk (default: 0.05)
8. n_rep: numeric, number of repetitions (default: 100)
9. H: numeric, project duration (years) (default: 50)
10. D: numeric, discount rate (default: 0.03)
11. warmup: numeric, warm-up period (years) (default: 5)
12. postproject_ratio: numeric, ratio of post-project release compared to during-project additionality accumulation rate (default: 2)
13. scc_df: numeric data frame, containing two columns, year and value (default: scc, a data frame which should be loaded and prepared)

## Required input data

The dataset that one needs to prepare to run a simulation for a specific real-life project is the yearly project/counterfactual carbon fluxes in the project.
It should be saved as a csv file in the "project_input_data" directory; the name of the file will be the "sites" variable.
The file must contain the following columns:
1. year: numeric
2. var: character, either "project", "counterfactual" or "additionality"
3. val: numeric, total carbon flux from the last year to this year (Mg CO2e)
4. n_sim: numeric, index of repetition
5. started: boolean, whether the year is larger than project start (t0)

## Other scripts

extend_scc.R and create_input_csv.R are maintenance scripts used to prepare input data.
ViewSnapshot.r is a function that selects one particular repetition and view its time series, but is not currently used for anything.