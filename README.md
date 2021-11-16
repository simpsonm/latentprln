# Interpolating Population Distributions using Public-use Data: An Application to Income Segregation using American Community Survey Data
This Github repository is a companion to the paper given by the title above. It contains code and instructions for reproducing all of the key results in the paper, including instructions for using the code to fit the models discussed in the paper to other data. 

# Reproduction Instructions
We include instructions to reproduce all analyses performed in the paper. This can be broken down into the following major steps:

* Download and clean the data
	* Create the synthetic population
	    * Create the synthetic datasets
		    * Run the simulation study
			    * Compile the results of the simulation study
	* Fit models to the ACS PUMA-level data
	   	* Exploratory analysis
	    * Compile the results of fitting models to ACS PUMA-level data
	* Fit models to the ACS Metro-level data and construct the income
      segregation indicies
		* Fit EIV regressions using income segregation indicies

Steps that are indented below a previous step depend on that step. So
every step below **Download and clean the data** depends on the that step being performed first, but once the data has been cleaned it does not matter whether **Exploratory Analysis** is performed first.

The following subsections will explain in detail how to perform each of these steps. Each subsection also lists all dependencies needed to complete it, including software dependencies and dependencies on previous steps, and contains an estimate for how long it will take to run that step. Most of these steps can be run on a moderately powerful laptop or desktop, though one step requires a powerful compute node.

## Download and clean the data
Most of the data used in the paper is downloaded from the Census Bureau using the tidycensus R package, and cleaned using the tidyverse R package. This is done in three files:
* `/code/acsclean.R`: Download and clean the data for the PUMA-level example in Section 3.2
* `/code/metroclean.R`: Download and clean the data for the metro-level example in Section 4.
* `/code/metroraceclean.R`: Download and clean the data for the metro-level by race example in Section 4.
Each script will save several files in the `/code/data` director. The metro scripts take some time to run, be patient. To use these script, you need three things:
* A Census API key, which may be obtained [here](https://api.census.gov/data/key_signup.html)
* A file that maps Census tracts to PUMAs, available [here](https://www2.census.gov/geo/docs/maps-data/data/rel/2010_Census_Tract_to_2010_PUMA.txt) from the Census Bureau, but also included in this repository: `/code/data/tract2puma2010.txt`.
* A file that maps Census tracts to metro areas, available [here](https://www.census.gov/geographies/reference-files/time-series/demo/metro-micro/delineation-files.html) from the Census Bureah, but also included in this repository: `/code/data/metro_deliniation_4-2018.csv`

Dependencies:
* The file `/code/clean_functions.R`.
* `R` version 4.1.0 and the following packages (and all dependencies of these packages):
   * `tidyverse` version 1.3.1
   * `tidycensus` version 1.1
Newer versions of these packages should work.

Time to complete: 5 minutes for `acsclean.R`. Several hours for each of the other scripts.

## Create the synthetic population
To complete this step, the `R` script `/code/simulation/sim_setup.R` needs to be run. This script creates the synthetic population used in the simulation study. The synthetic population is saved as `/code/simulation/population.RData`. Our R script will not reproduce our population in `R` version 4.1, porbably due to the change to the random number generator in version 3.6. Normally, the line `RNGkind(sample.kind = "Rounding")` would solve the problem, but it does not for our script. This is probably due to how packages we rely on use the random number generator. For this reason as well as convenience we include the `population.RData` file in this Github repository. Additionally, this script creates the plots in Figure A.2. 

Dependencies:
* previous steps
    * **Download and clean the data**
* Download the 2014 Missouri household PUMS from [here](https://www2.census.gov/programs-surveys/acs/data/pums/2014/5-Year/csv_hmo.zip), and place in `/code/data/MO/`. **This file is not included in the repository because it is too large.** Additionally the tidycensus R function `get_pums` can download PUMS files, though we had trouble using it to download this specific file.
* Download the 2014 Missouri shapefiles from [here](https://www.census.gov/cgi-bin/geo/shapefiles/index.php), and place in `/code/data/MO/shapefile/`. We have included these files in the repository for convenience, and because the names of these files may be different than what is expected by our scripts.
* `R` version 4.1.0 and the following packages (and all dependencies of these packages):
   * `tidyverse` version 1.3.1
   * `tidycensus` version 1.1
   * `sf` version 0.9.8
   * `sp` version 1.4-5
   * `maptools` version 1.1-1
   * `reldist` version 1.6-6
   * `gridExtra` version 2.3
   * `ggthemes` version 4.2.4
   * `broom` version 0.7.7
   
The above versions are sufficient to run the script, but to replicate `population.RData`, R version < 3.6 and corresponding packages are needed. See the apragraph at the beginning of this section. Unfortunately, we do not know the precise versions of the other packages necessary to completely reproduce the population, but we include the population in the repository for convenience.

Time to complete: 5 minutes.

## Create the synthetic datasets
To complete this step, the `R` script `/code/simulation/sim_create_sample.R` needs to be run. This script repeatedly does the following 1000 times:
1. Sample from the synthetic population using a stratified random sample.
2. Create synthetic tract-level estimates and their associated standard errors.
3. Subsample the initial sample to create a synthetic PUMS.
The results of each iteration of this process are compiled and stored in the list `samples.RData` which is saved in this directory. This file is too large to store on Github, so it is not included in this repository. Additionally, the datasets are processed into the correct form for the Stan model files, and stored in the list `standats.RData`. This file is also too large to store on Github, so it is not included in this repository. 

Dependencies:
* previous steps
    * **Download the data**
    * **Clean the data**
	* **Create the synthetic population**
* `R` version 4.1.0 and the following packages packages (and all dependencies of these packages):
   * `rstan` version 2.21.2
   * `sp` version 1.4.5
   * `maptools` version 1.1-1
   * `reldist` version 1.6-6
   * `survey` version 4.0
   * `Hmisc` version 4.5-0
   * `dplyr` version 1.0.6

Time to complete: about 30 minutes.

## Run the simulation study
**Warning:** This step is very computationally intensive and is intended to use multiple cores.

The code used to run this section is split into three files that are intended to be run simultaneously on a 12 core machine, each using 4 cores: `fit1.R`, `fit2.R`, and `fit3.R`, all in `/code/simulation`. The outputted files are `stat.out1.RData`, `stat.out2.RData`, and `stat.out3.RData`, all of which are too large to be included in the Github repository.

Dependencies:
* previous steps
    * **Download and clean the data**
	* **Create the synthetic population**
	* **Create the synthetic datasets**
* `R` version 4.10 and the following packages packages (and all dependencies of these packages):
   * `rstan` version 2.21.2
   * `reldist` version 1.6-6
   * `parallel` base package
   * `MASS` version 7.3-54
   * `data.table` version 1.14.0
   * `abind` version 1.4-5
   * `tidyverse` version 1.3.1

Time to complete: about 2 weeks on a moderately powerful desktop computer with 6 physical cores (12 logical cores) - running `fit1.R`, `fit2.R`, and `fit3.R` simultaneously.

## Compile the results of the simulation study
This step compiles consists of several parts: 1) fit the orignal PRLN to the simulation study, then 2) compile the results of the sumulation study. The `R` script `/code/simulation/fitprln.R` fits the original PRLN to each sample from the population and saves results into `prln.out.RData`. This file is not included in this repository. The `R` script `/code/simulation/postprocess.R` then compiles the results and creates tables 1, 2, and 3. Tables are created by `xtable` and are manually copy/pasted into the appropriate `.tex` file with minor editing.

Dependencies for `fitprln.R`:
* previous steps
    * **Download and clean the data**
	* **Create the synthetic population**
	* **Create the synthetic datasets**
	* **Run the simulation study**
* `R` version 4.1 (no packages needed)

Time to complete: about 5 minutes.

Dependencies for `postprocess.R`:
* previous steps
    * **Download and clean the data**
	* **Create the synthetic population**
	* **Create the synthetic datasets**
	* **Run the simulation study**
	* Run `fitprln.R`, earlier this step.	
* `R` version 4.1 and the following packages packages (and all dependencies of these packages):
   * `sp` version 1.4-5
   * `maptools` version 1.1-1
   * `Hmisc` version 4.5-0
   * `survey` version 4.0
   * `reldist` version 1.6-6
   * `dplyr` version 1.0.6
   * `MASS` version 7.3-54
   * `data.table` version 1.14.0
   * `rstan` version 2.21.2
   * `abind` version 1.4-5
   * `tidyverse` version 1.3.1
   * `parallel` base package
   * `xtable` version 1.8-4
   * `reshape2` verseion 1.4.4
   
Time to complete: about 5 minutes.

## Exploratory analysis
To complete this step, only one `R` script needs to be run: `/code/acspuma/explore.R`. This creates several of the tables and figures in Appendix A. Tables are created by `xtable` and are manually copy/pasted into the appropriate `.tex` file with minor editing. The plots are saved to `/doc`. These plots are included in the Github repository for convenience.

Dependencies:
* previous steps
    * **Download and clean the data**
* `R` version 4.1.0 and the following packages packages (and all dependencies of these packages):
   * `tidyverse` version 1.3.1
   * `sf` version 0.9.8
   * `ggthemes` version 4.2.4
   * `xtable` version 1.8-4

Time to complete: 1 minute.


## Fit models to the ACS PUMA-level data

The following scripts fit the tract-level models to each PUMA:
1. `/code/acspuma/fitco.R`
2. `/code/acspuma/fitil.R`
3. `/code/acspuma/fitmo.R`
4. `/code/acspuma/fitmt.R`
5. `/code/acspuma/fitny.R`
Each model is fit using 4 cores, one per Markov chain. The MCMC draws
are saved separately for each each PUMA in the folder
`/code/acspuma`. These files are too large to include in the Github
repository. Each model takes a about 5-10 minutes to fit.

Dependencies:
* previous steps
    * **Download and clean the data**
* `R` version 4.1 and the following packages packages (and all dependencies of these packages):
   * `rstan` version 2.21.2
   * `tidyverse` version 1.3.1

Time to complete: < 1 hour. On a machine with many cores, the scripts can be run simultaneously to save time.

## Compile the results of fitting models to ACS PUMA-level data

This section compiles the results of fitting all of the models to the ACS data, fits the original Pareto Linear procedure to the ACS data, then puts it all into the appropriate tables. Tables D.1-D.5 are created by `xtable` and are manually copy/pasted into the appropriate `.tex` file with minor editing. The script `/code/acs_puma/postprocess.R` accomplishes all of this.

Dependencies:
* previous steps
    * **Download and clean the data**
	* **Fit models to the ACS PUMA-level data**
* `R` version 4.1 and the following packages packages (and all dependencies of these packages):
   * `rstan` version 2.21.2
   * `tidyverse` version 1.3.1
   * `sf` version 0.9.8
   * `xtable` version 1.8-4	 
   * `sp` version 1.4-5
   * `maptools` version 1.1-1
   * `reldist` version 1.6-6
   * `survey` version 4.0
   * `Hmisc` version 4.5-0
   * `parallel` base package
   * `MASS` version 7.3-54
   * `data.table` version 1.14.0
   * `abind` version 1.4-5
   * `reshape2` verseion 1.4.4

Time to complete: 2-3 hours. On a non-windows machine much can be parallelized using the `ncore` parameter to save time.

## Fit models to the ACS Metro-level data and construct the income segregation indicies
This section fits all of the tract-level models to each tract for each metro area, then constructs the information theory and divergence indicies for the analysis in Section 4. This is performace for 1) all households, 2) black households only, and 3) white households only, for each tract in each metro area. This process is very memory and computation intensive. We divided this into 12 distinct jobs that could be run in parallel, each job with its own compute node with 28 cores and 8GB of memory per core. The shell scripts to schedule the jobs are included in the repository, though the details of how to set these up will depend on the particulars of the server you have access to:
* `code/metrorace/metro_fit_2018_black_1.sh`
* `code/metrorace/metro_fit_2018_black_2.sh`
* `code/metrorace/metro_fit_2018_black_3.sh`
* `code/metrorace/metro_fit_2018_black_4.sh`
* `code/metrorace/metro_fit_2018_combined_1.sh`
* `code/metrorace/metro_fit_2018_combined_2.sh`
* `code/metrorace/metro_fit_2018_combined_3.sh`
* `code/metrorace/metro_fit_2018_combined_4.sh`
* `code/metrorace/metro_fit_2018_white_1.sh`
* `code/metrorace/metro_fit_2018_white_2.sh`
* `code/metrorace/metro_fit_2018_white_3.sh`
* `code/metrorace/metro_fit_2018_white_4.sh`
Each job corresponds to a particular `R` script:
* `code/metrorace/metro_fit_2018_black_1.R`
* `code/metrorace/metro_fit_2018_black_2.R`
* `code/metrorace/metro_fit_2018_black_3.R`
* `code/metrorace/metro_fit_2018_black_4.R`
* `code/metrorace/metro_fit_2018_combined_1.R`
* `code/metrorace/metro_fit_2018_combined_2.R`
* `code/metrorace/metro_fit_2018_combined_3.R`
* `code/metrorace/metro_fit_2018_combined_4.R`
* `code/metrorace/metro_fit_2018_white_1.R`
* `code/metrorace/metro_fit_2018_white_2.R`
* `code/metrorace/metro_fit_2018_white_3.R`
* `code/metrorace/metro_fit_2018_white_4.R`
Each script runs a single chain for each tract in each metro area for each group of households. Then all the chains for the tracts in a particular metro area for a particular group of households are jointly used to contruct the information theory and divergence indicies along with their standard errors for each iteration of the chains. These are saved in the following files
* `code/metrorace/metro_kls_2018_black_1.RData`
* `code/metrorace/metro_kls_2018_black_2.RData`
* `code/metrorace/metro_kls_2018_black_3.RData`
* `code/metrorace/metro_kls_2018_black_4.RData`
* `code/metrorace/metro_kls_2018_combined_1.RData`
* `code/metrorace/metro_kls_2018_combined_2.RData`
* `code/metrorace/metro_kls_2018_combined_3.RData`
* `code/metrorace/metro_kls_2018_combined_4.RData`
* `code/metrorace/metro_kls_2018_white_1.RData`
* `code/metrorace/metro_kls_2018_white_2.RData`
* `code/metrorace/metro_kls_2018_white_3.RData`
* `code/metrorace/metro_kls_2018_white_4.RData`
These files are too large to include in the repository.

Dependencies:
* Note the computational requirements in the description above
* previous steps
    * **Download and clean the data**
* `R` version 4.1 and the following packages packages (and all dependencies of these packages):
   * `rstan` version 2.21.2
   * `tidyverse` version 1.3.1
   * `foreach` version 1.5.1
   * `doParallel` version 1.016

Time to complete: Each job takes up to 7 days to complete on a powerful compute node with 28 cores and >= 8GB memory per core, though the jobs can be run in parallel.

## EIV regressions using income segregation indicies
This section, 1) creates figures H.4, H.5, and H.6, then 2) estimates gini coefficients using latent PRLN for black households only and white households only for each tract in each metro area, then 3) fits the EIV regressions from section 4 for all households, black households only, and white households only, and finally 4) uses the results to construct tables H.1, H.2, H.3, H.4, H.5, and H.6, and tables 4 and 5. Fitting the EIV regressions is memory intensive, requiring >=32 GB of memory.

Dependencies:
* 32 GB of memory
* previous steps
    * **Download and clean the data**
	* **Fit models to the ACS Metro-level data and construct the income segregation indicies**
* `R` version 4.1 and the following packages packages (and all dependencies of these packages):
   * `rstan` version 2.21.2
   * `tidyverse` version 1.3.1
   * `sf` version 0.9.8
   * `xtable` version 1.8-4	 
   * `reldist` version 1.6-6

Time to complete: About 24 hours and a moderately powerful desktop computer (with 32 GB of memory).
