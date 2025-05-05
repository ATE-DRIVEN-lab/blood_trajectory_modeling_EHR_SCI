# Modeling of trajectories of routine blood values as dynamic biomarkers in spinal cord injury
Repo that reproduce the modeling of routine blood trajectories in SCI from EHR data study performed by the Health.data DRIVEN lab at the School of Public Health Sciences, University of Waterloo. Code created by Drs. Marzieh Mussavi Rizi and Abel Torres Espin.

**Pre-Print**
> Modeling trajectories of routine blood tests as dynamic biomarkers for outcome in spinal cord injury
Marzieh Mussavi Rizi, Daniel Fernandez, John LK Kramer, Rajiv Saigal, Anthony M. DiGiorgio, Michael S. Beattie, Adam R Ferguson, Nikos Kyritsis, Abel Torres-Espin, TRACK-SCI investigators
medRxiv 2025.01.20.25320728; doi: https://doi.org/10.1101/2025.01.20.2532072

**Peer-reviewed publication**
>TBD

# Notes on reproducing this work

## Datasets
We do not provide the datasets directly, and users of this code will need to download the data. Three datasets are used in this work: MIMIC-III version 1.4, MIMIC-IV version 1.0, and a subset of the [TRACK-SCI cohort study](https://spinalcordinjury.ucsf.edu/content/track-sci-0). The .Rmd script contains the necessary code to prepare the data for analysis.

### MIMIC
Data has been download from [PhysioNet](https://physionet.org/). Both MIMIC databases (DB) are relational DB structured in tables. Documentation about the DB schema can be found [here](https://mimic.mit.edu/docs/).

Note that data access need Data Use Agreement with PhysioNet. No data is provided in this document or repository. The code would not run without the data!

### TRACK-SCI dataset
The necessary TRACK-SCI data can be downloaded from the Open Data Commons for Spinal Cord Injury (SCI) [here](https://doi.org/10.34945/F5PK6X). If you use the data, please cite:
>Mussavi Rizi, M., Saigal, R., DiGiorgio, A. M., Ferguson, A. R., Beattie, M. S., Kyritsis, N., Torres Espin, A.. 2025. Blood laboratory values from 137 de-identified TRACK-SCI participants from routine collected real-world data. Open Data Commons for Spinal Cord Injury. ODC-SCI:1345. doi: 10.34945/F5PK6X

### SAPSII
Part of this work uses SAPS II values for both MIMIC datasets. If you want to reproduce our work using this code, you will need to calculate it first, and save it in a `mimic_SC_saps.csv` file that contains four columns: subject_id = subject identifier;	hadm_id = hospital admision identifier;	icustay_id = ICU stay identifier;	sapsii = calculated SAPS II.

For MIMIC-III, we compute SAPS II scores for the selected cohort using SQL code publicly available on GitHub. (https://github.com/MIT-LCP/mimic-code/blob/main/mimic-iii/concepts/severityscores/sapsii.sql). For MIMIC-IV, we used the equivalent script (https://github.com/MIT-LCP/mimic-code/blob/main/mimic-iv/concepts/score/sapsii.sql).
 
## Dependencies
The code should run with the following environment. Further information can be found in the .Rmd file.

"R version 4.4.1 (2024-06-14 ucrt)", "tidyverse 2.0.0", "data.table 1.17.0", "stringr 1.5.1", "DT 0.33", "gtsummary 2.2.0", "lcmm 1.9.4", "caret 7.0-1", "yardstick 1.3.2", "patchwork 1.3.0", "parallel (base with R 4.4.1)"

## Running the code
Some sections of the .Rmd script are not evaluated during knitting (rendering) due to their computational overhead. We have provided intermediate files containing the necessary objects and biproducts of the code, including the final trajectory models to facilitate reproducibility. By cloning this repo, you should be able to reproduce our results without having to re-fit all the models, but you can do so too. `To reproduce this work, you will need to run the code on the IDE by chunk`. 

>Running the full script from scratch will override some of the provided files and it can take hours to days to complete.
