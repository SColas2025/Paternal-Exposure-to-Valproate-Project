# Paternal-Exposure-to-Valproate-Project
GitHub repository containing SAS codes for generating numerical results and statistical estimates presented in the publication: **"Paternal Valproate Use and Neurodevelopmental Disorder and Congenital Malformation Risk in Offspring"**. JAMA Netw Open. Published November 10, 2025. doi:10.1001/jamanetworkopen.2025.42581. 

In this study, we evaluated the association between paternal exposure to valproate and neurodevelopmental disorder or congenital malformation risk in children.

### Input data

The analyses performed by the SAS code in this repository are based on a dataset described in the file "01_data_description.docx". This input dataset used in the SAS-code is created separately for each country (Sweden, Denmark, Norway) by combining data from multiple National registers within that country. 

### SAS-version 
SAS 9.4

### Document overview

#### `01_data_description.docx`
This file contains the description of the dataset used for the analysis,
including a detailed variable list. 

#### `Assoc_NDD_maternal`, `Assoc_NDD_offspring`, `Assoc_NDD_paternal`
These SAS codes produce numbers and estimate for the analyses of the Association between maternal/offspring/paternal characteristics and NDD (no related Table in the manuscript / supplementary material)

#### `CoxCrudeEffectEst_NDD`
This SAS code produces numbers and estimates for the Effect Estimate Crude analyses (Table 2) 

#### `CoxPS_EffectEst_KMeansAdj_NDD`
This SAS code produces numbers and estimates for the Effect estimates in the different Clusters of paternal exposition (no related Table in the manuscript / supplementary)

#### `CoxPS_EffectEst_NDD`
This SAS code produces numbers and estimates for the Effect Estimate Adjusted analyses (Table 2)

#### `PS_weight_logistic`, `PSLogisticModelEstimates`, `PSBalanceAssessment_NDD`
These SAS codes produce numbers and estimate for the Variable Estimates from Logistic Regression Propensity Score Model (eTable 11) and estimates for the PS Model Balance analyses (eFigures 3, 4 and 5 and eTable 12)

### Availability of data and materials

Data are not accessible outside of the specific authorization granted to the investigator.

### Contact

Sandrine Colas `Sandrine.Colas@sanofi.com`

project link: https://github.com/SColas2025/Paternal-Exposure-to-Valproate-Project













