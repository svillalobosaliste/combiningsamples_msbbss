# master_thesis
Master's thesis

Sofía Villalobos Aliste (6060714)

Research repository of thesis manuscript produced between September 2021 and May 2022. 

This repository contains all the information regarding the master thesis and internship project at the Centraal Bureau voor de Statstiek "Combining probability and non-probability samples for estimation", of the program Methodology and Statistics for the Behavioural, Biomedical, and Social Sciences, Utrecht University year 2022.

It is included here:

A -The approval files from the Ethical Review Board of the Faculty of Social and Behavioural Sciences of Utrecht University for the simulation study and the use of existing data.

B - The thesis manuscript.

C - The results are produced in R (version 4.1.2), and saved in excel files.

D - The code with which the results were produced.

The tables included in the thesis manuscript are constructed with the information of the outputs shown in this repository and can be replicated as follows:

1_create_function: This script creates the necessary functions to generate the simulated samples and compute three different estimators for each probability, non-probability, and combined sample.

2_storage: This script runs the simulation files and storage the mean of the ARMSE per domain, and the mean absolute bias per domain. After creating the functions, running this script takes around 7 days.

3_results_marmse: This script produces outputs for tables 3 to 7, and Appendices B, C, and D.

4_results_bias: This script produces outputs for tables 8 to 10 and Appendix E.

5_data_cbs: This script reproduces the outputs for tables 11 and 12 using data from CBS. The data from CBS is in a secure environment of CBS and should be required by email directly to them.

6_rnorm: This script includes an alternative approach mentioned at the end of section 3.1 of the thesis manuscript. This script is self-sufficient and includes the functions to generate the samples and compute the estimators, the storage of the results of the MARMSE, and a summary table.

