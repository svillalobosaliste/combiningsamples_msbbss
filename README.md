# master_thesis
Masther's thesis repository contains all the scripts for simulations and analysis of the data
belonging to the master thesis and internship project at the Centraal Bureau voor de 
Statstiek "Combining probability and non-probability samples for estimation", 
of the program Methodology and Statistics for the Behavioural, Biomedical, and Social Sciences, 
Utrecht University year 2022.


The outputs from the thesis_manuscript can be replicated as follows: 

1_create_function: This scripts creates the necessary functions to generate the
simulated samples and compute the estimators for each 
probability, non-probability, and combined samples.

2_storage: This script reads the simulation files and storage the results of the functions.
Then reads the Armse of each simulation and compute the mean per domain.

3_results_marmse: This scripts produces outputs for tables 3 to 7, and Appendices
B, C, and D.

4_results_bias: This scripts produces outputs for tables 9 to 10 and Appendix E.

5_data_cbs: This script reproduces outputs for tables 11 and 12 using data from
CBS. The data from CBS is in the secure environment of CBS and should be required by
email directly to them.
