Master's thesis

Sofía Villalobos Aliste (6060714)
s.f.villalobosaliste@students.uu.nl

Research repository of thesis manuscript produced between September 2021 and May 2022. 

This repository contains all the information regarding the master thesis and internship project at the
Centraal Bureau Voor de Statstiek "Combining probability and non-probability samples for estimation", 
of the program Methodology and Statistics for the Behavioural, Biomedical, and Social Sciences, 
Utrecht University year 2022.

In this study, we seek to find a method that combines a probability and a non-probability sample to reduce the 
mean squared error (MSE) of an estimator of proportions in comparison to the MSEs of the separate estimators 
from the probability and the non-probability sample. We performed simulation studies where two different 
methods of modelling the bias of an estimator based on the non-probability sample were tested, and the MSE 
was evaluated. We then applied these methods to a real dataset from Statistics Netherlands (CBS) to show 
that the MSEs can indeed be reduced. 

#######################
Three folders are included here:

"Pdf_files" folder which includes: 

- The approval files from the Ethical Review Board of the Faculty of Social and Behavioural Sciences of 
Utrecht University for the simulation study and the use of existing data.

- The thesis manuscript.

"Outputs" folder which include:

- The results produced in R (version 4.1.2) saved in excel and rds files divided into three folders:

  - "Alternative approach": This output is one summary table of an alternative approach that is mentioned at the end 
of the section 3.1 of the thesis manuscript, but its results are not included in it.

  - "Bias": These outputs were used to construct tables 8 to 10 and Appendix E of the thesis manuscript.

  - "MARMSE": These outputs were used to construct tables 3 to 7, and Appendices B, C, and D of the thesis manuscript.

- And 

  - Template: This is an excel file that is used as a template for tables 3 to 10 and Appendices B, C, D, and E
of the thesis manuscript.

"Scripts" folder which include:

- The code with which the results were produced. 

- 6 RDS files to replicate the outputs using the simulation results. The files needed are noted
at the beginning of the scripts.

The tables included in the thesis manuscript are constructed with the information of the outputs saved in the 
file "Outputs" and can be replicated by running the scripts in the following order:

1_create_function: This script creates the necessary functions to generate the simulated samples and compute 
three different estimators for the probability, non-probability, and combined sample.

2_storage: This script runs the simulation files and storages 12 elements of the simulation on a list 
(description on the script). Then it computes the mean of the ARMSE per domain, and the mean absolute 
bias per domain. Running this script takes around 7 days.

3_results_marmse: This script produces the outputs found in the folder "Outputs" -> "MARMSE". The results can be 
replicated by using the simulation results cond_eq.rds and cond_uneq.rds.

4_results_bias: This script produces the outputs found in the folder "Outputs" -> "Bias". The results can be 
replicated by using the simulation results bias_eq_1.rds, bias_eq_2.rds, bias_uneq_1,rds, and bias_uneq_2.rds

5_data_cbs: This script produces the outputs for tables 11 and 12 of the thesis manuscript using data from CBS. 
The data from CBS is in a secure environment of CBS and should be required by email directly to them.*

6_rnorm: This script includes an alternative approach mentioned at the end of section 3.1 of the thesis manuscript. 
This script is self-sufficient and includes the functions to generate the samples and compute the estimators, 
the storage of the results of the MARMSE, and a summary table.

#######################
*to get access visit: https://www.cbs.nl/en-gb/our-services/customised-services-microdata/microdata-conducting-your-own-research or email microdata@cbs.nl for more information.
#######################
Postprocessing:
To obtain tables 3 to 10 of the thesis manuscript, and Appendices B, C, D, and E, all the output tables with 
dimensions of 16x16 are plugged into the template of the excel file "Template" found in the folder "Outputs".
