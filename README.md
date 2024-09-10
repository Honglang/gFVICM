# gFVICM: Generalized Functional Varying Index Coefficient Model for Dynamic Gene-environment Interactions

Code for a paper on Generalized functional varying-index coefficient model for dynamic synergistic gene-environment interactions with binary longitudinal traits

* Authors: Jingyi Zhang, Honglang Wang, Yuehua Cui.
* In case of questions or comments please contact cuiy@msu.edu!

* The code was written/evaluated in R with the following software versions: R version 4.3.1 (2023-06-16 ucrt)

* This folder contains the following data and files that can be used to reproduce all analysis and figures of the manuscript. It contains three subfolders containing the following files:


### ./Real_Data_Analysis/:
    * functions_real_data.R
    An R script contains all of the functions used for the real data analysis in the manuscript. 
    
    * case_study.R
    An R script that contains the code of the real data analysis reported in the paper. This is the code
    that yields exactly the results (figures, tables, and in-text information) as reported in the manuscript 
    when being applied to the original data, which was used for analysis in the manuscript but can not be published 
    online due to data privacy reasons. Interested readers can contact the corresponding author to obtain access to 
    the original data. Theo riginal data has been made available to RR editors strictly for the purpose of the
    reproducibility audit in a separate file.
    
### ./Fig_Ex_With_Intermediate_Data/:
    * fig_plotting.R
    An R script that was used to generate one typical figure about confidence bands of the curve estimation of
    the manuscript.
    * All other six .txt files are intermediate results for the figure plotting.
    
### ./Simulation/
    * functions.R
    An R script contains all of the functions used for the estimation and testing problems in the manuscript.

    * main.R
    An R script that performs the simulations about estimation and testing reported in the manuscript. Simulation was
    performed in parallel on a linux server, but should also work under windows. The random number generation seeds 
    were 1 to 500.

  
