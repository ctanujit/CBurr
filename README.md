# Learning Patterns from Biological Networks: A Compounded Burr Probability Model

In this study, we consider Biological network data sets that are publicly available at https://networkrepository.com/bio.php

Usage of the repository for the reproducibility of the paper (https://arxiv.org/abs/2407.04465): 

1. In this repository, we present examples with 1 biological data set, namely bio-SC-HT. This data set can be directly imported to R statistical software for the analysis of the degree distribution of biological network data. The ".csv" file of these data sets contains the frequency degree for bio-SC-HT network. 

2. The "model_name.R" file contains the implementation of different distributions, including our proposed CBurr distribution. Other files like chi_square test and test_statistics (.R files) provide the accuracy metrics and statistical test implementations. 

3. Once the implementation is done, the predicted outputs of all the probability models are restored in output_bio-SC-HT_CBurr.csv (for our method). This file is further used for the computation of different metrics for finding the predictive accuracy of several models in the manuscript. Using the outputs of output_dataname_model.csv file, we obtain the graphs (Plots of degree distributions). 
