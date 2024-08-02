# Learning Patterns from Biological Networks: A Compounded Burr Probability Model

In this study, we consider Biological network data sets that are publicly available at https://networkrepository.com/bio.php

Usage of the repository for the reproducibility of the paper: 

1. In this repository, we present examples with 1 biological data set, namely #bio-SC-HT. These data set can be directly imported to R statistical software for the analysis of degree distribution of networks data. For example, the repository contains the following data sets (one data example from 10 different domains) for presentation: ego-Twitter(In) (Social Networks), cit-HepTh(In) (Citation Networks), ca-CondMat (Collaboration Networks), Google(In) (Web graphs), Yeast-PPIN (Biological Networks), amazon0601(In) (Product co-purchasing Networks), sx-mathoverflow(In) (Temporal Networks), Email-Enron (Communication Networks), Wiki-Topcats (Ground-truth Networks), and Human25890-session1 (Brain Networks). The ".csv" file of these data sets contain the frequency-degree for these real-world complex networks. 

2. The "models.R" file contains the implementation of popularly-used degree distributions, namely Lomax, power-law, power-law with cutoff, Log-normal, and Exponential distributions. Furthermore, the file "models.R" also contains the implementations of our proposed "Generalized Lomax" family of distributions, namely GLM Type-I, GLM Type-II, GLM Type-III and GLM Type-IV models. The decsirptions of all these models are provided in the manuscript titled "Searching for a new probability distribution for modeling non-scale-free heavy-tailed real-world networks". 

3. Once the implementation is done, the predicted outputs of Lomax, power-law, power-law with cutoff, Log-normal, Exponential, GLM Type-I, GLM Type-II, GLM Type-III and GLM Type-IV models are restored in dataname_output.csv file. For example, in case of ego-Twitter(In).csv data set; all the predicted values based on different probability models are presented in ego-Twitter(In)_outputs.csv file. This file is further used for the computation of different metrics for finding predictive accuracy of several models in the manuscript. 

4. Using the outputs of dataname_output.csv file, we obtain the graphs (Plots of degree distributions along with different pro
