############# Test Statistics Values ##########################

rm(list = setdiff(ls(), lsf.str()))

library("entropy")
#library('DescTools')
library('entropy')
library(nnet)
library(rpart)
#library(neuralnet)
#library(ISLR)
#library(caret)
#library(xts)
#library(zoo)
#library(xlsx)
library(Metrics)
library(MASS)
library(CORElearn)
library(pROC)


data<-read.csv('compared_output_bio-SC-HT_CBurr.csv')


  actual_freq<-data$Actual
  #estimated_freq<-data$LOMAX
  #estimated_freq<-data$barr

  
  estimated_freq<-data$CBurr
  
  
  #estimated_freq<-data$Estimated

  #estimated_freq<-data$Pareto_Estimated1
  #estimated_freq<-data$Lognormal_Estimated
  #estimated_freq<-data$Poisson_Estimated
  #estimated_freq<-data$Drift_Power_Law_Estimated
  #estimated_freq<-data$PowerLaw_Cutoff_Estimated
  #estimated_freq<-data$Exponential_Estimated

  
  X=length(actual_freq)
  
  
  ### REGRESSION STATISTICS  #############
  
  print('MAPE')
  print(smape(actual_freq,estimated_freq))
  #print(mape(actual_freq,estimated_freq))
  print('RMSE')
  print(rmse(actual_freq,estimated_freq))
  print('MAE')
  print(mae(actual_freq,estimated_freq))
  
  
  ############################


  kldivergence<-KL.plugin(actual_freq,estimated_freq)
  
  print('kldivergence')
  
  print(kldivergence)
  
 
  probdistrbn_actual_freq=actual_freq/sum(actual_freq)
  probdistrbn_estimated_freq=estimated_freq/sum(estimated_freq)

  
  




