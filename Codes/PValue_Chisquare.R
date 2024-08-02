
######### Bootstrap Chi-Square Value ################

rm(list = ls())



data<-read.csv('compared_output_bio-SC-HT_superm1.csv')


y1=data$Degree
actual_freq<-data$Actual


#estimated_freq<-data$DGPL
#estimated_freq<-data$LOMAX
#estimated_freq<-data$barr_type5

estimated_freq<-data$CBurr

#estimated_freq<-data$Estimated1
#estimated_freq<-data$Estimated2
#estimated_freq<-data$PowerLaw_Estimated
#estimated_freq<-data$PowerLaw_Estimated1
#estimated_freq<-data$PowerLaw_Estimated2
#estimated_freq<-data$PowerLaw_Estimated4
#estimated_freq<-data$Pareto_Estimated
#estimated_freq<-data$Pareto_Estimated2
#estimated_freq<-data$Lognormal_Estimated
#estimated_freq<-data$PowerLaw_Cutoff_Estimated
#estimated_freq<-data$Exponential_Estimated



actual_chisquare_value<-sum(((actual_freq-estimated_freq)^2)/estimated_freq)


array_of_synthetic_chisquare_value<-rep(0,10000)



for(i in 1:10000)
{
print(i)
synthetic_data=sample(1:max(y1),sum(actual_freq),prob=actual_freq,rep=T)

frequ_table_synthetic_data<-as.data.frame(table(synthetic_data))

unique_degree<-as.numeric(as.character(frequ_table_synthetic_data$synthetic_data))
unique_freq<-as.numeric(as.character(frequ_table_synthetic_data$Freq))


y<-1:max(y1)
synthetic_freq<-rep(0,max(y1))
synthetic_freq[unique_degree]<-unique_freq


synthetic_chisquare_value<-sum(((synthetic_freq-estimated_freq)^2)/estimated_freq)


array_of_synthetic_chisquare_value[i]<-synthetic_chisquare_value

#print(synthetic_chisquare_value)

}

hist(array_of_synthetic_chisquare_value)

print(mean(array_of_synthetic_chisquare_value>actual_chisquare_value))



