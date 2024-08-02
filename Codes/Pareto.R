
################# Pareto Distribution #################

rm(list = ls())



data<-read.csv('final_degree_frequency_bio-SC-HT.csv')





y1=data$C1
freq1<-data$C2

fulldata<-numeric(max(y1))

for(k in 1:max(y1))
{
  if(is.element(k,y1)==TRUE)
  {
    index=match(k,y1)
    fulldata[k]=freq1[index]
  }
  else
    fulldata[k]=0
}

y<-1:max(y1)
freq<-fulldata

var1<-sum(freq)
#x_min=1

estxmin1<-1

n<-sum(freq)

########### Estiamtion of Parameter ######################

estalpha1<-n/(t(freq)%*%(log(y)-log(estxmin1)))

######## Pareto Distribution Probability Function ##################

probability.fun<- function(gamma,datafile){

estxmin<-gamma[1]

estalpha<-gamma[2]

x<-datafile


large_N=600000

c1=1/(sum(1/((1:large_N)^(estalpha+1)))*(estalpha*(estxmin^estalpha)))

term1<-estalpha

term2<-estxmin^estalpha

term3<-x^(estalpha+1)


term4<-term2/term3

#term5<-term1*term4
term5<-term1*term4*c1


return(term5)

}


outfunction<-probability.fun(c(estxmin1,estalpha1),1:max(y1))

finaloutput<-outfunction*var1

print(sum(outfunction))

write.csv(finaloutput,'output_bio_SC-Ht_pareto.csv')

#print(finaloutput)
