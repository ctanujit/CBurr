############### Power-Law distribution #####################


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

x_min=1

########### Estimation of Parameters #############################

## continuous case
estalpha<-1+var1*((freq%*%log(y/x_min))^(-1))
#altestaplpha<-1+var1*((sum(log(rep(y1,freq1))))^(-1))


## discrete case
#estalpha<-1+var1*((freq%*%log(y/(x_min-0.5)))^(-1))

#estalpha<-2
#estalpha<-1.5

##################### Power-Law distribution Probability Function #####################

probability.fun<- function(gamma,datafile){
  
  estalpha1<-gamma
  
  x<-datafile
  
  large_N=500
  
  c1<-1/sum(1/((1:large_N)^(estalpha1)))
  
  h1<-c1/x^(estalpha1)
  
  #h1<-((estalpha1-1)/x_min)*((x/x_min)^(-estalpha1))
    
  return(h1)
  
}


outfunction<-probability.fun(estalpha,1:max(y1))

finaloutput<-outfunction*var1

#print(finaloutput)

print(sum(outfunction))

print('=============')

#print('actual')
#print(freq)


#print('predicted')
#print(finaloutput)


write.csv(finaloutput,'output_bio_SC-HT_poweralw.csv')


