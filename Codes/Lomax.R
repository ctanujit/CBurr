############### LOMAX Distribution ########################

rm(list = ls())

library(optimx)


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

############# LOMAX Distribution LIkelihood Function #########################

LOMAX.lik<-function(vector1,freq){
  alpha=vector1[1]
  lamda=vector1[2]
  
  
  n<-sum(freq)
  
  total_value=0
  
  large_N=500000
  
  
  for(i in 1:large_N)
  {
    
    
    term1=alpha*(lamda^(alpha))
    #print(term1)
    term2=(i+lamda)^(-(alpha+1))
    #print(term2)
   
    total_value = total_value + term1*term2
    #print(total_value)
    
    #print("---------")
  }
  
  C=1/total_value
  
  print(C)
  
  likterm1=n*log(C)
  
  likterm2=n*log(alpha)
  
  likterm3=n*alpha*log(lamda)
  
  likterm4=(alpha+1)*(t(freq)%*%(log(y+lamda)))
  
  final_likterm<-sum(c(likterm1,likterm2,likterm3,-likterm4))
  
  return(-final_likterm)
  
}

output11<-optim(c(1,1),LOMAX.lik, freq = freq)

print(output11)


################# LOMAX Distribution Probability Function ##########################


probability.fun_lomax<- function(vector3,datafile){
  
  estalpha=vector3[1]
  estlamda=vector3[2]

  x=datafile
  
  sum2=0
  large_N=500000
  
  for(i in 1:large_N)
  {
    
    
    pdf_term1=estalpha*(estlamda^(estalpha))
    #print(term1)
    pdf_term2=(i+estlamda)^(-(estalpha+1))
    
    sum2 = sum2 + pdf_term1*pdf_term2
    #print(sum2)
    
    #print("---------")
  }
  
  #print('%%%%%%%%%')
  
  C2=1/sum2
  #print(sum2)
  print(C2)
  
  pdf_term3=estalpha*(estlamda^(estalpha))
  #print(term1)
  pdf_term4=(x+estlamda)^(-(estalpha+1))
  
  final_pdf_term=C2*pdf_term3*pdf_term4
  
  return(final_pdf_term)
  
  #return(probability_array)
  
}


#outfun<-probability.fun_lomax(c(2,1),y)
outfun<-probability.fun_lomax(c(output11$par[1],output11$par[2]),y)
print('probability_sum')
print(sum(outfun))
plot(y,outfun)


finaloutput<-outfun*var1
#print(finaloutput)


write.csv(finaloutput,'output_bio_SC-HT_lomax_1.csv')




