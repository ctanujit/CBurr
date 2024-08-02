
########## Burr Exponentiated Distribution ###################

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


######### Burr Exponentiated  Likelihood Function ##################

barr_exponiented.lik<-function(vector1,freq){
  
  alpha=vector1[1]
  lambda=vector1[2]
  parak=vector1[3]
  parac=vector1[4]
  
  n<-sum(freq)
  
  total_value=0
  
  large_N=500000
  
  
  for(i in 1:large_N)
  {
    
    
    term1=alpha*parak*parac*(lambda^(-parac))*(i^(parac-1))
    
    term2=(1+(i/lambda)^parac)^(-parak-1)
    
    term3=(1-(1+(i/lambda)^parac)^(-parak))^(alpha-1)
    
    
    total_value = total_value + term1*term2*term3
    #print(total_value)
    
    #print("---------")
  }
  
  C=1/total_value
  
  print(C)
  
  likterm1=n*log(C)
  
  likterm2=n*log(alpha*parak*parac)
  
  likterm3=n*parac*log(lambda)
  
  likterm4=(parac-1)*(t(freq)%*%log(y))
  
  likterm5=(parak+1)*(t(freq)%*%(log(1+(y/lambda)^parac)))
  
  likterm6=(alpha-1)*(t(freq)%*%(log(1-(1+(y/lambda)^parac)^(-parak))))
  
  
  final_likterm<-sum(c(likterm1,likterm2,-likterm3,likterm4,-likterm5,likterm6))
  
  return(-final_likterm)
  
}

output11<-optim(c(1,1,1,1),barr_exponiented.lik, freq = freq)

print(output11)


##################Burr Exponentiated  probability Function #########################


probability.fun_barr_exponiented<- function(vector3,datafile){
  
  estalpha=vector3[1]
  estlambda=vector3[2]
  estparak=vector3[3]
  estparac=vector3[4]
  
  x=datafile
  
  sum2=0
  large_N=500000
  
  for(i in 1:large_N)
  {
    
    pdf_term1=estalpha*estparak*estparac*(estlambda^(-estparac))*(i^(estparac-1))
    
    pdf_term2=(1+(i/estlambda)^estparac)^(-estparak-1)
    
    pdf_term3=(1-(1+(i/estlambda)^estparac)^(-estparak))^(estalpha-1)
    
    
    sum2 = sum2 + pdf_term1*pdf_term2*pdf_term3
    #print(sum2)
    
    #print("---------")
  }
  
  #print('%%%%%%%%%')
  
  C2=1/sum2
  #print(sum2)
  print(C2)
  
  
  pdf_term4=estalpha*estparak*estparac*(estlambda^(-estparac))*(x^(estparac-1))
  
  pdf_term5=(1+(x/estlambda)^estparac)^(-estparak-1)
  
  pdf_term6=(1-(1+(x/estlambda)^estparac)^(-estparak))^(estalpha-1)
  
  final_pdf_term=C2*pdf_term4*pdf_term5*pdf_term6
  
  return(final_pdf_term)
  
  #return(probability_array)
  
}


#outfun<-probability.fun_lomax(c(1000,1),y)
outfun<-probability.fun_barr_exponiented(c(output11$par[1],output11$par[2],output11$par[3],output11$par[4]),y)
print('probability_sum')
print(sum(outfun))
plot(y,outfun)


finaloutput<-outfun*var1
#print(finaloutput)


write.csv(finaloutput,'output_bio-SC-HT_barr_exponiented.csv')




