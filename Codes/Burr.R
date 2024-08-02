
########### Burr Distribution #############

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

######### Generate Burr Likelihood function #############

barr.lik<-function(vector1,freq){
  
  parama=vector1[1]
  paramb=vector1[2]
  paramc=vector1[3]
  
  
  n<-sum(freq)
  
  total_value=0
  
  large_N=500000
  
  
  for(i in 1:large_N)
  {
    
    
    term1=(parama^(-paramc))*paramb*paramc
    
    term2=i^(paramc-1)
    
    term3=(1+((i/parama)^(paramc)))^(-paramb-1)
    

    total_value = total_value + term1*term2*term3
    #print(total_value)
    
    #print("---------")
  }
  
  C=1/total_value
  
  print(C)
  
  likterm1=n*log(C)
  
  likterm2=n*log(paramc)
  
  likterm3=n*log(paramb)
  
  likterm4=n*paramc*log(parama)
  
  likterm5=(paramc-1)*(t(freq)%*%log(y))
  
  likterm6=(paramb+1)*(t(freq)%*%(log(1+((y/parama)^(paramc)))))

  final_likterm<-sum(c(likterm1,likterm2,likterm3,-likterm4,likterm5,-likterm6))
  
  return(-final_likterm)
  
}

output11<-optim(c(1,1,1),barr.lik, freq = freq)

print(output11)


########################## Burr Probability Function ########################################


probability.fun_barr<- function(vector3,datafile){
  
  estparama=vector3[1]
  estparamb=vector3[2]
  estparamc=vector3[3]
  
  x=datafile
  
  sum2=0
  large_N=500000
  
  for(i in 1:large_N)
  {
    
    pdf_term1=(estparama^(-estparamc))*estparamb*estparamc
    
    pdf_term2=i^(estparamc-1)
    
    pdf_term3=(1+((i/estparama)^(estparamc)))^(-estparamb-1)
    
  
    sum2 = sum2 + pdf_term1*pdf_term2*pdf_term3
    #print(sum2)
    
    #print("---------")
  }
  
  #print('%%%%%%%%%')
  
  C2=1/sum2
  #print(sum2)
  print(C2)
  
  
  pdf_term4=(estparama^(-estparamc))*estparamb*estparamc
  
  pdf_term5=x^(estparamc-1)
  
  pdf_term6=(1+((x/estparama)^(estparamc)))^(-estparamb-1)
  

  final_pdf_term=C2*pdf_term4*pdf_term5*pdf_term6
  
  return(final_pdf_term)
  
  #return(probability_array)
  
}


#outfun<-probability.fun_lomax(c(1000,1),y)
outfun<-probability.fun_barr(c(output11$par[1],output11$par[2],output11$par[3]),y)
print('probability_sum')
print(sum(outfun))
plot(y,outfun)


finaloutput<-outfun*var1
#print(finaloutput)


write.csv(finaloutput,'output_bio-SC-HT_barr.csv')




