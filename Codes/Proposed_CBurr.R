################ Proposed CBurr Distribution #################################

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


################ Proposed CBurr Distribution Likelihood Function #######################

CBurr.lik<-function(vector1,freq){
  
  parama=vector1[1]
  paramb=vector1[2]
  paramc=vector1[3]
  theta=vector1[4]
  
  
  n<-sum(freq)
  
  total_value=0
  
  large_N=500000
  
  
  for(i in 1:large_N)
  {
    
    
    term1=(parama^(-paramc))*paramb*paramc
    
    term2=i^(paramc-1)
    
    term3=(1+((i/parama)^(paramc)))^(-paramb-1)
    
    term4=exp(-theta*(1-(1+(i/parama)^paramc)^(-paramb)))
    
    term5=1+theta*((1+(i/parama)^paramc)^(-paramb))

    total_value = total_value + term1*term2*term3*term4*term5
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
  
  likterm7=theta*(t(freq)%*%(1-(1+(y/parama)^paramc)^(-paramb)))
  
  likterm8=t(freq)%*%(log(1+theta*((1+(y/parama)^paramc)^(-paramb))))

  final_likterm<-sum(c(likterm1,likterm2,likterm3,-likterm4,likterm5,-likterm6,-likterm7,likterm8))
  
  return(-final_likterm)
  
}

output11<-optim(c(1,1,3.3,0.3),CBurr.lik, freq = freq)

## best parameters found : (1,1,3.3,0.3)/(1,1,2,1)
##

print(output11)


#########################Proposed CBurr Distribution Probability Function ##################


probability.fun_CBurr<- function(vector3,datafile){
  
  estparama=vector3[1]
  estparamb=vector3[2]
  estparamc=vector3[3]
  esttheta=vector3[4]
  
  x=datafile
  
  sum2=0
  large_N=500000
  
  for(i in 1:large_N)
  {
    
    pdf_term1=(estparama^(-estparamc))*estparamb*estparamc
    
    pdf_term2=i^(estparamc-1)
    
    pdf_term3=(1+((i/estparama)^(estparamc)))^(-estparamb-1)
    
    pdf_term4=exp(-esttheta*(1-(1+(i/estparama)^estparamc)^(-estparamb)))
    
    pdf_term5=1+esttheta*((1+(i/estparama)^estparamc)^(-estparamb))
    
    
  
    sum2 = sum2 + pdf_term1*pdf_term2*pdf_term3*pdf_term4*pdf_term5
    #print(sum2)
    
    #print("---------")
  }
  
  #print('%%%%%%%%%')
  
  C2=1/sum2
  #print(sum2)
  print(C2)
  
  
  pdf_term6=(estparama^(-estparamc))*estparamb*estparamc
  
  pdf_term7=x^(estparamc-1)
  
  pdf_term8=(1+((x/estparama)^(estparamc)))^(-estparamb-1)
  
  pdf_term9=exp(-esttheta*(1-(1+(x/estparama)^estparamc)^(-estparamb)))
  
  pdf_term10=1+esttheta*((1+(x/estparama)^estparamc)^(-estparamb))
  

  final_pdf_term=C2*pdf_term6*pdf_term7*pdf_term8*pdf_term9*pdf_term10
  
  return(final_pdf_term)
  
  #return(probability_array)
  
}


#outfun<-probability.fun_lomax(c(1000,1),y)
outfun<-probability.fun_CBurr(c(output11$par[1],output11$par[2],output11$par[3],output11$par[4]),y)
print('probability_sum')
print(sum(outfun))
plot(y,outfun)


finaloutput<-outfun*var1
#print(finaloutput)


write.csv(finaloutput,'output_bio-SC-HT_CBurr.csv')




