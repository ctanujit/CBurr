
####### Exponential Distribution #################

rm(list = ls())



data<-read.csv('final_degree_frequency_skitter1.csv')


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

############ Exponential Distribution Likelihood Function ##############


exponential.lik<-function(theta,U){
  lambda<-theta
  xmin<-1
  n<-sum(freq)

  sum=0

  large_N=500000


  for(i in 1:large_N)
  {
    sum = sum + (exp(-lambda*i)*(lambda))
  }

  C=1/sum

  t1<-n*log(C)

  t2<-n*log(lambda)

  t3<-lambda*(t(freq)%*%y)

  t4<- sum(c(t1,t2,-t3))

  return(-t4)

}


output<-optim(1,exponential.lik,U=freq)

print(output)

############### Estimate Parameters ####################

estlambda=(sum(freq))/(t(freq)%*%y)

################# Exponential Distribution Probability Function #############

probability.fun<- function(gamma,datafile){
  
  estlambda<-gamma
  
  xmin1<-1
  
  x<-datafile
  
  sum1=0
  
  large_N=500000
  
  
  for(i in 1:large_N)
  {
    sum1 = sum1 + (exp(-estlambda*i)*(estlambda))
  }
  
  C1=1/sum1
  
  h2<-((estlambda)*(exp(-x*estlambda)))
  
  h3<-C1*h2
  
  return(h3)
  
}

#outfunction<-probability.fun(output$par,y)
outfunction<-probability.fun(estlambda,y)


#print(outfunction)

finaloutput<-outfunction*var1

print(finaloutput)

write.csv(finaloutput,'output_skitter1_exponential.csv')

