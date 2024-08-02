
############### Poisson Distribution ########################

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

############## Poisson Distribution Likelihood Function ###################


      # poisson.lik<-function(para1,freq,y){
      #                                     
      #                                     lambda<-para1
      #                                     
      #                                     print(lambda)
      #                                     
      #                                     b1<-log(lambda)*(t(y)%*%freq)
      #                                     
      #                                     print('b1')
      #                                     print(b1)
      #                                     
      #                                     b2=t(freq)%*%log(lfactorial(y+1))
      # 
      #                                     print(b2)
      #                                     
      #                                     sum=0
      #                                     
      #                                     large_N=25000
      #                                     
      #                                     
      #                                     for(i in 1:large_N)
      #                                     {
      #                                       sum=sum+lambda^i/lfactorial(i)
      #                                       
      #                                     }
      #                                     
      #                                     print('sum')
      #                                     print(sum)
      #                                     
      #                                     C=exp(lambda)/sum
      #                                     
      #                                     print(C)
      #                                     
      #                                     
      #                                     if(C==0)
      #                                     {
      #                                       b3=0
      #                                     }
      #                                     else
      #                                     {
      #                                       b3=sum(freq)*log(C)
      #                                     }
      #                                    
      #                                    
      #                                     print('b3')
      #                                     print(b3)
      #                                     
      #                                     b4=sum(freq)*lambda
      #                                     
      #                                     print('b4')
      #                                     print(b4)
      #                                     
      #                                     b5= (b3-b4+b1-b2)
      #                                     
      #                                     print('b5')
      #                                     print(b5)
      #                                     
      #                                     return(-b5)
      #                                }
      # output=optim(1,poisson.lik,freq=freq,y=y)
      # lambdahat=output$par

############ Estimation of parameter ##################

estlambda=(t(freq)%*%y)/(sum(freq))
print('estlambda')
print(estlambda)

################ Poisson Distribution Probability Function ##################

        probability.fun<- function(gamma,datafile){
        
        estlambda<-gamma
        
        xmin1<-1
        
        x<-datafile
        
        sum1=0
        
        large_N=150
        
        
        for(i in 1:large_N)
        {
          sum1=sum1+estlambda^i/factorial(i)
          
        } 
        
        C1=1/sum1
        
        h2<-((estlambda^x)/(factorial(x)))
        
        h3<-C1*h2
        
        return(h3)
        
      }
      
    
      #outfunction<-probability.fun(output$par,y)
      outfunction<-probability.fun(estlambda,y)
      
      print(sum(outfunction[1:150]))
      
      finaloutput<-outfunction*var1
      
      #print(finaloutput)
      
      write.csv(finaloutput,'output_bio-SC-HT_poisson.csv')
      
      
   
