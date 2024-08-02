
################ Log-normal Distribution ##########################


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
n=var1


########  estimated parameters  Log-normal Distribution ###

estmu=(t(freq)%*%log(y))/n

estsigma=sqrt((t(freq)%*%((log(y)-estmu)^2))/n)


########## OR / Log-normal Distribution Likelihood Function #################

lognormal.lik<-function(theta,U){
mu<-theta[1]
sigma<-theta[2]
xmin<-1
n<-sum(freq)

#print(xmin)

y1<-((log(xmin)-mu)/(sqrt(2)*sigma))

#print(y1)

integrand <- function(x){exp(-x^2)}

t1<-integrate(integrand, lower = y1, upper = Inf)

# t11<-(2/sqrt(pi))

t2<-n*log(sqrt(2/pi))

#print(sigma)

t3<- n*log(sigma)

t4<- t(y)%*%freq

t5<- (t(freq)%*%(log(y)-mu)^2)/(2*sigma^2)

t6<- n*log(t1$value)

t7<- sum(c(t2,-t3,-t4,-t5,-t6))

return(-t7)

}

output<-optim(c(1,2),lognormal.lik,U=freq)

print(output)

############### Log-normal Distribution Probability Function #################

probability.fun<- function(gamma,datafile){

estmu1<-gamma[1]

estsigma1<-gamma[2]

xmin1<-1

x<-datafile


# y2<-((log(xmin1)-estmu)/(sqrt(2)*estsigma))
# 
# integrand1 <- function(x1){exp(-x1^2)}
# 
# h1<-integrate(integrand1, lower = y2, upper = Inf)
# 
# 
# h11<-(2/sqrt(pi))
# 
# h2<-(1/(h11*h1$value))

#print(h2)


h3<-1/(sqrt(2*pi)*(estsigma))

#print(h3)

large_N=500
c1<-1/(h3*sum(exp(-((log(1:large_N)-estmu1)/(sqrt(2)*estsigma1))^2)/(1:large_N)))



h4<-x

#print(h4)

h5<-exp(-((log(x)-estmu1)/(sqrt(2)*estsigma1))^2)

h6<-h5/h4

#print(h6)

#h7<-h3*h6
h7<-h3*h6*c1

#print(h7)

return(h7)

}

outfunction<-probability.fun(c(estmu,estsigma),1:max(y1))

print(sum(outfunction))

finaloutput<-outfunction*var1

#print(finaloutput)

write.csv(finaloutput,'output_bio_SC-HT_lognormal.csv')


