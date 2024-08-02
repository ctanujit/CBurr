############### Power-Law Cutoff Distribution #########################


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


############### Power-Law Cutoff Distribution Likelihood Function ####################


powerlawcutof.lik<-function(vector1,freq){
alpha<-vector1[1]
lamda<-vector1[2]
xmin<-1

n<-sum(freq)
#print(alpha)
#print(lamda)

##############################

# t1<-xmin
# 
# #integrand<-function(z){exp(-lamda*z)}
# integrand<-function(z){exp(-2*z)}
# 
# #integrand<-function(z){(exp(-lamda*z))*(z^(-alpha))}
# 
# k1<-integrate(integrand, lower=t1, upper=Inf)
# 
# C=1/k1$value
# 
# print(C)

##############################

sum=0

large_N=500000

for(i in 1:large_N)
{
  sum = sum + ((1/(i^alpha))*(1/exp(lamda*i)))
}

C=1/sum

##############################

k2<-n*log(C)


k4<-alpha*(t(freq)%*%log(y))

k5<-lamda*(t(y)%*%freq)

k6<-sum(c(k2,-k4,-k5))

#print(k6)

return(-k6)

}


output11<-optim(c(1,2),powerlawcutof.lik,freq=freq)

print(output11)

################### Power-Law Cutoff Distribution Probability Function #################

probability.fun<- function(vector2,datafile){

estalpha<-vector2[1]

estlamda<-vector2[2]

xmin1<-1

x<-datafile

######################################

# t11<-xmin1
# 
# #integrand<-function(z){exp(-estlamda*z)}
# integrand<-function(z){exp(-2*z)}
# 
# #integrand<-function(z){(exp(-estlambda*z))*(z^(-estalpha))}
# 
# k11<-integrate(integrand, lower=t11, upper=Inf)
# 
# C1=1/k11$value

########################################

sum1=0

large_N=500000

for(i in 1:large_N)
{
  sum1 = sum1 + ((1/(i^estalpha))*(1/exp(estlamda*i)))
}

C1=1/sum1

#####################################

j3<-x^(-estalpha)

j4<-exp(-estlamda*x)

j5<-(C1*j3*j4)

return(j5)

}


outfun<-probability.fun(c(output11$par[1],output11$par[2]),y)
print('probability sum')
print(sum(outfun))
finaloutput<-outfun*var1

#print(finaloutput)


write.csv(finaloutput,'output_bio-SC-HT_powerlaw_cutoff.csv')




