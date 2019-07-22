dat<-read.csv('datayigan.csv')
dat<-t(dat$num)/10000
#定义log-likelihood函数
LL<-function(params,data)
{
  t1<-dnorm(data,params[2],params[3])
  t2<-dnorm(data,params[4],params[5])
  f<-params[1]*t1+(1-params[1])*t2
  ll<-sum(log(f))
  return(-ll)
}

hist(dat,freq=F)
lines(density(dat))
geyser.res<-nlminb(c(0.5,9,10,10,10),LL,data=dat,
                   lower=c(0.0001,-Inf,0.0001,-Inf,0.0001),
                   upper=c(0.9999,Inf,Inf,Inf,Inf))
#LL是被最小化的函数。
geyser.res$par
X<-seq(8,12,length=1000)
p<-geyser.res$par[1]
mu1<-geyser.res$par[2]
sig1<-geyser.res$par[3]
mu2<-geyser.res$par[4]
sig2<-geyser.res$par[5]

f<-p*dnorm(X,mu1,sig1)+(1-p)*dnorm(X,mu2,sig2)
#作出数据的直方图
hist(dat,probability=T,col=0,ylab="Density",
     ylim=c(0,1),xlab="Eruption waiting times")
lines(X,f)


library(lhs)
n=1000
A<-randomLHS(n,8)
A[,7:8]<-A[,7:8]*0.38+0.0002

A[1,]<-c(0.4458,0.030645,0.023041,0.021969,0.082479,0.085013)

R=A[,1]*A[,2]/A[,5]+(1-A[,1])*A[,3]/(A[,6]+A[,4])
R[1]

library(randomForest)
D<-matrix(0,n,9)
colnames(D)<-c('R_0','q','lambd_i','lambd_e','mu','beta_i','beta_e','i_0','e_0')

D[,1]<-R
for (i in 1:6) {
  D[,i+1]<-A[,i]
}

Forest<-randomForest(R_0~.,D,mtry=4,importance=T,proximity=T)
print(importance(Forest)) 
E<-importance(Forest)

write.table(D,"./Dyg.dat",row.names =F,col.names =F)
write.table(E,"./Eyg.dat",row.names =F,col.names =F)




