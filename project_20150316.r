setwd("C:/Users/User/Dropbox/PAPERS/projects/brownian")

rm(list=ls())		# Clear objects
cransite='http://cran.cnr.Berkeley.edu'
pckgs<-c("moments","MASS","pscl","mnormt")
for (pck in pckgs){
	if(!pck%in%rownames(installed.packages())){
		install.packages(pck,repos=cransite)}	
}
	# Install packages if not installed
lapply(pckgs,require,character.only=T)
seed<-(7079305)	# Fix the random seed



### 0. Functions
r1<-function(x) round(x,digits=1)
r2<-function(x) round(x,digits=2)
r3<-function(x) round(x,digits=3)
r4<-function(x) round(x,digits=4)
stat95<-function(x,marg){
	Mean<-apply(x,marg,mean)
	SE<-apply(x,marg,FUN=function(u) sqrt(var(u)))
	qnt<-apply(x,marg,function(y) quantile(y,c(.025,.25,.5,.75,.975)))
	stats<-cbind(Mean,SE,t(qnt))
	return(stats)
}
preset<-function() par(mfrow=c(1,1))



### 1. Data
snp0<-snp00<-read.csv("SnP500_09-15_d.csv")
snp0$dt<-as.Date(snp0[,1],"%m/%d/%Y")
snp0$x<-snp0$Close; snp0<-snp0[,-(1:2)]
snp0$year<-as.numeric(format(snp0$dt,"%Y"))
snp_15<-snp0[snp0$year==2015,]
snp_dat<-snp0[snp0$year<2015,]

# Save plot to PDF (instead of EPS)
pdf("./output/raw.pdf", width = 8, height = 4.5)
par(mfrow = c(1, 2))
plot(snp_dat$dt, snp_dat$x, type = "l", xlab = "Date", ylab = "Closing Price (X)")
plot(snp_dat$dt, snp_dat$y, type = "l", xlab = "Date", ylab = "Log-Closing Price (Y)")
dev.off()
preset()



### 2. GBM - MCMC
mcmc1<-function(y=y,t=t,nsim=ns0,nburnin=nb0){
	nchain<-nburnin+nsim; set.seed(seed)
	y.d<-y[-length(y)]-y[-1]
	#t.d<-as.numeric(t[-length(t)]-t[-1])/365.25
	t.d<-rep(1,length(y.d))/252
	n<-length(y.d)
	
	# Initialization
	th.ls<-sig2.ls<-rep(NA,nchain)
	a0<-2; b0<-.001; sig0<-100
	th.ls[1]<-0; sig2<-sig2.ls[1]<-10^(-3)
	
	# Update theta
	samp_th<-function(sig2=sig2){
		sig2.th<-1/(1/sig0+sum(t.d)/sig2)
		mu.th<-sig2.th*(sum(y.d)/sig2)
		rnorm(1,mu.th,sqrt(sig2.th))
	}
	
	# Update sigma^2
	samp_sig2<-function(th=th){
		a.sig2<-a0+n/2
		b.sig2<-b0+sum((y.d-th*t.d)^2/t.d)/2
		rigamma(1,a.sig2,b.sig2)
	}
	
	# MCMC
	for (k in 2:nchain){
		th.ls[k]<-th<-samp_th(sig2=sig2)
		sig2.ls[k]<-sig2<-samp_sig2(th=th)
		if(k%%500==0) print(k)
	}
	
	# Transform back
	mu.ls<-th.ls+sig2.ls/2
	sig.ls<-sqrt(sig2.ls)
	
	return(list(th=th.ls,sig2=sig2.ls,mu=mu.ls,sig=sig.ls))
}

nb0<-1000; ns0<-5000
y<-snp_dat$y; t<-snp_dat$dt

out1<-mcmc1(y=y,t=t,nsim=ns0,nburnin=nb0)

pdf("./output/par1.pdf", width = 9, height = 3.25)
par(mfrow = c(1, 3))
hist(out1$mu[-(1:nb0)], main = "", prob = T, xlab = expression(mu))
lines(density(out1$mu[-(1:nb0)]))
hist(out1$sig[-(1:nb0)], main = "", prob = T, xlab = expression(sigma))
lines(density(out1$sig[-(1:nb0)]))
pacf(out1$mu[-(1:nb0)], main = expression("PACF of "*mu))
dev.off()
preset()

stats1.mu<-stat95(as.matrix(out1$mu[-(1:nb0)]),2)[,c(1,2,3,5,7)]
stats1.sig<-stat95(as.matrix(out1$sig[-(1:nb0)]),2)[,c(1,2,3,5,7)]
rbind(stats1.mu,stats1.sig)
stats1.sig[2]/stats1.sig[1]



# MLE
y.d<-y[-length(y)]-y[-1]
#t.d<-as.numeric(t[-length(t)]-t[-1])/365.25
t.d<-rep(1,length(y.d))/252
th.hat<-sum(y.d)/sum(t.d); th.hat
sig2.hat<-mean((y.d-th.hat*t.d)^2/t.d)
th.hat+sig2.hat/2; sqrt(sig2.hat)



# Fitted value
y0<-tail(y,1); n<-length(y)-1; t.d<-rep(1,n)/252
th<-out1$th[-(1:nb0)]; sig<-out1$sig[-(1:nb0)]
set.seed(seed+1)
Bt_i<-function(i) sig[i]*sqrt(t.d)*rnorm(n)
diffuse<-t(sapply(1:ns0,Bt_i))
tht_i<-function(i) th[i]*t.d
drift<-t(sapply(1:ns0,tht_i))
x.pred0<-exp(y0+t(apply(drift+diffuse,1,cumsum)))
x.pred0.ci<-apply(x.pred0,2,FUN=function(x) quantile(x,c(0.05,0.95)))
x.pred0.mean<-apply(x.pred0,2,mean)

## 1.2 Prediction
y0<-head(y,1); n.new<-dim(as.matrix(snp_15))[1]
t.d.new<-rep(1,n.new)/252
set.seed(seed+1)
Bt_i<-function(i) sig[i]*sqrt(t.d.new)*rnorm(n.new)
diffuse<-t(sapply(1:ns0,Bt_i))
tht_i<-function(i) th[i]*t.d.new
drift<-t(sapply(1:ns0,tht_i))
x.pred<-exp(y0+t(apply(drift+diffuse,1,cumsum)))
x.pred.ci<-apply(x.pred,2,FUN=function(x) quantile(x,c(0.05,0.95)))
x.pred.mean<-apply(x.pred,2,mean)

pdf("./output/pred1.pdf", width = 8, height = 4.5)
par(mfrow=c(1,2))
# Fitted
plot(snp_dat$dt,snp_dat$x,type="l",xlab="Date",ylab="Closing Price (X)",
     ylim=c(500,5000))
lines(rev(snp_dat$dt),c(tail(snp_dat$x,1),x.pred0.mean),
      lty=4,col="blue",lwd=2)
lines(rev(snp_dat$dt),c(tail(snp_dat$x,1),x.pred0.ci[1,]),
      lty=2,col="blue",lwd=2)
lines(rev(snp_dat$dt),c(tail(snp_dat$x,1),x.pred0.ci[2,]),
      lty=2,col="blue",lwd=2)
# Pred
plot(head(snp0$dt,100),head(snp0$x,100),type="l",xlab="Date",
     ylim=c(1800,2400),ylab="Closing Price (X)",xaxt="n")
# axis.dt<-as.Date(c("01Sep2014","01Oct2014","01Nov2014",
#                    "01Dec2014","01Jan2015","01Feb2015","01Mar2015"),
#                  format="%d%b%Y")
# axis.Date(1,at=axis.dt,labels=format(axis.dt,"%b-%Y"))
for (i in 1:6){
     lines(c(snp_dat$dt[1],rev(snp_15$dt)),
           c(snp_dat$x[1],x.pred[i,]),col="grey65")
}
lines(c(snp_dat$dt[1],rev(snp_15$dt)),lty=4,
      c(snp_dat$x[1],x.pred.mean),col="blue",lwd=2)
lines(c(snp_dat$dt[1],rev(snp_15$dt)),lty=2,
      c(snp_dat$x[1],x.pred.ci[1,]),col="blue",lwd=2)
lines(c(snp_dat$dt[1],rev(snp_15$dt)),lty=2,
      c(snp_dat$x[1],x.pred.ci[2,]),col="blue",lwd=2)
lines(c(snp_15$dt,snp_dat$dt[1]),c(snp_15$x,snp_dat$x[1]),
      col="red",lwd=2)
dev.off()
preset()




### 2. GBM w/ Jump - MCMC
mcmc2<-function(y=y,nsim=ns0,nburnin=nb0){
	nchain<-nburnin+nsim; set.seed(seed)
	y.d<-rev(y)[-1]-rev(y)[-length(y)]; n<-length(y.d)
	t.d<-rep(1,n)/252
	
	# Initialization
	th.ls<-sig2.ls<-lamt.ls<-muz.ls<-sigz2.ls<-rep(NA,nchain)
	jt.ls<-zt.ls<-matrix(NA,nchain,n)
	a0<-2; b0<-.001; a.lam0<-b.lam0<-1; sig0<-100
	th<-th.ls[1]<-muz<-muz.ls[1]<-0
	sig2<-sig2.ls[1]<-lamt<-lamt.ls[1]<-sigz2<-sigz2.ls[1]<-10^(-3)
	jt<-jt.ls[1,]<-zt<-zt.ls[1,]<-0
	
	# Update zt
	samp_zt<-function(th=th,sig2=sig2,jt=jt,muz=muz,sigz2=sigz2){
		sig2.zt<-1/(1/sigz2+jt^2/sig2/t.d)
		mu.zt<-sig2.zt*(muz/sigz2+jt*(y.d-th*t.d)/sig2/t.d)
		rnorm(n,mu.zt,sqrt(sig2.zt))
	}
	
	# Update jt
	samp_jt<-function(th=th,sig2=sig2,lamt=lamt,zt=zt){
		p1<-lamt*dnorm(y.d,th*t.d+zt,sqrt(sig2*t.d))
		p0<-(1-lamt)*dnorm(y.d,th*t.d,sqrt(sig2*t.d))
		rbinom(n,1,prob=p1/(p1+p0))
	}
	
	# Update lambda*dt
	samp_lamt<-function(jt=jt){
		rbeta(1,a.lam0+sum(jt),b.lam0+sum(1-jt))
	}
	
	# Update theta
	samp_th<-function(sig2=sig2,jt=jt,zt=zt){
		sig2.th<-1/(1/sig0+sum(t.d)/sig2)
		mu.th<-sig2.th*(sum(y.d-jt*zt)/sig2)
		rnorm(1,mu.th,sqrt(sig2.th))
	}
	
	# Update sigma^2
	samp_sig2<-function(th=th,jt=jt,zt=zt){
		a.sig2<-a0+n/2
		b.sig2<-b0+sum((y.d-th*t.d-jt*zt)^2/t.d)/2
		rigamma(1,a.sig2,b.sig2)
	}
	
	# Update mu.z
	samp_muz<-function(sigz2=sigz2,jt=jt,zt=zt){
		#sig2.muz<-1/(1/sig0+sum(jt)/sigz2)
		#mu.muz<-sig2.muz*sum(jt*zt)/sigz2
		sig2.muz<-1/(1/sig0+n/sigz2)
		mu.muz<-sig2.muz*sum(zt)/sigz2	
		rnorm(1,mu.muz,sqrt(sig2.muz))
	}
	
	# Update sigz^2
	samp_sigz2<-function(muz=muz,jt=jt,zt=zt){
		#a.sigz2<-a0+sum(jt)/2
		#b.sigz2<-b0+sum((jt*(zt-muz))^2)/2
		a.sigz2<-a0+n/2
		b.sigz2<-b0+sum((zt-muz)^2)/2
		rigamma(1,a.sigz2,b.sigz2)
	}
	
	# MCMC
	for (k in 2:nchain){
		zt.ls[k,]<-zt<-samp_zt(th=th,sig2=sig2,jt=jt,
						muz=muz,sigz2=sigz2)
		jt.ls[k,]<-jt<-samp_jt(th=th,sig2=sig2,lamt=lamt,zt=zt)
		lamt.ls[k]<-lamt<-samp_lamt(jt=jt)
		
		th.ls[k]<-th<-samp_th(sig2=sig2,jt=jt,zt=zt)
		sig2.ls[k]<-sig2<-samp_sig2(th=th,jt=jt,zt=zt)
		
		muz.ls[k]<-muz<-samp_muz(sigz2=sigz2,jt=jt,zt=zt)
		sigz2.ls[k]<-sigz2<-samp_sigz2(muz=muz,jt=jt,zt=zt)
		
		if(k%%500==0) print(k)
	}
	
	# Transform back
	mu.ls<-th.ls+sig2.ls/2
	sig.ls<-sqrt(sig2.ls)
	lam.ls<-lamt.ls/(1/252)
	
	return(list(th=th.ls,sig2=sig2.ls,mu=mu.ls,sig=sig.ls,
			lamt=lamt.ls,lam=lam.ls,jt=jt.ls,zt=zt.ls,
			muz=muz.ls,sigz2=sigz2.ls))
}

nb0<-5000; ns0<-5000
y<-snp_dat$y; t<-snp_dat$dt

out2<-mcmc2(y=y,nsim=ns0,nburnin=nb0)

pdf("./output/par2.pdf",width=9,height=3.25)
par(mfrow=c(1,3))
hist(out2$mu[-(1:nb0)],main="",prob=T,xlab=expression(mu))
lines(density(out2$mu[-(1:nb0)]))
hist(out2$sig[-(1:nb0)],main="",prob=T,xlab=expression(sigma))
lines(density(out2$sig[-(1:nb0)]))
pacf(out2$mu[-(1:nb0)],main=expression("PACF of "*mu))
dev.off()
preset()

stats2.mu<-stat95(as.matrix(out2$mu[-(1:nb0)]),2)[,c(1,2,3,5,7)]
stats2.sig<-stat95(as.matrix(out2$sig[-(1:nb0)]),2)[,c(1,2,3,5,7)]
stats2.muz<-stat95(as.matrix(out2$muz[-(1:nb0)]),2)[,c(1,2,3,5,7)]
stats2.sigz<-stat95(as.matrix(sqrt(out2$sigz2)[-(1:nb0)]),2)[
			,c(1,2,3,5,7)]
stats2.lamt<-stat95(as.matrix(out2$lamt[-(1:nb0)]),2)[,c(1,2,3,5,7)]
rbind(stats2.mu,stats2.sig,stats2.muz,stats2.sigz,stats2.lamt)
stats2.sig[2]/stats2.sig[1]

# Fitted value
y0<-tail(y,1); n<-length(y)-1; t.d<-rep(1,n)/252
th<-out2$th[-(1:nb0)]; sig<-out2$sig[-(1:nb0)]
muz<-out2$muz[-(1:nb0)]; sigz<-sqrt(out2$sigz2[-(1:nb0)])
jt<-out2$jt[-(1:nb0),]; zt<-out2$zt[-(1:nb0),]
lamt<-out2$lamt[-(1:nb0)]; set.seed(seed+1)
Bt_i<-function(i) sig[i]*sqrt(t.d)*rnorm(n)
diffuse<-t(sapply(1:ns0,Bt_i))
tht_i<-function(i) th[i]*t.d
drift<-t(sapply(1:ns0,tht_i))
jz_i<-function(i) jt[i,]*zt[i,]
jump<-t(sapply(1:ns0,jz_i))
x.pred0j<-exp(y0+t(apply(drift+diffuse+jump,1,cumsum)))
x.pred0j.ci<-apply(x.pred0j,2,FUN=function(x) quantile(x,c(0.05,0.95)))
x.pred0j.mean<-apply(x.pred0j,2,mean)

## 2.2 Prediction
y0<-head(y,1); n.new<-dim(as.matrix(snp_15))[1]
t.d.new<-rep(1,n.new)/252
set.seed(seed+1)
Bt_i<-function(i) sig[i]*sqrt(t.d.new)*rnorm(n.new)
diffuse<-t(sapply(1:ns0,Bt_i))
tht_i<-function(i) th[i]*t.d.new
drift<-t(sapply(1:ns0,tht_i))
jz_i<-function(i) rbinom(n.new,1,lamt[i])*rnorm(n.new,muz[i],sigz[i])
jump<-t(sapply(1:ns0,jz_i))
x.predj<-exp(y0+t(apply(drift+diffuse+jump,1,cumsum)))
x.predj.ci<-apply(x.predj,2,FUN=function(x) quantile(x,c(0.05,0.95)))
x.predj.mean<-apply(x.predj,2,mean)


# Plots
pdf("./output/pred2.pdf", width=8,height=4.5)
par(mfrow=c(1,2))
# Fitted
plot(snp_dat$dt,snp_dat$x,type="l",xlab="Date",ylab="Closing Price (X)",
		ylim=c(500,5000))
lines(rev(snp_dat$dt),c(tail(snp_dat$x,1),x.pred0j.mean),
		lty=1,col="blue",lwd=1)
lines(snp_dat$dt,snp_dat$x)
lines(rev(snp_dat$dt),c(tail(snp_dat$x,1),x.pred0j.ci[1,]),
		lty=6,col="blue",lwd=1)
lines(rev(snp_dat$dt),c(tail(snp_dat$x,1),x.pred0j.ci[2,]),
		lty=6,col="blue",lwd=1)
# Pred
plot(head(snp0$dt,100),head(snp0$x,100),type="l",xlab="Date",
	ylim=c(1800,2400),ylab="Closing Price (X)",xaxt="n")
# axis.dt<-as.Date(c("01sep2014","01oct2014","01nov2014",
# 			"01dec2014","01jan2015","01feb2015","01mar2015"),
# 			format="%d%b%Y")
# axis.Date(1,at=axis.dt,labels=format(axis.dt,"%b-%Y"))
for (i in 1:6){
	lines(c(snp_dat$dt[1],rev(snp_15$dt)),
			c(snp_dat$x[1],x.predj[i,]),col="grey65")
}
#lines(c(snp_dat$dt[1],rev(snp_15$dt)),
#			c(snp_dat$x[1],x.predj[2,]),col="slategrey",lwd=2)
lines(c(snp_dat$dt[1],rev(snp_15$dt)),lty=4,
			c(snp_dat$x[1],x.predj.mean),col="blue",lwd=2)
lines(c(snp_dat$dt[1],rev(snp_15$dt)),lty=2,
			c(snp_dat$x[1],x.predj.ci[1,]),col="blue",lwd=2)
lines(c(snp_dat$dt[1],rev(snp_15$dt)),lty=2,
			c(snp_dat$x[1],x.predj.ci[2,]),col="blue",lwd=2)
lines(c(snp_15$dt,snp_dat$dt[1]),c(snp_15$x,snp_dat$x[1]),
		col="red",lwd=2)
dev.off()
preset()
