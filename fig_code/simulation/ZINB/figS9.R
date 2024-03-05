source("ZINB_sim.R")
library(ggplot2)
library(ggsci)


M=function(k,para.list){
  Z1=para.list[[1]][,k]
  d=para.list[[2]]
  C1=para.list[[3]][,k]
  W=para.list[[4]]
  X=para.list[[5]]
  #logistic
  logis=data.frame(Z1,W)
  glm1=summary(glm(Z1~.-1,family=stats::quasibinomial(link = "logit"),control=list(maxit=100),data=logis))
  alpha=glm1[["coefficients"]][,1]
  alpha_p=glm1[["coefficients"]][,4]
  #weight NB
  NB=data.frame(C1,X)
  suppressWarnings(glmnb <- try(glm.nb(C1~.-1+offset(d), data = NB,link = log,weights=1-Z1), silent = TRUE))
  if (inherits(glmnb, "try-error")) {
    glmmnb <- glmmTMB(reformulate(termlabels = colnames(NB)[-c(1,2)], response='C1'),offset=d,se = TRUE,data =  NB,family =nbinom2(link = "log"),ziformula = ~0,weights=1-Z1)
    beta=summary(glmmnb)[["coefficients"]][["cond"]][,1]
    beta_sd=summary(glmmnb)[["coefficients"]][["cond"]][,2]
    beta_p=summary(glmmnb)[["coefficients"]][["cond"]][,4]
    phi=summary(glmmnb)[["sigma"]]
  }else if(any(is.na(glmnb[["coefficients"]]))) {
    glmmnb <- glmmTMB(reformulate(termlabels = colnames(NB)[-c(1,2)], response='C1'),offset=d,se = TRUE,data =  NB,family =nbinom2(link = "log"),ziformula = ~0,weights=1-Z1)
    beta=summary(glmmnb)[["coefficients"]][["cond"]][,1]
    beta_sd=summary(glmmnb)[["coefficients"]][["cond"]][,2]
    beta_p=summary(glmmnb)[["coefficients"]][["cond"]][,4]
    phi=summary(glmmnb)[["sigma"]]
  }else{
    beta=summary(glmnb)[["coefficients"]][,1]
    beta_sd=summary(glmnb)[["coefficients"]][,2]
    beta_p=summary(glmnb)[["coefficients"]][,4]
    phi=summary(glmnb)[["theta"]]
  }
  M_k=list(alpha,alpha_p,beta,beta_sd,beta_p,phi)
  return(M_k)
}

mbDecoda_pho=function(count,x,id,signif=0.05,maxit = 100,reltol = 1e-5,adjust="BH"){
  count=as.matrix(count)
  n=length(count[,1])
  K=length(count[1,])
  group=x
  tax_struc=c(which(colSums(count[group==unique(group)[1],])==0),which(colSums(count[group==unique(group)[2],])==0))
  tax_keep=which(colSums(count)==0)
  W=matrix(1,n,1)
  X=cbind(rep(1,n),group)
  p1=length(W[1,])
  p2=length(X[1,])
  #######initial parameters########
  d=rep(0,n)
  phi=rep(10,K)
  alpha=matrix(0,p1,K)
  beta=matrix(0,p2,K)
  Z=matrix(0,n,K)
  diff=1
  diff1=1
  Q=-Inf
  iter=0
  #######EM#########
  while(diff>reltol && iter <= maxit){
    ########E step#########
    iter=iter+1
    eta=1/(1+exp(-W%*%alpha))
    Y=exp(d+X%*%beta)
    PHI=kronecker(t(phi),rep(1,n))
    Z=eta/(eta+(1-eta)*PHI/((Y+PHI)^PHI))
    Z[count!=0]=0
    t=1:K
    para.list=list(Z,d,count,W,X)
    ########M step 1#########
    M_out=sapply(t,M,para.list=para.list)
    row.names(M_out)=c("alpha","alpha_p","beta","beta_sd","beta_p","phi")
    alpha=matrix(unlist(M_out[1,]),nrow=p1)
    beta=matrix(unlist(M_out[3,]),nrow=p2)
    phi=unlist(M_out[6,])
    Y=exp(d+X%*%beta)
    PHI=kronecker(t(phi),rep(1,n))
    ########M step 2#########,count,X,Z,beta,phi
    d_fn=function(d){
      l2_neg=-sum((1-Z)*dnbinom(x=count, size=PHI, mu=Y, log = T))
      return(l2_neg)
    }
    d_grad=function(d){
      d_grad=rowSums((1-Z)*(count-Y)*PHI/(Y+PHI))
      return(d_grad)
    }
    opti=optim(par=d, fn=d_fn, gr =d_grad,method = "BFGS", hessian = FALSE)
    d=opti[["par"]]
    l2=-d_fn(d)
    eta=1/(1+exp(-W%*%alpha))
    l1=sum(Z*log(eta)+(1-Z)*log(1-eta))
    diff=l1+l2-Q
    Q=l1+l2
  }
  beta_sd=matrix(unlist(M_out[4,]),nrow=p2)
  delta0=beta[2,]
  sd0=beta_sd[2,]
  p=sort(delta0)
  interval.out=matrix(NA,10,5)
  colnames(interval.out)=c("interval.prop","bias_inter","bias_sd_inter","power","fdr")
  bias_true=mean(delta0[-id])
  bias_sd_true=sd(delta0[-id])
  interval.out[,1]=seq(0.1,1,by=0.1)
  for (j in 1:10) {
    interval.prop=interval.out[j,1]
    pp=function(p,i){p[i+floor(K*interval.prop)]-p[i]}
    if(j==10){
      bias=mean(p)
      sd_bias=sd(p)
    }else{
      a=which.min(sapply(1:(K-floor(K*interval.prop)), pp,p=p))[1]
      bias=mean(p[a:(a+floor(K*interval.prop))])
      sd_bias=sd(p[a:(a+floor(K*interval.prop))])
    }
    delta=delta0-bias
    sd=sqrt(sd0^2+(sd_bias)^2)
    p.val = sapply(delta/sd, function(x) 2*pnorm(abs(x), mean = 0, sd = 1, lower.tail = F))
    p.val[tax_struc]=0  #these taxon are present in only one group, which lead to weighted nb variable dissatisfaction rank.
    p.val[tax_keep]=1 #this taxon are considered as inexistent
    q.val=p.adjust(p.val, method =adjust)
    id_ZINBA=which(q.val<0.05)
    interval.out[j,2]=bias
    interval.out[j,3]=sd_bias
    interval.out[j,4]=length(intersect(id_ZINBA,id))/length(id)
    interval.out[j,5]=1-length(intersect(id_ZINBA,id))/length(id_ZINBA)
  }
  
  # Prepare outputs
  out=list(interval.out,bias_true,bias_sd_true)
  names(out)=c("interval.out","bias.true","sd.true")
  return(out)
}



library(doSNOW)
K=100
n1=25
n2=25
n.sim=100
sig.prob = c(0.05,0.1,0.2,0.5)
bias = "small"
sim.seed=matrix(1:(n.sim*4),n.sim,4)
confounder = F
ncores=80


DATA=list()
for (i in 1:4) {
  sig.prob.i=sig.prob[i]
  for (j in sim.seed[,i]) {
    set.seed(j)
    DATA[[j]]=ZINB_sim(seed=j,K=K,n1=n1,n2=n2,p=sig.prob.i,bias =bias,zi=0.3,confounder =confounder)
  }
}


cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
simlist=foreach(i = DATA) %dopar% {
  library(glmmTMB)
  adjust="BH"
  data=i
  id=which(data[["diff.taxa"]]==1)
  group=data[["grp"]]
  count=data[["count"]]
  suppressWarnings(out <- try(mbDecoda_pho(count=count,x=group,id=id)
                              , silent = TRUE))
  out
}
stopCluster(cl)

f1 <- function(list) {
  bias.e=list[["interval.out"]][,2]
  bias.t=list[["bias.true"]]
  bias.e-bias.t
}

f2 <- function(list) {
  sd.e=list[["interval.out"]][,3]
  sd.t=list[["sd.true"]]
  sd.e-sd.t
}


bias.bias=foreach(i=simlist, .combine=cbind) %do% f1(i)
bias.sd=foreach(i=simlist, .combine=cbind) %do% f2(i)
power=foreach(i=simlist, .combine=cbind) %do% i[["interval.out"]][,4]
fdr=foreach(i=simlist, .combine=cbind) %do% i[["interval.out"]][,5]
fdr[is.na(fdr)]=0
f3<- function(i) {
  ind=matrix(1:400,100,4)
  res=rbind(cbind(seq(0.1,1,by=0.1),rowMeans(abs(bias.bias[,ind[,i]]))),
            cbind(seq(0.1,1,by=0.1),rowMeans(abs(bias.sd[,ind[,i]]))),
            cbind(seq(0.1,1,by=0.1),rowMeans(power[,ind[,i]])),
            cbind(seq(0.1,1,by=0.1),rowMeans(fdr[,ind[,i]])))
  res=as.data.frame(res)
  class=rep(c("estimation bias of c","estimation bias of sd(c)","power","FDR"),each=10)
  res=data.frame(res,class)
  res
}

out=as.data.frame(foreach(i=1:4, .combine=rbind) %do% f3(i))
p=rep(c(0.05,0.1,0.2,0.5),each=40)
out=data.frame(out,p)
colnames(out)=c("interval.prop","value","class","p")
out$class=factor(out$class, levels =c("estimation bias of c","estimation bias of sd(c)","power","FDR"))
out$p=factor(out$p, levels =c("0.05","0.1","0.2","0.5"))
out$interval.prop=factor(out$interval.prop, levels =as.character(seq(0.1,1,by=0.1)))


out1=out[out$class%in%c("power","FDR"),]
out1=out1[out1$interval.prop%in%as.character(seq(0.1,0.9,by=0.1)),]
p1=ggplot(out1, aes(x=interval.prop, y=value, group=p))+ theme_bw()
p1 +  geom_point(aes(col=p),pch=7,size=1) +  
  geom_line(aes(col=p)) +
  facet_wrap(~class,scales="free")+ 
  scale_color_npg()+
  labs( x = expression(rho),y = NULL, color = expression(pi))#theme(legend.position = "top")



