library(glmmTMB)

#######ZINBLDA mixture of logistic and NB EM##########
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

mbDecoda=function(count,x,Gamma=NULL,W=NULL,interval.prob=0.5,signif=0.05,d_init=F, prev.cut = 0.05, maxit = 100,reltol = 1e-5,adjust="BH"){
  count=as.matrix(count)
  n=length(count[,1])
  K=length(count[1,])
  group=x
  if(length(unique(group))==2){
    tax_struc=c(which(colSums(count[group==unique(group)[1],])==0),which(colSums(count[group==unique(group)[2],])==0))
  }else{
    tax_struc=integer(0)
  }
  tax_keep=which(colMeans(count!=0)<=prev.cut) #filtering rare taxa
  if(is.null(W)){
    W=matrix(1,n,1)
  }else{
    W=cbind(rep(1,n),W)
  }
  if(is.null(Gamma)){
    X=cbind(rep(1,n),group)
  }else{
    X=cbind(rep(1,n),group,Gamma)
  }
  p1=length(W[1,])
  p2=length(X[1,])
  
  #######initial parameters########
  if(d_init==T){
    d=log(rowSums(count))
    d[is.na(d)]=1
  }else{d=rep(0,n)}
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
      suppressWarnings(l2_neg <- try(-sum((1-Z)*dnbinom(x=count, size=PHI, mu=Y, log = T)), silent = TRUE))
      if (inherits(l2_neg, "try-error")) {
        l2_neg <- -sum((1-Z)*dnbinom(x=count, size=PHI, prob = PHI/(PHI+Y), log = T))
        return(l2_neg)
      }else{
        l2_neg <- -sum((1-Z)*dnbinom(x=count, size=PHI, mu=Y, log = T))
        return(l2_neg)
      }
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
  
  #estimation
  alpha_p=matrix(unlist(M_out[2,]),nrow=p1)
  beta_p=matrix(unlist(M_out[5,]),nrow=p2)
  beta_sd=matrix(unlist(M_out[4,]),nrow=p2)
  phi=unlist(M_out[6,])
  delta0=beta[2,]
  sd0=beta_sd[2,]
  
  #bias correction
  p=sort(delta0)
  pp=function(p,i){p[i+floor(K*interval.prob)]-p[i]}
  a=which.min(sapply(1:(K-floor(K*interval.prob)), pp,p=p))[1]
  bias=mean(p[a:(a+floor(K*interval.prob))])
  delta=delta0-bias
  sd_bias=sd(p[a:(a+floor(K*interval.prob))])
  sd=sqrt(sd0^2+(sd_bias)^2)
  
  #inference
  p.val = sapply(delta/sd, function(x) 2*pnorm(abs(x), mean = 0, sd = 1, lower.tail = F))
  p.val[is.na(p.val)]=1 
  p.val[tax_struc]=0  #these taxa are considered as structure zero, who present in only one group.
  p.val[tax_keep]=NA #these taxa are considered as rare taxa
  q.val=p.adjust(p.val, method =adjust)
  reject=c(q.val<signif)
  if(is.null(Gamma)){
    beta1=beta[1,]
    beta_sd1=beta_sd[1,]
    beta_p1=beta_p[1,]
  }else{
    beta1=beta[-2,]
    beta_sd1=beta_sd[-2,]
    beta_p1=beta_p[-2,]
  }
  
  # Prepare outputs
  ZINB=list(alpha,alpha_p,phi)
  names(ZINB)=c("alpha","alpha_p","phi")
  
  Confounders=list(beta1,beta_sd1,beta_p1)
  names(Confounders)=c("beta","beta_sd","beta_p")
  
  test=data.frame(delta,sd,p.val,q.val,reject)
  rownames(test)=colnames(count)
  Bias_correct=list(delta0,bias,sd_bias,test)
  names(Bias_correct)=c("delta_bias","bias","sd_bias","DAA")
  
  out=list(ZINB,Confounders,Bias_correct)
  names(out)=c("ZINB","Confounders","Bias_correct")
  return(out)
}
