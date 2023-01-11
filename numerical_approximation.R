WeightSumChisq = function(weight,df=rep(1,length(weight))){
  # Calculate the df and the noncentral parameter of an approximated noncentral chisq distribution for the weighted sum of independent chisq distributions.
  # Option weight: a vector of weights for each chisq.
  # Option df: a vector of df for each chisq.
  # Return the df(l) and the noncentrality (\delta) of the approximation noncentral chisq.
  # Reference: Liu (2009). Match the skewness.
  
  # Sum of weights to different powers and terms related to criterions
  S1 = sum(weight)    # equal to c1[1]
  S2 = sum(weight^2)  # equal to c1[2]
  S3 = sum(weight^3)
  S4 = sum(weight^4)  
  tmp_1 = S4/(S2^2)
  tmp_2 = (S3^2)/(S2^3)
  
  # Calculate a=sqrt(l+2\delta), \delta, and l
  if(tmp_1<tmp_2){
    a = 1/(S3/(S2^(3/2))-sqrt(tmp_2-tmp_1))
    delta = a^3*S3/(S2^(3/2))-a^2
    l = a^2-2*delta
  }
  else{
    a = (S2^(3/2))/S3    
    delta = 0
    l = a^2
  }
  #print(c(tmp_1,tmp_2,a,delta,l))
  return(c(mean.Q=S1,SD.Q=sqrt(2*S2),df=l,noncentral=delta,mean.X=l+delta,SD.X=sqrt(2*(l+2*delta))))
}

WeightSumChisq_mod = function(weight,df=rep(1,length(weight))){
  # Calculate the df and the noncentral parameter of an approximated noncentral chisq distribution for the weighted sum of independent chisq distributions.
  # Option weight: a vector of weights for each chisq.
  # Option df: a vector of df for each chisq.
  # Return the df(l) and the noncentrality (\delta) of the approximation noncentral chisq.
  # Reference: Liu (2009). Modified version with matching the kurtosis.
  
  # Sum of weights to different powers and terms related to criterions
  S1 = sum(weight)    
  S2 = sum(weight^2)
  S3 = sum(weight^3)
  S4 = sum(weight^4)  
  tmp_1 = S4/(S2^2)
  tmp_2 = (S3^2)/(S2^3)
  
  # Calculate a=sqrt(l+2\delta), \delta, and l
  if(tmp_1<tmp_2){
    a = 1/(S3/(S2^(3/2))-sqrt(tmp_2-tmp_1))
    delta = a^3*S3/(S2^(3/2))-a^2
    l = a^2-2*delta
  }
  else{
    a = sqrt(1/tmp_1)
    delta = 0
    l = a^2
  }
  #print(c(tmp_1,tmp_2,a,delta,l))
  return(c(mean.Q=S1,SD.Q = sqrt(2*S2),df=l,noncentral=delta,mean.X=l+delta,SD.X=sqrt(2*(l+2*delta))))
}

quan_rho_calculator = function(grid=seq(0,1,by=0.05),df_f,lambda_r,min_pvalue,acc=5e-10,lim=1e4,max_core=4){
  grid_inner = grid[2:(length(grid)-1)]
  numCores = parallel::detectCores()
  
  quan_rho_grid_inner = mclapply(grid_inner,FUN=function(grid_point){
    diff.prob.davies = function(q,tailprob,lambda,df=rep(1,length(lambda)),delta=rep(0,length(lambda)),lim=1e4,acc=5e-10){
      tailprob-CompQuadForm::davies(q,lambda=lambda,lim=lim,acc=acc)$Qq
    }
    quan_tmp = uniroot(diff.prob.davies,lower=0,upper=1e7,tailprob=min_pvalue,lambda=c(grid_point*rep(1,df_f),(1-grid_point)*lambda_r),acc=acc,lim=lim)$root
  },
  mc.cores=min(numCores,max_core))
  quan_rho_grid_inner = as.numeric(unlist(quan_rho_grid_inner))
  names(quan_rho_grid_inner) = as.character(grid_inner)
  return(list(grid=grid,quan_rho_grid_inner=quan_rho_grid_inner))
}

pvalue_minimizePValue_davies = function(u_f,u_r,df_f,lambda_r,min_pvalue=NULL,quan_grid=NULL,acc=5e-10,lim=1e4,ChisqApp="3M",max_core=4){
  if(is.null(min_pvalue)){
    p_rho = function(rho,df_f,lambda_r,u_f,u_r,acc=5e-10){
      pvalue = CompQuadForm::davies(q=rho*u_f+(1-rho)*u_r,
                      lambda=c(rho*rep(1,df_f),(1-rho)*lambda_r),acc=acc,lim=1/acc)$Qq
      if(pvalue<=acc) pvalue = pvalue.liu(Q=rho*u_f+(1-rho)*u_r,
                                          weight=c(rho*rep(1,df_f),(1-rho)*lambda_r),
                                          method=ChisqApp)$pvalue
      return(pvalue)
    }
    min_pvalue_inside01 = optimize(p_rho,interval=c(0,1),maximum=FALSE,df_f=df_f,lambda_r=lambda_r,u_f=u_f,u_r=u_r,acc=acc)$objective
    min_pvalue = min(c(min_pvalue_inside01,
                       p_rho(0,df_f=df_f,lambda_r=lambda_r,u_f=u_f,u_r=u_r,acc=acc),
                       p_rho(1,df_f=df_f,lambda_r=lambda_r,u_f=u_f,u_r=u_r,acc=acc)))
  }
  

  integrand = function(a,min_pvalue,df_f,lambda_r,acc=5e-10,grid=NULL,quan_rho_grid_inner=NULL,max_core=4){
    quan_rho = function(rho,df_f,lambda_r,p_rho_min,u_f){
      if(rho==1) stop("rho is not allowed to be 1.")
      else{
        diff.prob.davies = function(q,tailprob,lambda,df=rep(1,length(lambda)),delta=rep(0,length(lambda)),lim=1e4,acc=5e-10){
          tailprob-CompQuadForm::davies(q,lambda=lambda,lim=lim,acc=acc)$Qq
        }
        quan_tmp = uniroot(diff.prob.davies,lower=0,upper=1e7,tailprob=p_rho_min,lambda=c(rho*rep(1,df_f),(1-rho)*lambda_r),acc=acc,lim=lim)$root
        names(quan_tmp) = c()
        quan_tmp_u_r = (quan_tmp - rho*u_f)/(1-rho)
        return(quan_tmp_u_r)
      }
    }
    
    optimize_interval = c(0,1)
    
    numCores = parallel::detectCores()
    res = mclapply(a,FUN=function(a.tmp){
      if(!is.null(grid) & !is.null(quan_rho_grid_inner)){
        grid_inner = as.numeric(names(quan_rho_grid_inner))
        ### refne interval for optimization by running the following over a grid (no 0 and 1)
        grid_inner_index = order((quan_rho_grid_inner - grid_inner*a.tmp)/(1-grid_inner))[1]
        optimize_interval = c(grid[grid_inner_index],grid[grid_inner_index+2])
      }
      min_quantile = optimize(quan_rho,interval=optimize_interval,df_f=df_f,lambda_r=lambda_r,p_rho_min=min_pvalue,u_f=a.tmp)$objective
      tail.prob = 0
      if(class(try(CompQuadForm::davies(q=min_quantile,lambda=lambda_r,acc=acc,lim=1/acc)))!="try-error") tail.prob = CompQuadForm::davies(q=min_quantile,lambda=lambda_r,acc=acc,lim=1/acc)$Qq
      if(tail.prob<=acc) tail.prob = pvalue.liu(Q=min_quantile,weight=lambda_r,method=ChisqApp)$pvalue
      (1-tail.prob) * dchisq(a.tmp,df=df_f)
    },
    mc.cores=min(numCores,max_core))
    res = as.numeric(unlist(res))
    return(res)
  }
  
  pvalue = 1 - integrate(integrand,lower=0,upper=qchisq(min_pvalue,df=df_f,lower.tail=FALSE),min_pvalue=min_pvalue,df_f=df_f,lambda_r=lambda_r,acc=acc,grid=quan_grid[[1]],quan_rho_grid_inner=quan_grid[[2]],max_core=max_core)$value
  if(pvalue<0) pvalue = abs(pvalue)
  return(pvalue)
}

pvalue_minimizePValue_liu = function(u_f,u_r,df_f,lambda_r,min_pvalue=NULL,acc=5e-10,lim=1e4,ChisqApp="4M"){
  if(is.null(min_pvalue)){
    p_rho = function(rho,df_f,lambda_r,u_f,u_r,acc=5e-10){
      pvalue = CompQuadForm::davies(q=rho*u_f+(1-rho)*u_r,
                      lambda=c(rho*rep(1,df_f),(1-rho)*lambda_r),acc=acc,lim=1/acc)$Qq
      if(pvalue<=acc) pvalue = pvalue.liu(Q=rho*u_f+(1-rho)*u_r,
                                          weight=c(rho*rep(1,df_f),(1-rho)*lambda_r),
                                          method=ChisqApp)$pvalue
    }
    min_pvalue_inside01 = optimize(p_rho,interval=c(0,1),maximum=FALSE,df_f=df_f,lambda_r=lambda_r,u_f=u_f,u_r=u_r,acc=acc)$objective
    min_pvalue = min(c(min_pvalue_inside01,pvalue.r,pvalue.f))
  }
  
  integrand = function(a,min_pvalue,df_f,lambda_r,acc=5e-10){
    quan_rho = function(rho,df_f,lambda_r,p_rho_min,u_f){
      if(rho==1) stop("rho is not allowed to be 1.")
      else{
        if(ChisqApp=="4M") parameters = WeightSumChisq_mod(weight=c(rho*rep(1,df_f),(1-rho)*lambda_r),df=rep(1,length(weight)))
        if(ChisqApp=="3M") parameters = WeightSumChisq(weight=c(rho*rep(1,df_f),(1-rho)*lambda_r),df=rep(1,length(weight)))
        # mean.Q, SD.Q, df, noncentral, mean.X, SD.X
        quan_tmp = (qchisq(p_rho_min,df=parameters["df"],ncp=parameters["noncentral"],lower.tail=FALSE) - parameters["mean.X"]) / parameters["SD.X"] * parameters["SD.Q"] + parameters["mean.Q"] # Note: use upper tail for quantile calcualtion to avoid Inf when p_rho_min is extremely small.
        names(quan_tmp) = c()
        quan_tmp_u_r = (quan_tmp - rho*u_f)/(1-rho)
        return(quan_tmp_u_r)
      }
    }
    res = sapply(a,function(a.tmp){
      # min_quantile = optimize(quan_rho,interval=c(0,1-1e-3),df_f=df_f,lambda_r=lambda_r,p_rho_min=min_pvalue,u_f=a.tmp)$objective
      min_quantile_tmp = optimize(quan_rho,interval=c(0,1),df_f=df_f,lambda_r=lambda_r,p_rho_min=min_pvalue,u_f=a.tmp)$objective # try 2017-04-06
      min_quantile = min(quan_rho(rho=0,df_f=df_f,lambda_r=lambda_r,p_rho_min=min_pvalue,u_f=a.tmp),
                         min_quantile_tmp) # try 2017-04-06
      tail.prob = 0
      if(class(try(CompQuadForm::davies(q=min_quantile,lambda=lambda_r,acc=acc,lim=1/acc)))!="try-error") tail.prob = CompQuadForm::davies(q=min_quantile,lambda=lambda_r,acc=acc,lim=1/acc)$Qq
      if(tail.prob<=acc) tail.prob = pvalue.liu(Q=min_quantile,weight=lambda_r,method=ChisqApp)$pvalue
      (1-tail.prob) * dchisq(a.tmp,df=df_f)
    })
    return(res)
  }
  pvalue = 1 - integrate(integrand,lower=0,upper=qchisq(min_pvalue,df=df_f,lower.tail=FALSE),min_pvalue=min_pvalue,df_f=df_f,lambda_r=lambda_r,acc=acc)$value
  return(pvalue)
}

pvalue.liu = function(Q,weight,method="3M",df=rep(1,length(weight))){
  if(method=="4M") para = WeightSumChisq_mod(weight=weight,df=df)
  else para = WeightSumChisq(weight=weight,df=df)
  mean.Q = para[1]
  SD.Q = para[2]
  mean.X = para[5]
  SD.X = para[6]
  df.app = para[3]
  noncentral.app = para[4]
  Q.new = (Q-mean.Q)/SD.Q * SD.X + mean.X
  pvalue = pchisq(Q.new,df=df.app,ncp=noncentral.app,lower.tail=FALSE)
  return(list(pvalue=unname(pvalue),df=unname(df.app),ncp=unname(noncentral.app)))
}

Weight_ScrnStat = function(data,d=0,p,method){
  
  n = nrow(data)
  y = data[,1]
  e = data[,2]
  if(d!=0) {
    G = data[,(3+d):(3+d+p-1)]
    X = data[,3:(3+d-1)]
  }
  else G = data[,3:(3+p-1)]
  
  if(length(setdiff(c(0,1),e))==0)
  {
    cor.func = function(g,e){
      if(d!=0) margin.fit = summary(glm(e~g+X,family="binomial"))$coef
      else margin.fit = summary(glm(e~g,family="binomial"))$coef
      if("g" %in% rownames(margin.fit)) res = margin.fit[2,c(1,4,3)]
      else res = c(0,1,0)
      return(res)
    }
  }
  else
  {
    cor.func = function(g,e){
      if(d!=0) margin.fit = summary(lm(e~g+X))$coef
      else margin.fit = summary(lm(e~g))$coef
      if("g" %in% rownames(margin.fit)) res = margin.fit[2,c(1,4,3)]
      else res = c(0,1,0)
      return(res)
    }
  }
  cor.wt = apply(G,2,cor.func,e)
  
  marg.func = function(g,d){
    if(d!=0) summary(glm(y~g+X,family="binomial"))$coef[2,c(1,4,3)]
    else summary(glm(y~g,family="binomial"))$coef[2,c(1,4,3)]
  }
  marg.wt = apply(G,2,marg.func,d)
  
  if(method=="ScrnStat1"){
    wt = cor.wt
    wt[,cor.wt[2,]>=marg.wt[2,]] = marg.wt[,cor.wt[2,]>=marg.wt[2,]]
    wt.return = wt[3,]
  }
  else if(method=="ScrnStat2"){
    wt.sign = sign(cor.wt[3,])
    wt.return = wt.sign* (cor.wt[3,]^2+marg.wt[3,]^2)
  }
  else stop("method needs to be either ScrnStat1 or ScrnStat2")
  
  return(wt.return)
}
