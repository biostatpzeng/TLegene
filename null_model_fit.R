NullRegressionLargeData = function(Data,n,p,type="Continuous",null.cov,nullvar=NULL){

  if(is.null(null.cov)) MainCov = rep(1,n)
  else MainCov = cbind(rep(1,n),as.matrix(Data[,null.cov]))
  ## Above adjustment is to accommodate extra unpenalized covariates. No scaling is used for unpenalized covariates.

  if(type=="Continuous"){
    if(is.null(null.cov)) nullfit = lm(Data[,1]~1)
    else nullfit = lm(Data[,1]~1+as.matrix(Data[,null.cov]))
  }
  else{
    if(is.null(null.cov)) nullfit = glm(Data[,1]~1, family="binomial")
    else nullfit = glm(Data[,1]~1+as.matrix(Data[,null.cov]), family="binomial")
  }
  halpha = coef(nullfit)
  hresid = residuals(nullfit,type="response")
  hfit = fitted(nullfit)
  if(is.null(nullvar)){
    if(type=="Continuous"){
      hsigma2 = sum(residuals(nullfit)^2)/(n-length(null.cov)-1)
      hfit2 = rep(hsigma2,length(hfit))
    }
    else hfit2 = (hfit*(1-hfit))
  }
  else{
    if(type=="Continuous"){
      hsigma2 = nullvar
      hfit2 = rep(nullvar,length(hfit))
    }
    else hfit2 = nullvar
  }

  MainCov.scale = MainCov
  if(!is.null(null.cov)){
    MainCov.null.scale = scale(MainCov[,2:(length(null.cov)+2-1)])
    MainCov.scale[,2:(length(null.cov)+2-1)] = MainCov.null.scale
  }

  hfit2Sqrt = sqrt(hfit2)
  if(is.matrix(MainCov.scale)){
    MainCovStd.scale = matrix(rep(hfit2Sqrt,times=ncol(MainCov.scale)),nrow=nrow(MainCov.scale),ncol=ncol(MainCov.scale)) * MainCov.scale
  }
  else{
    MainCovStd.scale = hfit2Sqrt * MainCov.scale
  }
  MainCovSq.scale = t(MainCovStd.scale) %*% MainCovStd.scale
  Omega = MainCovSq.scale

  gc()

  if(type=="Continuous") return(list(MainCov=MainCov.scale,Omega=Omega,resid=hresid,fit2=unique(hfit2),alpha=halpha))
  else return(list(MainCov=MainCov.scale,Omega=Omega,resid=hresid,fit2=hfit2,alpha=halpha))
}
