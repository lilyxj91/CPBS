

#' Title inverse logit
#'
#' @param x
#'
#' @return
#' @export
#'
expit = function(x){ exp(x)/(1+exp(x)) }



#' Title weight.fun: calculate weight w_{ij}
#'
#' @param theta theta
#' @param beta beta
#' @param labda lambda
#' @param time time
#' @param yi Yi
#' @param xi Xi
#'
#' @return wi
#' @export
#'
#'
weight.fun<-function(theta,beta,labda,time,yi,xi)
{
  wi<-rep(0,length(time))
  pi<-rep(0,length(time))
  Hi<-pweibull (time, theta[1],theta[2])/pweibull (20, theta[1],theta[2])

  temp2<-cumsum(labda)
  pi<-labda*exp(xi[1]*beta[1]+xi[2]*beta[2])*exp(-temp2*exp( xi[1]*beta[1]+xi[2]*beta[2]))

  ui<-sum(pi*Hi)
  wi<- 1/ui*(1-Hi)*pi
  return(wi)
}

#' Title Ey
#'
#' @param beta
#' @param labda
#' @param alpha
#' @param time
#' @param yi
#' @param xi
#'
#' @return
#' @export
#'
Ey = function(beta,labda,alpha,time,yi,xi){
  zz1=c(1,xi)
  Pi_M1= expit(zz1%*%alpha)

  temp1<-sum(labda*(time<=yi))
  Si = exp(-temp1*exp( xi[1]*beta[1]+xi[2]*beta[2]))

  num = Pi_M1*Si
  denom = 1-Pi_M1+num
  return(num/denom)
}

#' Title loglikelihood of alpha
#'
#' @param alpha
#' @param data
#' @param weight
#' @param omega
#'
#' @return
#' @export
#'
loglik.alpha = function(alpha,data,weight,omega){

  mu = cbind(rep(1,length(data[,1])),data[,3:4])%*% alpha
  pi= expit(mu)
  temp1=matrix(pi,nr=dim(weight)[1],nc=dim(weight)[2])

  temp2 = sum(omega*weight*log(temp1))
  temp3 = sum((1-omega[,1])*log(1-pi))
  return(-sum(temp2+temp3))
}

#' Title loglikelihood of theta
#'
#' @param theta
#' @param betaold
#' @param labdaold
#' @param thetaold
#' @param time
#' @param data
#' @param omega
#'
#' @return
#' @export
#'
loglik.theta <-function(theta, betaold,labdaold,thetaold,time,data,omega)
{
  temp2<-cumsum(labdaold)
  hi = dweibull(time,thetaold[1],thetaold[2])/pweibull (20, thetaold[1],thetaold[2])
  Hi = pweibull(time,thetaold[1],thetaold[2])/pweibull (20, thetaold[1],thetaold[2])
  s1=0

  for (i in 1:length(data[,1])) {
    xi = data[i,3:4]
    pi<-labdaold*exp(xi[1]*betaold[1]+xi[2]*betaold[2])*exp(-temp2*exp(xi[1]*betaold[1]+xi[2]*betaold[2]))
    ui<-sum(pi*Hi)
    cdf = cumsum(pi)

    sfun = stepfun(time,c(cdf,1))

    integrand = function(x,thetaold,theta){
      sfun(x)*dweibull(x,thetaold[1],thetaold[2])/pweibull (20, thetaold[1],thetaold[2])*
        (dweibull(x,theta[1],theta[2],log=TRUE)-pweibull (20, theta[1],theta[2],log=TRUE))
    }
    a=0;b=20
    out=vector("list")
    out$nodes=c(-0.9739065, -0.8650634, -0.6794096, -0.4333954, -0.1488743,  0.1488743,  0.4333954,  0.6794096,  0.8650634,  0.9739065)
    out$weights=c(0.06667134, 0.14945135, 0.21908636, 0.26926672, 0.29552422, 0.29552422, 0.26926672, 0.21908636, 0.14945135, 0.06667134)

    try1 = 0.5*(b-a)*sum(out$weights*integrand(0.5*((b-a)*out$nodes+a+b),thetaold=thetaold,theta=theta))

    s1 = s1+omega[i,1]*(dweibull(data[i,5],theta[1],theta[2],log=TRUE)-
                          pweibull (20, theta[1],theta[2],log=TRUE)  + try1/ui)
  }
  return(-sum(s1))
}


#' Title proposed EM algorithm
#'
#' @param tol A small number, default 1e-4
#' @param maxit max iteration, default 1000
#' @param data dataset
#'
#' @return a list
#' @import survival
#' @export
#'
cpbs_em <- function(tol=1e-4, maxit=1000,data){
  ## jump at both Y=1 and Y unobserved
  index=rep(0,dim(data)[1])
  for ( i in 1:(dim(data)[1])){
    if(is.na(data[i,6])){index[i]=1} else if(data[i,6]==1)
    {index[i]=1} else{
      index[i]=0
    }
  }

  junk=data[index==1,]

  tt<-sort(unique(junk[,1]));
  len<-length(data[,1])
  lensp<-length(tt)

  alphaold = c(1,1,1)
  thetaold = c(1,2)

  #data.1=data[complete.cases(data),]
  #tempcox1<-coxph( Surv(data.1[,1], data.1[,2]) ~ data.1[,3]+data.1[,4])
  #betaes3 <-tempcox1$coef
  #betaold<-betaes3
  betaold<-  c(-0.5,1)

  ttd<-rep(tt[1], lensp)
  for(iii in 2:lensp)
  {
    ttd[iii]<-tt[iii]-tt[iii-1]
  }

  labdaold<-ttd/mean(data[,1])

  weight1 =  omega = matrix(0,len,lensp)
  addon.num = addon.denom =  matrix(0,len,lensp)
  for ( i in 1:len){
    addon.num[i,]= (data[i,1]==tt)*data[i,2]
    addon.denom[i,]= (data[i,1]==tt)
  }
  count<-1
  #cat("iter = ", count, "beta = ", betaold,"alpha=",alphaold,"theta=",thetaold,"\n")

  while (count <= maxit)  {
    count <- count + 1

    for(ii in 1:len)
    {temp<-weight.fun(thetaold,betaold,labdaold, tt,data[ii,1],data[ii,3:4])
    omega[ii,]=ifelse(is.na(data[ii,6]),Ey(betaold,labdaold,alphaold,tt,data[ii,1],data[ii,3:4]),data[ii,6])
    weight1[ii,]<-temp
    }
    sumweight<-apply(omega*(weight1+addon.num),2,sum)

    alphanew = optim(alphaold,loglik.alpha,data=data,weight=weight1+addon.denom,omega=omega)$par
    #alphanew = alphaold

    ### only for Y_i==1 and Y_i unobserved  group
    data.1=data[index==1,]
    cov1<-c(data.1[,1],rep(tt,len))
    cov2<-c(data.1[,2],rep(1,len*lensp))
    cov3<-c(data.1[,3],rep(data[,3], each=lensp))
    cov4<-c(data.1[,4],rep(data[,4], each=lensp))

    cov5<-c(omega[index==1,1],as.vector(t(omega*weight1)))
    cov6<-ifelse(cov5>0,cov5,min(cov5[cov5>0]))
    tempcox<-coxph( Surv(cov1, cov2) ~ cov3+cov4, weights=cov6)
    betanew <- tempcox$coef

    weightn<-(omega*(weight1+addon.denom))*exp(data[,3]*betanew[1]+data[,4]*betanew[2])
    sumweightn<-apply(weightn,2,sum)

    csumweightn<-rep(sum(sumweightn), lensp)
    #csumweightn1<-rep(sum(sumweightn), lensp)
    for(kk in 2:lensp)
    { csumweightn[kk]<-csumweightn[kk-1]-sumweightn[kk-1]}
    #{csumweightn1[kk]<-sum(sumweightn)-cumsum(sumweightn)[kk-1]}

    labdanew<-sumweight/csumweightn

    thetanew = optim(thetaold,loglik.theta,betaold=betanew,labdaold=labdanew,thetaold=thetaold,
                     time=tt,data=data,omega=omega)$par

    #cat("iter = ", count, "beta = ", round(betanew,3),"alpha=",round(alphanew,3),"theta=",round(thetanew,3),"\n")

    if (max(c(abs(betanew-betaold),abs(alphanew-alphaold),abs(thetanew-thetaold))) < tol ){  # calculate the l1 norm
      cat("\nSuccessfully Converged\n")
      return(list(beta = betanew,labda=labdanew,alpha=alphanew,theta=thetanew,conv=1,
                  weight=weight1,addon.denom=addon.denom,addon.num=addon.num,omega=omega,tt=tt))
    } else {
      alphaold <- alphanew
      betaold <- betanew
      labdaold <- labdanew
      thetaold = thetanew
    }
  }
  print("Convergence Failed")
  return(list(beta = betanew,labda=labdanew,alpha=alphanew,theta=thetanew,conv=0,
              weight=weight1,addon.denom=addon.denom,addon.num=addon.num,omega=omega,tt=tt))
}
