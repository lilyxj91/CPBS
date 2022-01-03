#' Title
#'
#' @param theta
#' @param A
#'
#' @return
#' @export
#'
theta.function = function(theta,A){
  result = dweibull(A,theta[1],theta[2],log=TRUE)
  return(-sum(result))
}



#' Title naive method of ignoring the unique data structure.
#'
#' @param IND
#' @param X
#' @param A
#' @param Y_SAB
#'
#' @return
#' @export
#'
naive.method = function(T,IND,X,A,Y_SAB){
  mylogit = glm(Y_SAB[!is.na(Y_SAB)] ~ X[!is.na(Y_SAB),], family = "binomial")
  alpha = coef(mylogit)
  alpha.se = summary(mylogit)$coefficients[,2]

  theta.output = optim(c(2,mean(A)),theta.function,A=A,hessian=TRUE)
  theta = theta.output$par
  theta.se = sqrt(diag(solve(theta.output$hessian)))

  tempcox = coxph( Surv(A[Y_SAB==1],T[Y_SAB==1], IND[Y_SAB==1]) ~  X[Y_SAB==1,])
  beta = summary(tempcox)$coefficients[,1]
  beta.se = summary(tempcox)$coefficients[,3]

  return(list(alpha=alpha,alpha.se=alpha.se, beta = beta, beta.se=beta.se, theta=theta, theta.se=theta.se))
}

