# Likang Xu (s2295871)

# In applied statistics, smoothing plays an important role in model regression. 
# The project is about smoothing with basis expansions and penalized regression 
# using R functions given x, y data. In detail, we have the model: 
# y=X*beta+epsilon, where X_{ij} = b_j(x_i). b_j(x) are basis functions and in 
# this study, evenly spaced B-spline basis functions are used. To avoid 
# overfitting, a smoothing penalty is imposed to decrease the differences of
# beta_js corresponding to b_js from one to the next, hence the model is estimated
# by penalized least squares:
#      beta_hat =argmin_beta ||y-X*beta||^2+lambda*beta^T*D^T*D*beta.
# The value of lambda is chosen when the generalized cross validation criterion
# (GCV) is minimized. Smooth functions estimated by the above equation are called
# P-splines.

# This project has five functions to fit P-splines to x, y data. The function 
# fit.pspline (x,y,k=20,sp,bord=3,pord=2) would fit a spline smoother given six 
# inputs (explained below), return a list contains beta_hat (coef), fitted values
# (fitted), residual variance (sig2), and some other values. The function 
# pspline(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100) takes 7 arguments and 
# ould return the best fit spline smoother selected from a sequence of lambdas by
# GCV. The function print.pspline(m) could report some details of the model fit 
# and silently return a list containing some elements. The forth function, 
# predict.pspline(m,x,se=TRUE), would make predictions given the model fit and 
# new x data, as well as calculating the corresponding standard errors. Finally,
# plot.pspline(m) could generate 3 plots, a plot of the original x, y data with
# the estimated smooth function and approximate 95% credible intervals, a plot 
# of the model residuals against fitted values, and a qqplot of the residuals, 
# then silently return a list. 



# fit.pspline arguments:
# x, y – the vectors of data to smooth, 
# k – the number of basis functions
# sp – the smoothing parameter (lambda)
# bord – the B-spline order to use
# pord – the order of difference to use in the penalty

# Given the equation above, using QR decomposition X=QR, and eigen decomposition
# U*Lambda*U^T=R^{-T}*D^T*D*R^{-1}: 
# beta_hat=R^{-1}*U*(I+lambda*Lambda)^{-1}*U^T*Q^T*y, 
# fitted=X*beta_hat, 
# kappa=tr{( I+lambda*Lambda)^{-1}}, 
# sig2=||y-fitted||^2/(n-kappa), 
# gcv=sig2/(n-kappa). 
# This function firstly constructs the B-spline basis matrix X using x, k, bord,
# and the D matrix given k and pord. Next, defining Q, R, U and Lamda. 
# A=(I+lambda*Lambda)^{-1} is defined for computational efficiency. Then, beta_hat,
# fitted, sig2, gcv, r2, the covariance matrix (V), and the standard errors of 
# fitted values (stde) is calculated. Finally, this function returns a list 
# containing beta_hat (coef), fitted values (fitted), residual variance (sig2), 
# the generalized cross validation criterion (gcv), the effective degrees of 
# freedom (edf), bord, pord, number of coefficients (k), R square (r2), the 
# covariance matrix (V), the standard errors of fitted values (se), x and y. 

fit.pspline <- function(x,y,k=20,sp,bord=3,pord=2){
  # construct the X
  dk <- diff(range(x))/(k-bord)
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  # form D
  D <- diff(diag(k),differences=pord)
  # QR decomposition of X
  qrX = qr(X)
  Q <- qr.Q(qrX); R <- qr.R(qrX)
  R_inv <- solve(R) # the inverse of R
  # eigen decomposition, define U and Lambda
  ed <- eigen(t(R_inv) %*% (t(D) %*% (D %*% R_inv))) 
  U <- ed$vectors; Lamda <- diag(ed$values)
  
  # calculation of all the needed values:
  # beta_hat and fitted values
  A <- solve(diag(k) + sp*Lamda)
  beta_hat <- solve(R, U %*% (A %*% (t(U) %*% (t(Q) %*% y))))
  fit_value <- X %*% beta_hat 
  
  # sig2 and gcv
  kappa <- sum(diag(A))
  n <- nrow(X)
  sig2 <- sum((y-fit_value)^2)/(n-kappa)
  gcv <- sig2/(n-kappa)
  
  # R square, covariance matrix, and standard errors
  r2 <- 1-(n-1)*sig2/sum((y-mean(y))^2)
  V <- solve(R, U %*% (A %*% (t(U) %*% t(R_inv)))) * sig2
  stde <- rowSums(X * (X %*% V))^.5
  
  # return the list
  list(coef = beta_hat, fitted = fit_value, sig2 = sig2, gcv = gcv, edf = kappa,
       bord = bord, pord = pord, coefficients = k, r2 = r2, V = V, se = stde, 
       x = x, y = y)
}



# pspline arguments:
# logsp – the interval of the smoothing parameter in log scale
# ngrid – the number of lambdas to try
# other arguments are the same as fit.pspline
# This function aims to find the best model fit given a sequence of smoothing 
# parameters. Given logsp and ngrid, it firstly forms a sequence sp that contains
# all the lambdas (if logsp is a number, sp would be exp(logsp), the function 
# would return the model fit given this sp). Then, using a loop to fit the pspline
# models for each lambda and store the gcv. Finally select the largest lambda with
# the minmum gcv, and return the model fit (the same list as fit.pspline)

pspline <- function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100){
  # judge whether logsp is a number
  if (length(logsp)!=1){
    # form the sequence of log(lambda)
    logsp <- seq(logsp[1], logsp[2], length = ngrid)
  }
  sp <- exp(logsp) # the lambda(s)
  
  gcv <- 0*sp # used to store each gcv of each lambda
  for (i in 1:length(sp)) {
    gcv[i] <- fit.pspline(x,y, sp = sp[i])$gcv
  }
  # find the largest lambda with the minmum gcv and return the model fit
  opt_gcv <- max(which(gcv == min(gcv)))
  fit.pspline(x,y,sp = sp[opt_gcv])
}



# print.pspline argument: m is just the list returned by pspline or fit.pspline
# This function prints some details of the model fit m, and silently returns a 
# list contains gcv, edf, and r2.

print.pspline <- function(m){
  cat('Order ',{m$bord},' p-spline with order ',{m$pord},
      ' penalty\nEffective degrees of freedom: ',{m$edf},
      '    Coefficients: ',{m$coefficients},'\nresidual std dev: ',
      {sqrt(m$sig2)},'    r-squared: ',{m$r2},'   GCV: ',{m$gcv},sep = '')
  
  invisible(list(gcv = m$gcv, edf = m$edf, r2 = m$r2))
}


# predict.pspline arguments:
# m is the same as above
# x is the new data (must within the range of the original data)
# se decides wither to return the standard errors
# This function firstly judges whether x is within the range of the original 
# data, if the answer is yes, form X (Xp) using the parameters in m, then predict
# the fitted value (pred_value) and calculate the corresponding standard error 
# (stde); otherwise the values would not be calculated. Finally, if se=T, return
# a list of fit (pred_value) and se (stde); otherwise return a list of fit.

predict.pspline <- function(m,x,se=TRUE){
  # judge whether x is within the range of the original data
  if (min(x) >= range(m$x)[1] & max(x) <= range(m$x)[2]){
    # the basis matrix Xp
    k = m$coefficients; bord = m$bord
    dk <- diff(range(x))/(k-bord) 
    knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
    Xp <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
    
    pred_value <- Xp %*% m$coef # the fitted value
    stde <- rowSums(Xp * (Xp %*% m$V))^.5 # corresponding standard error
    
    # judge whether se is needed
    if (se) {
      list(fit = pred_value, se = stde)
    } else {
      list(fit = pred_value)
    }
  } else {
    cat('New x values should within the range of the original data')
  }
}



# plot.pspline has the same argument as print.pspline
# It firstly calculates the lower and upper confidence limits (ll, ul), residuals,
# given m, then plot the original data. The smooth function and approximate 95% 
# credible intervals are plotted using a new sequence of data (ordered). After 
# this, the model residuals against fitted values is plotted, along with a line 
# on the x-axis. The third plot is the qqplot of residuals. Finally, the function
# returns a list of ll, ul and x silently.

plot.pspline <- function(m){
  # lower and upper confidence limits
  ll <- m$fitted-1.96*sqrt(m$se)
  ul <- m$fitted+1.96*sqrt(m$se)
  
  Residuals <- m$y-m$fitted
  
  # plot the original data x, y
  plot(m$x, m$y, main = 'Original Data & Smooth Function & 95% CI', pch=20,
       xlab = 'x', ylab = 'y', ylim = c(floor(min(m$y))-floor(min(m$y))%%10-10,
       ceiling(max(m$y))-ceiling(max(m$y))%%10+20))
  
  # plot smooth function and CI
  # generate ordered x
  x_1 <- seq(min(m$x), max(m$x), length = 200)
  pred <- predict.pspline(m, x_1) 
  y_1 <- pred$fit # ordered fitted values
  lines(x_1, y_1, col = "blue") # smooth function
  lines(x_1, y_1-1.96*pred$se, col = "red", lty = 2) # ll
  lines(x_1, y_1+1.96*pred$se, col = "red", lty = 2) # ul

  # plot the model residuals against fitted values
  plot(m$fitted, Residuals, main = 'Residuals VS Fitted Values',
       xlab = 'Fitted Values')
  # add a line on x-axis
  x_axis <- seq(min(m$fitted), max(m$fitted), length = 100)
  lines(x_axis, x_axis*0, col = 'grey', lty = 2)
  
  # qqplot
  qqnorm(Residuals)
  qqline(Residuals, col = 2, lty = 2)
  
  # return the list silently
  invisible(list(ll = ll, ul = ul, x = m$x))
}

