QRNLMM = function(y,x,groups,initial,exprNL,covar=NA,p=0.5,
                  precision=0.0001,MaxIter=500,M=20,cp=0.25,
                  beta=NA,sigma=NA,Psi=NA,
                  show.convergence=TRUE,CI=95,
                  verbose=TRUE)
{
  if(any(is.na(groups)==TRUE)) stop("There are some NA's values in groups")
  if(length(y) != length(groups)) stop("groups does not match with  the provided data. (length(y) != length(groups))")
  
  resexp = validate_str(exprNL)
  
  if(nchar(resexp)>0)
  {
    cat('\n')
    cat('Some defined variables or maths expressions not recognized:\n')
    cat('\n')
    cat(resexp,'\n')
    cat('\n')
    cat('* For the NL function just "x", "covar", "fixed" and "random" can be defined.\n')
    cat('\n')
    cat('* For derivating the deriv R function recognizes the arithmetic operators +, -, *, / and ^, and the single-variable functions exp, log, sin, cos, tan, sinh, cosh, sqrt, pnorm, dnorm, asin, acos, atan, gamma, lgamma, digamma and trigamma, as well as psigamma for one or two arguments (but derivative only with respect to the first).')
    stop(paste("Expression/s \"",resexp,"\" do/es not defined. More details above.",sep=""))
  }
  
  nj   = c(as.data.frame(table(groups))[,2])
  dqnc = countall(exprNL)
  d    = dqnc[1]
  q    = dqnc[2]
  nc   = dqnc[3]
  
  if(all(is.na(covar)==FALSE)){
    if(any(is.na(covar)==TRUE)) stop("There are some NA's values in covar")
    covar = as.matrix(covar)
    
    if(nc != dim(covar)[2]){stop("The number of declared covariates in exprNL must coincide with the column number of covar.")}
  }
  
  if(length(p)==1)
  {
    ## Verify error at parameters specification
    # Validacion de las dimensiones de x, groups y covar (by Gemini)
    if (any(is.na(covar))) {
      # Code to execute if any element of covar is NA
      if (length(groups) != nrow(as.matrix(x))) {
        stop("x variable does not have the same number of rows as groups.")
      }
    }
    
    #No data
    if( (length(x) == 0) | (length(y) == 0)) stop("All parameters must be provided.")
    
    #Validating if exists NA's
    if(sum(y[is.na(y)==TRUE]) > 0) stop("There are some NA's values in y")
    if(sum(y[is.na(x)==TRUE]) > 0) stop("There are some NA's values in x")
    if(any(is.na(covar)==TRUE) && all(is.na(covar))==FALSE) stop("There are some NA's values in covar")
    
    #Validating dims data set
    if(ncol(as.matrix(y)) > 1) stop("y must have just one column")
    if( length(y) != sum(nj) ) stop("nj does not match with  the provided data. (length(y) != sum(nj))")
    if( length(y) != nrow(as.matrix(x)) ) stop("x variable does not have the same number of lines than y")
    if(length(y) != nrow(as.matrix(covar)) && is.na(covar)==FALSE) stop("covar variable does not have the same number of lines than y")
    
    
    if(is.character(exprNL)==FALSE && is.expression(exprNL)==FALSE) stop("exprNL must be of class expression or character.")
    
    if(length(initial) != d) stop("The vector of initial parameter must have dimensions equal to the number of fixed effects declared in exprNL.")
    if(is.numeric(initial) == FALSE) stop("The vector of initial parameter must be of class numeric.")
    
    #Validating supports
    if(p > 1 | p < 0) stop("p must be a real number in (0,1)")
    if(precision <= 0) stop("precision must be a positive value (suggested to be small)")
    if(MaxIter <= 0 |MaxIter%%1!=0) stop("MaxIter must be a positive integer value")
    if(M <= 0 |M%%1!=0) stop("M must be a positive integer value >= 10")
    if(cp > 1 | cp < 0) stop("cp must be a real number in [0,1]")
    if(is.logical(show.convergence) == FALSE) stop("show.convergence must be TRUE or FALSE.")
    
    #Matrix column labels
    namesz <- ('b1')
    if(q>1){
      for(i in 2:q){namesz <- c(namesz, paste("b",i,sep=""))}
    }
    
    # No data
    if (!all(is.na(beta))) {
      if (length(beta) != d) {
        stop("beta must have dimensions equal to the number of fixed effects declared in exprNL.")
      }
    }
    
    if (!is.na(sigma)) {
      if (sigma <= 0 || length(sigma) != 1) {
        stop("sigma must be a positive real number")
      }
    }
    
    
    if (!all(is.na(Psi))) {
      if (ncol(as.matrix(Psi)) != q | nrow(as.matrix(Psi)) != q) {
        stop("Psi must be a square matrix of dims equal to the number of random effects declared in exprNL.")
      }
      if (det(Psi) <= 0) {
        stop("Psi must be a square symmetrical real positive definite matrix.")
      }
    }
    
    exprNL0 = exprNL
    inter  = paste("function(x,fixed,random,covar=NA){resp = ",exprNL0,";return(resp)}",sep="")
    nlmodel0 = nlmodel = eval(parse(text = inter))
    
    if(nc > 0){
      exprNL = gsub('covar\\[','covar\\[,',as.character(exprNL))
      inter = paste("function(x,fixed,random,covar){covar=as.matrix(covar);resp = ",
                    exprNL,";return(resp)}",sep="")
      nlmodel = eval(parse(text = inter))
    }
    
    # Initial values
    
    # Case 1:  beta and sigma are both NA
    if (all(is.na(beta)) && is.na(sigma)) {
      OPT <- optim(par = initial, fn = minbeta, y = y, x = x, covar = covar, p = p, q = q, nlmodel = nlmodel)
      beta <- OPT$par
      sigmae <- (1 / length(y)) * OPT$value
    }
    
    # Case 2: beta is NA, sigma is not NA
    if (all(is.na(beta)) && !is.na(sigma)) {
      OPT <- optim(par = initial, fn = minbeta, y = y, x = x, covar = covar, p = p, q = q, nlmodel = nlmodel)
      beta <- OPT$par
    }
    
    # Case 3: beta is not NA, sigma is NA
    if (all(!is.na(beta)) && is.na(sigma)) {
      OPT <- optim(par = beta, fn = minbeta, y = y, x = x, covar = covar, p = p, q = q, nlmodel = nlmodel)
      sigmae <- (1 / length(y)) * OPT$value
    }
    
    # Default value for Psi if it is NA
    if (all(is.na(Psi))) {
      Psi <- diag(q)
    }
    
    #Running the algorithm
    out <- QSAEM_NL(y = y,x = x,nj = nj,initial = initial,exprNL = exprNL0,covar = covar,p = p,precision = precision,M=M,pc=cp,MaxIter=MaxIter,beta = beta,sigmae = sigmae,D=Psi,nlmodel=nlmodel0,d=d,q=q)
    
    if(verbose){
      cat('\n')
      cat('---------------------------------------------------\n')
      cat('Quantile Regression for Nonlinear Mixed Model\n')
      cat('---------------------------------------------------\n')
      cat('\n')
      cat("Quantile =",p)
      cat('\n')
      cat("Subjects =",length(nj),";",'Observations =',sum(nj),
          ifelse(sum(nj==nj[1])==length(nj),'; Balanced =',""),
          ifelse(sum(nj==nj[1])==length(nj),nj[1],""))
      cat('\n')
      cat('\n')
      cat('- Nonlinear function \n')
      cat('\n')
      cat('nlmodel = \n')
      cat(as.character(inter))
      cat('\n')
      cat("return(resp)}")
      cat('\n')
      cat('\n')
      cat('-----------\n')
      cat('Estimates\n')
      cat('-----------\n')
      cat('\n')
      cat('- Fixed effects \n')
      cat('\n')
      print(round(out$res$table,5))
      cat('\n')
      cat('sigma =',round(out$res$sigmae,5),'\n')
      cat('\n')
      cat('Random effects \n')
      cat('\n')
      cat('i) Weights \n')
      print(head(round(out$res$weights,5)))
      cat('\n')
      cat('ii) Variance-Covariance Matrix \n')
      dimnames(out$res$D) <- list(namesz,namesz)
      print(round(out$res$D,5))
      cat('\n')
      cat('------------------------\n')
      cat('Model selection criteria\n')
      cat('------------------------\n')
      cat('\n')
      critFin <- c(out$res$loglik, out$res$AIC, out$res$BIC, out$res$HQ)
      critFin <- round(t(as.matrix(critFin)),digits=3)
      dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","HQ"))
      print(critFin)
      cat('\n')
      cat('-------\n')
      cat('Details\n')
      cat('-------\n')
      cat('\n')
      cat("Convergence reached? =",(out$res$iter < MaxIter))
      cat('\n')
      cat('Iterations =',out$res$iter,"/",MaxIter)
      cat('\n')
      cat('Criteria =',round(out$res$criterio,5))
      cat('\n')
      cat('MC sample =',M)
      cat('\n')  
      cat('Cut point =',cp)
      cat('\n')
      cat("Processing time =",out$res$time,units(out$res$time))
    }
    
    if(show.convergence == TRUE)
    {
      cpl = cp*MaxIter
      ndiag  = (q*(1+q)/2)
      npar   = d+1+ndiag
      
      labels = list()
      for(i in 1:d){labels[[i]] = bquote(beta[.(i)])}
      labels[[d+1]] = bquote(sigma)
      for(i in 1:ndiag){labels[[i+d+1]] = bquote(psi[.(i)])}
      
      par(mar=c(4, 4.5, 1, 0.5))
      op <- par(mfrow=c(ifelse(npar%%3==0,npar%/%3,(npar%/%3)+1),3))
      
      for(i in 1:npar)
      {
        plot.ts(out$conv$teta[i,],xlab="Iteration",ylab=labels[[i]])
        abline(v=cpl,lty=2)
      }
    }
    
    par(mfrow=c(1,1))
    par(mar= c(5, 4, 4, 2) + 0.1)
    
    fitted.values = rep(NA,sum(nj))
    
    if(nc == 0){
      for (j in 1:length(nj)){
        pos = (sum(nj[1:j-1])+1):(sum(nj[1:j]))
        rand = as.matrix(out$res$weights)[j,]
        fitted.values[pos] = nlmodel(x = x[pos],
                                     fixed = out$res$beta,
                                     random = rand)
      }
    }else{
      covar = as.matrix(covar)
      for (j in 1:length(nj)){
        pos = (sum(nj[1:j-1])+1):(sum(nj[1:j]))
        rand = as.matrix(out$res$weights)[j,]
        fitted.values[pos] = nlmodel(x = x[pos],
                                     fixed = out$res$beta,
                                     random = rand,
                                     covar = covar[pos,,drop = FALSE])
      }
    }
    
    res      = list(p = p,
                    iter = out$res$iter,
                    criteria = out$res$criterio,
                    nlmodel = nlmodel,
                    beta = out$res$beta,
                    weights = out$res$weights,
                    sigma= out$res$sigmae,
                    Psi = out$res$D,
                    SE=out$res$EP,
                    table = out$res$table,
                    loglik=out$res$loglik,
                    AIC=out$res$AIC,
                    BIC=out$res$BIC,
                    HQ=out$res$HQ,
                    fitted.values = fitted.values,
                    residuals = fitted.values - y,
                    time = out$res$time)
    
    
    par(mfrow=c(1,1))
    par(mar= c(5, 4, 4, 2) + 0.1)
    obj.out = list(conv=out$conv,res = res)
    class(obj.out)  =  "QRNLMM"
    return(obj.out)  
  }
  else
  {
    p = sort(unique(p))
    obj.out  = vector("list", length(p))
    
    
    ## Verify error at parameters specification
    
    #No data
    if( (length(x) == 0) | (length(y) == 0)) stop("All parameters must be provided.")
    
    #Validating if exists NA's
    if(sum(y[is.na(y)==TRUE]) > 0) stop("There are some NA's values in y")
    if(sum(y[is.na(x)==TRUE]) > 0) stop("There are some NA's values in x")
    if(any(is.na(covar)==TRUE) && all(is.na(covar))==FALSE) stop("There are some NA's values in covar")
    
    #Validating dims data set
    if(ncol(as.matrix(y)) > 1) stop("y must have just one column")
    if( length(y) != sum(nj) ) stop("nj does not match with  the provided data. (length(y) != sum(nj))")
    if( length(y) != nrow(as.matrix(x)) ) stop("x variable does not have the same number of lines than y")
    if(length(y) != nrow(as.matrix(covar)) && is.na(covar)==FALSE) stop("covar variable does not have the same number of lines than y")
    
    
    if(is.character(exprNL)==FALSE && is.expression(exprNL)==FALSE) stop("exprNL must be of class expression or character.")
    
    if(length(initial) != d) stop("The vector of initial parameter must have dimensions equal to the number of fixed effects declared in exprNL.")
    if(is.numeric(initial) == FALSE) stop("The vector of initial parameter must be of class numeric.")
    
    #Validating supports
    if (!all(p > 0 & p < 1)) stop("p vector must contain real values in (0,1)")
    if(precision <= 0) stop("precision must be a positive value (suggested to be small)")
    if(MaxIter <= 0 |MaxIter%%1!=0) stop("MaxIter must be a positive integer value")
    if(M <= 0 |M%%1!=0) stop("M must be a positive integer value >= 10")
    if(cp > 1 | cp < 0) stop("cp must be a real number in [0,1]")
    if(is.logical(show.convergence) == FALSE) stop("show.convergence must be TRUE or FALSE.")
    
    #Matrix column labels
    namesz <- ('b1')
    if(q>1){
      for(i in 2:q){namesz <- c(namesz, paste("b",i,sep=""))}
    }
    
    #pb2 = tkProgressBar(title = "QRNLMM for several quantiles",
    #                    min = 0,max = length(p), width = 300)

    cat("\n")
    pb2 <- progress_bar$new(
      format = ":what [:bar] :percent eta: :eta \n",
      total = length(p),
      clear = TRUE,
      width= 60,
      show_after = 0)
    
    #pb2$tick(len = 0,tokens = list(what = "QRNLMM: Preparing"))
    
    
    # No data
    if (!all(is.na(beta))) { # Corrected: Check if *all* are not NA
      if (length(beta) != d) {
        stop("beta must have dimensions equal to the number of fixed effects declared in exprNL.")
      }
    }
    
    if (!is.na(sigma)) { # Corrected:  scalar check is OK
      if (sigma <= 0 || length(sigma) != 1) {
        stop("sigma must be a positive real number")
      }
    }
    
    # Load required libraries (No changes needed here, but good practice to have)
    #  (Assumed this is handled elsewhere, but if needed, libraries are loaded with library())
    
    if (!all(is.na(Psi))) { # Corrected: Check if *all* are not NA
      if (ncol(as.matrix(Psi)) != q | nrow(as.matrix(Psi)) != q) {
        stop("Psi must be a square matrix of dims equal to the number of random effects declared in exprNL.")
      }
      if (det(Psi) <= 0) {
        stop("Psi must be a square symmetrical real positive definite matrix.")
      }
    }
    
    
    exprNL0 = exprNL
    inter  = paste("function(x,fixed,random,covar=NA){resp = ",exprNL0,";return(resp)}",sep="")
    nlmodel0 = nlmodel = eval(parse(text = inter))
    
    if(nc > 0){
      exprNL = gsub('covar\\[','covar\\[,',as.character(exprNL))
      inter = paste("function(x,fixed,random,covar){covar = as.matrix(covar);resp = ",
                    exprNL,";return(resp)}",sep="")
      nlmodel = eval(parse(text = inter))
    }
    
    #nlmodel = eval(parse(text = paste("function(x,fixed,random,covar=NA){resp = ",exprNL,";return(resp)}",sep="")))
    
    # Initial values
    
    # Case 1: beta and sigma are both NA
    if (all(is.na(beta)) && is.na(sigma)) {
      OPT <- optim(par = initial, fn = minbeta, y = y, x = x, covar = covar, p = p[1], q = q, nlmodel = nlmodel)
      beta <- OPT$par
      sigmae <- (1 / length(y)) * OPT$value
    }
    
    # Case 2: beta is NA, sigma is not NA
    if (all(is.na(beta)) && !is.na(sigma)) {
      OPT <- optim(par = initial, fn = minbeta, y = y, x = x, covar = covar, p = p[1], q = q, nlmodel = nlmodel)
      beta <- OPT$par
    }
    
    # Case 3: beta is not NA, sigma is NA
    if (all(!is.na(beta)) && is.na(sigma)) {
      OPT <- optim(par = beta, fn = minbeta, y = y, x = x, covar = covar, p = p[1], q = q, nlmodel = nlmodel)
      sigmae <- (1 / length(y)) * OPT$value
    }
    
    # Default value for Psi if it is NA
    if (all(is.na(Psi))) {
      Psi <- diag(q)
    }
    
    for(k in 1:length(p))
    {
      #setTkProgressBar(pb2, k-1, label=paste("Running quantile ",p[k],"   -   ",k-1,"/",length(p),sep = ""))
      
      #Running the algorithm
      
      pb2$tick(k-1,tokens = list(what = "QRNLMM: Total progress  "))
      
      out <- QSAEM_NL(y = y,x = x,nj = nj,initial = initial,exprNL = exprNL0,covar = covar,p = p[k],precision = precision,M=M,pc=cp,MaxIter=MaxIter,beta = beta,sigmae = sigmae,D=Psi,nlmodel=nlmodel0,d=d,q=q)
      
      if(verbose){
        cat('\n')
        cat('---------------------------------------------------\n')
        cat('Quantile Regression for Nonlinear Mixed Model\n')
        cat('---------------------------------------------------\n')
        cat('\n')
        cat("Quantile =",p[k])
        cat('\n')
        cat("Subjects =",length(nj),";",'Observations =',sum(nj),
            ifelse(sum(nj==nj[1])==length(nj),'; Balanced =',""),
            ifelse(sum(nj==nj[1])==length(nj),nj[1],""))
        cat('\n')
        cat('\n')
        cat('- Nonlinear function \n')
        cat('\n')
        cat('nlmodel = \n')
        cat(as.character(inter))
        cat('\n')
        cat("return(resp)}")
        cat('\n')
        cat('\n')
        cat('-----------\n')
        cat('Estimates\n')
        cat('-----------\n')
        cat('\n')
        cat('- Fixed effects \n')
        cat('\n')
        print(round(out$res$table,5))
        cat('\n')
        cat('sigma =',round(out$res$sigmae,5),'\n')
        cat('\n')
        cat('Random effects \n')
        cat('\n')
        cat('i) Weights \n')
        print(head(round(out$res$weights,5)))
        cat('\n')
        cat('ii) Variance-Covariance Matrix \n')
        dimnames(out$res$D) <- list(namesz,namesz)
        print(round(out$res$D,5))
        cat('\n')
        cat('------------------------\n')
        cat('Model selection criteria\n')
        cat('------------------------\n')
        cat('\n')
        critFin <- c(out$res$loglik, out$res$AIC, out$res$BIC, out$res$HQ)
        critFin <- round(t(as.matrix(critFin)),digits=3)
        dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","HQ"))
        print(critFin)
        cat('\n')
        cat('-------\n')
        cat('Details\n')
        cat('-------\n')
        cat('\n')
        cat("Convergence reached? =",(out$res$iter < MaxIter))
        cat('\n')
        cat('Iterations =',out$res$iter,"/",MaxIter)
        cat('\n')
        cat('Criteria =',round(out$res$criterio,5))
        cat('\n')
        cat('MC sample =',M)
        cat('\n')  
        cat('Cut point =',cp)
        cat('\n')
        cat("Processing time =",out$res$time,units(out$res$time))
        cat('\n')
        cat('\n')
      }
      
      fitted.values = rep(NA,sum(nj))
      
      if(nc == 0){
        for (j in 1:length(nj)){
          pos = (sum(nj[1:j-1])+1):(sum(nj[1:j]))
          rand = as.matrix(out$res$weights)[j,]
          fitted.values[pos] = nlmodel(x = x[pos],
                                       fixed = out$res$beta,
                                       random = rand)
        }
      }else{
        covar = as.matrix(covar)
        for (j in 1:length(nj)){
          pos = (sum(nj[1:j-1])+1):(sum(nj[1:j]))
          rand = as.matrix(out$res$weights)[j,]
          fitted.values[pos] = nlmodel(x = x[pos],
                                       fixed = out$res$beta,
                                       random = rand,
                                       covar = covar[pos,,drop = FALSE])
        }
      }
      
      res      = list(p = p[k],
                      iter = out$res$iter,
                      criteria = out$res$criterio,
                      nlmodel = nlmodel,
                      beta = out$res$beta,
                      weights = out$res$weights,
                      sigma= out$res$sigmae,
                      Psi = out$res$D,
                      SE=out$res$EP,
                      table = out$res$table,
                      loglik=out$res$loglik,
                      AIC=out$res$AIC,
                      BIC=out$res$BIC,
                      HQ=out$res$HQ,
                      fitted.values = fitted.values,
                      residuals = fitted.values - y,
                      time = out$res$time)
      
      obj.outk = list(conv=out$conv,res = res)
      obj.out[[k]] = obj.outk
      
      beta   = out$res$beta
      sigmae = out$res$sigmae
      Psi    = out$res$D
    }
    
    pb2$terminate()
    
    #close(pb2)
    
    par(mfrow=c(1,1))
    betas = eps = matrix(NA,length(p),d+1)
    
    for (i in 1:length(p))
    {
      j = p[i]
      betas[i,] = rbind(obj.out[[i]]$res$beta,obj.out[[i]]$res$sigma)
      eps[i,] = obj.out[[i]]$res$SE[1:(d+1)]
    }
    
    LIMSUP = t(betas + qnorm(1-((1-(CI/100))/2))*eps)
    LIMINF = t(betas - qnorm(1-((1-(CI/100))/2))*eps)
    labels = list()
    for(i in 1:d){labels[[i]] = bquote(beta[.(i)])}
    labels[[d+1]] = bquote(sigma)
    par(mar=c(4, 4.5, 1, 0.5))
    op <- par(mfrow=c(ifelse((d+1)%%2==0,(d+1)%/%2,((d+1)%/%2)+1),2),oma=c(0,0,2,0))
    
    for(i in 1:(d+1)){
      
      if(length(p)<4)
      {
        ymin = min(betas[,i],LIMSUP[i,],LIMINF[i,])
        ymax = max(betas[,i],LIMSUP[i,],LIMINF[i,])
        
        plot(p,betas[,i],ylim=c(ymin-2*mean(eps[,i]),ymax+2*mean(eps[,i])),xaxt='n', type='n',xlab='quantiles',ylab=labels[[i]])
        axis(side=1, at=p)
        polygon(c(p,rev(p)),c(LIMSUP[i,],rev(LIMINF[i,])), col = "gray50", border = NA)
        lines(p,betas[,i])
        lines(p,LIMSUP[i,])
        lines(p,LIMINF[i,])
      }
      else
      {
        smoothingSpline = smooth.spline(p, betas[,i], spar=0.1)
        smoothingSplineU = smooth.spline(p, betas[,i]+(qnorm(1-((1-(CI/100))/2)))*eps[,i], spar=0.1)
        smoothingSplineL = smooth.spline(p, betas[,i]-(qnorm(1-((1-(CI/100))/2)))*eps[,i], spar=0.1)
        plot(p, betas[,i], type='n',xaxt='n',xlab='quantiles',lwd=2,ylim=c(min(smoothingSplineL$y)-2*mean(eps[,i]),max(smoothingSplineU$y)+2*mean(eps[,i])),ylab=labels[[i]])
        axis(side=1, at=p)
        
        #create filled polygon in between the lines
        polygon(c(smoothingSplineL$x,rev(smoothingSplineU$x)),c(smoothingSplineU$y,rev(smoothingSplineL$y)), col = "gray50", border = NA)
        
        #plot lines for high and low range
        lines(p, betas[,i], type='l',lwd=1)
        lines(smoothingSplineU,lwd=1)
        lines(smoothingSplineL,lwd=1)
      }
    }
    title("Point estimative and 95% CI for model parameters", outer=TRUE)
    
    
    if(show.convergence=="TRUE")
    {
      cpl = cp*MaxIter
      ndiag  = (q*(1+q)/2)
      npar   = d+1+ndiag
      
      labels = list()
      for(i in 1:d){labels[[i]] = bquote(beta[.(i)])}
      labels[[d+1]] = bquote(sigma)
      for(i in 1:ndiag){labels[[i+d+1]] = bquote(psi[.(i)])}
      
      for(k in 1:length(p))
      {
        cat('\n')
        cat('-------------------------------\n')
        cat("Press [ENTER] to see next plot:")
        line <- readline()
        
        par(mar=c(4, 4.5, 1, 0.5))
        op <- par(mfrow=c(ifelse(npar%%3==0,npar%/%3,(npar%/%3)+1),3),oma=c(0,0,2,0))
        
        for(i in 1:npar)
        {
          plot.ts(obj.out[[k]]$conv$teta[i,],xlab="Iteration",ylab=labels[[i]])
          abline(v=cpl,lty=2)
        }
        title(paste("Convergence plots for quantile",p[k]), outer=TRUE)
      }
    }
    par(mfrow=c(1,1))
    par(mar= c(5, 4, 4, 2) + 0.1)
    class(obj.out)  =  "QRNLMM"
    return(obj.out)
  }
  par(mfrow=c(1,1))
  par(mar= c(5, 4, 4, 2) + 0.1)
}
