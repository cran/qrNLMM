"predict.QRNLMM" = function(object,x = NULL,groups = NULL,covar = NULL,y = NULL,MC=1000,...){
  
  thin = 5
  
  if(class(object) != "QRNLMM")
    stop("The object must be an object of the class QRNLMM.")
  
  if(is.null(x) & is.null(groups) & is.null(covar) & is.null(y)){
    
    message("NOTE: No newdata provided. Returning fitted values for original dataset.")
    
    if(is.null(object$res$beta)){
      
      # several quantiles
      
      nquant = length(object)
      p = as.numeric(lapply(object, function (x) x$res$p))
      new.values = matrix(unlist(lapply(object, function (x) x$res$fitted.values)),
                          ncol = nquant,
                          byrow = FALSE)
      
      new.values = as.data.frame(new.values)
      colnames(new.values) = p
      
    }else{
      
      # one quantile
      
      nquant = 1
      
      new.values = data.frame(new.values = object$res$fitted.values)
      colnames(new.values) = object$res$p
      
    }
    
    return(new.values)
    
  }
  
  # Second case
  # no 'y' provided
  
  if(is.null(x) | is.null(groups))
    stop("At least both 'x' and 'groups' must be provided.")
  
  if(any(is.na(x)))
    stop("There are some NA's values in x")
  
  if(any(is.na(groups)==TRUE)) 
    stop("There are some NA's values in groups")
  
  if(nrow(as.matrix(x)) != length(groups)) 
    stop("groups does not match with  the provided data in x. (nrow(x) != length(groups))")
  
  if(is.null(object$res$beta)){
    # many quantiles
    obj = object[[1]]
    p = as.numeric(lapply(object, function (x) x$res$p))
    
  }else{
    #one quantile
    obj = object
    p   = object$res$p
  }
  
  dqnc = countall2(
    gsub(" ",
         "",
         paste0(
           deparse(
             obj$res$nlmodel),
           collapse = ""),
         fixed = TRUE)
  )
  
  d    = dqnc[1]
  q    = dqnc[2]
  nc   = dqnc[3]
  
  if(nc > 0){
    
    # are there covariates?
    
    if(is.null(covar)){
      stop("Covariates must be provided for the fitted model.")
    }
    
    if(any(is.na(covar)))
      stop("There are some NA's values in covar")
    
    if(nrow(as.matrix(covar)) != length(groups)) 
      stop("'groups' does not match with  the provided data in 'covar'. (nrow(covar) != length(groups))")
    
  }
  
  
  if(is.null(y)){
    message("NOTE: response 'y' not provided. Population curves will be computed.")
  }else{
    
    if(nrow(as.matrix(y)) != length(groups)) 
      stop("NOTE: response y does not match with the provided data in groups. (length(y) != length(groups))")
    
    message("NOTE: response 'y' provided. Individual curves will be computed.")
  }
  
  groups = as.numeric(groups)
  
  #d,q,nc
  
  n = dim(obj$res$weights)[1] #subjects
  nquant = length(p)
  nj   = c(as.data.frame(table(groups))[,2])
  
  rand = rep(0,q)
  
  new.values = matrix(NA,sum(nj),nquant)
  
  # no covariates
  
  pb <- progress_bar$new(
    format = "  Predicting... [:bar] :percent eta: :eta",
    total = length(nj)*nquant, clear = TRUE, width= 60, show_after = 0)
  
  pb$tick(0)
  #Sys.sleep(0.2)
  
  if(nc == 0){
    
    # no covariates
    # one quantile
    
    if(nquant == 1){
      
      for (j in 1:length(nj)){
        
        pos = (sum(nj[1:j-1])+1):(sum(nj[1:j]))
        
        if(is.null(y)){
          
          
          # sds = apply(
          #   X = object$res$weights,
          #   MARGIN = 2,
          #   FUN = function(x) return(
          #     fitdistrplus::fitdist(
          #       x,
          #       distr = "norm",
          #       fix.arg = list(mean = 0))$estimate
          #   )
          # )
          # 
          # rand = qnorm(p = p,sd = sds)
          
          rand = apply(object$res$weights,2,quantile,probs = p)
          
          new.values[pos,1] = object$res$nlmodel(x = x[pos],
                                                 fixed = object$res$beta,
                                                 random = rand)
          
          pb$tick()
          
        }else{
          
          chutebi0 = chutebi = rep(0,q)
          chutebibi0 = chutebibi = object$res$Psi # bibi^T
          
          bmetro = matrix(MHbi3(j=j,M=MC,
                                y1 = y[pos],x1 = x[pos],
                                cov1 = covar[pos,,drop = FALSE],
                                bi=chutebi0,
                                bibi=chutebibi0,
                                d=d,q=q,p=p,nj=nj,
                                beta=object$res$beta,
                                sigmae=object$res$sigma,
                                D=object$res$Psi,
                                nlmodel=object$res$nlmodel),q,MC)
          
          rand_j = apply(bmetro[,seq(1,MC,by = thin)],1,mean)
          
          new.values[pos,1] = object$res$nlmodel(x = x[pos],
                                                 fixed = object$res$beta,
                                                 random = rand_j)
          
          pb$tick()
          
        }
      }
      
      #no covariates
      # several quantiles
      
    }else{
      
      for (j in 1:length(nj)){ 
        
        pos = (sum(nj[1:j-1])+1):(sum(nj[1:j]))
        for(k in 1:nquant){
          
          if(is.null(y)){
            
            # sds = apply(
            #   X = object[[k]]$res$weights,
            #   MARGIN = 2,
            #   FUN = function(x) return(
            #     fitdistrplus::fitdist(
            #       x,
            #       distr = "norm",
            #       fix.arg = list(mean = 0))$estimate
            #   )
            # )
            # 
            # rand_k = qnorm(p = p[k],sd = sds)
            
            rand_k = apply(object[[k]]$res$weights,2,quantile,probs = p[k])
            
            new.values[pos,1] = obj$res$nlmodel(x = x[pos],
                                                fixed = object[[k]]$res$beta,
                                                random = rand_k)
            
            pb$tick()
            
          }else{
            
            chutebi0 = chutebi = rep(0,q)
            chutebibi0 = chutebibi = object[[k]]$res$Psi/n # bibi^T
            
            bmetro = matrix(MHbi3(j=j,M=MC,
                                  y1 = y[pos],x1 = x[pos],
                                  cov1 = covar[pos,,drop = FALSE],
                                  bi=chutebi0,
                                  bibi=chutebibi0,
                                  d=d,q=q,p=p[k],nj=nj,
                                  beta=object[[k]]$res$beta,
                                  sigmae=object[[k]]$res$sigma,
                                  D=object[[k]]$res$Psi,
                                  nlmodel=object[[k]]$res$nlmodel),q,MC)
            
            rand_j = apply(bmetro[,seq(1,MC,by = thin)],1,mean)
            
            new.values[pos,k] = obj$res$nlmodel(x = x[pos],
                                                fixed = object[[k]]$res$beta,
                                                random = rand_j)
            
            pb$tick()
            
          }
        }
      }
    } 
  }else{
    
    # with covariates
    
    # one quantile
    
    covar = as.matrix(covar)
    
    # one quantile
    
    if(nquant == 1){
      for (j in 1:length(nj)){
        
        pos = (sum(nj[1:j-1])+1):(sum(nj[1:j]))
        
        if(is.null(y)){
          
          
          # sds = apply(
          #   X = object$res$weights,
          #   MARGIN = 2,
          #   FUN = function(x) return(
          #     fitdistrplus::fitdist(
          #       x,
          #       distr = "norm",
          #       fix.arg = list(mean = 0))$estimate
          #   )
          # )
          # 
          # rand = qnorm(p = p,sd = sds)
          
          rand = apply(object$res$weights,2,quantile,probs = p)
          
          new.values[pos,1] = object$res$nlmodel(x = x[pos],
                                                 fixed = object$res$beta,
                                                 random = rand,
                                                 covar=covar[pos,,drop=FALSE])
          
          pb$tick()
          
        }else{
          
          
          chutebi0 = chutebi = rep(0,q)
          chutebibi0 = chutebibi = object$res$Psi # bibi^T
          
          bmetro = matrix(MHbi3(j=j,M=MC,
                                y1 = y[pos],x1 = x[pos],
                                cov1 = covar[pos,,drop = FALSE],
                                bi=chutebi0,
                                bibi=chutebibi0,
                                d=d,q=q,p=p,nj=nj,
                                beta=object$res$beta,
                                sigmae=object$res$sigma,
                                D=object$res$Psi,
                                nlmodel=object$res$nlmodel),q,MC)
          
          rand_j = apply(bmetro[,seq(1,MC,by = thin)],1,mean)
          
          pb$tick()
          
          
          new.values[pos,1] = object$res$nlmodel(x = x[pos],
                                                 fixed = object$res$beta,
                                                 random = rand_j,
                                                 covar=covar[pos,,drop=FALSE])
          
        }
      }
      
      # several quantiles
      
    }else{
      for (j in 1:length(nj)){ 
        
        pos = (sum(nj[1:j-1])+1):(sum(nj[1:j]))
        
        for(k in 1:nquant){
          
          if(is.null(y)){
            
            # sds = apply(
            #   X = object[[k]]$res$weights,
            #   MARGIN = 2,
            #   FUN = function(x) return(
            #     fitdistrplus::fitdist(
            #       x,
            #       distr = "norm",
            #       fix.arg = list(mean = 0))$estimate
            #   )
            # )
            # 
            # rand_k = qnorm(p = p[k],sd = sds)
            
            rand_k = apply(object[[k]]$res$weights,2,quantile,probs = p[k])
            
            new.values[pos,k] = obj$res$nlmodel(x = x[pos],
                                                fixed = object[[k]]$res$beta,
                                                random = rand_k,
                                                covar=covar[pos,,drop=FALSE])
            
            pb$tick()
            
          }else{
            
            chutebi0 = chutebi = rep(0,q)
            chutebibi0 = chutebibi = object[[k]]$res$Psi/n # bibi^T
            
            bmetro = matrix(MHbi3(j=j,M=MC,
                                  y1 = y[pos],x1 = x[pos],
                                  cov1 = covar[pos,,drop = FALSE],
                                  bi=chutebi0,
                                  bibi=chutebibi0,
                                  d=d,q=q,p=p[k],nj=nj,
                                  beta=object[[k]]$res$beta,
                                  sigmae=object[[k]]$res$sigma,
                                  D=object[[k]]$res$Psi,
                                  nlmodel=obj$res$nlmodel),q,MC)
            
            rand_j = apply(bmetro[,seq(1,MC,by = thin)],1,mean)
            
            
            pb$tick()
            
            
            new.values[pos,k] = obj$res$nlmodel(x = x[pos],
                                                fixed = object[[k]]$res$beta,
                                                random = rand_j,
                                                covar=covar[pos,,drop=FALSE])
          }
        }
      }
    }
  }
  
  new.values = as.data.frame(new.values)
  colnames(new.values) = p
  return(new.values)
  
}
