
# library(mvtnorm)
# library(psych)
# library(tcltk)
# library(ald)
# library(lqr)
# library(psych)
# library(progress)

sqrtm <- function(A)
{
  if(length(A)==1)
    Asqrt=sqrt(A)
  else{
    sva <- svd(A)
    if (min(sva$d)>=0){
      Asqrt <- sva$u%*%diag(sqrt(sva$d))%*%t(sva$v)  # svd e decomposi??o espectral
    }else{
      stop("Matrix square root is not defined")
    }
  }
  return(as.matrix(Asqrt))
}

##########################################################################################
#PAUSE BETWEEN PLOTS
##########################################################################################

readkey <- function()
{
  cat ("Press [enter] to continue")
  line <- readline()
}

####################################################################
####################################################################

minbeta = function(y,x,fixed,covar,p,q,nlmodel)
{
  dif = y - nlmodel(x,fixed,rep(0,q),covar)
  return((sum(p*dif[dif>0]) - sum((1-p)*dif[dif<0])))
}

#######################################################################################
#######################################################################################

densbiv=function(j,bi,beta,sigmae,D,p,y1,x1,cov1,q=q,nj=nj,nlmodel=nlmodel)
{
  bi = matrix(bi,nrow = q,ncol = 1)
  dprod = 1
  for(k in 1:nj[j])
  {
    dprod = dprod * (dALD(y=y1[k],mu=nlmodel(x1[k],beta,bi,cov1[k,]),sigma=sigmae,p=p))
  }
  dens=as.numeric(dprod*dmvnorm(x=as.numeric(bi),mean=matrix(rep(0,q),q,1),sigma=D))
  return(dens)
}



####################################################################
####################################################################


densbiv3 =function(j,bi,beta,sigmae,D,p,y1,x1,cov1,q=q,nj=nj,nlmodel=nlmodel)
{
  bi = matrix(bi,nrow = q,ncol = 1)
  dprod = 1
  for(k in 1:nj[j])
  {
    dprod = dprod * (dALD(y=y1[k],mu=nlmodel(x1[k],beta,bi,cov1[k,,drop = FALSE]),sigma=sigmae,p=p))
  }
  dens=as.numeric(dprod*dmvnorm(x=as.numeric(bi),mean=matrix(rep(0,q),q,1),sigma=D))
  return(dens)
}

####################################################################
####################################################################

MHbi3 = function(j,M,y1,x1,cov1,bi,bibi,d,q,p,nj,beta=beta,sigmae=sigmae,D=D,nlmodel=nlmodel)
{  
  E_bi = bi
  V_bi = bibi - E_bi%*%t(E_bi)
  GEN = matrix(NA,nrow=q,ncol=(M+1))
  count = 1
  GEN[,1]=rmvnorm(n = 1,mean = E_bi,sigma=V_bi)
  
  while(count <= M)
  {
    cand = rmvnorm(n=1,mean=as.vector(GEN[,count]),sigma=V_bi)
    
    c1 = densbiv3(j=j,bi=cand,beta=beta,sigmae=sigmae,D=D,p=p,y1,x1,cov1,q,nj,nlmodel=nlmodel)*dmvnorm(x=GEN[,count],mean=as.vector(cand),sigma=V_bi)
    c2 = densbiv3(j=j,bi=GEN[,count],beta=beta,sigmae=sigmae,D=D,p=p,y1,x1,cov1,q,nj,nlmodel=nlmodel)*dmvnorm(x=as.vector(cand),mean=as.vector(GEN[,count]),sigma=V_bi)
    alfa = c1/c2
    
    if(is.nan(alfa)+0==1) {alfa=0.0001}
    if (runif(1) < min(alfa, 1))
    {
      count = count + 1
      GEN[,count] = cand
    }
  }
  return(GEN[,2:(M+1)])
}


####################################################################
####################################################################

MHbi2 = function(j,M,y1,x1,cov1,bi,bibi,d,q,p,nj,beta=beta,sigmae=sigmae,D=D,nlmodel=nlmodel)
{  
  E_bi = bi
  V_bi = bibi - E_bi%*%t(E_bi)
  GEN = matrix(NA,nrow=q,ncol=(M+1))
  count = 1
  GEN[,1]=rmvnorm(n = 1,mean = E_bi,sigma=V_bi)
  
  while(count <= M)
  {
    cand = rmvnorm(n=1,mean=as.vector(GEN[,count]),sigma=V_bi)
    
    c1 = densbiv(j=j,bi=cand,beta=beta,sigmae=sigmae,D=D,p=p,y1,x1,cov1,q,nj,nlmodel=nlmodel)*dmvnorm(x=GEN[,count],mean=as.vector(cand),sigma=V_bi)
    c2 = densbiv(j=j,bi=GEN[,count],beta=beta,sigmae=sigmae,D=D,p=p,y1,x1,cov1,q,nj,nlmodel=nlmodel)*dmvnorm(x=as.vector(cand),mean=as.vector(GEN[,count]),sigma=V_bi)
    alfa = c1/c2
    
    if(is.nan(alfa)+0==1) {alfa=0.0001}
    if (runif(1) < min(alfa, 1))
    {
      count = count + 1
      GEN[,count] = cand
    }
  }
  return(GEN[,2:(M+1)])
}

####################################################################
####################################################################

MElim = function(q)
{
  Dtest = matrix(data = 0,ncol = q,nrow = q)
  Dtest[upper.tri(Dtest, diag = T)]=1
  vDtest = as.vector(Dtest)
  Elim = matrix(data = 0,nrow = (q*(1+q)/2),ncol = q^2)
  count=1
  for(i in 1:q^2)
  {
    if(vDtest[i]==1)
    {
      Elim[count,i]=1
      count = count +1
    }
  }
  return(Elim)
}

####################################################################
####################################################################

####################################################################
####################################################################

logveroIS = function(beta,sigmae,D,y,x,covar,nj,bi,bibi,MIS,n,d,q,p,n.covar,nlmodel=nlmodel)
{
  logvero = 0
  bi = matrix(data = bi,nrow = n,ncol = q)
  
  for(j in 1:n)
  {
    y1  = y[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]#APPROVED
    x1  = matrix(x[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=1)#APPROVED
    if(n.covar>0){cov1 = matrix(covar[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=n.covar)}
    
    if(q==1){bibij=bibi[j]}else{bibij=bibi[j,,]}
    
    E_bi = bi[j,]
    V_bi = bibij - E_bi%*%t(E_bi)
    Bgen = rmvnorm(n = MIS,mean = E_bi,sigma=matrix(V_bi,q,q))
    sum = 0
    
    for(l in 1:MIS)
    {
      multden = 1
      for(k in 1:nj[j])
      {
        multden = multden * (dALD(y=y1[k],mu=nlmodel(x1[k],beta,Bgen[l,],cov1[k,]),sigma=sigmae,p=p))
      }
      multt = multden*(dmvnorm(x = Bgen[l,],mean=rep(0,q),sigma=D)/dmvnorm(x = Bgen[l,],mean = E_bi,sigma=V_bi))
      sum = sum + multt
    }
    prom    = sum/MIS
    logvero = logvero + log(prom)
  }
  return(logvero)
}

####################################################################
####################################################################

logveroMC = function(beta,sigmae,D,y,x,covar,nj,MC,n,d,q,p,n.covar,nlmodel=nlmodel)
{
  logvero = 0
  
  for(j in 1:n)
  {
    y1  = y[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]#APPROVED
    x1  = matrix(x[(sum(nj[1:j-1])+1):(sum(nj[1:j]))],ncol=1)#APPROVED
    if(n.covar>0){cov1 = matrix(covar[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=n.covar)}
    
    Bgen = rmvnorm(n = MC,mean=rep(0,q),sigma=matrix(D,q,q))
    sum = 0
    
    for(l in 1:MC)
    {
      multden = 1
      for(k in 1:nj[j])
      {
        multden = multden * (dALD(y=y1[k],mu=nlmodel(x1[k],beta,Bgen[l,],cov1[k,]),sigma=sigmae,p=p))
      }
      sum = sum + multden
    }
    prom    = sum/MC
    logvero = logvero + log(prom)
  }
  return(logvero)
}

###############################################################################
###############################################################################

countCharOccurrences <- function(char, s) {
  s2 <- gsub(char,"",s)
  return (nchar(s) - nchar(s2))
}

countCharOccurrences2 <- function(char, s) {
  ll = nchar(char)
  s2 <- gsub(char,"",s,fixed = TRUE)
  return ((nchar(s) - nchar(s2))/ll)
}

countall = function(s){
  frecF = rep(NA,10)
  frecR = rep(NA,10)
  frecC = rep(NA,10)
  for(i in 1:10){
    strF = paste0("fixed[",i,"]")
    strR = paste0("random[",i,"]")
    strC = paste0("covar[",i,"]")
    frecF[i] = countCharOccurrences2(strF,s)
    frecR[i] = countCharOccurrences2(strR,s)
    frecC[i] = countCharOccurrences2(strC,s)
  }
  return(c(sum(frecF > 0),sum(frecR > 0),sum(frecC > 0)))
}



countall2 = function(s){
  frecF = rep(NA,10)
  frecR = rep(NA,10)
  frecC = rep(NA,10)
  for(i in 1:10){
    strF = paste0("fixed[",i,"]")
    strR = paste0("random[",i,"]")
    strC = paste0("covar[,",i,"]")
    frecF[i] = countCharOccurrences2(strF,s)
    frecR[i] = countCharOccurrences2(strR,s)
    frecC[i] = countCharOccurrences2(strC,s)
  }
  return(c(sum(frecF > 0),sum(frecR > 0),sum(frecC > 0)))
}


#countfixedrandom(exprNL)

###############################################################################
###############################################################################

validate_str = function(prueba)
{
  prueba = gsub("\n","   ",prueba)
  prueba = gsub("[[:punct:]]","   ",gsub("[|[[:digit:]|]]","    ",prueba))
  prueba = gsub("\\d","   ",prueba)
  prueba = paste("  ",prueba,"  ",sep = "")
  functions = c(" fixed "," random "," covar "," exp "," x "," log "," sinh "," cosh "," sqrt "," pnorm "," dnorm "," asin "," acos "," atan "," lgamma "," digamma "," trigamma "," psigamma "," sin "," cos "," tan "," gamma ")
  for(i in 1:length(functions))
  {
    prueba = gsub(functions[i],"  ",prueba)
    
  }
  prueba = gsub("[[:blank:]]"," ",prueba)
  prueba = gsub("^ *|(?<= ) | *$","", prueba, perl=T)
  return(prueba)
}

###############################################################################
###############################################################################


group.plot = function(x, y, groups, ...) {
  if (length(y) != length(groups))
    stop("groups does not match with the provided data. (length(y) != length(groups))")
  if (length(x) != length(groups))
    stop("groups does not match with the provided data. (length(x) != length(groups))")
  if (length(y) != length(x))
    stop("the provided data does not match. (length(y) != length(x))")
  
  # Calcular nj DESPUeS de filtrar los datos
  nj <- as.numeric(table(as.character(groups)))
  
  Ymatrix = Xmatrix = matrix(data = NA, nrow = length(nj), ncol = max(nj))
  
  for (j in 1:length(nj)) {
    # Corregir el calculo de los indices para que coincidan con los datos filtrados
    elementos_grupo_j <- which(groups == unique(groups)[j])
    Xmatrix[j, 1:nj[j]] <- x[elementos_grupo_j]
    Ymatrix[j, 1:nj[j]] <- y[elementos_grupo_j]
  }
  matplot(x = t(Xmatrix), y = t(Ymatrix), ...)
}


#group.plot(x = Time,y = weight,groups = Plot)


group.lines = function(x, y, groups, ...) {
  if (length(y) != length(groups))
    stop("groups does not match with the provided data. (length(y) != length(groups))")
  if (length(x) != length(groups))
    stop("groups does not match with the provided data. (length(x) != length(groups))")
  if (length(y) != length(x))
    stop("the provided data does not match. (length(y) != length(x))")
  
  # Calculate nj AFTER filtering the data
  nj <- as.numeric(table(as.character(groups)))
  
  Ymatrix = Xmatrix = matrix(data = NA, nrow = length(nj), ncol = max(nj))
  
  for (j in 1:length(nj)) {
    # Correct the index calculation to match the filtered data
    elementos_grupo_j <- which(groups == unique(groups)[j])
    Xmatrix[j, 1:nj[j]] <- x[elementos_grupo_j]
    Ymatrix[j, 1:nj[j]] <- y[elementos_grupo_j]
  }
  matlines(x = t(Xmatrix), y = t(Ymatrix), ...)
}

#group.lines(x = Time,y = weight,groups = Plot)


group.points = function(x, y, groups, ...) {
  if (length(y) != length(groups))
    stop("groups does not match with the provided data. (length(y) != length(groups))")
  if (length(x) != length(groups))
    stop("groups does not match with the provided data. (length(x) != length(groups))")
  if (length(y) != length(x))
    stop("the provided data does not match. (length(y) != length(x))")
  
  # Calculate nj AFTER filtering the data
  nj <- as.numeric(table(as.character(groups)))
  
  Ymatrix = Xmatrix = matrix(data = NA, nrow = length(nj), ncol = max(nj))
  
  for (j in 1:length(nj)) {
    # Correct the index calculation to match the filtered data
    elementos_grupo_j <- which(groups == unique(groups)[j])
    Xmatrix[j, 1:nj[j]] <- x[elementos_grupo_j]
    Ymatrix[j, 1:nj[j]] <- y[elementos_grupo_j]
  }
  matpoints(x = t(Xmatrix), y = t(Ymatrix), ...)
}


#group.points(x = Time,y = weight,groups = Plot)