
# functions for reparametrisations of the NRM

makeD = function (a, first, last)
{
  k = length(a)
  D = matrix(0, k, k)
  tel = 1
  for (i in 1:length(first)) {
    for (j in 1:(last[i] - first[i] + 1)) {
      if (j == 1) {
        D[tel, tel] = -1/a[tel]
      }
      else {
        D[tel, tel - 1] = -1/(a[tel - 1] - a[tel])
        D[tel, tel] = 1/(a[tel - 1] - a[tel])
      }
      tel = tel + 1
    }
  }
  return(D)
}
# 1pl vs 2pl (input): plot(as.double(f1$em$b[-1,]),-as.double(f2$em$b[-1,]*f2$em$a[-1,]))

# the dexter/oplm parametrisation, without renormalisation
to_dexter = function(a,logb,ncat,item_id, H=NULL, fixed_items=NULL, ref_group=-1)
{
  n=1:length(ncat)
  any_fixed = !(is.null(fixed_items) || all(fixed_items==0))

  # without zero's
  a = as.integer(unlist(mapply(function(i,k){ a[2:k,i] },n,ncat)))
  logb = as.double(unlist(mapply(function(i,k){ logb[2:k,i] },n,ncat)))
  last = cumsum(ncat-1L)
  first = c(1,last[-length(last)]+1)

  cov.beta = NULL
  cov.all = NULL

  DD = makeD(a, first, last)

  beta = DD %*% logb

  items = tibble(item_id=rep(item_id,ncat-1), item_score = a, beta=drop(beta))
  SE_pop = NULL
  if (!is.null(H))
  {
    cov.all = solve(H)

    if(any_fixed)
    {
      par_indx = which(rep(fixed_items, ncat-1L) == 0)
      k = length(par_indx)
      DD = DD[par_indx,par_indx]
      cov.beta = DD %*% cov.all[1:k,1:k] %*% t(DD)
      items$SE_beta = NA_real_
      items$SE_beta[par_indx] = sqrt(-diag(cov.beta))
    } else
    {
      k = length(logb)
      cov.beta = DD %*% cov.all[1:k,1:k] %*% t(DD)
      items$SE_beta = sqrt(-diag(cov.beta))
    }
    w = (k+1):nrow(cov.all)
    SE_pop = sqrt(-diag(cov.all[w,w,drop=FALSE]))
    if(ref_group==1)
    {
      SE_pop = c(NA,SE_pop)
    } else if(ref_group>1)
    {
      w = ref_group*2-1
      SE_pop = c(SE_pop[1:(w-1)],NA,SE_pop[w:length(SE_pop)])
    }
  }

  list(items=items, SE_pop=SE_pop,
       cov.beta=cov.beta, cov.all=cov.all)
}



# inverse of to_dexter
# a=matrix(0:3,4,1)
# logb=matrix(c(1,runif(3,0,3)),4,1)
#
# items = dexterMML:::to_dexter(a,logb,4,'bla')$items
#
# dexterMML:::from_dexter(a=items$item_score, beta=items$beta)
# for single item
from_dexter = function(a,beta)
{
  DD = makeD(a, 1, length(a))
  solve(DD,beta)
}



start_2pl = function(a, ncat, icatg, ref_group, item_id, fixed_param=NULL)
{
  nit = ncol(a)
  A = rep(1,nit)
  fixed_items = rep(0L,nit)
  fixb = matrix(0, nrow(a),ncol(a))
  
  if(!is.null(fixed_param))
  {
    fixed_param = tibble(item_id=item_id, index=1:nit) %>%
      inner_join(fixed_param, by='item_id', suffix = c('','.ignore')) %>%
      arrange(.data$index, .data$item_score)
    
    fpar = split(fixed_param,fixed_param$index)
    
    for(x in fpar)
    {
      i = x$index[1]
      if(all(a[2:ncat[i],i] == x$item_score))
      {
        fixed_items[i] = 1L
        fixb[2:ncat[i],i] = x$beta
        A[i] = x$alpha[1]
      } else
      {
        stop("Not implemented: mismatch between fixed parameters and data vs item scores")
      }
    }
    ref_group = -1L
    A[fixed_items==0L] = exp(mean(log(A[fixed_items==1L])))
  }
  b = start_beta(a,ncat,icatg,ref_group,fixed_items, fixb)
  list(A=A, b=b, fixed_items=fixed_items, ref_group=ref_group)
}

start_1pl = function(a, ncat, icatg, ref_group, item_id, fixed_param=NULL)
{
  nit = ncol(a)
  fixed_items = rep(0L,nit)
  fixb = matrix(0, nrow(a),ncol(a))
  
  if(!is.null(fixed_param))
  {
    fixed_param = tibble(item_id=item_id, index=1:nit) %>%
      inner_join(fixed_param, by='item_id', suffix = c('','.ignore')) %>%
      arrange(.data$index, .data$item_score)
    
    fpar = split(fixed_param,fixed_param$index)
    
    for(x in fpar)
    {
      i = x$index[1]
      if(all(a[2:ncat[i],i] == x$item_score))
      {
        fixed_items[i] = 1L
        fixb[2:ncat[i],i] = x$beta
      } else
      {
        stop("Not implemented: mismatch between fixed parameters and data vs item scores")
      }
    }
    ref_group = -1L
  }
  
  b = start_beta(a,ncat,icatg,ref_group,fixed_items, fixb)
  
  for(i in 1:ncol(b))
    b[2:ncat[i],i] = from_dexter(a[2:ncat[i],i], b[2:ncat[i],i])

  list(b=b, fixed_items=fixed_items, ref_group=ref_group)
}
  



# simplify pars, arranges by items, 
# uses a different (beta) parametrisation that is used by ability functions
simple_pars = function(parms, items=NULL)
{
  if(inherits(parms,'data.frame'))
  {
    df = parms
    if(!is.null(items))
    {
      if(!all(items %in% df$item_id))
        stop('not all items present in your data have parameters')
      df = filter(df,.data$item_id %in% items)
      df$item_id = as.integer(factor(df$item_id, levels=items))
    } else
    {
      df$item_id = as.integer(factor(df$item_id))
    }
    
    df = arrange(df,.data$item_id,.data$item_score)
    
    
    ncat = max(table(df$item_id))
    max_score = max(df$item_score)
    nit = n_distinct(df$item_id)
    out = list(icat = matrix(0L, max_score+1,nit), ncat = integer(nit),
             imax = integer(nit), b = matrix(0,ncat+1,nit), a = matrix(0L,ncat+1,nit))
    
    if('alpha' %in% colnames(df))
    {
      out$model='2PL'
      out$A = distinct(df,.data$item_id,.keep_all=TRUE) %>%
        arrange(.data$item_id) %>%
        pull(.data$alpha)
    } else
    {
      out$model='1PL'
      out$A = rep(1,ncol(out$a))
    }
    
    out$icat[1,] = 1L
    items = split(df,df$item_id)
    
    for(i in seq_along(items))
    {
      itm = items[[i]]
      b = if(out$model=='1PL'){-from_dexter(itm$item_score, itm$beta)/itm$item_score}else{itm$beta}
      n = nrow(itm)
      
      out$b[2:(n+1),i] = b
      out$a[2:(n+1),i] = itm$item_score
      out$icat[itm$item_score+1L,i] = 1L
      out$ncat[i] = n + 1L
      out$imax[i] = max(itm$item_score)
    }
    return(out)
  }
  if(inherits(parms,'parms_mml'))
  {
    if(parms$model == '1PL')
    {
      b = -parms$em$b/parms$em$a
      b[1,] = 0
      A = rep(1,ncol(b))
    } else
    {
      b = parms$em$b 
      A = parms$em$A
    }
    
    if(is.null(items))
    {
      return(list(A=A,a=parms$em$a,b=b,
                  model=parms$model, icat=parms$pre$icat,
                  imax=parms$pre$imax,ncat=parms$pre$ncat))
    }
    if(!all(items %in% parms$item_id))
      stop('not all items present in your data have parameters')
    w = match(items,parms$item_id)
    return(list(A=A[w],a=parms$em$a[,w],b=b[,w],
                model=parms$model, icat=parms$pre$icat[,w],
                imax=parms$pre$imax[w], ncat=parms$pre$ncat[w]))
  }
  stop("expected a parms_mml object")
}


