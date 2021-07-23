
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
# a=items$item_score; beta=items$beta
# for single item
from_dexter = function(a,beta)
{
  DD = makeD(a, 1, length(a))
  solve(DD,beta)
}

# to do: test for items with >2 categories
start.1pl = function(a,icat,ncat)
{
  b=matrix(0,nrow(a),ncol(a))
  nc = sqrt(1.702)
  for(i in seq_along(ncat))
  {
    icat_i = icat[a[1:ncat[i],i]+1L,i]
    b[2:ncat[i],i] = nc*(log(icat_i[2:ncat[i]]/icat_i[1]) - a[1:(ncat[i]-1)])
  }
  b
}

start.2pl = function(a,icat,ncat)
{
  b=matrix(0,nrow(a),ncol(a))
  nc = -sqrt(1.702)
  for(i in seq_along(ncat))
  {
    icat_i = icat[a[1:ncat[i],i]+1L,i]
    b[2:ncat[i],i] = nc*(log(icat_i[2:ncat[i]]/icat_i[1]))
  }
  b
}

#simplify pars, arranges by items
simple_pars = function(parms, items=NULL)
{
  if(inherits(parms,'data.frame'))
  {
    warning("data.frame parameters in dexterMML is experimental")
    df = parms
    if(!is.null(items))
    {
      if(!all(items %in% df$item_id))
        stop('not all items present in your data have parameters')
      df = filter(df,.data$item_id %in% items)
    }
    df$item_id = as.integer(factor(df$item_id, levels=items))
    
    df = arrange(df,.data$item_id,.data$item_score)
    
    
    ncat = max(df$item_score)
    nit = n_distinct(df$item_id)
    out = list(icat = matrix(0L, ncat+1,nit), ncat = integer(nit),
             imax = integer(nit), b = matrix(0,ncat+1,nit), a = matrix(0L,ncat+1,nit))
    
    if('alpha' %in% colnames(df))
    {
      out$A = distinct(df,.data$item_id,.keep_all=TRUE) %>%
        arrange(.data$item_id) %>%
        pull(.data$alpha)
      out$model='2PL'
    } else
    {
      out$model='1PL'
    }
    
    out$icat[1,] = 1L
    itm = split(df,df$item_id)
    
    for(i in seq_along(itm))
    {
      out$b[2:(nrow(itm[[i]])+1),i] = itm[[i]]$beta
      out$a[2:(nrow(itm[[i]])+1),i] = itm[[i]]$item_score
      out$icat[itm[[i]]$item_score+1L,i] = 1L
      out$ncat[i] = nrow(itm[[i]]) + 1L
      out$imax[i] = max(itm[[i]]$item_score)
    }
    return(out)
  }
  if(inherits(parms,'parms_mml'))
  {
    if(is.null(items))
    {
      return(list(A=parms$em$A,a=parms$em$a,b=parms$em$b,
                  model=parms$model, icat=parms$pre$icat,
                  imax=parms$pre$imax,ncat=parms$pre$ncat))
    }
    if(!all(items %in% parms$item_id))
      stop('not all items present in your data have parameters')
    w = match(parms$item_id, items)
    w = w[!is.na(w)]
    return(list(A=parms$em$A[w],a=parms$em$a[,w],b=parms$em$b[,w],
                model=parms$model, icat=parms$pre$icat[,w],
                imax=parms$pre$imax[w], ncat=parms$pre$ncat[w]))
  }
  stop("expected a parms_mml object")
}


