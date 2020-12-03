
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

# the dexter/oplm parametrisation, without renormalisation
to_dexter = function(a,logb,ncat,item_id, H=NULL, fixed_items=NULL, ref_group=-1)
{
  n=1:length(ncat)
  any_fixed = !(is.null(fixed_items) || all(fixed_items==0))

  # without zero's
  a = as.integer(mapply(function(i,k){ a[2:k,i] },n,ncat))
  logb = as.double(mapply(function(i,k){ logb[2:k,i] },n,ncat))
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
      DD = DD[par_indx, par_indx]
      k = length(par_indx)
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

beta_matrix = function(beta,ncat)
{
  b=matrix(0,max(ncat),length(ncat))
  x=1L
  for(i in seq_along(ncat))
  {
    b[2:ncat[i],i] = beta[x:(x+ncat[i]-2L)]
    x=x+ncat[i]-1L
  }
  b
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
    b[2:ncat[i],i] = nc*(log(icat[2:ncat[i],i]/icat[1,i]) - a[1:(ncat[i]-1)])
  }
  b
}






