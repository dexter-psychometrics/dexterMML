
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


to_dexter = function(a,logb,ncat,item_id, H=NULL)
{
  n=1:length(ncat)

  # without zero's
  a = as.integer(mapply(function(i,k){ a[2:k,i] },n,ncat))
  logb = as.double(mapply(function(i,k){ logb[2:k,i] },n,ncat))
  last = cumsum(ncat-1L)
  first = c(1,last[-length(last)]+1)

  cov.beta = NULL
  cov.all = NULL

  DD = makeD(a, first, last)

  beta = DD %*% logb

  k = length(logb)
  if (!is.null(H))
  {
    cov.all = solve(H)
    cov.beta = DD %*% cov.all[1:k,1:k] %*% t(DD)
  }


  list(items=tibble(item_id=rep(item_id,ncat-1), item_score = a, beta=drop(beta)),
       cov.beta=cov.beta,cov.all=cov.all)
}


