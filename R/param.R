


to_dexter = function(a,b,ncat,item_id, H=NULL)
{
  n=1:length(ncat)

  a = as.integer(mapply(function(i,k){ a[1:k,i] },n,ncat))
  b = as.double(mapply(function(i,k){ b[1:k,i] },n,ncat))
  last = cumsum(ncat)
  first = c(1,last[-length(last)]+1)



  res=dexter.toOPLM(a,b,first,last,H)

  list(items=tibble(item_id=rep(item_id,ncat-1), item_score = res$a, beta=res$beta),
       acov=res$cov.beta)
}


dexter.toOPLM = function (a, b, first, last, H = NULL, fixed_b = NULL)
{
  b_rn = b
  a_org = a

  if (!is.null(fixed_b))
    fixed_b = fixed_b[-first]
  tmp = remove_zero(a, b, first, last)
  b = tmp$b
  a = tmp$a
  first = tmp$first
  last = tmp$last
  logb = log(b)
  cov.beta = NULL


    DD = makeD(a, first, last)
    if (is.null(fixed_b)) {
      beta = DD %*% logb
      b_rn = b_rn * exp(mean(beta) * a_org)
      k = length(b)
      CC = matrix(-1/k, k, k)
      diag(CC) = (k - 1)/k
      beta = CC %*% beta
      if (!is.null(H)) {
        A = CC %*% DD
        cov.beta = solve(H)
        cov.beta = A %*% cov.beta %*% t(A)
      }
    }
    else {
      beta = DD %*% logb
      if (!is.null(H)) {
        fixed_set = which(!is.na(fixed_b))
        cov.beta = solve(H[-fixed_set, -fixed_set])
        cov.beta = DD[, -fixed_set] %*% cov.beta %*%
          t(DD[, -fixed_set])
      }
    }

  return(list(beta = beta, cov.beta = cov.beta, a = a, b_renorm = b_rn,
              first = first, last = last))
}

remove_zero = function (a, b, first, last)
{
  if (is.matrix(b)) {
    b = b[, -first]
    if (is.null(dim(b)))
      b = as.matrix(t(b))
  }
  if (is.vector(b))
    b = b[-first]
  a = a[-first]
  new_first = first
  new_last = last - 1L
  for (i in 2:length(first)) {
    ncat = last[i] - first[i]
    new_first[i] = new_last[i - 1] + 1L
    new_last[i] = new_first[i] + ncat - 1L
  }
  return(list(a = a, b = b, first = new_first, last = new_last))
}


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
