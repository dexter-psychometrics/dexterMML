


to_dexter = function(a,b,ncat,item_id)
{
  n=1:length(ncat)

  a = as.integer(mapply(function(i,k){ a[1:k,i] },n,ncat))
  b = as.double(mapply(function(i,k){ b[1:k,i] },n,ncat))
  last = cumsum(ncat)
  first = c(1,last[-length(last)]+1)
  res=dexter:::toOPLM(a,b,first,last)
  tibble(item_id=rep(item_id,ncat-1), item_score = res$a, beta=res$beta)
}
