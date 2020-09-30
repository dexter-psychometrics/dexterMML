



#' 2pl dichotomous items, 1 population
#'
est = function(dat, discriminations=TRUE)
{
  mode(dat) = 'integer'
  pre = lapply(mat_pre(dat), drop)

  if(any(pre$imax < 1))
     stop('found items with maximum score 0')

  # assume for now max>1 means poly and not teacher weighted

  if(max(pre$imax) == 1L)
  {
    # dichotomous case
    if(any(pre$isum==pre$inp))
      stop('found items without 0 scores')

    #starting values
    # prox algorithm for dichotomous
    prox = prox_dich(
      pre$isum, pre$psum,
      pre$inp, pre$pni, pre$icnp, pre$pcni,
      pre$ip, pre$pi)

    theta = drop(prox$theta)
    beta = drop(prox$beta)

    # normalize
    m_ = weighted.mean(theta, pre$pni)
    theta = theta - m_
    s_ = sqrt(weighted.mean(theta^2,pre$pni))
    theta = theta/s_
    beta = (beta-m_)/s_

    a=NULL
    if(discriminations)
    {
      j = dexterMML:::start_lr(theta, pre$ip, pre$ix,
                   pre$inp, pre$icnp,
                   beta)
      a = drop(j$alpha)
      beta = drop(j$beta)
    }
    # so far
    return(list(start = list(a=a,beta=beta,theta=theta),
                pre = pre))


  }
  else
  {
    # see https://web.archive.org/web/20190719030511/https://www.rasch.org/rmt/rmt84k.htm
    stop('not started poly yet')
  }

}


test_it = function()
{
  library(dexter)
  library(dplyr)
  library(dexterMML)

  items = tibble(item_id=sprintf('i%03i',1:60),item_score=sample(1:4,60,replace=T),beta=rnorm(60))

  theta = rnorm(10000)

  dat = r_score(items)(theta)
  dat[dat>1]=1L

  dat[1:5000,1:20]=NA_integer_
  dat[5001:10000,41:60]=NA_integer_

  test = est(dat)
  plot(test$start$a,items$item_score)
  plot(test$start$beta,items$beta)
  plot(test$start$theta,theta)


}
