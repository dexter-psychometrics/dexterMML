



#' 2pl dichotomous items, 1 population
#'
est = function(dat, group = NULL,se=FALSE)
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

    if(is.null(group))
    {

      # normalize
      m_ = weighted.mean(theta, pre$pni)
      theta = theta - m_
      s_ = sqrt(weighted.mean(theta^2,pre$pni))
      theta = theta/s_
      beta = (beta-m_)/s_
    }
    else
    {
      # normalize
      group = as.factor(group)
      lev = levels(group)
      group = as.integer(group)
      group_n = as.integer(table(group))
      ref_group = which.max(group_n)
      m_ = weighted.mean(theta[group==ref_group], pre$pni[group==ref_group])
      theta = theta - m_
      s_ = sqrt(weighted.mean(theta[group==ref_group]^2,pre$pni[group==ref_group]))
      theta = theta/s_
      beta = (beta-m_)/s_

      split_theta = split(theta,group)
      start_mu = sapply(split_theta, mean)
      start_var = sapply(split(theta,group), var)
      group = group -1L
    }

    # starting values for a
    j = start_lr(theta, pre$ip, pre$ix,
                   pre$inp, pre$icnp,
                   beta)
    a = drop(j$alpha)
    beta = drop(j$beta)
    #a=rep(1,length(beta))
    theta_grid = seq(-6,6,.6)


    if(is.null(group))
    {
      em = estimate_2pl_dich(a, beta, pre$pni, pre$pcni, pre$pi, pre$px,
                             theta_grid)
    } else
    {
      em = estimate_2pl_dich_multigroup(a, beta, pre$pni, pre$pcni, pre$pi, pre$px,
                                         theta_grid, start_mu, start_var, group_n, group)
    }
    items = tibble(item_id=colnames(dat),a=drop(em$a),b=drop(em$b))
    pop = tibble(group=lev,mu=drop(em$mu),sigma=drop(em$sd))

    if(se)
    {
      J = oakes(em$a, em$b, pre$pni, pre$pcni, pre$pi, pre$px,
                                       theta_grid, em$mu, em$sd, group_n, group)
      # to do: figure out hessian for mu,sigma
      J=J[1:ncol(em$obs),1:ncol(em$obs)]
      hess = em$obs+(J+t(J))/2
      SE=sqrt(-diag(solve(hess)))
      items$se_a=SE[1:ncol(dat)]
      items$se_b=SE[(1+ncol(dat)):(2*ncol(dat))]

    }
    # so far
    return(list(start = list(a=a,beta=beta,theta=theta),
                items=items,pop=pop,theta=em$thetabar,ll=em$LL,niter=em$niter,
                pre = pre))
  }
  else
  {
    # see https://web.archive.org/web/20190719030511/https://www.rasch.org/rmt/rmt84k.htm
    stop('not started poly yet')
  }

}
