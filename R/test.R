



#' 2pl dichotomous items, 1 population
#'
est = function(dat, group = NULL, model= c('1PL','2PL'), se=FALSE)
{
  model = match.arg(model)
  
  mode(dat) = 'integer'
  max_score = max(dat,na.rm=TRUE)
  
  pre = lapply(mat_pre(dat, max_score), drop)

  if(any(pre$imax < 1))
     stop('found items with maximum score 0')
  if(any(pre$icat[1,]==0))
    stop('found items without a zero score')
  
  if(is.null(group))
  {
    has_groups=FALSE
    group = integer(nrow(dat))
    ref_group = 0L
    group_n = nrow(dat)
    group_id = 'population'
  } else
  {
    has_groups = TRUE
    group = as.factor(group)
    group_id = levels(group)
    group = as.integer(group) -1L
    group_n = as.integer(table(group))
    ref_group = which.max(group_n) -1L
  }

  theta_grid = seq(-6,6,.6)

  # estimation
  
  if(model == '2PL' && all(pre$ncat==2))
  {
    # dichotomous case
    # if item score other than 0/1 to do: recode

    #starting values
    # prox algorithm for dichotomous
    # for adaptive data the prox algo may be very wrong, just have to test
    prox = prox_dich(
      pre$isum, pre$psum,
      pre$inp, pre$pni, pre$icnp, pre$pcni,
      pre$ip, pre$pi)

    theta = drop(prox$theta)
    beta = drop(prox$beta)

    #normalize
    m_ = weighted.mean(theta[group==ref_group], pre$pni[group==ref_group])
    theta = theta - m_
    s_ = sqrt(weighted.mean(theta[group==ref_group]^2,pre$pni[group==ref_group]))
    theta = theta/s_
    beta = (beta-m_)/s_

    split_theta = split(theta,group)
    start_mu = sapply(split_theta, mean)
    start_var = sapply(split(theta,group), var)

    # starting values for a
    j = start_lr(theta, pre$ip, pre$ix,
                   pre$inp, pre$icnp,
                   beta)
    a = drop(j$alpha)
    beta = drop(j$beta)

    em = estimate_2pl_dich_multigroup(a, beta, pre$pni, pre$pcni, pre$pi, pre$px,
                                         theta_grid, start_mu, start_var, group_n, group, ref_group)

    items = tibble(item_id=colnames(dat),a=drop(em$a),b=drop(em$b))
    pop = tibble(group=group_id,mu=drop(em$mu),sigma=drop(em$sd))

    if(se)
    {
      nit = nrow(items)
      npop = nrow(pop)

      res = Oakes_2pl_dich(items$a, items$b, em$r0, em$r1,
                           pre$pni, pre$pcni, pre$pi, pre$px, theta_grid,
                           pop$mu, pop$sigma, group_n, group,ref_group)

      SE = sqrt(-diag(solve(res$H)))

      items$se_a=SE[seq(1,2*nit,2)]
      items$se_b=SE[seq(2,2*nit,2)]
      if(has_groups)
      {
        pop$se_mu = NA_real_
        pop$se_sigma = NA_real_
        pop$se_mu[-(ref_group+1)] = SE[seq(2*nit+1,length(SE),2)]
        pop$se_sigma[-(ref_group+1)] = SE[seq(2*nit+2,length(SE),2)]
      }
    }
    # so far
    return(list(start = list(a=a,beta=beta,theta=theta),pre = pre,
                items=items, pop=pop,theta=em$thetabar,ll=em$LL,niter=em$niter))
  }
  else if(model=='1PL')
  {
    # nominal response model
    
    # this changes the respons vectors px and ix in pre
    a = categorize(pre$inp, pre$pni, pre$icnp, pre$pcni,pre$ip, pre$pi,
                       pre$icat, pre$ncat, pre$ix, pre$px)

    # prox is een lelijk gedoetje voor poly, even gelaten
    # see https://web.archive.org/web/20190719030511/https://www.rasch.org/rmt/rmt84k.htm
    
    b = matrix(1,nrow(a),ncol(a))
    
    mu = rep(0, length(group_n))
    sigma = rep(1, length(group_n))
    
    em = estimate_nrm(a, b, pre$ncat,
                   pre$pni, pre$pcni, pre$pi, pre$px, 
                   theta_grid, mu, sigma, group_n, group, ref_group)

    return(em);
    
    
  }
  stop('2pl poly nog niet gedaan')
}
