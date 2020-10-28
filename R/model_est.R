



#' Estimate a model using MML
#'
#' Estimate a one or two parameter model using Marginal Maximum Likelihood
#'
#' @param dataSrc a matrix, long format data.frame or a dexter database
#' @param predicate logical predicate to filter dataSrc, has no effect when dataSrc is a matrix
#' @param group if dataSrc is a matrix then a vector of length nrows, otherwise one or more person
#' properties together grouping people. See details.
#' @param model 1PL or 2PL, see details.
#' @param se should standard errors be determined. For large datasets with many items this can take some time. Set
#' to false to save time.
#'
#' @return a list of things, to do: organize return value
#'
#' @details
#' In a 1PL item difficulties on a logistic scale are estimated in a marginal Nominal Response Model.
#' In a 2PL model items also get a discrimination. This
#' can be interpreted as a noise factor in the item measurement analogous to item test
#' correlations in classical test theory. Both the 1PL and 2PL model can handle polytomous data (to do: 2PL poly not finished yet) and respect the item
#' scores. Missing categories (e.g. an item scored 0,1,4) are allowed.
#'
#' Specifying grouping variables for test takers is very important in MML estimation. Failure to include
#' relevant grouping can seriously bias parameter and subsequent ability estimation.
#'
#' Note that MML estimation requires extra assumptions about the population distribution compared to CML.
#' Consequently there is rarely a good reason to use MML estimation for an 1PL since it is an exponential family
#' model and can be better estimated with CML. Only in case of adaptive data (where CML is not allowed) should you
#' use MML to estimate a 1PL.
#'
#' A 2PL cannnot be estimated using CML, except in case of complete data (see the interaction model in dexter).
#' So for 2PL models MML is usually the method of choice.
#'
#'
est = function(dataSrc, predicate=NULL, group = NULL, model= c('1PL','2PL','old_2PLd'), se=TRUE)
{
  model = match.arg(model)

  # prepare data from possibly different sources
  # to do: also accept mst db
  if(inherits(dataSrc, 'matrix'))
  {
    dat = dataSrc
    mode(dat) = 'integer'
    if(!is.null(group) && length(group) != nrow(dat))
      stop(sprintf("Length of group (%i) is not equal to number of rows in data (%i)",
                   length(group),nrow(dat)))
  } else
  {
    qtpredicate = eval(substitute(quote(predicate)))
    env = caller_env()
    dat = get_resp_matrix(dataSrc,qtpredicate,env)
    # to do: possibility to handle groups in resp_matrix should become part of dexter
    if(!is.null(group))
    {
      if(!is.character(group))
        stop("Group should be a character variable")

      if(inherits(dataSrc, "DBIConnection"))
      {
        # to do: allow booklet_id??
        g = dbGetQuery(dataSrc,sprintf("SELECT person_id, %s FROM dxPersons;",
                                  paste0(group, collapse=',')))
      } else if(inherits(dataSrc, 'data.frame'))
      {
        g = distinct(dataSrc, .data$person_id, .keep_all=TRUE)
      }
      g = g %>%
        mutate(person_id = factor(.data$person_id, levels=rownames(dat))) %>%
        filter(!is.na(.data$person_id)) %>%
        arrange(as.integer(.data$person_id))

      if(length(group)== 1) group = g[[group]]
      else
      {
        nc = c(sapply(group[length(group)-1], function(x) max(nchar(g[[x]]))),0L)
        fmt = paste0("%-",nc,"s",collapse='_')
        g = as.list(g[,group])
        g$fmt = fmt
        group = do.call(sprintf, g)
      }
    }
  }

  ## datasrc preparation done

  ## Pre and groups
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

  if(model == 'old_2PLd' ) # keep for a while for testing and example starting values
  {
    if(!all(pre$ncat==2))
      stop("old 2pl routine is only for dichotomous")
    # dichotomous case
    # if item score other than 0/1 to do: recode?

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
    start_sigma = sapply(split(theta,group), sd)

    # starting values for a
    j = start_lr(theta, pre$ip, pre$ix,
                   pre$inp, pre$icnp,
                   beta)
    a = drop(j$alpha)
    beta = drop(j$beta)

    em = estimate_2pl_dich_multigroup(a, beta, pre$pni, pre$pcni, pre$pi, pre$px,
                                         theta_grid, start_mu, start_sigma, group_n, group, ref_group)

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
                       pre$icat, pre$imax,max(pre$ncat), pre$ix, pre$px)

    # prox is een lelijk gedoetje voor poly, even gelaten
    # see https://web.archive.org/web/20190719030511/https://www.rasch.org/rmt/rmt84k.htm


    b = apply(pre$icat,2, function(x) 1-log(2*x/lag(x)))
    b[1,] = 0

    mu = rep(0, length(group_n))
    sigma = rep(1, length(group_n))

    em = estimate_nrm(a, b, pre$ncat,
                   pre$pni, pre$pcni, pre$pi, pre$px,
                   theta_grid, mu, sigma, group_n, group, ref_group)

    if(se)
    {
      design = design_matrices(pre$pni, pre$pcni, pre$pi, group, ncol(dat), length(group_n))
      res = Oakes_nrm(a, em$b, pre$ncat, em$r,
                      pre$pni, pre$pcni, pre$pi, pre$px,
                      theta_grid, em$mu, em$sd, group_n, group,
                      design$items, design$groups, ref_group)

      # the Jacobian does not seem wholly senang but I cannot find a mistake in the code
      # maybe it should be done on a rerun of the estep with final parameters?

      ipar = sum(pre$ncat-1)
      dx = to_dexter(em$a,em$b,pre$ncat,colnames(dat),res$H)
      items = dx$items
      items$SE_beta = sqrt(diag(dx$cov.beta))

      pop = tibble(group=group_id,mu=drop(em$mu),sd=drop(em$sd))
      s = sqrt(-diag(dx$cov.all)[-(1:ipar)])
      r=ref_group
      if(r==0)
        s=c(NA,s)
      else
        s = c(s[1:(2*r)],NA,s[(2*r+1):length(s)])

      pop$SE_mu = s[seq(1,length(s),2)]
      pop$SE_sd = s[seq(2,length(s),2)]
      return(list(items=items,pop=pop,em=em,pre=pre))
    }

    return(list(items=to_dexter(em$a,exp(em$b),pre$ncat,colnames(dat))$items,em=em,pre=pre));


  } else
  {
    # this changes the respons vectors px and ix in pre
    a = categorize(pre$inp, pre$pni, pre$icnp, pre$pcni,pre$ip, pre$pi,
                   pre$icat, pre$imax,max(pre$ncat), pre$ix, pre$px)

    b = apply(pre$icat,2, function(x) 1-log(2*x/lag(x)))
    b[1,] = 0
    A=rep(1,ncol(dat))

    mu = rep(0, length(group_n))
    sigma = rep(1, length(group_n))

    em = estimate_poly2(a, A, b, pre$ncat,
                        pre$pni, pre$pcni, pre$pi, pre$px,
                        theta_grid, mu, sigma, group_n, group, ref_group)

    items = tibble(item_id = rep(colnames(dat),pre$ncat-1),
                   alpha = rep(em$A,pre$ncat-1L),
                   item_score = as.integer(a[-1,]),
                   beta = as.double(em$b[-1,]))

    pop = tibble(group=group_id,mu=drop(em$mu),sd=drop(em$sd))
    if(se)
    {
      design = design_matrices(pre$pni, pre$pcni, pre$pi, group, ncol(dat), length(group_n))
      res = Oakes_poly2(a, em$A, em$b, pre$ncat, em$r,
                      pre$pni, pre$pcni, pre$pi, pre$px,
                      theta_grid, em$mu, em$sd, group_n, group,
                      design$items, design$groups, ref_group)


      SE = sqrt(-diag(solve(res$H)))
      ipar = sum(pre$ncat)
      first = cumsum(lag(pre$ncat,default=1L))
      items$SE_alpha = rep(SE[first],pre$ncat-1L)
      items$SE_beta = (SE[1:ipar])[-first]
      if(has_groups)
      {
        s = SE[-(1:ipar)]
        r=ref_group
        if(r==0)
          s=c(NA,s)
        else
          s = c(s[1:(2*r)],NA,s[(2*r+1):length(s)])
        pop$SE_mu = s[seq(1,length(s),2)]
        pop$SE_sd = s[seq(2,length(s),2)]
      }
    }
    return(list(items=items,pop=pop,em=em,pre=pre))
  }

}
