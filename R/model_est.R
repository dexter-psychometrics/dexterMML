
### to do: test hessian poly2 with prior, +/- correct????

#' DexterMML: MML addition to dexter
#'
#' For use in cases where CML is not possible, such as an adaptive test
#' or a multi stage test with unknown routing rules (it happens).
#'
#' In all other cases using CML (i.e. dexter or dexterMST) is strongly preferrable.
#'
#' DexterMML distinguishes itself from other MML R-packages by:
#'
#' \itemize{
#' \item including far fewer models and options
#' \item being considerably faster
#' \item support for the dexter data(base) structure
#' }
#'
"_PACKAGE"

em_gridsep = list(`1PL`=c(.6,.3,.2),
                  `2PL`=c(.3,.2,.1))

#' Estimate a model using MML
#'
#' Estimate a one or two parameter model using Marginal Maximum Likelihood
#'
#' @param dataSrc a matrix, long format data.frame or a dexter database
#' @param predicate logical predicate to filter dataSrc, has no effect when dataSrc is a matrix
#' @param group if dataSrc is a matrix then a vector of length nrows, otherwise one or more person
#' properties together grouping people. See details.
#' @param model 1PL or 2PL, see details.
#' @param fixed_param data.frame with columns: item_id, item_score, beta and, if model is 2PL, also alpha.
#' @param se should standard errors be determined. For large datasets with many items this can take some time. Set
#' to false to save time.
#' @param priorA if the estimation does not converge or gives extreme results, usually in an adaptive test or with too few
#' observations for some items in your data, you can attempt to use a prior to improve the results. Choice of
#' lognormal or normal
#' @param priorA_mu first moment of prior distribution on discrimination parameters.
#' @param priorA_sigma second moment of prior distribution on discrimination parameters.
#'
#' @return a list of things, to do: organize return value
#'
#' @details
#' In a 1PL item difficulties are estimated according to the Nominal Response Model.
#' In a 2PL, the items also get a discrimination. This
#' can be interpreted as a noise factor in the item measurement analogous to item test
#' correlations in classical test theory. Both the 1PL and 2PL model can handle polytomous data and respect the item
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
#' A 2PL cannnot be estimated using CML, except in the case of complete data (see the interaction model in dexter).
#' So for 2PL models MML is usually the method of choice.
#'
#'
est = function(dataSrc, predicate=NULL, group = NULL, model= c('1PL','2PL'),
               fixed_param=NULL, se=TRUE,
               priorA = c('none','lognormal','normal'),
               priorA_mu = ifelse(priorA=='lognormal',0,1),
               priorA_sigma = ifelse(priorA=='lognormal',0.5,0.2))
{
  model = match.arg(model)
  priorA = match.arg(priorA)

  force(priorA_mu)
  if(priorA_sigma <=0)
    stop("priorA_sigma must be larger than 0")

  priorA = switch(priorA, lognormal=1L, normal=2L, 0L)

  max_em_iterations = 500L
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()

  data = get_mml_data(dataSrc,qtpredicate,env,group)

  dat = data$dat
  group = data$group

  ## Pre and groups
  max_score = max(dat,na.rm=TRUE)

  pre = lapply(mat_pre(dat, max_score), drop)

  if(any(pre$imax < 1))
  {
    cat('Items with maximum score 0:\n')
    print(colnames(dat)[pre$imax < 1])
    stop('Some items have a maximum score of 0, model cannot be calibrated')
  }
  if(any(pre$icat[1,]==0))
  {
    cat('Items without a 0 score:\n')
    print(colnames(dat)[pre$icat[1,]==0])
    stop('Some items have no 0 score category, model cannot be calibrated')
  }

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



  # estimation
  design = design_matrices(pre$pni, pre$pcni, pre$pi, group, ncol(dat), length(group_n))


  # this needs to be done somewhat smarter to save a few itr
  b = apply(pre$icat,2, function(x) 1-log(2*x/lag(x)))
  b[1,] = 0

  if(model=='1PL')
  {
    # nominal response model
    theta_grid = seq(-6,6,.6)

    # this changes the respons vectors px and ix in pre
    a = categorize(pre$inp, pre$pni, pre$icnp, pre$pcni,pre$ip, pre$pi,
                       pre$icat, pre$imax,max(pre$ncat), pre$ix, pre$px)

    fixed_items = rep(0L,ncol(dat))
    if(!is.null(fixed_param))
    {
      fixed_param = tibble(item_id=colnames(dat), .indx.=1:ncol(dat)) %>%
        inner_join(fixed_param, by='item_id') %>%
        arrange(.data$.indx., .data$item_score)

      if(nrow(b)==2)
      {
        b[2,] = b[2,]- mean(b[2,]) + mean(fixed_param$beta)
      }

      lapply(split(fixed_param, fixed_param$item_id), function(ipar)
      {
        i = ipar$.indx.[1]
        if(all(a[2:pre$ncat[i],i] == ipar$item_score))
        {
          b[2:pre$ncat[i],i] <<- from_dexter(ipar$item_score, ipar$beta)
        } else
        {
          stop("Not implemented: mismatch between fixed parameters and data vs item scores")
        }
      })
      fixed_items[unique(fixed_param$.indx.)] = 1L
      ref_group = -1L
    }
    check_connected(design, fixed_items)

    mu = rep(0, length(group_n))
    sigma = rep(1, length(group_n))


    em = estimate_nrm(a, b, pre$ncat,
                     pre$pni, pre$pcni, pre$pi, pre$px,
                     theta_grid, mu, sigma, group_n, group, fixed_items, ref_group,
                     max_iter=max_em_iterations)
    if(se)
    {
      if(!is.null(fixed_param))
      {
        w = which(fixed_items==1L)
        design$items[w,] = 0L
        design$items[,w] = 0L
        design$groups[,w] = 0L
      }

      res = Oakes_nrm(a, em$b, pre$ncat, em$r,
                      pre$pni, pre$pcni, pre$pi, pre$px,
                      theta_grid, em$mu, em$sd, group_n, group,
                      design$items, design$groups, fixed_items, ref_group)
      # needs some rearranging after this
      dx = to_dexter(em$a,em$b,pre$ncat,colnames(dat),res$H, fixed_items,ref_group+1L)
      pop = tibble(group=group_id,mu=drop(em$mu),sd=drop(em$sd),
                   SE_mu=dx$SE_pop[seq(1,nrow(em$mu)*2,2)], SE_sigma=dx$SE_pop[seq(2,nrow(em$mu)*2,2)])

      return(list(items=dx$items,pop=pop,em=em,pre=pre,model=model))
    }
    pre$a=a
    pop=tibble(group=group_id,mu=drop(em$mu),sd=drop(em$sd))
    return(list(items=to_dexter(em$a,em$b,pre$ncat,colnames(dat))$items,
                pop=pop,em=em,pre=pre,model=model))


  } else
  {
    # this changes the response vectors px and ix in pre
    a = categorize(pre$inp, pre$pni, pre$icnp, pre$pcni,pre$ip, pre$pi,
                   pre$icat, pre$imax,max(pre$ncat), pre$ix, pre$px)

    A=rep(1,ncol(dat))

    fixed_items = rep(0L,ncol(dat))
    if(!is.null(fixed_param))
    {
      fixed_param = tibble(item_id=colnames(dat), .indx.=1:ncol(dat)) %>%
        inner_join(fixed_param, by='item_id') %>%
        arrange(.data$.indx., .data$item_score)

      if(nrow(b)==2)
      {
        b[2,] = b[2,]- mean(b[2,]) + mean(fixed_param$beta)
      }

      lapply(split(fixed_param, fixed_param$item_id), function(ipar)
      {
        i = ipar$.indx.[1]
        if(all(a[2:pre$ncat[i],i] == ipar$item_score))
        {
          b[2:pre$ncat[i],i] <<- ipar$beta
          A[i] <<- ipar$alpha[1]
        } else
        {
          stop("Not implemented: mismatch between fixed parameters and data vs item scores")
        }
      })
      fixed_items[unique(fixed_param$.indx.)] = 1L
      A[fixed_items==0L] = mean(A[fixed_items==1L])
      ref_group = -1L
    }
    check_connected(design, fixed_items)
    mu = rep(0, length(group_n))
    sigma = rep(1, length(group_n))
    for(iter in 1:3)
    {
      theta_grid = seq(-6,6,em_gridsep[['2PL']][iter])
      em = estimate_poly2(a, A, b, pre$ncat,
                          pre$pni, pre$pcni, pre$pi, pre$px,
                          theta_grid, mu, sigma, group_n, group, fixed_items, ref_group,
                          A_prior=as.integer(priorA), A_mu=priorA_mu, A_sigma=priorA_sigma,
                          max_iter=max_em_iterations)
      if(em$err==0) break
      if(iter==3)
      {
        # this is a little harsh, check whether max change is acceptably close.
        warning("estimates do not converge, results cannot be trusted")
      }
      A=em$A; b=em$b; mu=em$mu; sigma=em$sd
    }

    items = tibble(item_id = rep(colnames(dat),pre$ncat-1),
                   alpha = rep(em$A,pre$ncat-1L),
                   item_score = as.integer(a[-1,]),
                   beta = as.double(em$b[-1,]))

    pop = tibble(group=group_id,mu=drop(em$mu),sd=drop(em$sd))
    if(se)
    {
      if(!is.null(fixed_param))
      {
        w = which(fixed_items==1L)
        design$items[w,] = 0L
        design$items[,w] = 0L
        design$groups[,w] = 0L
      }

      res = Oakes_poly2(a, em$A, em$b, pre$ncat, em$r,
                      pre$pni, pre$pcni, pre$pi, pre$px,
                      theta_grid, em$mu, em$sd, group_n, group,
                      design$items, design$groups, fixed_items,ref_group,
                      A_prior=as.integer(priorA), A_mu=priorA_mu, A_sigma=priorA_sigma)

      SE = sqrt(-diag(solve(res$H)))
      items$SE_alpha = NA_real_
      items$SE_beta = NA_real_
      i=1; px=1
      for(ix in 1:ncol(dat))
      {
        k = pre$ncat[ix]
        if(fixed_items[ix]==0L)
        {
          items$SE_alpha[i:(i+k-2)] = SE[px]
          items$SE_beta[i:(i+k-2)] = SE[(px+1):(px+k-1)]
          px=px+k
        }
        i = i+k-1
      }
      pop$SE_mu = NA_real_
      pop$SE_sd = NA_real_
      for(g in 1:length(group_n))
      {
        if(g != (ref_group+1))
        {
          pop$SE_mu[g] = SE[px]
          pop$SE_sd[g] = SE[px+1L]
          px=px+2L
        }
        i=i+1L
      }
    }
    pre$a=a
    return(list(items=items,pop=pop,em=em,pre=pre,model=model))
  }

}
