
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

bitflag = function(flag,i=1:32) as.integer(intToBits(flag))[i]==1L

em_report = function(em)
{
  if(em$err>0)
  {
    flags = bitflag(em$err,1:3)
    msg = c("minimization error occurred","decreasing likelihood","maximum iterations reached")[flags]
    msg = paste0(paste0(msg,collapse=' and '),'.')
    precision = max(em$maxdif_A, em$maxdif_b)

    if(precision < .001)
    {
      message("The EM solution has a lower accuracy (~",round(precision,5),"). Reasons: ",msg,
              " See the section 'troubleshooting' in the help documentation for possible solutions.")
    } else if(precision < .01)
    {
      warning("The EM solution has low accuracy (~",round(precision,4),"). Reasons: ",msg,
              " See the section 'troubleshooting' in the help documentation for possible solutions.",call.=FALSE)
    } else
    {
      warning("The EM algorithm did not converge. Reasons: ",msg,
              " See the section 'troubleshooting' in the help documentation for possible solutions.",call.=FALSE)
    }
  }
}


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
#' @return an object of type parms_mml, see \code{\link{coef.parms_mml}}.
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
#' @section troubleshooting:
#' The EM algorithm tries to converge on a solution up to a precision of 0.0001. It usually succeeds.
#' If it is not successful
#' a message or warning is given (dependent on the severity of the situation). A message can usually be ignored if
#' you are happy with a slightly lower precision. In case of a warning, the (less precise) results are
#' still returned to facilitate identification of the problem but you should generally not trust the results very much.
#'
#' The following possible solutions can be tried in such a case:
#'
#' \describe{
#' \item{omit fixed parameters}{If you use fixed parameters, try to calibrate without fixed parameters
#' first and plot the results against your fixed parameters. If these do not fall approximately on
#' a straight line, you might need to omit some of the fixed parameters that show the most misfit.}
#' \item{linear subsets 1PL}{If your testdata includes adaptive or random tests, but it contains
#' a significant subset of items that was administered only as linear or multi stage tests, you can use the functions fit_enorm
#' or fit_enorm_mst in dexter and dexterMST respectively to fit a CML solution on this subset. Next you
#' can fit the complete dataset with dexterMML while fixing the parameters of the linearly administered items.}
#' \item{priors in 2PL}{If the results of the calibration are extreme (e.g. parameters with absolute values >50) it
#' might be necessary to use a prior distribution on the discrimination parameters.
#' This may happen with adaptive test data.}
#' }
#'
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

  pgw = getOption("width")
  theta_grid = seq(-6,6,.3)
  max_em_iterations = 800L
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
  
  # this changes the respons vectors px and ix in pre
  a = categorize(pre$pni, pre$pcni, pre$pi,
                       pre$icat, pre$imax,max(pre$ncat), pre$px)

  if(se) cat("(1/2) Parameter estimation\n")

  if(model=='1PL')
  {
    # nominal response model
    b=start.1pl(a,pre$icat,pre$ncat)

    
    fixed_items = rep(0L,ncol(dat))
    if(!is.null(fixed_param))
    {
      fixed_param = tibble(item_id=colnames(dat), index=1:ncol(dat)) %>%
        inner_join(fixed_param, by='item_id',suffix=c('','.y')) %>%
        arrange(.data$index, .data$item_score)

      fpar = split(fixed_param,fixed_param$index)
      shift = 0
      for(x in fpar)
      {
        i = x$index[1]
        if(all(a[2:pre$ncat[i],i] == x$item_score))
        {
          shift = shift - log(sum(exp(b[2:pre$ncat[i],i])))
          b[2:pre$ncat[i],i] = from_dexter(x$item_score, x$beta)
          shift = shift + log(sum(exp(b[2:pre$ncat[i],i])))
        } else
        {
          stop("Not implemented: mismatch between fixed parameters and data vs item scores")
        }
      }

      fixed_items[unique(fixed_param$index)] = 1L
      shift = shift / sum(fixed_items)
      b[2:nrow(b),fixed_items==0L] = b[2:nrow(b),fixed_items==0L] + shift
      
      ref_group = -1L
    }
    check_connected(design, fixed_items)

    mu = rep(0, length(group_n))
    sigma = rep(1, length(group_n))

    em = estimate_nrm(a, b, pre$ncat,
                        pre$pni, pre$pcni, pre$pi, pre$px,
                        theta_grid, mu, sigma, group_n, group, fixed_items, ref_group,
                        max_iter=max_em_iterations,pgw=pgw)

    em_report(em)
    if(em$err != 0L && em$maxdif_b>.001) se=FALSE

    if(se)
    {
      cat("(2/2) Computing standard errors\n")
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
                      design$items, design$groups, fixed_items, ref_group,pgw=pgw)
      # needs some rearranging after this
      dx = to_dexter(em$a,em$b,pre$ncat,colnames(dat),res$H, fixed_items,ref_group+1L)
      pop = tibble(group=group_id,mu=drop(em$mu),sd=drop(em$sd),
                   SE_mu=dx$SE_pop[seq(1,nrow(em$mu)*2,2)], SE_sigma=dx$SE_pop[seq(2,nrow(em$mu)*2,2)])
      out = list(items=dx$items,pop=pop,em=em,pre=pre,model=model)
    } else
    {
      pop=tibble(group=group_id,mu=drop(em$mu),sd=drop(em$sd))
      out = list(items=to_dexter(em$a,em$b,pre$ncat,colnames(dat))$items,
                  pop=pop,em=em,pre=pre,model=model)
    }
  } else
  {
    A=rep(1,ncol(dat))

    b=start.2pl(a,pre$icat,pre$ncat)

    fixed_items = rep(0L,ncol(dat))
    if(!is.null(fixed_param))
    {
      fixed_param = tibble(item_id=colnames(dat), index=1:ncol(dat)) %>%
        inner_join(fixed_param, by='item_id',suffix=c('','.y')) %>%
        arrange(.data$index, .data$item_score)
      
      fpar = split(fixed_param,fixed_param$index)
      shift = 0
      for(x in fpar)
      {
        i = x$index[1]
        if(all(a[2:pre$ncat[i],i] == x$item_score))
        {
          shift = shift - log(sum(exp(b[2:pre$ncat[i],i])))
          b[2:pre$ncat[i],i] = x$beta
          shift = shift + log(sum(exp(b[2:pre$ncat[i],i])))
          A[i] = x$alpha[1]
        } else
        {
          stop("Not implemented: mismatch between fixed parameters and data vs item scores")
        }
      }
      
      fixed_items[unique(fixed_param$index)] = 1L
      shift = shift / sum(fixed_items)
      b[2:nrow(b),fixed_items==0L] = b[2:nrow(b),fixed_items==0L] + shift
      A[fixed_items==0L] = mean(A[fixed_items==1L])
      ref_group = -1L
    }
    check_connected(design, fixed_items)
    mu = rep(0, length(group_n))
    sigma = rep(1, length(group_n))

    em = estimate_poly2(a, A, b, pre$ncat,
                          pre$pni, pre$pcni, pre$pi, pre$px,
                          theta_grid, mu, sigma, group_n, group, fixed_items, ref_group,
                          A_prior=as.integer(priorA), A_mu=priorA_mu, A_sigma=priorA_sigma,
                          max_iter=max_em_iterations,pgw=pgw)

    em_report(em)
    if(em$err != 0L && max(em$maxdif_b,em$maxdif_A)>.001) se=FALSE

    items = tibble(item_id = rep(colnames(dat),pre$ncat-1),
                   alpha = rep(em$A,pre$ncat-1L),
                   item_score = as.integer(a[-1,]),
                   beta = as.double(em$b[-1,]))

    pop = tibble(group=group_id,mu=drop(em$mu),sd=drop(em$sd))
    if(se)
    {
      cat("(2/2) Computing standard errors\n")
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
                      A_prior=as.integer(priorA), A_mu=priorA_mu, A_sigma=priorA_sigma,pgw=pgw)

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
    out = list(items=items,pop=pop,em=em,pre=pre,model=model)
  }
  out$theta_grid = theta_grid
  out$item_id=colnames(dat)
  out$em$a=a
  class(out) = append('parms_mml',class(out))
  out
}


#' Extract information from MML fit object
#' 
#' @param object object returned by est
#' @param what information to extract
#' @param ... ignored
#' 
#' 
coef.parms_mml = function(object, what=c('items','populations','likelihood'),...)
{
  what=match.arg(what)
  if(what=='items')
    return(object$items)
  if(what=='populations')
    return(object$pop)
  if(what=='likelihood')
    return(c("log likelihood"=object$em$LL))
}


print.parms_mml = function(x,...)
{
  m = if(any(x$items$item_score) > 1) ' polytomous ' else ' dichotomous '

  p = paste0( 'MML parameters for',m,x$model,
              '\nitems: ', ncol(x$em$b),
              '\npersons: ', length(x$pre$pni),
              '\niterations:',x$em$niter,
              '\nUse coef() or coefficients() to extract the item parameters.\n')

  cat(p)
  invisible(x)
}

merge_arglists = function(args, default = NULL, override = NULL)
{
  if(!is.null(default))
    args = modifyList(default, args)
  
  if(!is.null(override))
    args = modifyList(args, override)
  
  args
}

#' Plot for fitted MML models
#' 
#' The plot shows 'fit' by comparing the expected score based on the model (grey line)
#' with the average scores based on the data (black line with dots) for groups of students
#' with similar estimated ability.
#' 
#' @param x object produced by fit_enorm
#' @param items item_id's of items to plot, if NULL, one plot for each item is made
#' @param nbins number of ability groups
#' @param ci confidence interval for the error bars, between 0 and 1. 
#' Default = 0.95 for a 95\% confidence interval
#' @param ... further arguments to plot
#' 
#' @method plot parms_mml
#' 
plot.parms_mml = function(x,items=NULL,nbins=5,ci=.95,...)
{
  parms=x
  if(parms$model=='1PL')
  {
    A=rep(1,ncol(parms$em$b))
    b=-parms$em$b/parms$em$a
  } else
  {
    A=parms$em$A
    b=parms$em$b
  }
  if(is.null(items))
    items = parms$item_id
  ii = match(items,parms$item_id)
  if(anyNA(ii))
    stop(paste('Items:', paste(items[is.na(ii)],collapse=', '),'not found.'),call.=FALSE)
  
  qnt = abs(qnorm((1-ci)/2))
  cmin = function(p, n) pmax(0, p - qnt * sqrt(p*(1-p)/n))
  cmax = function(p, n) pmin(1, p + qnt * sqrt(p*(1-p)/n))
  
  user.args = list(...)
  default.args = list(bty='l',xlab = expression(theta), ylab='score',main='$item_id',col='grey80')
  
  for(i in ii)
  {
    max_score = parms$pre$imax[i]
    x = plot_data(parms$pre$pcni, parms$pre$pi, parms$pre$px, parms$pre$inp,parms$em$thetabar, 
                  parms$em$a,i-1L) %>%
      mutate(bin=ntile(.data$theta,nbins)) %>%
      group_by(.data$bin) %>%
      summarise(m=mean(.data$theta),obs=mean(.data$item_score),n=n()) %>%
      ungroup() %>%
      mutate(expected = E_score(.data$m, A, parms$em$a, b, i-1L, parms$pre$ncat)) %>%
      mutate(conf_min = max_score * cmin(.data$expected/max_score, .data$n),
             conf_max = max_score * cmax(.data$expected/max_score, .data$n)) %>%
      mutate(outlier = .data$obs < .data$conf_min | .data$obs > .data$conf_max)
    
    
    plot.args = merge_arglists(user.args,
                              default=default.args,
                               override=list(x = x$m,y = x$expected, type="l",
                                             ylim=c(0,parms$pre$imax[i])))
    plot.args$main = gsub('$item_id',parms$item_id[i],plot.args$main,fixed=TRUE)
    do.call(plot,plot.args)
    
    suppressWarnings({
      arrows(x$m, x$conf_min, 
             x$m, x$conf_max, 
             length=0.05, angle=90, code=3, col=plot.args$col)})
    
    points(x$m,x$obs,bg=x$outlier*2,pch=21)
    lines(x$m,x$obs)
  }
}