

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
    flags = bitflag(em$err,1:4)
    msg = c("minimization error occurred","decreasing likelihood","maximum iterations reached",
            "alpha parameters near zero")[flags]
    msg = paste0(paste0(msg,collapse=' and '),'.')
    precision = max(em$maxdif_A, em$maxdif_b)
    if(flags[4])
    {
      warning("Some discrimination parameters are too close to zero, the model is unidentified")
    }

    if(precision < .001)
    {
      message("The EM solution has a lower accuracy (~",round(precision,5),"). Reasons: ",msg,
              " See the section 'troubleshooting' in the help documentation for possible solutions.")
    } else
    {
      warning("The EM algorithm did not converge. Reasons: ",msg,
              " See the section 'troubleshooting' in the help documentation for possible solutions.",call.=FALSE)
    }
  }
}

cal_settings = list(
  theta_grid = seq(-6,6,length.out=41),
  max_em_iterations = 800L)


#' Fit a marginal model
#'
#' Estimate a one or two parameter model using Marginal Maximum Likelihood
#'
#' @param dataSrc a matrix, long format data.frame or a dexter database
#' @param predicate logical predicate to filter dataSrc, has no effect when dataSrc is a matrix
#' @param group if dataSrc is a matrix then a vector of length nrows, otherwise the names of one or more person
#' properties together grouping people. See details.
#' @param fixed_param data.frame with columns: item_id, item_score, beta and, if model is 2PL, also alpha.
#' @param se should standard errors be determined. For large datasets with many items this can take some time. Set
#' to false to save time.
#' @param prior_alpha if the estimation does not converge or gives extreme results, usually in an adaptive test or with too few
#' observations for some items in your data, you can attempt to use a prior to improve the results. Choice of
#' lognormal or normal
#' @param prior_alpha_mu first moment of prior distribution on discrimination parameters.
#' @param prior_alpha_sigma second moment of prior distribution on discrimination parameters.
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
#' In the 2pl the probability of response \ifelse{html}{\out{a<sub>ij</sub>}}{\eqn{{a}_{i,j}}} for person p (that is a response to item i, category j with score \ifelse{html}{\out{a<sub>ij</sub>}}{\deqn{{a}_{i,j}}}) is:
#' \ifelse{html}{\out{<table style="display: inline;"><tr><td style="text-align: center;border-bottom: 1px solid black;">e<sup>(-&beta;<sub>ij</sub>+&theta;<sub>p</sub>)&alpha;<sub>i</sub>a<sub>ij</sub></sup></td></tr><tr><td style="text-align: center;"><table style="display: inline;"><tr><td rowspan="2">&sum;</td><td style="font-size:50\%;">m</td><td rowspan="2">e<sup>(-&beta;<sub>ik</sub>+&theta;<sub>p</sub>)&alpha;<sub>i</sub>a<sub>ik</sub></sup></td></tr><tr><td style="font-size:50\%;">k=1</td></tr></table></td></tr></table>}}{\deqn{\frac{e^{\left(- {\beta}_{i,j} + {\theta}_{p}\right) {a}_{i,j} {\alpha}_{i}}}{\sum_{k=1}^{m} e^{\left(- {\beta}_{i,k} + {\theta}_{p}\right) {a}_{i,k} {\alpha}_{i}}}}}       
#'
#' @section CML or MML:                                
#'
#' MML estimation requires extra assumptions about the population distribution compared to CML.
#' Consequently there is rarely a good reason to use MML estimation for an 1PL since it is an exponential family
#' model and can be better estimated with CML. Only in case of adaptive data (where CML is not allowed) should you
#' use MML to estimate a 1PL.
#'
#' A 2PL cannnot be estimated using CML, except in the case of complete data (see the interaction model in dexter).
#' So for 2PL models MML is usually the method of choice.
#'
#' Note that correctly specifying grouping variables for test takers is very important in MML estimation. Failure to include
#' relevant grouping can seriously bias parameter and subsequent ability estimation.
#'
#' @section Troubleshooting:
#' The EM algorithm tries to converge on a solution up to a precision of 0.0001. It usually succeeds.
#' If it is not successful
#' a message or warning is given (dependent on the severity of the situation). The (less precise) results are
#' still returned to facilitate identification of the problem but you should generally not trust the results very much.
#'
#' The following possible solutions can be tried in such a case:
#'
#' \describe{
#' \item{omit items with too few observations}{Items may have an insufficient number of observations. This is especially problematic
#' for a 2PL. There is no clear cut lower bound independent of the purpose of the estimation. However, items with fewer
#' than 100 observations will often cause technical problems in the estimation.}
#' \item{omit problematic items}{For a 2PL, for items that have an alpha parameter very near zero, the beta parameter is undefined
#' and will not converge. You can either look for key errors, set a prior on the alpha parameter or omit these items (based on the results of the failed calibration)}
#' \item{omit fixed parameters}{If you use fixed parameters, try to calibrate without fixed parameters
#' first and plot the results against your fixed parameters. If these do not fall approximately on
#' a straight line, you might need to omit some of the fixed parameters that show the most misfit.}
#' \item{use priors in 2PL}{If the results of the calibration are extreme or otherwise unbelievable (e.g. parameters with absolute values >50) it
#' may be necessary to use a prior distribution on the discrimination parameters.
#' This may happen with adaptive test data.}
#' }
#'
fit_1pl = function(dataSrc, predicate=NULL, group = NULL, 
                   fixed_param=NULL, se=TRUE)
{
  env = caller_env()
  qtpredicate = eval(substitute(quote(predicate)))
  res = est(dataSrc, qtpredicate, env,group=group,model='1PL',
            fixed_param=fixed_param,se=se) 
  
  res
}

#' @rdname fit_1pl
fit_2pl = function(dataSrc, predicate=NULL, group = NULL, 
                   fixed_param=NULL, se=TRUE,
                   prior_alpha = c('none','lognormal','normal'),
                   prior_alpha_mu = ifelse(prior_alpha=='lognormal',0,1),
                   prior_alpha_sigma = ifelse(prior_alpha=='lognormal',0.5,0.2))
{
  prior_alpha = match.arg(prior_alpha)
  force(prior_alpha_mu)
  if(prior_alpha_sigma <=0)
    stop("prior_alpha_sigma must be larger than 0")
  priorA = switch(prior_alpha, lognormal=1L, normal=2L, 0L)
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  
  est(dataSrc, qtpredicate, env,group=group,model='2PL',
       fixed_param=fixed_param,se=se,
       priorA=priorA, priorA_mu=prior_alpha_mu, priorA_sigma=prior_alpha_sigma)
}

est = function(dataSrc, qtpredicate=NULL, env=NULL, group = NULL, model= c('1PL','2PL'),
               fixed_param=NULL, se=TRUE,
               priorA = 0L,
               priorA_mu = 0,
               priorA_sigma = 0.2)
{
  model = match.arg(model)
  
  pgw = progress_width()
  theta_grid = cal_settings$theta_grid
  max_em_iterations = cal_settings$max_em_iterations
  
  data = mml_pre(dataSrc,qtpredicate,env,group)
  
  if(nrow(data$persons) < length(data$item_id)+1)
    stop("You need more persons than items in your dataset to estimate an IRT model",call.=FALSE)
  
  pre = data$pre
  nit = length(data$item_id)
  
  if(any(pre$imax < 1))
  {
    cat('Items with maximum score 0:\n')
    print(pre$item_id[pre$imax < 1])
    stop('Some items have a maximum score of 0, model cannot be calibrated',call.=FALSE)
  }
  if(any(pre$icat[1,]==0))
  {
    cat('Items without a 0 score:\n')
    print(pre$item_id[pre$icat[1,]==0])
    stop('Some items have no 0 score category, model cannot be calibrated',call.=FALSE)
  }
  
  ref_group = which.max(data$groups$group_n) - 1L
  
  if(any(data$groups$group_n<10L))
  {
    # an arbitrary number
    # 10 is ridiculously small but at least reasonably safe against divide by 0
    stop("Each subgroup needs to contain at least 10 persons.")
  }
  
  # estimation
  design = design_matrices(pre$pni, pre$pcni, pre$pi, data$persons$c_group_nbr, length(data$item_id), nrow(data$groups))
  
  mu = rep(0, nrow(data$groups))
  sigma = rep(1, nrow(data$groups))
  
  # this changes the respons vectors px and ix in pre
  a = categorize(pre$pni, pre$pcni, pre$icnp, pre$pi,
                 pre$icat, pre$imax,max(pre$ncat), pre$px, pre$ix)
  
  if(se && pgw>0) cat("(1/2) Parameter estimation\n")
  
  if(model=='1PL')
  {
    # nominal response model
    b=start.1pl(a,pre$icat,pre$ncat)
    
    
    fixed_items = rep(0L,nit)
    if(!is.null(fixed_param))
    {
      fixed_param = tibble(item_id=data$item_id, index=1:nit) %>%
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
    
    em = estimate_nrm(a, b, pre$ncat,
                      pre$pni, pre$pcni, pre$pi, pre$px,
                      theta_grid, mu, sigma, data$groups$group_n, data$persons$c_group_nbr, fixed_items, ref_group,
                      max_iter=max_em_iterations,pgw=pgw)
    
    em_report(em)
    pop = select(data$groups, -.data$c_group_nbr)
    pop$mu = drop(em$mu)
    pop$sd = drop(em$sd)
    
    if(em$err != 0L && em$maxdif_b>.001) se=FALSE
    
    if(se)
    {
      if(pgw>0) cat("(2/2) Computing standard errors\n")
      if(!is.null(fixed_param))
      {
        w = which(fixed_items==1L)
        design$items[w,] = 0L
        design$items[,w] = 0L
        design$groups[,w] = 0L
      }
      
      hess = full_hessian_nrm(a,  em$b, pre$ncat, theta_grid, fixed_items,
                              pre$ix, pre$pni, pre$pcni, pre$pi, pre$px, data$persons$c_group_nbr, data$groups$group_n,
                              pre$ip,pre$inp, pre$icnp,
                              em$mu, em$sd, ref_group,design$items,design$groups,
                              prog_width=pgw)
      
      dx = to_dexter(em$a,em$b,pre$ncat,data$item_id,H=hess, fixed_items,ref_group+1L)
      pop$SE_mu=dx$SE_pop[seq(1,nrow(em$mu)*2,2)]
      pop$SE_sigma=dx$SE_pop[seq(2,nrow(em$mu)*2,2)]
      
      out = list(items=dx$items,pop=pop,em=em,pre=pre,hess=hess)
    } else
    {
      out = list(items=to_dexter(em$a,em$b,pre$ncat,data$item_id)$items,
                 pop=pop,em=em,pre=pre)
    }
  } else
  {
    hess=NULL
    A=rep(1,nit)
    
    b=start.2pl(a,pre$icat,pre$ncat)
    
    fixed_items = rep(0L,nit)
    if(!is.null(fixed_param))
    {
      fixed_param = tibble(item_id=data$item_id, index=1:nit) %>%
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
    
    em = estimate_pl2(a, A, b, pre$ncat,
                      pre$pni, pre$pcni, pre$pi, pre$px,
                      theta_grid, mu, sigma, data$groups$group_n, data$persons$c_group_nbr, fixed_items, 
                      pre$ip, pre$inp, pre$icnp,
                      ref_group,
                      A_prior=as.integer(priorA), A_mu=priorA_mu, A_sigma=priorA_sigma,
                      use_m2=150L,max_iter=max_em_iterations,pgw=pgw)
    
    
    
    em_report(em)
    if(em$err != 0L && max(em$maxdif_b,em$maxdif_A)>.001) se=FALSE
    
    items = tibble(item_id = rep(data$item_id,pre$ncat-1),
                   alpha = rep(em$A,pre$ncat-1L),
                   item_score = as.integer(unlist(mapply(function(i,k){ a[2:k,i] },1:nit,pre$ncat))),
                   beta = as.double(unlist(mapply(function(i,k){ em$b[2:k,i] },1:nit,pre$ncat))))
    
    pop = select(data$groups, -.data$c_group_nbr)
    pop$mu = drop(em$mu)
    pop$sd = drop(em$sd)
    
    if(se)
    {
      if(pgw>0) cat("(2/2) Computing standard errors\n")
      if(!is.null(fixed_param))
      {
        w = which(fixed_items==1L)
        design$items[w,] = 0L
        design$items[,w] = 0L
        design$groups[,w] = 0L
      }
      
      hess = full_hessian_2pl(a, em$A, em$b, pre$ncat, theta_grid, fixed_items,
                              pre$ix, pre$pni, pre$pcni, pre$pi, pre$px, data$persons$c_group_nbr, data$groups$group_n,
                              pre$ip,pre$inp, pre$icnp,
                              em$mu, em$sd, ref_group,design$items,design$groups,
                              A_prior=as.integer(priorA), A_mu=priorA_mu, A_sigma=priorA_sigma,
                              prog_width=pgw)
      
      SE = sqrt(-diag(solve(hess)))
      items$SE_alpha = NA_real_
      items$SE_beta = NA_real_
      
      i=1; px=1
      for(ix in 1:nit)
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
      
      for(g in 1:nrow(pop))
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
    out = list(items=items,pop=pop,em=em,pre=pre,hess=hess,
               prior=list(A_prior=as.integer(priorA), A_mu=priorA_mu, A_sigma=priorA_sigma))
  }
  out$theta_grid = theta_grid
  out$item_id = data$item_id
  out$pre$fixed_items = fixed_items
  out$pre$group_n = data$groups$group_n
  out$pre$ref_group = ref_group
  out$pre$group = data$persons$c_group_nbr
  out$model = model
  out$em$a = a
  class(out) = append('parms_mml',class(out))
  out
}

#' Extract information from MML fit object
#' 
#' @param object object returned by fit_1pl or fit_2pl
#' @param what information to extract
#' @param ... ignored
#' 
coef.parms_mml = function(object, what=c('items','populations'),...)
{
  what=match.arg(what)
  if(what=='items')
    return(object$items)
  if(what=='populations')
    return(object$pop)
}


print.parms_mml = function(x,...)
{
  m = if(any(x$items$item_score > 1)) ' polytomous ' else ' dichotomous '

  p = paste0( 'MML parameters for',m,x$model,
              '\nitems: ', ncol(x$em$b),
              '\npersons: ', length(x$pre$pni),
              '\niterations: ',x$em$niter,
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

# plot looks weird, check, at least sorted is not necessary

#' Plot for fitted MML models
#' 
#' The plot shows 'fit' by comparing the expected score based on the model (grey line)
#' with the average scores based on the data (black line with dots) for groups of students
#' with similar estimated ability.
#' 
#' @param x object produced by fit_1pl or fit_2pl
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
    indx = (parms$pre$icnp[i]+1L):parms$pre$icnp[i+1]
    x = tibble(item_score=parms$em$a[parms$pre$ix[indx]+1L,i], 
               theta=parms$em$thetabar[1L+parms$pre$ip[indx]]) %>%
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
    
    arw = filter(x,.data$conf_min<.data$conf_max)
    arrows(arw$m, arw$conf_min, 
             arw$m, arw$conf_max, 
             length=0.05, angle=90, code=3, col=plot.args$col)

    points(x$m,x$obs,bg=x$outlier*2,pch=21)
    lines(x$m,x$obs)
  }
}

logLik.parms_mml = function(object,...)
{
  nodes = list(...)$nodes
  if(is.null(nodes))
  {
    ll = object$em$LL
  } else
  {
    #internal
    nodes = sort(nodes)
    pre=object$pre
    e=object
    if(e$model=='1PL')
    {
      ll = loglikelihood_nrm(e$em$a, e$em$b, pre$ncat,
                             pre$pcni, pre$pi, pre$px,
                             nodes, e$em$mu, e$em$sd, pre$group)
    } else
    {
     
      ll = loglikelihood_2pl(e$em$a, e$em$A, e$em$b, pre$ncat,
                          pre$pni, pre$pcni, pre$pi, pre$px,
                          nodes, e$em$mu, e$em$sd, pre$group)
    }
  }
  c("log likelihood" = ll)
}