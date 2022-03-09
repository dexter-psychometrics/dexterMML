

#pars is pars_mml object or data.frame
#brings 1pl to 2pl parametrisation
abl_pre = function(dataSrc, pars, qtpredicate, env, group=NULL)
{
  data = mml_pre(dataSrc,qtpredicate,env,group)
  
  pars = simple_pars(pars, data$item_id)
  
  pre = data$pre
  if(any(pre$imax>pars$imax))
  {
    message("items with scores in your data for which no parameters are available")
    print(colnames(dataSrc)[pre$imax>pars$imax])
    stop("mismatch between parameters and data")
  }
  
  data_a = categorize(pre$pni, pre$pcni, pre$icnp, pre$pi,
                      pars$icat, pars$imax, max(pars$ncat), pre$px, pre$ix)
  
  #have to test this, it might be wrong
  mismatch = sapply(1:ncol(data_a), function(i){length(setdiff(data_a[,i],pars$a[,i]))>0})
  
  if(any(mismatch))
  {
    message("items with scores in your data for which no parameters are available")
    print(data$item_id[mismatch])
    stop("mismatch between parameters and data")
  }
  
  
  list(persons=data$persons, group=data$groups, a=data_a ,A=pars$A, b=pars$b, ncat=pars$ncat,
       pni=pre$pni, pcni=pre$pcni, pi=pre$pi, px=pre$px,
       model=pars$model)
}


#' Abilities for a 1 and 2pl
#'
#' note: for a 1PL this function calls dexter::ability
#'
#' @param dataSrc	a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param parms	object produced by function fit_1pl or fit_2pl or possibly a data.frame of parameters
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param method Maximum Likelihood (MLE), Weighted Likelihood (WLE)
#' @param unweight whether to use the weighted score or the unweighted score. Has no effect for 1PL, see details
#'
#' @returns data.frame with variables person_id and theta and se
#'
#' @details
#' for a 2pl you have the option to use the weighted or unweighted score. The weighted score gives ML estimates
#' for the regular 2pl. The unweighted score is an adaptation where ability is computed conditional on the unweighted
#' sumscore of the respondent. This means people with the same unweighted score get the same ability estimate.
#' 
#' When using a data.frame of parameters, be sure that you use the correct parametrisation. See \code{\link{fit_1pl}} for details. 
#'
ability.mml = function(dataSrc, parms, predicate=NULL, method=c('MLE','WLE'), unweight=FALSE)
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  method=match.arg(method)
  
  pre = abl_pre(dataSrc, parms, qtpredicate=qtpredicate, env=env)
  
  
  if(pre$model=='1PL')
    unweight=TRUE
  
  res = theta_2pl(pre$a, pre$A, pre$b, pre$ncat,
                  pre$pni, pre$pcni, pre$pi, pre$px,
                  WLE=(method=='WLE'), USE_A = (!unweight))
  
  se = se_theta_2pl(pre$a, pre$A, pre$b, pre$ncat,
                    pre$pni, pre$pcni, pre$pi, res$theta,
                    WLE=(method=='WLE'), USE_A = (!unweight))
  
  if(!res$success)
  {
    warning("ability estimates for some some persons did not converge.",call.=FALSE)
    res$theta[res$convergence!=0L] = NA_real_
  }
  
  tibble(person_id=pre$persons$person_id, theta=drop(res$theta), se=drop(se))
}