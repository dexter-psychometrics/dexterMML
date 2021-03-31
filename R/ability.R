
#dataSrc is matrix
#pars is pars_mml object or data.frame
#brings 1pl to 2pl parametrisation
abl_pre = function(dataSrc, pars)
{
  pars = simple_pars(pars, colnames(dataSrc))
  # this is all still kludgy
  if(pars$model=='1PL' && !inherits(pars,'data.frame'))
  {
    pars$b = -pars$b/pars$a
    pars$b[1,] = 0
    pars$A = rep(1,ncol(dataSrc))
  }
  
  max_score = max(dataSrc, na.rm=TRUE)
  pre = lapply(mat_pre(dataSrc, max_score), drop)
  if(any(pre$imax>pars$imax))
  {
    message("items with scores in your data for which no parameters are available")
    print(colnames(dataSrc)[pre$imax>pars$imax])
    stop("mismatch between parameters and data")
  }
 
  data_a = categorize(pre$pni, pre$pcni, pre$icnp, pre$pi,
                      pars$icat, pars$imax, max(pars$ncat), pre$px, pre$ix)
  
  #have to test this, it might be wrong
  mismatch = sapply(1:ncol(dataSrc), function(i)
  {
    length(setdiff(data_a[,i],pars$a[,i]))>0  
  })
  
  if(any(mismatch))
  {
    message("items with scores in your data for which no parameters are available")
    print(colnames(dataSrc)[mismatch])
    stop("mismatch between parameters and data")
  }
  
  pid = rownames(dataSrc)
  if(is.null(pid))
    pid = 1:nrow(dataSrc)
  
  list(pid=pid, a=data_a ,A=pars$A, b=pars$b, ncat=pars$ncat,
       pni=pre$pni, pcni=pre$pcni, pi=pre$pi, px=pre$px,
       model=pars$model)
}


#' Abilities for a 1 and 2pl
#'
#' note: for a 1PL this function calls dexter::ability
#'
#' @param dataSrc	a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param parms	object produced by function fit_marginal
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param method Maximum Likelihood (MLE), Weighted Likelihood (WLE)
#' @param unweight whether to use the weighted score or the unweighted score. Has no effect for 1PL, see details
#'
#' @returns data.frame with variables person_id and theta
#'
#' @details
#' for a 2pl you have the option to use the weighted or unweighted score. The weighted score gives ML estimates
#' for the regular 2pl. The unweighted score is an adaptation where ability is computed conditional on the unweighted
#' sumscore of the repondent. This means people with the same unweighted score get the same ability estimate.,
#'
ability.mml = function(dataSrc, parms, predicate=NULL, method=c('MLE','WLE'), unweight=FALSE)
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  method=match.arg(method)

  if(!inherits(dataSrc, 'matrix'))
    dataSrc = get_resp_matrix(dataSrc,qtpredicate,env=env) # a bit wasteful
  
  pre = abl_pre(dataSrc,parms)
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

  tibble(person_id=pre$pid, theta=drop(res$theta), se=drop(se))
}
