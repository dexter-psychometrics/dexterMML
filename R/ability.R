


#' Abilities for a 1 and 2pl
#'
#' note: for a 1PL this function calls dexter::ability
#'
#' @param dataSrc	a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param parms	object produced by function est
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param method Maximum Likelihood (MLE), Weighted Likelihood (WLE)
#' @param unweight whether to use the weighted score or the unweighted score. Has no effect for 1PL, see details
#'
#' @returns data.frame with vraibels person_id and theta
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

  if(parms$model=='1PL')
  {
    rsp = get_resp_data(dataSrc, qtpredicate=qtpredicate, env=env)
    return(ability(rsp,method=method,parms=parms$items))

  } else if(parms$model=='2PL')
  {
    if(!inherits(dataSrc, 'matrix'))
      dataSrc = get_resp_matrix(dataSrc,qtpredicate,env=env) # a bit wasteful

    max_score = max(dataSrc,na.rm=TRUE)
    pre = lapply(mat_pre(dataSrc, max_score), drop)

    # to do: this needs protections and checks
    data_a = categorize(pre$inp, pre$pni, pre$icnp, pre$pcni,pre$ip, pre$pi,
                        parms$pre$icat, parms$pre$imax,max(parms$pre$ncat), pre$ix, pre$px)

    pid = rownames(dataSrc)
    if(is.null(pid))
      pid = 1:nrow(dataSrc)

    theta = theta_2pl(data_a, parms$em$A, parms$em$b, parms$pre$ncat,
                        pre$pni, pre$pcni, pre$pi, pre$px,
                        WLE=(method=='WLE'), USE_A = (!unweight))

        return(tibble(person_id=pid,theta=drop(theta)))
  }

}
