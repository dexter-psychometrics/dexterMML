

ability.mml = function(dataSrc, parms, predicate=NULL, method=c('MLE','WLE'), unweight=FALSE)
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  method=match.arg(method)

  if(parms$model=='1PL')
  {
    rsp = get_resp_data(db, qtpredicate=qtpredicate, env=env)
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

    if(unweight)
    {
      theta = theta_2plu(data_a, parms$em$A, parms$em$b, parms$pre$ncat,
                      pre$pni, pre$pcni, pre$pi, pre$px, (method=='WLE'))
    } else
    {
      theta = theta_2pl(data_a, parms$em$A, parms$em$b, parms$pre$ncat,
                        pre$pni, pre$pcni, pre$pi, pre$px, (method=='WLE'))
    }
    return(tibble(person_id=pid,theta=theta))
  }

}
