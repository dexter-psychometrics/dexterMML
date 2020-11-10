

#' palusible values
#'
#' @param dataSrc	a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param parms	object produced by function est
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param covariates character vector indicating discrete person properties or, in case dataSrc is a matrix,
#' a vector of groups
#' @param npv number of plausible values to draw per person
#'
#' @returns data.frame
#'
plausible_values.mml = function(dataSrc, parms, predicate=NULL, covariates=NULL, npv=1)
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  extra_columns = NULL
  if(!is.matrix(dataSrc))
    extra_columns = covariates

  respdat = get_resp_data(dataSrc,qtpredicate, env=env, extra_columns=extra_columns,summarised=TRUE,
                      parms_check=parms$items[,c('item_id','item_score')])

  rsp = respdat$x
  design = respdat$design
  rm(respdat)
  nbk = nlevels(design$booklet_id)

  if(is.matrix(dataSrc) && !is.null(covariates))
  {
    extra_columns = 'population'
    rsp$population = covariates # to do: need join probably
  }

  if(is.null(extra_columns))
  {
    rsp$pop__ = 0L
    npop = 1L
  } else
  {
    ii=-1L
    ind=function(){ii <<- ii+1L; ii}

    rsp = group_by(rsp, across(extra_columns)) %>%
      mutate(pop__ = ind()) %>%
      ungroup()

    npop = ii + 1L
  }

  rsp = arrange(rsp, .data$booklet_id, .data$pop__, .data$booklet_score)

  scoretab = count(rsp, .data$booklet_id, .data$pop__, .data$booklet_score)

  scoretab$booklet_id = as.integer(scoretab$booklet_id) - 1L
  design$booklet_id = as.integer(design$booklet_id) - 1L
  design$item_id = as.integer(design$item_id) - 1L
  bnit = count(design, .data$booklet_id)$n


  pre = pre_scoretab(scoretab$booklet_id, scoretab$pop__, scoretab$booklet_score, scoretab$n,
                     design$booklet_id, design$item_id, parms$pre$a, parms$pre$ncat, nbk, npop)


  A = if(parms$model == '1PL') rep(1,ncol(parms$em$a)) else parms$em$A


  pv = plausible_values(pre$booklet_id, pre$pop, pre$pbn, cumsum(c(0L,pre$pbn)),
                        pre$pbnp,cumsum(c(0L,pre$pbnp)),
                        pre$scoretab, pre$popn,
                        design$item_id, bnit, cumsum(c(0L,bnit)),
                        A, parms$pre$a, parms$em$b, parms$pre$ncat,
                        as.integer(npv))

  colnames(pv) = sprintf("PV%i",1:ncol(pv))
  cbind(rsp[,c('person_id',extra_columns)], as.data.frame(pv))
}
