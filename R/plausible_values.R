

#' plausible values
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
plausible_values.mml_old = function(dataSrc, parms, predicate=NULL, covariates=NULL, npv=1)
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


  pv = plausible_values_c(pre$booklet_id, pre$pop, pre$pbn, cumsum(c(0L,pre$pbn)),
                        pre$pbnp,cumsum(c(0L,pre$pbnp)),
                        pre$scoretab, pre$popn,
                        design$item_id, bnit, cumsum(c(0L,bnit)),
                        A, parms$pre$a, parms$em$b, parms$pre$ncat,
                        as.integer(npv))

  colnames(pv) = sprintf("PV%i",1:ncol(pv))
  cbind(rsp[,c('person_id',extra_columns)], as.data.frame(pv))
}

#' plausible values
#'
#'This function uses the SVE algorithm (Marsman,...;) to produce a dependent sample of plausible values
#'
#' @param dataSrc	a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param parms	object produced by function est
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param covariates character vector indicating discrete person properties or, in case dataSrc is a matrix,
#' a vector of groups (integer or character)
#' @param npv number of plausible values to draw per person
#'
#' @returns data.frame
#'
plausible_values.mml = function(dataSrc, parms, predicate=NULL, covariates=NULL, npv=1)
{
  # this part copied form model est, should become function, to do: clean up
  group=covariates
  if(inherits(dataSrc, 'matrix'))
  {
    dat = dataSrc
    mode(dat) = 'integer'
    if(is.null(colnames(dat)))
      colnames(dat) = sprintf("item%06i",1:ncol(dat))

    if(!is.null(group) && length(group) != nrow(dat))
      stop(sprintf("Length of group (%i) is not equal to number of rows in data (%i)",
                   length(group),nrow(dat)))

    person_id = if(is.null(rownames(dat))) sprintf("p%09i",1:nrow(dat)) else rownames(dat)
    out=tibble(person_id=person_id)
    if(!is.null(group))
      out$group=group

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

      out = g

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

  max_score = max(dat,na.rm=TRUE)

  pre = lapply(mat_pre(dat, max_score), drop)

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


  A = if(parms$model == '1PL') rep(1,ncol(parms$em$a)) else parms$em$A

  data_a = categorize(pre$inp, pre$pni, pre$icnp, pre$pcni,pre$ip, pre$pi,
                      parms$pre$icat, parms$pre$imax,max(parms$pre$ncat), pre$ix, pre$px)

  starting_values = theta_2pl(data_a, A, parms$em$b, parms$pre$ncat,
                                      pre$pni, pre$pcni, pre$pi, pre$px,
                                      WLE=TRUE, USE_A = TRUE)

  pv =  plausible_values_c2(A, data_a, parms$em$b, parms$pre$ncat,
                         pre$pni, pre$pcni, pre$pi, pre$px, group, group_n,
                         as.integer(npv), starting_values,
                         n_prior_updates=70L, thin=70L)

  colnames(pv) = sprintf("PV%i",1:ncol(pv))
  cbind(out, as.data.frame(pv))
}
