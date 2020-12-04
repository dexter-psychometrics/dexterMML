#' plausible values
#'
#' This function uses the Single Variable Exchange algorithm
#' (Murray, Ghahramani & MacKay, 2012; see Marsman, Maris, Bechger & Glas, 2014, for application to IRT models)
#' to produce a dependent sample of plausible values.
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
#' @details
#' This function uses the Single Variable Exchange algorithm
#' (Murray, Ghahramani & MacKay,2012; see Marsman, Maris, Bechger & Glas, 2014, for application to IRT models)
#' to produce a dependent sample of plausible values. The settings are such that the autocorrelation is reasonably low
#' for tests with up to 40 items. To further lower the autocorrelation, e.g. for longer tests, the user can
#' simply draw more plausible values than needed (see argument npv) and apply additional thinning.
#'
#' @references
#' Marsman, M., Maris, G., Bechger, T. M., and Glas, C.A.C. (2014) Composition algorithms for Conditional Distributions.
#' PhD thesis
#' Murray, I., Ghahramani, Z., & MacKay, D. (2012). MCMC for doubly-intractable distributions. via: https://arxiv.org/abs/1206.6848
#'
plausible_values.mml = function(dataSrc, parms, predicate=NULL, covariates=NULL, npv=1)
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()

  data = get_mml_data(dataSrc,qtpredicate,env,covariates)

  dat = data$dat
  group = data$group
  out = data$persons
  group=covariates

  max_score = max(dat,na.rm=TRUE)

  pre = lapply(mat_pre(dat, max_score), drop)

  if(is.null(group))
  {
    group = integer(nrow(dat))
    group_n = nrow(dat)
  } else
  {
    group = as.factor(group)
    group = as.integer(group) -1L
    group_n = as.integer(table(group))
  }

  if(parms$model=='1PL')
  {
    A=rep(1,ncol(parms$em$b))
    b=-parms$em$b/parms$em$a
  } else
  {
    A=parms$em$A
    b=parms$em$b
  }

  data_a = categorize(pre$pni, pre$pcni, pre$pi,
                      parms$pre$icat, parms$pre$imax,max(parms$pre$ncat), pre$px)

  starting_values = theta_2pl(data_a, A, b, parms$pre$ncat,
                                      pre$pni, pre$pcni, pre$pi, pre$px,
                                      WLE=TRUE, USE_A = TRUE)

  pv =  plausible_values_c(A, data_a, b, parms$pre$ncat,
                         pre$pni, pre$pcni, pre$pi, pre$px, group, group_n,
                         as.integer(npv), starting_values,
                         n_prior_updates=70L, thin=70L,pgw = getOption("width"))

  colnames(pv) = sprintf("PV%i",1:ncol(pv))
  cbind(out, as.data.frame(pv))
}
