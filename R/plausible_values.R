#' plausible values
#'
#' This function uses the Single Variable Exchange algorithm
#' (Murray, Ghahramani & MacKay, 2012; see Marsman, Maris, Bechger & Glas, 2014, for application to IRT models)
#' to produce a dependent sample of plausible values.
#'
#' @param dataSrc	a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param parms	object produced by function fit_marginal
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
  
  pre = abl_pre(dat, parms)

  # if WLE not entirely converged, no so much of a problem
  starting_values = theta_2pl(pre$a, pre$A, pre$b, pre$ncat,
                                      pre$pni, pre$pcni, pre$pi, pre$px,
                                      WLE=TRUE, USE_A = (pre$model!='1PL'))$theta

  pv =  plausible_values_c(pre$A, pre$a, pre$b, pre$ncat,
                         pre$pni, pre$pcni, pre$pi, pre$px, group, group_n,
                         as.integer(npv), starting_values,
                         n_prior_updates=70L, thin=70L,pgw = progress_width())

  colnames(pv) = sprintf("PV%i",1:ncol(pv))
  cbind(out, as.data.frame(pv))
}

#' Simulate data for a 2pl
#'
#' note: to simulate from a 1pl, you can omit the alpha column or set it to 1
#'
#' @param pars data.frame with columns alpha, item_score and beta
#' @param theta vector of abilities
#' 
#' @returns matrix with persons as rows and items as columns
#' 
sim_2pl = function(pars,theta)
{
  colnames(pars) = tolower(colnames(pars))
  if(!'alpha' %in% colnames(pars))
    pars$alpha=1
  
  pars = select(ungroup(pars),.data$item_id,.data$item_score,.data$alpha,.data$beta) %>%
    mutate(item_score=as.integer(.data$item_score)) %>%
    arrange(.data$item_id,.data$item_score)
  
  #check alpha equal?
  
  items = pars %>%
    group_by(.data$item_id) %>%
    summarise(ncat=n()+1,alpha=first(.data$alpha)) 
  
  ncat=items$ncat
  
  a=matrix(0L,max(items$ncat),nrow(items))
  b=matrix(0,max(items$ncat),nrow(items))
  
  pars = split(pars,pars$item_id)
  for(i in seq_along(pars))
  {
    a[2:ncat[i],i]=pars[[i]]$item_score
    b[2:ncat[i],i]=pars[[i]]$beta
  }
  
  dat = sim_2plc(a, items$alpha, b, ncat,theta)
  colnames(dat) = items$item_id
  dat
}



