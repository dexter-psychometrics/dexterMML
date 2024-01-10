#' plausible values
#'
#' This function uses the Single Variable Exchange algorithm
#' (Murray, Ghahramani & MacKay, 2012; see Marsman, Maris, Bechger & Glas, 2014, for application to IRT models)
#' to produce a dependent sample of plausible values.
#'
#' @param dataSrc	a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param parms	object produced by function fit_1pl or fit_2pl or possibly a data.frame of parameters
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
#' When using a data.frame of parameters, be sure that you use the correct parametrisation. See \code{\link{fit_1pl}} for details. 
#'
#' @references
#' Marsman, M., Maris, G., Bechger, T. M., and Glas, C.A.C. (2014) Composition algorithms for Conditional Distributions.
#' PhD thesis
#' 
#' Murray, I., Ghahramani, Z., & MacKay, D. (2012). MCMC for doubly-intractable distributions. via: https://arxiv.org/abs/1206.6848
#'
plausible_values.mml = function(dataSrc, parms, predicate=NULL, covariates=NULL, npv=1)
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  
  pre = abl_pre(dataSrc, parms, qtpredicate=qtpredicate, env=env, group=covariates)
  
  
  starting_values = rnorm(nrow(pre$persons),0,1)
  
  pv =  plausible_values_c(pre$A, pre$a, pre$b, pre$ncat,
                           pre$pni, pre$pcni, pre$pi, pre$px, pre$persons$c_group_nbr, pre$group$group_n,
                           as.integer(npv), starting_values,
                           n_prior_updates=80L, thin=60L,pgw = progress_width())
  
  colnames(pv) = sprintf("PV%i",1:ncol(pv))
  cbind(select(pre$persons,-"c_group_nbr"), as.data.frame(pv))
}


#' Simulate data for a 2pl
#'
#' note: to simulate from a 1pl, you can omit the alpha column or set it to 1
#'
#' @param pars data.frame with columns item_id, alpha, item_score and beta
#' @param theta vector of abilities
#' 
#' @returns matrix with persons as rows and items as columns
#' 
sim_2pl = function(pars,theta)
{
  colnames(pars) = tolower(colnames(pars))
  if(!'alpha' %in% colnames(pars))
    pars$alpha=1
  
  pars = select(ungroup(pars),"item_id","item_score","alpha","beta") |>
    mutate(item_score=as.integer(.data$item_score),
           item_id=factor(as.character(.data$item_id))) |>
    arrange(.data$item_id,.data$item_score) |>
    mutate(index=dense_rank(.data$item_id))
  
  #check alpha equal?
  
  items = pars |>
    group_by(.data$item_id) |>
    summarise(ncat=n()+1,alpha=first(.data$alpha)) |>
    ungroup() |>
    arrange(.data$item_id)
  
  ncat=items$ncat
  
  a=matrix(0L,max(items$ncat),nrow(items))
  b=matrix(0,max(items$ncat),nrow(items))
  
  for(itm in split(pars,pars$item_id))
  {
    i=itm$index[1]
    a[2:ncat[i],i]=itm$item_score
    b[2:ncat[i],i]=itm$beta
  }
  
  dat = sim_2plc(a, items$alpha, b, ncat,theta)
  colnames(dat) = as.character(items$item_id)
  dat
}


