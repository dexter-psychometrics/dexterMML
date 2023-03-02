
#experimental
# probably needs a prior on the guessing
fit_1plG = function(dataSrc, predicate=NULL, group = NULL)
{
  env = caller_env()
  qtpredicate = eval(substitute(quote(predicate)))
  
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
    print(data$item_id[pre$imax < 1])
    stop('Some items have a maximum score of 0, model cannot be calibrated',call.=FALSE)
  }
  if(any(pre$icat[1,]==0))
  {
    cat('Items without a 0 score:\n')
    print(data$item_id[pre$icat[1,]==0])
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
  design = design_matrices(pre$pni, pre$pcni, pre$pi, data$persons$c_group_nbr, nit, nrow(data$groups))
  
  mu = rep(0, nrow(data$groups))
  sigma = rep(1, nrow(data$groups))
  
  
  stopifnot(all(pre$imax==1))
  
  a = categorize(pre$pni, pre$pcni, pre$icnp, pre$pi,
                 pre$icat, pre$imax,max(pre$ncat), pre$px, pre$ix)
  
  start = start_1pl(a, pre$ncat, pre$icatg, ref_group, data$item_id)
  
  ref_group = start$ref_group
  fixed_items = start$fixed_items
  check_connected(design, fixed_items)
  
  em = estimate_1plG(-start$b[2,],
                    pre$pni, pre$pcni, pre$pi, pre$px,
                    theta_grid, mu, sigma, data$groups$group_n, data$persons$c_group_nbr, 
                    #fixed_items, 
                    ref_group,
                    max_iter=max_em_iterations)
  
  em$LL = loglikelihood_1plG_GH(em$b, em$lg, pre$pni, pre$pcni, pre$pi, pre$px, 
                                    quadpoints$nodes, quadpoints$weights, em$mu, em$sigma, data$persons$c_group_nbr)
  
  
  list(em=em,items=tibble(item_id=data$item_id,item_score=1,beta=em$b[,1],g=exp(em$lg)/(1+exp(em$lg))))
}

sim_1plG = function(theta,items)
{
  stopifnot(all(items$guess>=0))
  dat = sim_1plGc(items$beta,items$guess,theta) 
  colnames(dat) = items$item_id
  dat
}



fit_1plAG = function(dataSrc, predicate=NULL, group = NULL)
{
  env = caller_env()
  qtpredicate = eval(substitute(quote(predicate)))
  
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
    print(data$item_id[pre$imax < 1])
    stop('Some items have a maximum score of 0, model cannot be calibrated',call.=FALSE)
  }
  if(any(pre$icat[1,]==0))
  {
    cat('Items without a 0 score:\n')
    print(data$item_id[pre$icat[1,]==0])
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
  design = design_matrices(pre$pni, pre$pcni, pre$pi, data$persons$c_group_nbr, nit, nrow(data$groups))
  
  mu = rep(0, nrow(data$groups))
  sigma = rep(1, nrow(data$groups))
  
  
  stopifnot(all(pre$imax==1))
  
  a = categorize(pre$pni, pre$pcni, pre$icnp, pre$pi,
                 pre$icat, pre$imax,max(pre$ncat), pre$px, pre$ix)
  
  start = start_1pl(a, pre$ncat, pre$icatg, ref_group, data$item_id)
  
  ref_group = start$ref_group
  fixed_items = start$fixed_items
  check_connected(design, fixed_items)
  
  em = estimate_1plAG(-start$b[2,],
                     pre$pni, pre$pcni, pre$pi, pre$px,
                     theta_grid, mu, sigma, data$groups$group_n, data$persons$c_group_nbr, 
                     #fixed_items, 
                     ref_group,
                     max_iter=max_em_iterations)
  
  em$LL = loglikelihood_1plAG_GH(em$g, em$alpha,em$b, pre$pni, pre$pcni, pre$pi, pre$px, 
                                quadpoints$nodes, quadpoints$weights, em$mu, em$sigma, data$persons$c_group_nbr)
  
  
  list(em=em,items=tibble(item_id=data$item_id,item_score=1,beta=em$b[,1],guess=em$g[,1],alpha=em$alpha))
}

sim_1plAG = function(theta,items)
{
  #stopifnot(all(items$guess>=0))
  dat = sim_1plAGc(items$guess, items$alpha,items$beta,theta) 
  colnames(dat) = items$item_id
  dat
}


