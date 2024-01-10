
progress_width = function()
{
  if(!getOption("knitr.in.progress",FALSE) && getOption("dexter.progress", TRUE) && interactive())
  {
    getOption("width")
  } else {-1L}
}

num_hess = function(e,delta=1e-5)
{
  pre = e$pre
  if(e$model=='1PL')
  {
    num_hessian_nrm(e$em$a, e$em$b, pre$ncat, pre$pcni, pre$pi, pre$px,
                    e$theta_grid, e$em$mu, e$em$sigma, pre$group, delta)
  }
  else
  {
    num_hessian_2pl(e$em$a, e$em$A, e$em$b, pre$ncat,
                           pre$pni, pre$pcni, pre$pi, pre$px,
                           e$theta_grid, e$em$mu, e$em$sigma, pre$group, delta)
  }
}



# 
# grd_2pl = function(e)
# {
#   pre = e$pre
#   gradient_2pl(e$em$a, e$em$A, e$em$b, pre$ncat,pre$pcni, pre$pi, pre$px,
#                e$theta_grid, e$em$mu, e$em$sigma, pre$group, 
#                pre$ip, pre$inp, pre$icnp)
# }


# library(statmod)
# GH = gauss.quad.prob(160,'normal',mu=0,sigma=1)
# quadpoints = tibble(nodes=GH$nodes,weights=GH$weights) |>
#   filter(weights>1e-22) |>
#   arrange(weights) |>
#   mutate(weights=weights/sum(weights))
# 
# # quadpoints
# usethis::use_data(quadpoints, internal = TRUE)