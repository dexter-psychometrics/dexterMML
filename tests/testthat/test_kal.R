context('calibration checks')
library(dexter)

# to do: save sim data to remove chance component

test_that('1pl verb agg confirms to dexter',{
  db = start_new_project(verbAggrRules, ":memory:")
  add_booklet(db, verbAggrData, "agg")
  f=fit_enorm(db,item_id!='S3DoShout')
  p=coef(f)
  
  e = fit_1pl(db,item_id!='S3DoShout', fixed_param=p[1:2,])
  expect_lt(mean(abs(p$beta-coef(e)$beta)),.01)

  expect_gt(cor(p$SE_beta,coef(e)$SE_beta,use='complete.obs'),.99)

  if('standard_errors' %in% names(formals(dexter::ability)))
  {
    a.cml=ability(db,f,item_id!='S3DoShout',method='WLE',standard_errors=TRUE)
  } else
  {
    a.cml=ability(db,f,item_id!='S3DoShout',method='WLE')
  }
  
  a.mml=ability.mml(db,e,item_id!='S3DoShout',method='WLE')
  
  expect_lt(mean(abs(a.cml$theta-a.mml$theta)),.01)
  expect_lt(mean(abs(a.cml$se-a.mml$se)),.001)
  close_project(db)
})


test_that('2pl with simdata and missing cats',{
  set.seed(1234) 
  theta=c(rnorm(1000),rnorm(500,1,2),rnorm(500,-1,.5))
  group=rep(letters[1:3],c(1000,500,500))
  db = start_new_project(verbAggrRules, ":memory:")
  add_booklet(db, verbAggrData, "agg")
  f=fit_enorm(db)
  p=coef(f)
  close_project(db)
  p$item_score[c(2,4,6)] = p$item_score[c(2,4,6)] + 2L
  p$item_score[7:8] = p$item_score[7:8] *2L
  p$item_score[9:12] = p$item_score[9:12] + 1L
  p$alpha = rep(runif(24,.5,1.5),each=2)
  
  dat = sim_2pl(p,theta)
  dat[1:1000,seq(1,ncol(dat),2)] = NA_integer_
  dat[1001:1500,seq(2,ncol(dat),2)] = NA_integer_

  e=fit_2pl(dat,group=group)
  
  cf = inner_join(coef(e), p, by=c('item_id','item_score'), suffix=c('','.true')) %>%
    mutate(df_alpha = (alpha-alpha.true),
           df_beta = (beta-beta.true)) 
  
  expect_gt(cor(cf$alpha,cf$alpha.true),.9)
  expect_true(between(mean(cf$alpha-cf$alpha.true),-.1,.1))
  
  expect_gt(cor(cf$beta,cf$beta.true),.95)
  expect_true(between(mean(cf$beta-cf$beta.true),-.05,.05))

  abl=ability.mml(dat,e,method='WLE')
  expect_gt(cor(abl$theta,theta),.9)
 
})

test_that('2pl with fixed_parameters',{
  set.seed(1234)
  db = start_new_project(verbAggrRules, ":memory:")
  add_booklet(db, verbAggrData, "agg")
  f=fit_enorm(db)
  p=coef(f)
  close_project(db)
  p$alpha = rep(runif(24,.5,1.5),each=2)
  theta = rnorm(5000,-.5,2)

  dat = sim_2pl(p,theta)
  
  e = fit_2pl(dat, fixed_param=p[1:6,],se=FALSE)
  
  pop=coef(e,what='population')
  expect_lt(abs(pop$mean-mean(theta)),.05)
  expect_lt(abs(pop$sd-sd(theta)),.05)
  i=coef(e,what='items')
  
  expect_lt(mean(abs(i$beta-p$beta)),.05)
  expect_lt(mean(abs(i$alpha-p$alpha)),.05)

})