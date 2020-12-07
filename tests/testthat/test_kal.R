context('calibration checks')
library(dexter)

# to do: save sim data to remove chance component

test_that('1pl verb agg confirms to dexter',{
  db = start_new_project(verbAggrRules, ":memory:")
  add_booklet(db, verbAggrData, "agg")
  f=fit_enorm(db)
  p=coef(f)
  
  e = est(db, fixed_param=p[1:2,])
  expect_lt(mean(abs(p$beta-coef(e)$beta)),.02)

  expect_lt(mean(abs(p$SE_beta-coef(e)$SE_beta)[-(1:2)]),.01)
  
  a.cml=ability(db,f,method='WLE')
  a.mml=ability.mml(db,e,method='WLE')
  
  expect_lt(mean(abs(a.cml$theta-a.mml$theta)),.01)
  close_project(db)
})


test_that('2pl with simdata and missing cats',{
  theta=c(rnorm(1000),rnorm(500,1,2),rnorm(500,-1,.5))
  group=rep(letters[1:3],c(1000,500,500))
  db = start_new_project(verbAggrRules, ":memory:")
  add_booklet(db, verbAggrData, "agg")
  f=fit_enorm(db)
  p=coef(f)
  p$item_score[c(2,4,6)] = p$item_score[c(2,4,6)] + 2L
  p$item_score[7:8] = p$item_score[7:8] *2L
  p$item_score[9:12] = p$item_score[9:12] + 1L
  p$alpha = rep(runif(24,.5,1.5),each=2)
  
  dat = sim_2pl(p,theta)
  #dat[1:1000,seq(1,ncol(dat),2)] = NA_integer_
  #dat[1001:1500,seq(2,ncol(dat),2)] = NA_integer_

  e=est(dat,model='2PL',group=group,se=FALSE)
  
  expect_lt(max(abs(coef(e)$beta-p$beta)),.4)
  expect_lt(mean(abs(coef(e)$beta-p$beta)),.2)

  expect_lt(max(abs(coef(e)$alpha-p$alpha)),.3)
  expect_lt(mean(abs(coef(e)$alpha-p$alpha)),.15)

  abl=ability.mml(dat,e,method='WLE')
  expect_gt(cor(abl$theta,theta),.9)
 
})

test_that('2pl with fixed_parameters',{
  db = start_new_project(verbAggrRules, ":memory:")
  add_booklet(db, verbAggrData, "agg")
  f=fit_enorm(db)
  p=coef(f)
  p$alpha = rep(runif(24,.5,1.5),each=2)
  theta = rnorm(5000,-.5,2)

  dat = sim_2pl(p,theta)
  
  e = est(dat, fixed_param=p[1:6,],model='2PL',se=FALSE)
  
  pop=coef(e,what='population')
  expect_lt(abs(pop$mu-mean(theta)),.05)
  expect_lt(abs(pop$sd-sd(theta)),.1)
  i=coef(e,what='items')
  
  expect_lt(mean(abs(i$beta-p$beta)),.05)
  expect_lt(mean(abs(i$alpha-p$alpha)),.05)

})