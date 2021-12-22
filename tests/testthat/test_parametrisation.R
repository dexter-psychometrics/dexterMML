context('parametrisation issues')
library(dexter)
library(dplyr)

test_that('1pl verb agg confirms to dexter',{
  
  db = start_new_project(verbAggrRules, ":memory:")
  add_booklet(db, verbAggrData, "agg")
  
  touch_rules(db,tibble(item_id='S1DoScold',response=c(1,2), item_score=c(2,5)))
  
  f = fit_enorm(db)
  f1 = fit_1pl(db,se=FALSE)  
  
  expect_lt(mean(abs(coef(f)$beta + mean(coef(f1)$beta) - coef(f1)$beta)), 0.03, label='mml vs dexter 1pl calibration')
  
  # try dexter to be sure
  d1 = ability_tables(f, method='WLE')
  d2 = ability_tables(coef(f), method='WLE')
  
  expect_lt(max(abs(d1$theta - d2$theta)), 1e-10, label='double check dexter coef')
  
  # use mml pars for abl checks
  theta.dx = ability(db, coef(f1), method='WLE')
  theta.mml = ability.mml(db, f1, method='WLE')
  theta.mml_cf = ability.mml(db, coef(f1), method='WLE')
  
  tst1 = inner_join(theta.dx, theta.mml, by='person_id')
  expect_lt(max(abs(tst1$theta.x-tst1$theta.y)), 1e-10, label=' dexter vs mml theta')
  
  tst2 = inner_join(theta.mml, theta.mml_cf, by='person_id')
  expect_lt(max(abs(tst2$theta.x-tst2$theta.y)), 1e-10, label='theta mml object vs coef')
 

  # we would possibly test 2pl theta against some other program, this is just internal consistency
  
  f2 = fit_2pl(db,se=FALSE)
  theta.mml = ability.mml(db, f2, method='WLE')
  theta.mml_cf = ability.mml(db, coef(f2), method='WLE')
  
  tst = inner_join(theta.mml_cf, theta.mml, by='person_id')
  expect_lt(max(abs(tst$theta.x-tst$theta.y)), 1e-10, label=' dexter vs mml theta')
  
})


# fixed items al goed getest met parametrisatie?
