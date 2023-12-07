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
  
  rsp = get_responses(db) %>%
    group_by(person_id) %>%
    slice_sample(n=20) %>%
    ungroup()
    
  theta.mml = ability.mml(rsp, f2, method='WLE')
  
  # also test possible issues with sorting of capital letters
  
  rsp = rsp %>% mutate(item_id = gsub('S2','s2',item_id))
  coef_f2 = coef(f2) %>% mutate(item_id = gsub('S2','s2',item_id))
  theta.mml_cf = ability.mml(rsp, coef_f2, method='WLE')
  
  tst = inner_join(theta.mml_cf, theta.mml, by='person_id')
  expect_lt(max(abs(tst$theta.x-tst$theta.y)), 1e-10, label='coef vs parms object theta')
  
  close_project(db)
  
})


# fixed items al goed getest met parametrisatie?
