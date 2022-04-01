context('design checks')

library(dplyr)

test_that('Cpp design check functions correctly',{
  nit=50
  items = matrix(1L,nit,nit)
  items[41:nit,1:10] = 0L

  items_fixed = integer(nit)
  groups = matrix(1L,1,nit)

  expect_true(dexterMML:::check_connected_c(items,groups,items_fixed)==1L,label='connected==1')

  items[] = 0L
  items[1:15,1:15] = 1L
  items[16:nit,16:nit] = 1L
  expect_true(dexterMML:::check_connected_c(items,groups,items_fixed)==4L,label='connected via groups')
  groups = matrix(1L,2,nit)
  expect_true(dexterMML:::check_connected_c(items,groups,items_fixed)==4L,label='connected via groups')
  groups[1,1:15]=0L
  groups[2,16:nit]=0L
  expect_true(dexterMML:::check_connected_c(items,groups,items_fixed)==0L,label='not connected')

  #now it gets convoluted
  items[] = 0L
  items[1:10,1:10] = 1L
  items[11:20,11:20] = 1L
  items[21:30,21:30] = 1L
  items[31:nit,31:nit] = 1L
  groups[]=0L
  groups[1,1:20]=1L
  groups[2,21:nit]=1L
  expect_true(dexterMML:::check_connected_c(items,groups,items_fixed)==0L,label='not connected')
  items_fixed[c(1,30)]=1L
  expect_true(dexterMML:::check_connected_c(items,groups,items_fixed)==8L,label='tenuous==8')

  groups = matrix(1L,1,nit)
  expect_true(dexterMML:::check_connected_c(items,groups,items_fixed)==4L,label='connected by fixed==4')

})

test_that('design matrix is produced correctly',{
  
  # 49 persons, 50 items
  dat = tibble(person_id = rep(1:49,each=2), item_id=rep(1:50,each=2)[-c(1,100)], item_score=sample(0:3,98,replace=TRUE),
               g = rep(c('a','b'),c(50,48)))

  data = dexterMML:::mml_pre(dat,NULL,NULL,group='g')
  
  pre = data$pre
  nit = length(data$item_id)
  
  ds = dexterMML:::design_matrices(pre$pni, pre$pcni, pre$pi, data$persons$c_group_nbr, nit, nrow(data$groups))
  
  # ds$items should be a block design
  # diagonal and off diagonal 1 should be 1, the rest 0
  expect_true(all(diag(ds$items[-1,]) == 1) && all(diag(ds$items[,-1]) == 1) && all(diag(ds$items) == 1), 
              label='block diagonal design == 1')
  
  expect_equal(sum(ds$items),sum(diag(ds$items[-1,]))  + sum(diag(ds$items[,-1])) + sum(diag(ds$items)),
               label='off block diagonal design == 0')

  expect_true(all(ds$groups[1,unique(filter(dat,g=='a')$item_id)] ==1) && all(ds$groups[2,unique(filter(dat,g=='b')$item_id)] ==1), 
              label='design groups 1')

  expect_equal(sum(ds$groups),
               sum(ds$groups[1,unique(filter(dat,g=='a')$item_id)]) + sum(ds$groups[2,unique(filter(dat,g=='b')$item_id)]),
               label='design groups 0')
  
})
