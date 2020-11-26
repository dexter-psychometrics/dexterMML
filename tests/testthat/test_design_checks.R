context('design checks')


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
