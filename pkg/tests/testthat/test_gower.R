
context("Basic distance elements")

dL <- expand.grid(c(TRUE,FALSE),c(TRUE,FALSE))
test_that("distance between logicals",{
  expect_equal(gower_dist(data.frame(x=dL[,1]),data.frame(x=dL[,2])),c(0,1,1,NaN))
})

bands <- c("Grand Magus","Skull Fist")
dF <- expand.grid(bands,bands)
test_that("distance between factor variables",{
  expect_equal(gower_dist(data.frame(x=dF[,1]),data.frame(x=dF[,2])),c(0,1,1,0))
})

dN <- data.frame(x = as.numeric(1:4),y=as.numeric(c(1,1,2,3)))
test_that("distance between numerical variables",{
  expect_equal(gower_dist(data.frame(x=dN[,1]),data.frame(x=dN[,2])),c(0,1/3,1/3,1/3))
})

test_that("multivariate dataset",{
  dM1 <- data.frame(x=dL[,1],y=dF[,1],z=dN[,1])  
  dM2 <- data.frame(x=dL[,2],y=dF[,2],z=dN[,2])
  expect_equal(gower_dist(x=dM1,y=dM2), c(0,7/9,7/9,1/6))
  # check symmetry
  expect_equal(gower_dist(dM1,dM2),gower_dist(dM2,dM1))
  # not counting NA's in the denominator
  dM1[array(c(2,3,4,1,2,3),dim=c(3,2))] <- NA
  expect_equal(gower_dist(dM1,dM2), c(0,3/4,3/4,0))
})




