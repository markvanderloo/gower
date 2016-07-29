
context("Basic distance elements")

test_that("distance between logicals",{
  dL <- expand.grid(c(TRUE,FALSE),c(TRUE,FALSE))
  x = data.frame(x=dL[,1])
  y = data.frame(x=dL[,2])
  expect_equal(gower_dist(x = x,y = y),c(0,1,1,NaN))
})

test_that("distance between factor variables",{
  bands <- c("Grand Magus","Skull Fist")
  dF <- expand.grid(bands,bands)
  expect_equal(gower_dist(data.frame(x=dF[,1]),data.frame(x=dF[,2])),c(0,1,1,0))
})

test_that("distance between numerical variables",{
  dN <- data.frame(x = as.numeric(1:4),y=as.numeric(c(1,1,2,3)))
  expect_equal(gower_dist(data.frame(x=dN[,1]),data.frame(x=dN[,2])),c(0,1/3,1/3,1/3))
})

test_that("distance between character variables",{
  dC <- data.frame(x=letters[1:3],y=letters[3:1],stringsAsFactors=FALSE)
  expect_equal(gower_dist( data.frame(x=dC[,1]), data.frame(x=dC[,2])),c(1,0,1))

})

test_that("multivariate dataset",{
  bands <- c("Grand Magus","Skull Fist")
  dL <- expand.grid(c(TRUE,FALSE),c(TRUE,FALSE))
  dN <- data.frame(x = as.numeric(1:4),y=as.numeric(c(1,1,2,3)))
  dF <- expand.grid(bands,bands)
  dM1 <- data.frame(x=dL[,1],y=dF[,1],z=dN[,1])  
  dM2 <- data.frame(x=dL[,2],y=dF[,2],z=dN[,2])
  expect_equal(gower_dist(x=dM1,y=dM2), c(0,7/9,7/9,1/6))
  # check symmetry
  expect_equal(gower_dist(dM1,dM2),gower_dist(dM2,dM1))
  # not counting NA's in the denominator
  dM1[array(c(2,3,4,1,2,3),dim=c(3,2))] <- NA
  expect_equal(gower_dist(dM1,dM2), c(0,3/4,3/4,0))
})

test_that("recycling",{
  expect_equal(length(gower_dist(x=iris[1,],y=iris)), nrow(iris))
  expect_equal(length(gower_dist(x=iris,y=iris[1,])), nrow(iris))
  expect_equal(length(gower_dist(x=iris[1:3,],y=iris)), nrow(iris))
  expect_equal(length(gower_dist(x=iris,y=iris[1:3,])), nrow(iris))
})


test_that("exceptions",{
  expect_warning(gower_dist(
    x = data.frame(x=c(1.2,1.2,1.2))
    , y = data.frame(x=c(1.2,1.2,1.2))
  ))
  expect_warning(gower_dist(
    x = data.frame(x=c(1.2,1.2,1.2))
    , y = data.frame(x=c(1.2,1.2,1.3))
    , eps=0.2
  ))

  # should not be a warning, but the case with one columns NA is interesting.  
  # expect_warning(gower_dist(
  #   x = data.frame(x=rep(NA_integer_,3))
  #   , y = data.frame(x=c(NaN,-Inf,Inf))
  #   , eps=0.2
  # ))
  


  })
