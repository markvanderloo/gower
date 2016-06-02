
context("Basic distance elements")

test_that("distance between logicals",{
  d <- expand.grid(c(TRUE,FALSE),c(TRUE,FALSE))
  expect_equal(gower_dist(d[1],d[2]),c(0,1,1,NaN))
})

test_that("distance between factor variables",{
  bands <- c("Grand Magus","Skull Fist","Cathedral")
  d <- expand.grid(bands,bands[1:2])
  expect_equal(gower_dist(d[1],d[2]),c(0,1,1,1,0,1))
})

test_that("distance between numerical variables",{
  d <- data.frame(x = as.numeric(1:3),y=as.numeric(c(1,1,5)))
  expect_equal(gower_dist(d[1],d[2]),c(0,0.25,0.5))
})
