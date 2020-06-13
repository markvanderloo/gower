
if (at_home()){
  dat1 <- iris[1:10,]
  dat2 <- iris[6:15,]

  nthrd1 <- gower_dist(dat1, dat2, nthread=1)
  nthrd2 <- gower_dist(dat1, dat2, nthread=2)
  nthrd4 <- gower_dist(dat1, dat2, nthread=4)
  nthrd7 <- gower_dist(dat1, dat2, nthread=7)

  expect_equal(nthrd1, nthrd2)
  expect_equal(nthrd1, nthrd4)
  expect_equal(nthrd1, nthrd7)

}



