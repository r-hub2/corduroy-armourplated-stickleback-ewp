test_that("null model runthrough", {
  fit_null <- ewp_reg(eggs ~ 1, data = linnet)
  expect_s3_class(fit_null, "ewp")
  print(fit_null)
  summary(fit_null)
  #get estimates and compare against published estimates from Ridout & Besbeas (2004)
  coefs <- unname(coef(fit_null))
  expect_equal(exp(coefs[1]),4.88,tolerance=0.017)
  expect_equal(coefs[2],1.46,tolerance=0.056)
  expect_equal(coefs[3],2.36,tolerance=0.056)
  #test that predict and fitted are non-empty
  expect_false(any(is.na(fit_null$fitted)))
  expect_equal(length(fit_null$fitted), fit_null$n)
  expect_false(any(is.na(predict(fit_null))))
  expect_equal(length(predict(fit_null)), fit_null$n)
  })

test_that("NA handling and dropped levels", {
  set.seed(1234)
  linnet2 <- linnet
  linnet2$eggs[sample(nrow(linnet2), size = 100)] <- NA
  linnet2$arbfac <- factor(sample(letters[1:10], size = nrow(linnet2), replace = TRUE))
  linnet2 <- subset(linnet2, arbfac %in% letters[1:5], drop = FALSE)
  fit_na <- ewp_reg(eggs ~ arbfac, data = linnet2)
  expect_s3_class(fit_na, "ewp")
})

test_that("in-formula transforms", {
  set.seed(1234)
  fit_na <- ewp_reg(eggs ~ scale(cov2), data = linnet)
  expect_s3_class(fit_na, "ewp")
})

test_that("simulation method", {
  fite <- ewp_reg(eggs ~ cov1 + cov2, data = linnet)
  sime <- simulate(fite, nsim = 1)
  expect_s3_class(sime, "data.frame")
  expect_equal(ncol(sime), 1)
  sime3 <- simulate(fite, nsim = 3)
  expect_equal(ncol(sime3), 3)
})
