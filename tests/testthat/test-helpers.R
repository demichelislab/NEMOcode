test_that("norm works", {
  expect_equal(sum(norm_one(c(1,2,3,4,5))), 1)
})

test_that("clamp works", {
    clmp = .clamp_val(rnorm(100), 0, 1)
    expect_true(max(clmp)-min(clmp) <= 1)
})

test_that("beta sin works", {
    x = .5
    expect_true(x == invSi(siTr(x)))
})

test_that("test beta mval", {
    x = .5
    expect_true(x == .to_beta(.to_mval(x)))
})
