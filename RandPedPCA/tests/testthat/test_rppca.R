
pc <- rppca(pedLinv)

# test
test_that("center is always FALSE", {
  expect_false(pc$center)
})

test_that("scale is always FALSE", {
  expect_false(pc$scale)
})

test_that("vProp is not set by default", {
  expect_null(pc$vProp)
})



