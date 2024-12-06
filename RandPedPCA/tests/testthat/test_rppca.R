



# test
test_that("rppca on Linv", {
  expect_no_condition(
    pc <- rppca(pedLinv)
  )

  expect_warning(summary(pc))
  expect_false(pc$center)
  expect_false(pc$scale)
  expect_null(pc$varProps)

})

test_that("rppca on pedigree", {
  ped <- pedigree(sire  = pedMeta$fid,
                  dam   = pedMeta$mid,
                  label = pedMeta$id)
  expect_no_error(pc2 <- rppca(ped)) # this annoyingly throws a warning with an unhelpful suggestion
  expect_no_condition(summary(pc2))
  expect_false(pc2$center)
  expect_false(pc2$scale)
  expect_failure(expect_null(pc2$varProps))
})


test_that("rppca on Linv with totVar supplied", {
  expect_no_condition(pc3 <- rppca(pedLinv, totVar = 3549.629))
  expect_no_condition(summary(pc3))
  expect_false(pc3$center)
  expect_false(pc3$scale)
  expect_failure(expect_null(pc3$varProps))
})















