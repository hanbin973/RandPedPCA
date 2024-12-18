



# test
test_that("rppca on Linv", {
  expect_no_condition(
    pc <- rppca(pedLInv)
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

test_that("Comparing STD values", {
  expect_no_condition(
    pc <- rppca(pedLInv)
  )

  ped <- pedigree(sire  = pedMeta$fid,
                  dam   = pedMeta$mid,
                  label = pedMeta$id)
  expect_no_error(pc2 <- rppca(ped))
  expect_true(all(pc$sdev[1:2] - pc2$sdev[1:2] < 1e-10))
})



test_that("rppca on Linv with totVar supplied", {
  ped <- pedigree(sire  = pedMeta$fid,
                  dam   = pedMeta$mid,
                  label = pedMeta$id)
  tv <- sum(inbreeding(ped)+1)
  expect_no_condition(pc3 <- rppca(pedLInv, totVar = tv))
  expect_no_condition(summary(pc3))
  expect_false(pc3$center)
  expect_false(pc3$scale)
  expect_failure(expect_null(pc3$varProps))
})















