



# test

# Lin not centred

test_that("rppca on Linv", {
  expect_no_condition(
    pc <- rppca(pedLInv)
  )

  expect_warning(summary(pc))
  expect_false(pc$center)
  expect_false(pc$scale)
  expect_null(pc$varProps)

})

# Lin not centred

test_that("rppca on Linv with totVar", {
  expect_no_condition(
    pc <- rppca(pedLInv, totVar = 3521.534)
  )

  expect_no_condition(summary(pc))
  expect_true(dim(summary(pc)$importance)[1] == 3) # three rows
  expect_false(pc$center)
  expect_false(pc$scale)
  expect_vector(pc$varProps)

})

# Linv  centred
test_that("rppca on Linv centred", {
  expect_no_condition(
    pc <- rppca(pedLInv, center=T)
  )
  expect_warning(summary(pc))
  expect_true(pc$center)
  expect_false(pc$scale)
  expect_null(pc$varProps)

})

# Linv  centred
test_that("rppca on Linv centred with totVar", {
  expect_no_condition(
    pc <- rppca(pedLInv, center=T, totVar=2694.038)
  )
  expect_no_condition(summary(pc))
  expect_true(dim(summary(pc)$importance)[1] == 3) # three rows
  expect_true(pc$center)
  expect_false(pc$scale)
  expect_vector(pc$varProps)

})



test_that("rppca on pedigree", {
  ped <- pedigree(sire  = pedMeta$fid,
                  dam   = pedMeta$mid,
                  label = pedMeta$id)
  expect_no_error(pc2 <- rppca(ped))
  expect_no_condition(summary(pc2))
  expect_true(dim(summary(pc2)$importance)[1] == 3) # three rows
  expect_false(pc2$center)
  expect_false(pc2$scale)
  expect_vector(pc2$varProps)
})

test_that("rppca on pedigree with (redundant) totVar", {
  ped <- pedigree(sire  = pedMeta$fid,
                  dam   = pedMeta$mid,
                  label = pedMeta$id)
  expect_warning(pc2 <- rppca(ped, center=F, totVar=123))
  expect_true(dim(summary(pc2)$importance)[1] == 3) # three rows
  expect_no_condition(summary(pc2))
  expect_false(pc2$center)
  expect_false(pc2$scale)
  expect_vector(pc2$varProps)
})


test_that("rppca on pedigree, centred", {
  ped <- pedigree(sire  = pedMeta$fid,
                  dam   = pedMeta$mid,
                  label = pedMeta$id)
  expect_no_error(pc2 <- rppca(ped, center=T))
  expect_warning(summary(pc2)) # total variance unknown
  expect_true(pc2$center)
  expect_false(pc2$scale)
  expect_null(pc2$varProps)
})

test_that("rppca on pedigree, centered with totVar", {
  ped <- pedigree(sire  = pedMeta$fid,
                  dam   = pedMeta$mid,
                  label = pedMeta$id)
  expect_no_condition(pc2 <- rppca(ped, center=T, totVar=2694.038)) # this annoyingly throws a warning with an unhelpful suggestion
  expect_no_condition(summary(pc2))
  expect_true(dim(summary(pc2)$importance)[1] == 3) # three rows
  expect_true(pc2$center)
  expect_false(pc2$scale)
  expect_vector(pc2$varProps)
})


test_that("Comparing STD values between rppca on pedigree and L^-1 input", {
  expect_no_condition(
    pc <- rppca(pedLInv)
  )

  ped <- pedigree(sire  = pedMeta$fid,
                  dam   = pedMeta$mid,
                  label = pedMeta$id)
  expect_no_error(pc2 <- rppca(ped))
  expect_true(all(pc$sdev[1:2] - pc2$sdev[1:2] < 1e-10))
})















