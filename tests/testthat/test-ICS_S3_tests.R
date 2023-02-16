
# ics - ics2 --------------------------------------------------------------

test_that("ics - ics2 eigenvalues", {
  X <- iris[,1:4]
  expect_equal(ics(X, S1 = cov, S2 = cov4, stdKurt = FALSE)@gKurt,
               ics2(X, S1 = MeanCov, S2 = Mean3Cov4)@gKurt)
})



# ICS ---------------------------------------------------------------------

## Eigenvalues -------------------------------------------------------------

test_that("ics - ICS_S3 eigenvalues - S1 and S2 are functions", {
  X <- iris[,1:4]
  expect_equal(ics(X, S1 = cov, S2 = cov4, stdKurt = FALSE)@gKurt,
               ICS(X, S1 = ICS_cov, S2 = ICS_cov4)$gKurt)
})


test_that("ics - ICS_S3 eigenvalues - S1 and S2 are functions", {
  X <- iris[,1:4]
  expect_equal(ics(X, S1 = cov, S2 = cov4, stdKurt = FALSE)@gKurt,
               ICS(X, S1 = cov, S2 = ICS_cov4)$gKurt)
})

test_that("ics - ICS_S3 eigenvalues - S1 and S2 are functions", {
  X <- iris[,1:4]
  mcd <- function(x,...) robustbase::covMcd(x, ...)$cov
  expect_equal(ics(X, S1 = mcd, S2 = cov, S1args = list(nsamp = "deterministic", alpha = 0.5), stdKurt = FALSE)@gKurt,
               ICS(X, S1 = mcd, S2 = ICS_cov,
                   S1_args = list(nsamp = "deterministic", alpha = 0.5))$gKurt)
})

test_that("ics - ICS_S3 scores - S1 and S2 are functions - eigenvalues standardization", {
  X <- iris[,1:4]
  expect_equal(ics(X, S1 = cov, S2 = cov4, stdKurt = TRUE)@gKurt,
               ICS(X, S1 = ICS_cov, S2 = ICS_cov4, stdKurt = TRUE, center = FALSE)$gKurt)
})



test_that("ics2 - ICS_S3 eigenvalues - S1 and S2 are functions", {
  X <- iris[,1:4]
  expect_equal(ics2(X, S1 = MeanCov, S2 = Mean3Cov4)@gKurt,
               ICS(X, S1 = ICS_cov, S2 = ICS_cov4)$gKurt)
})



test_that("ics - ICS_S3 eigenvalues - S1 and S2 are functions", {
  X <- iris[,1:4]
  expect_equal(ICS(X, S1 = ICS_cov, S2 = ICS_cov4, withen = FALSE)$gKurt,
               ICS(X, S1 = ICS_cov, S2 = ICS_cov4, withen = TRUE)$gKurt)
})


test_that("ics - ICS_S3 eigenvalues - S1 and S2 are matrices", {
  X <- iris[,1:4]
  expect_equal(ics(X, S1 = cov(X), S2 = cov4(X), stdKurt = FALSE,
              stdB = "B")@gKurt,
              ICS(X, S1 = cov(X), S2 = cov4(X), standardize = "scores")$gKurt)
})
test_that("ics - ICS_S3 eigenvalues - S1 and S2 are matrices", {
  X <- iris[,1:4]
  expect_equal(ics(X, S1 = cov(X), S2 = cov4(X), stdKurt = FALSE)@gKurt,
               ICS(X, S1 = ICS_cov(X)$scatter,
                   S2 =  ICS_cov4(X)$scatter)$gKurt)
})


test_that("ics - ICS_S3 eigenvalues - S1 and S2 are matrices/ICS_scatter", {
  X <- iris[,1:4]
  expect_equal(ics(X, S1 = cov(X), S2 = cov4(X), stdKurt = FALSE)@gKurt,
               ICS(X, S1 = ICS_cov(X), S2 = ICS_cov4)$gKurt)
})

# TODO: not ok
# Aurore: what do we want to do here? should we put some warnings saying it is not possible?
test_that("ics - ICS_S3 eigenvalues -  S2 is ICS_scatter", {
  X <- iris[,1:4]
  expect_equal(ics(X, S1 = cov(X), S2 = cov4(X), stdKurt = FALSE)@gKurt,
               ICS(X, S1 = ICS_cov(X), S2 = ICS_cov4(X))$gKurt)
})


# Scores ------------------------------------------------------------------
test_that("ics - ICS_S3 scores - S1 and S2 are functions - scores standardization", {
  X <- iris[,1:4]
  expect_equal(ics(X, S1 = cov, S2 = cov4, stdKurt = FALSE, stdB = "Z")@Scores,
               ICS(X, S1 = ICS_cov, S2 = ICS_cov4, standardize = "scores", center = FALSE)$Scores)
})


test_that("ics - ICS_S3 scores - S1 and S2 are functions - eigenvalues standardization", {
  X <- iris[,1:4]
  expect_equal(ics(X, S1 = cov, S2 = cov4, stdKurt = FALSE, stdB = "B")@Scores,
               ICS(X, S1 = ICS_cov, S2 = ICS_cov4, standardize = "eigenvalues", center = FALSE)$Scores)
})


test_that("ics2 - ICS_S3 scores - S1 and S2 are functions - centering", {
  X <- iris[,1:4]
  expect_equal(ics2(X, S1 = MeanCov, S2 = Mean3Cov4)@Scores,
               ICS(X, S1 = ICS_cov, S2 = ICS_Mean3Cov4, standardize = "scores", center = TRUE)$Scores)
})

## QR ----------------------------------------------------------------------
test_that("ics - ICS_S3 eigenvalues - S1 and S2 are functions - QR", {
  X <- iris[,1:4]
  expect_equal(ics2(X, S1 = MeanCov, S2 = Mean3Cov4)@gKurt,
               ICS(X, S1 = ICS_cov, S2 =  ICS_covW, S2_args = list(alpha = 1, cf = 1/(ncol(X)+2)), QR = TRUE)$gKurt)
})


test_that("ics - ICS_S3 eigenvalues - S1 and S2 are functions - QR", {
  X <- iris[,1:4]
  expect_equal(ics2(X, S1 = MeanCov, S2 = Mean3Cov4)@gKurt,
               ICS(X, S1 = ICS_cov, S2 =  ICS_covW, S2_args = list(alpha = 1, cf = 1/(ncol(X)+2)), QR = FALSE)$gKurt)
})


test_that("ics - ICS_S3 scores - S1 and S2 are functions - QR ", {
  X <- iris[,1:4]
  expect_equal(ics(scale(X, center = TRUE, scale = FALSE), S1 = cov, S2 = cov4, stdB = "Z", stdKurt = FALSE)@Scores,
               ICS(X, S1 = ICS_cov, S2 =  ICS_covW, S2_args = list(alpha = 1, cf = 1/(ncol(X)+2)), QR = TRUE, standardize = "scores")$Scores)
})

test_that("ics2 - S1 and S2 are functions - QR", {
  X <- iris[,1:4]
  X_rank_deficient <- sweep(X, 2, c(10^(-12), 10^(-3), 1, 10^6), "*")
  expect_error(ics2( X_rank_deficient, S1 = MeanCov, S2 = Mean3Cov4)@gKurt)
})

test_that("ics2 - S1 and S2 are functions - QR", {
  X <- iris[,1:4]
  X_rank_deficient <- sweep(X, 2, c(10^(-12), 10^(-3), 1, 10^6), "*")
  expect_error( ICS( X_rank_deficient, S1 = ICS_cov, S2 =  ICS_covW, S2_args = list(alpha = 1, cf = 1/(ncol(X)+2)), QR = FALSE)$gKurt)
})

test_that("ics2 - S1 and S2 are functions - QR", {
  X <- iris[,1:4]
  X_rank_deficient <- sweep(X, 2, c(10^(-12), 10^(-3), 1, 10^6), "*")
  expect_equal(ics(X, S1 = cov, S2 = cov4, stdB = "Z", stdKurt = FALSE)@gKurt,
    ICS( X_rank_deficient, S1 = ICS_cov, S2 =  ICS_covW, S2_args = list(alpha = 1, cf = 1/(ncol(X)+2)), QR = TRUE)$gKurt)
})

