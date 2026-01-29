test_that("probabilities are valid", {

  alphabet <- LETTERS[1:5]
  x <- sample(alphabet, 50, replace = TRUE)

  order <- 3
  cm <- dynamicModel(x, N = order)
  ppm <- compute_ppm_table(cm, order, alphabet)

  expect_true(all(ppm$P > 0))
  expect_true(all(ppm$P <= 1))
})

test_that("resolved rows have exactly one order_used", {

  alphabet <- c("I","IV","V")
  x <- c("I","IV","V","I","IV","V")

  order <- 3
  cm <- dynamicModel(x, N = order)
  ppm <- compute_ppm_table(cm, order, alphabet)

  expect_true(all(!is.na(ppm$Order)))
})
