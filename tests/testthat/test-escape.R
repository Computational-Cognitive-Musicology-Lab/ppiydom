test_that("escape methods compute correct values", {

  C  <- 10
  t  <- 3
  t1 <- 2

  a <- escape_A(C, t, t1)
  expect_equal(a$denom, 11)
  expect_equal(a$esc, 1/11)
  expect_equal(a$subtract, 0)

  b <- escape_B(C, t, t1)
  expect_equal(b$denom, 10)
  expect_equal(b$esc, 3/10)
  expect_equal(b$subtract, 1)

  c <- escape_C(C, t, t1)
  expect_equal(c$denom, 13)
  expect_equal(c$esc, 3/13)
  expect_equal(c$subtract, 0)

  d <- escape_D(C, t, t1)
  expect_equal(d$denom, 10 + 0.5 * 3)
  expect_equal(d$esc, 3 / (10 + 0.5 * 3))
  expect_equal(d$subtract, 0.5)

  ax <- escape_AX(C, t, t1)
  expect_equal(ax$denom, 10 + 3 + 1)
  expect_equal(ax$esc, (2 + 1) / (10 + 3 + 1))
  expect_equal(ax$subtract, 0)
})
