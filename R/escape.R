# C: Total times this context has occurred
# t: Number of distinct events seen after this context
# t1: Number of events seen exactly once (singletons) after this context
# esc: escape probability
# denom: denominator for seen events
# subtract: optional count subtraction

escape_A <- function(C, t, t1) {
  list(
    denom = C + 1,
    esc = 1 / (C + 1),
    subtract = 0
  )
}

escape_B <- function(C, t, t1) {
  list(
    denom = C,
    esc = t / C,
    subtract = 1
  )
}

escape_C <- function(C, t, t1) {
  list(
    denom = C + t,
    esc = t / (C + t),
    subtract = 0
  )
}

escape_D <- function(C, t, t1) {
  list(
    denom = C + 0.5 * t,
    esc = t / (C + 0.5 * t),
    subtract = 0.5
  )
}

escape_AX <- function(C, t, t1) {
  list(
    denom = C + t + 1,
    esc = (t1 + 1) / (C + t + 1),
    subtract = 0
  )
}


discount_A <- function(C, t, t1) {
  list(
    lambda = C / (C + 1),
    k = 0
  )
}

discount_B <- function(C, t, t1) {
  list(
    lambda = C / (C + t),
    k = -1
  )
}

discount_C <- function(C, t, t1) {
  list(
    lambda = C / (C + t),
    k = 0
  )
}

discount_D <- function(C, t, t1) {
  list(
    lambda = C / (C + t/2),
    k = -1/2
  )
}

discount_AX <- function(C, t, t1) {
  list(
    lambda = C / (C + t1),
    k = 0
  )
}

