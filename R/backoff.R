library(data.table)


#' Compute Local Probabilities for a Single Order; Used by Backoff PPM
#'
#' @param dt Data.table from `dynamic_count_table_for_each_order` for a single order.
#'           Must have columns: index, context_id, Event, Ce, C, t, t1
#' @param escape_func Function to compute escape probability.
#'                    Signature: function(C, t, t1) returns list(denom, esc, subtract)
#'
#' @return Data.table with columns: index, Event, prob_local, esc
#' @export
compute_local_probs <- function(dt, escape_func, normalize = FALSE) {

  # Initialize vectors
  prob_local <- numeric(nrow(dt))
  esc        <- numeric(nrow(dt))

  # Compute escape and denominator per row
  escape_stats <- escape_func(dt$C, dt$t, dt$t1)

  denom <- escape_stats$denom
  esc[] <- escape_stats$esc
  subtract <- escape_stats$subtract

  # Numerator = Ce - subtract (cannot be negative)
  numer <- pmax(dt$Ce - subtract, 0)

  # Only assign probabilities for valid rows (C > 0 and numer > 0)
  valid_idx <- which(!is.na(denom) & denom > 0 & numer >= 0)
  prob_local[valid_idx] <- numer[valid_idx] / denom[valid_idx]

  out <- data.table(
    index = dt$index,
    Event = dt$Event,
    prob_local = prob_local,
    esc = esc
  )
  out
}


#' Compute PPM Probabilities with Backoff
#'
#' Vectorized backoff implementation of PPM.
#' For each timestep, uses highest available order with non-zero counts,
#' falling back to lower orders and finally to base probabilities.
#'
#' @param x Character vector of events
#' @param N Maximum order
#' @param alphabet Character vector of full alphabet
#' @param escape_func Escape function (e.g., `escape_C`).
#'
#' @return data.table with columns: index, Event, P, IC
#' @export
ppm_backoff <- function(x, N, alphabet, escape_func=escape_C) {

  T <- length(x)
  alpha_len <- length(alphabet)

  # 1. Dynamic count tables
  dt_orders <- dynamic_count_table_for_each_order(x, N = N, alphabet = alphabet)

  # 2. Base probability: 1 / (∣alphabet∣)
  base_prob <- rep(1 / alpha_len, T * alpha_len)
  print(base_prob)

  # 3. Compute PPM probabilities
  dt_final <- copy(dt_orders[[N + 1]])[, .(index, Event)]
  dt_final[, P := NA_real_]

  # Initialize leftover mass
  p_mass <- rep(1, nrow(dt_final))

  # Compute local probabilities for all orders in advance
  local_probs <- lapply(dt_orders, compute_local_probs,
                        escape_func = escape_func)

  # Backoff loop from highest to lowest order
  for (n in N:0) {
    dt_n <- local_probs[[n + 1]]

    # Rows with seen events (Ce > 0)
    seen <- dt_orders[[n + 1]]$Ce > 0

    dt_final[seen, P := p_mass[seen] * dt_n$prob_local[seen]]
    p_mass[seen] <- p_mass[seen] * (1 - dt_n$prob_local[seen])
  }

  # Assign remaining mass to base probabilities
  remaining <- is.na(dt_final$P)
  dt_final[remaining, P := p_mass[remaining] * base_prob[remaining]]

  dt_final[, IC := -log2(P)]

  # 4. Compute entropy per timestep
  Entropy <- compute_entropy(dt_final)

  # 5. Merge actual events with computed probabilities
  # select only the probability for the actual event
  # TODO: order
  dt_result <- dt_final[
    data.table(index = seq_len(T), Event = x),
    on = .(index, Event)
  ][
    , .(index, Event, P, IC, Entropy)
  ]

  dt_result
}

test <- c("A", "B", "A", "C", "A", "B", "A", "C", "A")
alphabet <- c("A", "B", "C")
max_order <- 3
result <- ppm_backoff(
  x = test,
  N = max_order,
  alphabet = alphabet,
  escape_func = escape_C
)

print(result)
