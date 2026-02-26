library(data.table)


#' Compute Discounted Probabilities for a Single Order; Used by Interpolation PPM
#'
#' @param dt Data.table from `dynamic_count_table_for_each_order` for a single order.
#'           Must have columns: index, context_id, Event, Ce, C, t, t1
#' @param discount_func Function to compute discounted probability mass.
#'                      Signature: function(C, t, t1) returns list(lambda, k)
#'
#' @return Data.table with columns: index, Event, probability
#' @export
compute_discounted_probs <- function(dt, discount_func) {
  dt[, {
    stats  <- discount_func(C[1], t[1], t1[1])
    lambda <- if (C[1] > 0) stats$lambda else 0

    list(lambda = lambda)
  }, by = .(index, context_id, Event)]
  #disc_stats <- discount_func(dt$C, dt$t, dt$t1)
  #lambda <- disc_stats$lambda

  #out <- data.table(
   # index = dt$index,
    #Event = dt$Event,
    #lambda = lambda
  #)

  # normalization
  #out[, lambda := {
   # s <- sum(lambda)
    #if (s > 0) lambda / s else lambda
  #}, by = index]
}

#' Compute PPM Probabilities with Interpolation
#'
#' Vectorized interpolated PPM.
#' Each order contributes to the final probability weighted by its discounted probability mass.
#'
#' @param x Character vector of events
#' @param N Maximum order
#' @param alphabet Character vector of full alphabet
#' @param discount_func Discount function (e.g., `discount_C`).
#'
#' @return data.table with columns: index, Event, P, IC
#' @export
ppm_interpolated <- function(x, N, alphabet, discount_func=discount_C) {
  T <- length(x)
  alpha_len <- length(alphabet)

  # 1. Dynamic count tables
  dt_orders <- dynamic_count_table_for_each_order(x, N = N, alphabet = alphabet)

  # 2. Base probabilities: 1 / (∣alphabet∣+1−tseen)
  base_prob <- numeric(T * alpha_len)
  seen_symbols <- character(0)
  for (t in seq_len(T)) {
    if (t > 1) seen_symbols <- unique(c(seen_symbols, x[t - 1]))
    denom <- length(alphabet) + 1 - length(seen_symbols)
    idx <- ((t - 1) * alpha_len + 1):(t * alpha_len)
    base_prob[idx] <- 1 / denom
  }

  # 3. Compute PPM probabilities
  # Initialize probability matrix: rows = timesteps × symbols
  dt_final <- copy(dt_orders[[N + 1]])[, .(index, Event)]
  P <- base_prob

  # Precompute discounted probabilities
  discount_probs <- lapply(dt_orders, compute_discounted_probs,
                           discount_func = discount_C)
  # Interpolation loop (from lowest to highest order)
  for (n in 0:N) {
    dt_n <- dt_orders[[n + 1]]
    lambda_dt <- discount_probs[[n + 1]]
    # Apply interpolation formula:
    # Interpolated P = lambda * (Ce / C if C > 0 else 0) + (1 - lambda) * P_lower
    # TODO: check what k does
    Ce_over_C <- ifelse(dt_n$C > 0,
                        dt_n$Ce / dt_n$C,
                        0)
    lambda <- lambda_dt$lambda

    P <- lambda * Ce_over_C + (1 - lambda) * P
  }

  dt_final[, P := P]
  # Normalization
  dt_final[, P := {
    s <- sum(P)
    if (s > 0) P / s else P
  }, by = index]
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

print("PPIDyOM Result:")
test <- c("A", "B", "A", "C", "A", "B", "A", "C", "A")
alphabet <- c("A", "B", "C")
max_order <- 3
result <- ppm_interpolated(
  x = test,
  alphabet = alphabet,
  N = max_order,
  discount_func = discount_C
)
print(result)


print("PPM Result:")
seq <- factor(x, levels = alphabet)
mod <- new_ppm_simple(
  order_bound = max_order,
  alphabet_levels = alphabet,
  escape = "c",
  shortest_deterministic = FALSE,
  exclusion = FALSE,
  update_exclusion = FALSE
)
res <- model_seq(mod, seq)
print(res)

