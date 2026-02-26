library(data.table)
library(ppm)
#' @import data.table

# ===========COUNT MATRIX===========

#' Create lagged matrix for N-grams
#'
#' @param x Vector of tokens.
#' @param N Maximum N-gram length.
#' @param ... Additional vectors of same length as \code{x}.
#' @return A \code{data.table} with lagged columns LagN..Lag0 and an index column.

lagMatrix <- function(x, N = 10, ...) {

  Ntab <- as.data.table(lapply(N:0, \(n) shift(x, n)))

  colnames(Ntab) <- paste0('Lag', N:0)

  Ntab <- cbind(Ntab, as.data.table(list(...)))

  Ntab[ , index := 1:nrow(Ntab)]

  Ntab[]
}

#' Build dynamic count matrix
#'
#' @param x Vector of tokens.
#' @param N Maximum N-gram length.
#' @param prior Optional prior model (countMatrix) to add counts to.
#' @param ... Additional vectors of same length as \code{x}.
#' @return A \code{data.table} with counts for each N-gram order.
dynamicModel <- function(x, N = 5, prior = NULL, ...) {

  countMatrix <- lagMatrix(x, N)

  # Lag0 is the next event; Lag0+ constitute the context
  lags <- setdiff(colnames(countMatrix), c('Lag0', 'index'))

  for (n in 0:N) {

    Ncolname <- paste0('N', n)

    if (n == 0) {
      groupby <- "Lag0"
    } else {
      groupby <- tail(lags, n)
    }

    countMatrix <- countMatrix[ , (Ncolname) := seq_along(index) - 1, by = groupby]


    #
    if (!is.null(prior)) {
      priorFinal <- prior[, list(PriorCount = max(get(Ncolname)) + 1), by = groupby] # + 1 because prior countMatrix doesn't count the LAST time of each thing

      countMatrix <- merge(countMatrix, priorFinal, by = groupby, all.x = TRUE)
      countMatrix[, (Ncolname) := get(Ncolname) + ifelse(is.na(PriorCount), 0, PriorCount)]
      countMatrix[, PriorCount := NULL]
      #
    }
  }

  setorder(countMatrix, index)

  countMatrix[]
}



# ========ESCAPE FUNCTIONS==========

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



# ========= PREDICT =============

#' Compute PPM probabilities for a given order
#'
#' Computes the probability of each event given its context at a specific order,
#' using a specified escape method and count type.
#'
#' @param countMatrix A \code{data.table} with lagged counts from \code{dynamicModel}.
#' @param N Integer. The order to compute probabilities for.
#' @param escape_method Function. An escape function that returns \code{denom}, \code{esc}, and \code{subtract}.
#' @param prob_count_type Character. Determines which counts are used to compute probabilities:
#'   \describe{
#'     \item{"C"}{Divident for probability of an event is the number of times this context has occurred, consistent with Harrison (2020)'s PPM.}
#'     \item{"Ce"}{Divident for probability of an event is the number of times this event has occurred after this context, consistent with Pearce (2005)'s paper.}
#'   }
#'
#' @return A \code{data.table} identical to \code{countMatrix} with added columns:
#'   \code{prob}, \code{esc}, \code{C}, \code{t}, \code{t1}, \code{event_count}.
#' @examples
#' x <- c("A", "B", "A", "C")
#' cm <- dynamicModel(x, N = 2)
#' compute_probs_for_n(cm, N = 2, escape_method = escape_C, prob_count_type = "Ce")

compute_probs_for_n <- function(
    countMatrix, N,
    escape_method = escape_C, prob_count_type = "Ce"
) {
  count_col <- paste0("N", N)
  context_cols <- paste0("Lag", N:1)
  probMatrix <- copy(countMatrix)

  # compute context
  stats <- probMatrix[
    ,
    {
      # number of occurrences for each event after the context
      occ <- ave(seq_along(Lag0), Lag0, FUN = seq_along)
      # If occurred once, singleton; twice, no longer singleton
      singleton_delta <- fifelse(occ == 1,  1,
                                 fifelse(occ == 2, -1, 0))

      .(
        index = index,
        # C: Context count; Total times this context has occurred
        C  = get(count_col),
        # t: Number of distinct events seen after this context
        t  = shift(cumsum(!duplicated(Lag0)), fill = 0),
        # t1: Number of singletons (events seen exactly once) after this context
        t1 = shift(cumsum(singleton_delta), fill = 0)
      )
    },
    by = context_cols
  ]
  setorder(stats, index)
  probMatrix[, `:=`(
    C = stats$C,
    t = stats$t,
    t1 = stats$t1
  )]
  # compute event counts for each row within its context
  probMatrix[
    ,
    event_count := seq_len(.N) - 1,
    by = c(context_cols, "Lag0")
  ]


  # compute event probability and escape probability
  probMatrix[, `:=`(prob = 0, esc = NA_real_)]

  # if the context C(N) is unseen, immediately back off to C(N-1);
  # No need to compute event/escape probabilities for this context.
  valid <- probMatrix$C > 0

  if (any(valid)) {
    esc_info <- escape_method(
      probMatrix$C[valid],
      probMatrix$t[valid],
      probMatrix$t1[valid]
    )
    if (prob_count_type == "C") {
      probMatrix[valid, probCount := C]
    } else {
      probMatrix[valid, probCount := event_count]
    }

    probMatrix[valid,
              `:=`(
                mle = pmax(probCount - esc_info$subtract, 0),
                prob = pmax(probCount - esc_info$subtract, 0) / esc_info$denom,
                esc  = esc_info$esc
              )
    ]
  }
  # MLE: Maximum Likelihood Estimates (the probability mass at each order)
  probMatrix[, mle := ifelse(C > 0, event_count / C, 0)]

  probMatrix[]
}

#' Compute PPM probabilities for a sequence
#'
#' Computes the probabilities of events in a sequence using a Prediction by Partial Matching (PPM) model.
#' Supports both **backoff** and **interpolated** smoothing, with flexible base distribution for unseen symbols.
#'
#' @param countMatrix A \code{data.table} containing lagged counts for each N-gram order, typically generated by \code{dynamicModel()}.
#' @param max_order Integer. Maximum N-gram order to use for probability estimation.
#' @param alphabet Character vector. The full set of possible symbols.
#' @param escape_method Function. Escape function to calculate escape probability and denominators. Should return a list with \code{denom}, \code{esc}, and \code{subtract}.
#' @param prob_count_type Character. Determines the numerator for event probabilities:
#'   \describe{
#'     \item{"C"}{Total number of times the context has occurred. Matches Harrison (2020)'s approach.}
#'     \item{"Ce"}{Number of times this event has occurred after the context. Matches Pearce (2005)'s approach.}
#'   }
#' @param smoothing Character. Smoothing method to use. Options:
#'   \describe{
#'     \item{"backoff"}{Backoff PPM: recursively back off to lower-order contexts when an event is unseen.}
#'     \item{"interpolated"}{Interpolated PPM: mix probabilities from higher and lower orders using escape weights.}
#'   }
#'   Default: \code{c("backoff", "interpolated")}, selects \code{"backoff"} if unspecified.
#' @param base_fallback Character. Base distribution to use for globally unseen symbols (order −1):
#'   \describe{
#'     \item{"uniform"}{Assign equal probability to all symbols in \code{alphabet}.}
#'     \item{"unseen"}{Assign probability only to symbols not yet observed: 1 / (|alphabet| + 1 − t), where \code{t} is the number of distinct symbols observed.}
#'   }
#'   Default: \code{c("uniform", "unseen")}.
#'
#' @return A \code{data.table} with columns:
#'   \describe{
#'     \item{index}{Row index.}
#'     \item{Event}{The current token (Lag0).}
#'     \item{Order}{The N-gram order used to resolve the event probability.}
#'     \item{P}{Probability of the event.}
#'     \item{IC}{Information content of the event, \(-\log_2 P\).}
#'   }
#'
#' @examples
#' # Example with a simple sequence
#' alphabet <- c("A", "B", "C")
#' x <- c("A", "B", "A", "C", "A")
#' cm <- dynamicModel(x, N = 2)
#' compute_ppm_table(
#'   countMatrix = cm,
#'   max_order = 2,
#'   alphabet = alphabet,
#'   escape_method = escape_C,
#'   prob_count_type = "Ce",
#'   smoothing = "backoff",
#'   base_fallback = "unseen"
#' )
#'
#' @export
compute_ppm_table <- function(countMatrix,
                              max_order,
                              alphabet,
                              escape_method = escape_C,
                              prob_count_type = "Ce",
                              smoothing = c("backoff", "interpolated"),
                              base_fallback = c("uniform", "unseen")) {

  smoothing <- match.arg(smoothing)
  base_fallback <- match.arg(base_fallback)

  if (smoothing == "backoff") {
    return(compute_ppm_backoff(
      countMatrix, max_order, alphabet,
      base_fallback,
      escape_method, prob_count_type
    ))
  }

  if (smoothing == "interpolated") {
    return(compute_ppm_interpolated(
      countMatrix, max_order, alphabet,
      base_fallback,
      escape_method, prob_count_type
    ))
  }
}

compute_ppm_backoff <- function(countMatrix,
                                max_order,
                                alphabet,
                                base_fallback,
                                escape_method,
                                prob_count_type) {
  dt <- copy(countMatrix)

  dt[, `:=`(
    P_final = NA_real_,
    order_used = NA_integer_,
    p_mass = 1,
    resolved = FALSE
  )]

  # compute all orders
  for (n in max_order:1) {
    probs_n <- compute_probs_for_n(dt, n, escape_method, prob_count_type)

    valid <- !dt$resolved & probs_n$C > 0
    event_seen <- valid & probs_n$prob > 0
    event_unseen <- valid & probs_n$prob <= 0

    # resolve events seen after context
    dt[event_seen, `:=`(
      P_final = p_mass * probs_n$prob[event_seen],
      order_used = n,
      resolved = TRUE
    )]

    # update escape mass for unresolved events
    dt[event_unseen, p_mass := p_mass * probs_n$esc[event_unseen]]
  }

  # Order -1 base fallback
  fallback <- !dt$resolved
  if (any(fallback)) {
    t_base <- length(unique(dt$Lag0[!is.na(dt$Lag0)]))

    if (base_fallback == "uniform") {
      base_prob <- 1 / length(alphabet)
    } else {
      base_prob <- 1 / (length(alphabet) + 1 - t_base)
    }

    dt[fallback, `:=`(
      P_final = p_mass * base_prob[fallback],
      order_used = -1,
      resolved = TRUE
    )]
  }


  # Information Content
  dt[, IC := -log2(P_final)]

  dt[, .(
    index,
    Event = Lag0,
    Order = order_used,
    P = P_final,
    IC = IC
  )]
}

compute_ppm_interpolated <- function(countMatrix,
                                     max_order,
                                     alphabet,
                                     base_fallback,
                                     escape_method,
                                     prob_count_type) {

  dt <- copy(countMatrix)

  # store per-order stats
  order_stats <- vector("list", max_order + 1)

  for (n in 1:max_order) {

    probs_n <- compute_probs_for_n(
      dt, n,
      escape_method = escape_method,
      prob_count_type = prob_count_type
    )

    order_stats[[n + 1]] <- probs_n[, .(
      index,
      C,
      esc,
      mle
    )]
  }
  print(order_stats)

  # base case (order -1)
  # number of distinct symbols observed in the sequence
  t_base <- length(unique(dt$Lag0[!is.na(dt$Lag0)]))
  print(t_base)

  if (base_fallback == "uniform") {
    base_prob <- rep(1 / length(alphabet), nrow(dt))
  } else {
    base_prob <- 1 / (length(alphabet) + 1 - t_base)
  }

  print(base_prob)
  P <- rep(base_prob, nrow(dt))

  for (n in 1:max_order) {

    stats_n <- order_stats[[n + 1]]

    has_context <- stats_n$C > 0

    P[has_context] <-
      stats_n$esc[has_context] * stats_n$mle[has_context] +
      (1 - stats_n$esc[has_context]) * P[has_context]
  }

  dt[, P_final := P]

  # normalize per time index
  dt[, P_final := P_final / sum(P_final), by = index]

  dt[, IC := -log2(P_final)]
  dt[, .(
    index,
    Event = Lag0,
    P = P_final,
    IC = IC
  )]
}




# =========TEST==========

# simple comparison with Harrison's PPM library

alphabet <- c("A", "B", "C")
x <- c("A", "B", "A", "C", "A", "B", "A", "C", "A")

max_order <- 3
cm <- dynamicModel(x, N = max_order)
cm
ppm <- compute_ppm_table(
  cm,
  max_order = max_order,
  alphabet = alphabet,
  escape_method = escape_C,
  prob_count_type = "C",
  smoothing = "interpolated",
  base_fallback = "unseen"
)
ppm

seq <- factor(x, levels = alphabet)
mod <- new_ppm_simple(
  order_bound = 3,
  alphabet_levels = alphabet,
  escape = "c",
  shortest_deterministic = FALSE,
  exclusion = FALSE,
  update_exclusion = FALSE
)
res <- model_seq(mod, seq)
res$distribution

# -------------------------


max_order <- 3
alphabet <- c('I', 'ii', 'III', 'IV', 'V', 'vi', 'vii')
test <- c('IV', 'I', 'I', 'IV', 'I', 'IV', 'I', 'V', 'I', 'IV', 'V')
countMatrix <- dynamicModel(test, N=max_order)
countMatrix
ppm_table <- compute_ppm_table(countMatrix, 3, alphabet, escape_method = escape_A)
ppm_table


