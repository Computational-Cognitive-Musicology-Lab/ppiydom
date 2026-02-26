library(data.table)
#' @import data.table

#' Generate Lagged N-gram Matrix
#'
#' Creates a lagged representation of a sequence for N-gram modeling.
#' @param x Character vector of symbols/events.
#' @param N Maximum N-gram order.
#' @return A `data.table` with columns LagN..Lag0 and index.
#' @export
lag_matrix <- function(x, N = 3) {
  dt <- data.table::as.data.table(lapply(N:0, function(n) data.table::shift(x, n)))
  data.table::setnames(dt, paste0("Lag", N:0))
  dt[, index := seq_len(.N)]
  dt
}

#' Build Dynamic Count Table
#'
#' Computes dynamic PPM counts for all symbols at all timesteps and orders.
#' Counts are cumulative up to (but excluding) the current timestep.
#'
#' @param x Character vector of symbols.
#' @param N Maximum order.
#' @param alphabet Character vector of full alphabet.
#'
#' @return A list of length N+1.
#' Each element is a data.table for each order with columns:
#' index, context_id, Event, Ce, C, t, t1
#' Event would enumerate across all symbols in the alphabet
#' @export
#'
#' TODO: add prior
dynamic_count_table_for_each_order <- function(x, N, alphabet) {

  T <- length(x)
  dt_lag <- lag_matrix(x, N)
  results <- vector("list", N + 1)

  for (n in 0:N) {
    context_cols <- if (n == 0) character(0) else paste0("Lag", n:1)

    # create context identifier
    if (n == 0) {
      dt_lag[, context_id := "ROOT"]
    } else {
      dt_lag[, context_id := do.call(paste, c(.SD, sep = "_")), .SDcols = context_cols]
    }

    # environment to store running counts per context
    counts <- new.env(hash = TRUE, parent = emptyenv())
    out <- vector("list", T)

    for (t in seq_len(T)) {

      ctx <- dt_lag$context_id[t]
      ctx_counts <- counts[[ctx]]
      if (is.null(ctx_counts)) ctx_counts <- integer(0)

      # Ce for all symbols in the alphabet (fill missing with 0)
      Ce_full <- setNames(integer(length(alphabet)), alphabet)
      if (length(ctx_counts) > 0) Ce_full[names(ctx_counts)] <- ctx_counts

      # Build row per symbol
      dt_step <- data.table(
        index      = t,
        context_id = ctx,
        Event      = alphabet,
        Ce         = as.integer(Ce_full),
        C          = sum(Ce_full),
        t          = sum(Ce_full > 0L),
        t1         = sum(Ce_full == 1L)
      )

      out[[t]] <- dt_step

      # update counts for this timestep AFTER storing
      sym <- x[t]
      if (is.null(counts[[ctx]])) {
        counts[[ctx]] <- setNames(1L, sym)
      } else if (sym %in% names(counts[[ctx]])) {
        counts[[ctx]][sym] <- counts[[ctx]][sym] + 1L
      } else {
        counts[[ctx]][sym] <- 1L
      }
    }

    results[[n + 1]] <- rbindlist(out)
  }

  results
}
