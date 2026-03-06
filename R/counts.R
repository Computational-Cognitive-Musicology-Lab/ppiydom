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

#' Compute Count Tables for STM and LTM with Optional Prior
#'
#' Generates count tables for STM (per-timestep) and/or LTM (per-context),
#' optionally incorporating previously accumulated LTM counts.
#'
#' @param x Character vector of symbols/events.
#' @param N Maximum N-gram order.
#' @param alphabet Character vector of all possible symbols.
#' @param model_type Character: one of `"stm"`, `"ltm"`, `"both"`.
#' @param prior Optional: previously accumulated LTM tables (list of length N+1).
#'   Used to initialize/accumulate counts for LTM or `both` type.
#' @return A list with elements depending on `model_type`:
#'   - `$stm`: list of length N+1, each a data.table of counts per timestep (if `model_type` includes `"stm"`).
#'   - `$ltm`: list of length N+1, each a data.table of counts per context (if `model_type` includes `"ltm"`).
#' @export
count_tables <- function(x, N, alphabet, model_type = c("stm","ltm","both"), prior = list()) {

  model_type <- match.arg(model_type)
  T <- length(x)
  dt_lag <- lag_matrix(x, N)

  stm_tables <- if(model_type %in% c("stm","both")) vector("list", N+1) else NULL
  ltm_tables <- if(model_type %in% c("ltm","both")) vector("list", N+1) else NULL


  for(n in 0:N) {
    context_cols <- if(n == 0) character(0) else paste0("Lag", n:1)
    # Context identifiers
    if(n==0) dt_lag[, context_id := "ROOT"]
    else dt_lag[, context_id := do.call(paste, c(.SD, sep="_")), .SDcols=context_cols]

    counts <- new.env(hash=TRUE, parent=emptyenv())

    # Initialize STM or LTM
    if(model_type %in% c("stm","both")) stm_list <- vector("list", T)
    ctx_ids <- unique(dt_lag$context_id)
    if(model_type %in% c("ltm","both")) ltm_list <- vector("list", length(ctx_ids))

    # Loop over timesteps for STM, or LTM accumulation
    for(t in seq_len(T)) {
      ctx <- dt_lag$context_id[t]

      # If STM, create table for the current timestep
      if(model_type %in% c("stm","both")) {
        ctx_counts <- counts[[ctx]]
        if(is.null(ctx_counts)) ctx_counts <- integer(0)

        # Ce for all symbols in the alphabet (fill missing with 0)
        Ce_full <- setNames(integer(length(alphabet)), alphabet)
        if (length(ctx_counts) > 0) Ce_full[names(ctx_counts)] <- ctx_counts

        stm_list[[t]] <- data.table(
          index = t,
          context_id = ctx,
          Event = alphabet,
          Ce = as.integer(Ce_full),
          C = sum(Ce_full),
          t = sum(Ce_full > 0L),
          t1 = sum(Ce_full == 1L)
        )
      }

      # Update symbol-by-context counts for either STM or LTM
      sym <- x[t]
      if(is.null(counts[[ctx]])) counts[[ctx]] <- setNames(1L, sym)
      else if(sym %in% names(counts[[ctx]])) counts[[ctx]][sym] <- counts[[ctx]][sym] + 1L
      else counts[[ctx]][sym] <- 1L
    }
    if(model_type %in% c("stm","both")) stm_tables[[n+1]] <- rbindlist(stm_list)

    # Loop over contexts for LTM
    if(model_type %in% c("ltm","both")) {
      # Add prior counts, if provided, to current sequence's counts
      if (length(prior) > 0) {
        prior_dt <- prior[[n+1]]
        if (!is.null(prior_dt) && nrow(prior_dt) > 0) {
          for (ctx in unique(prior_dt$context_id)) {
            prior_ctx <- prior_dt[context_id == ctx]
            prior_vec <- setNames(as.integer(prior_ctx$Ce), prior_ctx$Event)

            # Ensure counts[[ctx]] has all alphabet symbols
            if(is.null(counts[[ctx]])) {
              counts[[ctx]] <- setNames(integer(length(alphabet)), alphabet)
            } else {
              # Fill missing symbols with 0 to avoid NA
              missing_syms <- setdiff(names(prior_vec), names(counts[[ctx]]))
              if(length(missing_syms) > 0) counts[[ctx]][missing_syms] <- 0L
            }

            # safe addition
            for(sym in names(prior_vec)) {
              counts[[ctx]][sym] <- counts[[ctx]][sym] + prior_vec[sym]
            }
          }
        }
      }

      ltm_idx <- 1
      ctx_ids_prior <- if(length(prior) > 0 && !is.null(prior[[n+1]])) unique(prior[[n+1]]$context_id) else character(0)
      combined_ctx_ids <- unique(c(ctx_ids, ctx_ids_prior))

      for(ctx in combined_ctx_ids) {
        if(is.null(counts[[ctx]])) {
          counts[[ctx]] <- setNames(integer(length(alphabet)), alphabet)
        }
        Ce_full <- setNames(integer(length(alphabet)), alphabet)
        if(length(counts[[ctx]]) > 0) Ce_full[names(counts[[ctx]])] <- counts[[ctx]]
        ltm_list[[ltm_idx]] <- data.table(
          index = -1L,
          context_id = ctx,
          Event = alphabet,
          Ce = as.integer(Ce_full),
          C = sum(Ce_full),
          t = sum(Ce_full > 0L),
          t1 = sum(Ce_full == 1L)
        )
        ltm_idx <- ltm_idx + 1
      }
      ltm_tables[[n+1]] <- rbindlist(ltm_list)
    }
  }

  list(stm = stm_tables, ltm = ltm_tables)
}

#' Expand LTM Count Tables to Timestep-Aligned Order Tables
#'
#' Converts aggregated Long-Term Memory (LTM) count tables into
#' timestep-aligned tables compatible with PPM probability computation.
#'
#' LTM count tables contain counts aggregated over a training corpus
#' (typically with `index = -1`). For prediction on a new sequence `x`,
#' we must derive the counts associated with each timestep's context.
#'
#' The result matches the structure expected by PPM implementations
#' (`ppm_backoff`, `ppm_interpolated`, etc.).
#'
#' @param x Character vector of events.
#' @param N Maximum context order.
#' @param alphabet Character vector of the full symbol alphabet.
#' @param order_counts List of LTM count tables for orders `0..N`.
#'   Each element must be a `data.table` containing:
#'   `context_id`, `Event`, `Ce`, `C`, `t`, `t1`.
#'
#' @return A list of length `N + 1`.
#' Each element is a `data.table` with columns:
#' `index`, `context_id`, `Event`, `Ce`, `C`, `t`, `t1`.
#'
#' The tables contain one row per `(timestep, event)` pair and can be
#' directly passed to probability computation routines.
#'
#' @export
ltm_to_timestep_counts <- function(x, N, alphabet, order_counts) {

  dt_orders <- vector("list", N + 1)

  lag_dt <- lag_matrix(x, N)

  for (n in 0:N) {

    context_cols <- if (n == 0) character(0) else paste0("Lag", n:1)

    # Compute context identifiers
    if (n == 0) {
      lag_dt[, context_id := "ROOT"]
    } else {
      lag_dt[, context_id := do.call(paste, c(.SD, sep = "_")), .SDcols = context_cols]
    }

    # Expand contexts to include all alphabet symbols
    context_dt <- lag_dt[
      , .(Event = alphabet),
      by = .(index, context_id)
    ]

    ltm_dt <- order_counts[[n + 1]]

    # Join LTM counts
    # TODO: possible to optimize here?
    dt_n <- merge(
      context_dt,
      ltm_dt[, .(context_id, Event, Ce, C, t, t1)],
      by = c("context_id", "Event"),
      all.x = TRUE
    )

    # Fill unseen contexts with zero counts
    dt_n[is.na(C), `:=`(
      Ce = 0L,
      C = 0L,
      t = 0L,
      t1 = 0L
    )]

    dt_orders[[n + 1]] <- dt_n[, .(index, context_id, Event, Ce, C, t, t1)]
  }

  dt_orders
}

is_stm <- function(order_counts) {
  # Take the first order table as reference
  dt <- order_counts[[1]]
  # STM if any index > 0
  any(dt$index > 0)
}

