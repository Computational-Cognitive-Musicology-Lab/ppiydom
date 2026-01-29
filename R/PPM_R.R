#' @import data.table
library(data.table)

# ===========COUNT MATRIX===========

lagMatrix <- function(x, N = 10, ...) {
  # x is a sequence of token
  # N is the maximum N-gram length
  # ... can be additional vectors, all the same length as x

  Ntab <- as.data.table(lapply(N:0, \(n) shift(x, n)))

  colnames(Ntab) <- paste0('Lag', N:0)

  Ntab <- cbind(Ntab, as.data.table(list(...)))

  Ntab[ , index := 1:nrow(Ntab)]

  Ntab[]
}


dynamicModel <- function(x, N = 5, prior = NULL) {
  # x is a sequence of token
  # N is the maximum N-gram length
  # ... can be additional vectors, all the same length as x
  # prior is Null, or a previous model created by dynamicModel

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
compute_probs_for_n <- function(countMatrix, N, escape_method = escape_C) {
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

    probMatrix[valid,
              `:=`(
                prob = pmax(event_count - esc_info$subtract, 0) / esc_info$denom,
                esc  = esc_info$esc
              )
    ]
  }

  probMatrix[]
}

compute_ppm_table <- function(countMatrix,
                              max_order,
                              alphabet,
                              escape_method = escape_C) {
  dt <- copy(countMatrix)

  dt[, `:=`(
    P_final = NA_real_,
    order_used = NA_integer_,
    p_mass = 1,
    resolved = FALSE
  )]

  # compute all orders
  for (n in max_order:1) {
    probs_n <- compute_probs_for_n(dt, n, escape_method)
    # print(probs_n)

    valid <- !dt$resolved & probs_n$C >0
    event_seen <- valid & probs_n$prob > 0
    event_unseen <- valid & probs_n$prob <= 0

    # resolve event that are seen after the context with their event probability
    dt[event_seen, `:=`(
      P_final = p_mass * probs_n$prob[event_seen],
      order_used = n,
      resolved = TRUE
    )]

    # update escape mass for unseen events after the context (still unresolved)
    dt[event_unseen,
       p_mass := p_mass * probs_n$esc[event_unseen]
    ]

    # print(dt)
  }

  # 0-gram fall back
  fallback <- !dt$resolved
  dt[fallback, `:=`(
    P_final = p_mass * (1 / length(alphabet)),
    order_used = 0,
    resolved = TRUE
  )]

  # IC & Entropy
  dt[, `:=`(
    IC = -log2(P_final),
    Entropy = NA_real_
  )]

  dt[, c(
      .(index = index),
      mget(paste0("Lag", max_order:1)),
      .(
        Event = Lag0,
        Order = order_used,
        P = P_final,
        IC = IC,
        Entropy = Entropy
      )
    )]
}





# =========TEST==========

max_order <- 3
alphabet <- c('I', 'ii', 'III', 'IV', 'V', 'vi', 'vii')
test <- c('IV', 'I', 'I', 'IV', 'I', 'IV', 'I', 'V', 'I', 'IV', 'V')
countMatrix <- dynamicModel(test, N=max_order)
countMatrix
ppm_table <- compute_ppm_table(countMatrix, 3, alphabet, escape_method = escape_A)
ppm_table

