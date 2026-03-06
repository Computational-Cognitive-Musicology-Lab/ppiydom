library(data.table)


ppm <- setRefClass(
  "ppm",

  fields = list(
    N = "numeric",
    alphabet = "character",
    counts_ltm = "list"
  ),

  methods = list(

    #' Initialize a new PPM counter
    initialize = function(N, alphabet) {
      .self$N <- N
      .self$alphabet <- alphabet
      .self$counts_ltm <- list()
    },

    #' Train on a single sequence (incremental LTM update)
    #' @param x Character vector sequence to train
    train_sequence = function(x) {
      count_tables <- count_tables(
        x = x,
        N = .self$N,
        alphabet = .self$alphabet,
        model_type="ltm",
        prior = .self$counts_ltm
      )
      # Update LTM
      .self$counts_ltm <- count_tables$ltm
      invisible(NULL)
    },

    #' Predict IC and entropy for a sequence
    #' @param x Character vector sequence
    #' @param model_type Model type: "stm", "ltm", "both", "ltm+", "both+"
    #' @param ppm_type "interpolation" or "backoff"
    #' @param lambda escape or discount function (default = "C")
    #' @param b Bias parameter for relative-entropy weighting (used for + models)
    #' @return data.table with columns: index, Event, P, IC, Entropy
    predict_sequence = function(x, model_type = c("stm", "ltm", "both", "ltm+", "both+"),
                                ppm_type = c("interpolation", "backoff"),
                                lambda = "C",
                                b = 1) {
      model_type <- match.arg(model_type)
      ppm_type <- match.arg(ppm_type)
      T <- length(x)
      alpha_len <- length(.self$alphabet)

      # -------------------------
      # lambda function
      # -------------------------
      escape_func <- escape_functions[[lambda]]
      if(is.null(escape_func)) stop("Unknown escape lambda: ", lambda)
      discount_func <- discount_functions[[lambda]]
      if(is.null(discount_func)) stop("Unknown discount lambda: ", lambda)


      # ------------------------------------------------
      # Build count tables
      # ------------------------------------------------
      # Determine which counts to compute
      if (model_type == "stm") {
        count_type <- "stm"
      } else if (grepl("ltm", model_type) && !grepl("both", model_type)) {
        count_type <- "ltm"
      } else if (grepl("both", model_type)) {
        count_type <- "both"
      } else {
        stop("Invalid model_type")
      }

      counts <- count_tables(
        x = x,
        N = .self$N,
        alphabet = .self$alphabet,
        model_type = count_type,
        prior = .self$counts_ltm
      )

      # ------------------------------------------------
      # STM probabilities
      # ------------------------------------------------

      P_stm <- NULL
      if(model_type %in% c("stm","both","both+")) {

        P_stm <- if(ppm_type == "interpolation")
          ppm_interpolated(x, .self$N, .self$alphabet, counts$stm, discount_func=discount_func)
        else
          ppm_backoff(x, .self$N, .self$alphabet, counts$stm, escape_func=escape_func)

      }

      # ------------------------------------------------
      # LTM probabilities
      # ------------------------------------------------

      P_ltm <- NULL
      if(model_type %in% c("ltm","both","ltm+","both+")) {

        P_ltm <- if(ppm_type == "interpolation")
          ppm_interpolated(x, .self$N, .self$alphabet, counts$ltm, discount_func=discount_func)
        else
          ppm_backoff(x, .self$N, .self$alphabet, counts$ltm, escape_func=escape_func)
      }

      # ------------------------------------------------
      # Model selection
      # ------------------------------------------------

      if(model_type == "stm")
        result <- P_stm
      else if(model_type == "ltm" || model_type == "ltm+")
        result <- P_ltm
      else if(model_type %in% c("both","both+"))
        result <- combine_models(P_stm, P_ltm, b)

      # ------------------------------------------------
      # Online learning (+ models)
      # ------------------------------------------------

      if(model_type %in% c("ltm+","both+")) {
        .self$counts_ltm <- counts$ltm
      }

      result
    }
  )

)


combine_models <- function(p_stm, p_ltm, b=1) {

  dt <- merge(
    p_stm[, .(index, Event, P_stm=P, H_stm=Entropy)],
    p_ltm[, .(index, Event, P_ltm=P, H_ltm=Entropy)],
    by=c("index","Event")
  )

  dt[, w_stm := 2^(-b * H_stm)]
  dt[, w_ltm := 2^(-b * H_ltm)]

  dt[, norm := w_stm + w_ltm]

  dt[, w_stm := w_stm / norm]
  dt[, w_ltm := w_ltm / norm]

  dt[, P := w_stm * P_stm + w_ltm * P_ltm]

  dt[, IC := -log2(P)]

  dt[, .(index, Event, P, IC)]
}

