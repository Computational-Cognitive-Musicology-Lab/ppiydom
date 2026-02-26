test_that("scalability benchmark (informational only)", {

  skip_on_cran()
  skip_if_not(Sys.getenv("RUN_SCALABILITY_TESTS") == "true")

  library(ppm)

  # set.seed(123)

  alphabet <- letters
  max_order <- 3

  # sequence lengths to test
  seq_lens <- c(
    1e3,
    1e4,
    1e5,
    1e6
    # , 1e7
  )

  results <- data.frame(
    seq_len = seq_lens,
    ppidyom = NA_real_,
    ppm = NA_real_
  )

  cat("\n--- Scalability test ---\n")

  for (i in seq_along(seq_lens)) {
    seq_len <- seq_lens[i]

    x <- sample(alphabet, seq_len, replace = TRUE)
    seq <- factor(x, levels = alphabet)

    # ---- PPIDYOM ----
    t_ppidyom <- system.time({
      cm <- dynamicModel(x, N = max_order)
      tbl <- compute_ppm_table(cm, max_order, alphabet)
    })["elapsed"]

    rm(cm, tbl)
    gc(full = TRUE)

    # ---- PPM ----
    t_ppm <- system.time({
      mod <- new_ppm_simple(
        order_bound = max_order,
        alphabet_levels = alphabet
      )
      res <- model_seq(mod, seq)
    })["elapsed"]

    rm(mod, res)
    gc(full = TRUE)

    results$ppidyom[i] <- t_ppidyom
    results$ppm[i] <- t_ppm

    cat(
      sprintf(
        "n = %d | PPIDYOM: %.4f s | PPM: %.4f s\n",
        seq_len, t_ppidyom, t_ppm
      )
    )
  }

  # ---- Plot (linear scale, seconds) ----
  x_exp <- log10(results$seq_len)

  matplot(
    x_exp,
    results[, c("ppidyom", "ppm")],
    type = "b",
    pch = c(16, 17),
    lty = 1,
    xaxt = "n",
    xlab = expression("Sequence length (10"^k*")"),
    ylab = "Elapsed time (seconds)",
    main = "Scalability: PPIDYOM vs PPM"
  )

  axis(
    1,
    at = x_exp,
    labels = paste0("10^", x_exp)
  )

  legend(
    "topleft",
    legend = c("PPIDYOM", "PPM"),
    col = 1:2,
    pch = c(16, 17),
    lty = 1,
    bty = "n"
  )
})
