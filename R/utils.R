compute_entropy <- function(dt_dist) {
  dt_dist[, .(Entropy = -sum(P * log2(P))), by = index][, Entropy := Entropy]$Entropy
}



test <- c("A", "B", "A", "C", "A", "B", "A", "C", "A")
alphabet <- c("A", "B", "C")
max_order <- 3
dynamic_count_table_for_each_order(x, max_order, alphabet)

result <- ppm_compute(
  x = test,
  alphabet = alphabet,
  N = max_order,
  ppm_type = "interpolation",
  escape_func = escape_C,
  base_type = "unseen"
)

print(result)

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
