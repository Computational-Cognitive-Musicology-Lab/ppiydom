compute_entropy <- function(dt_dist) {
  dt_dist[, .(Entropy = -sum(P * log2(P))), by = index][, Entropy := Entropy]$Entropy
}

# ---------------------------
# Example sequence
# ---------------------------
train_seq1 <- c("A", "B", "A", "C", "A", "B", "A", "C", "A")
train_seq2 <- c("A", "B", "C", "A", "B", "C", "A", "B", "C")
test_seq  <- c("B", "A", "B", "C", "A")

alphabet <- c("A", "B", "C")
max_order <- 2

# ---------------------------
# Initialize PPM model
# ---------------------------
ppm_model <- ppm$new(N = max_order, alphabet = alphabet)

# ---------------------------
# Train model incrementally
# ---------------------------
ppm_model$train_sequence(train_seq1)
ppm_model$train_sequence(train_seq2)

# ---------------------------
# Predict on a test sequence
# ---------------------------
result <- ppm_model$predict_sequence(
  x = test_seq,
  model_type = "ltm+",   # use LTM with online update
  ppm_type = "backoff",  # or "interpolation"
  lambda = "C",
  b = 1
)

# ---------------------------
# View results
# ---------------------------
print(result)
