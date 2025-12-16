x <- sample(letters[1:7], 5e5, replace= T, prob = c(20,10,5,4,3,2,1))
# library(humdrumR)
# readHumdrum('~/Bridge/Research/Data/Humdrum/Kern/JSBach/371chorales/*.krn') -> chor

# with(chor, kern(Token, generic=T,simple=T)) -> x

# file <- chor$File

library(data.table)


d <- as.data.table(setNames(lapply(5:0, \(n) shift(x, n)), paste0('N', 5:0)))
d[ , I := seq_len(nrow(d))]

anachron <- function(n, data) {
  j <- paste0('N', n:0)

  y <- data

  y <- y[ , list(Count = length(I)), j]

  j <- head(j, -1)
  y[ , {
    total <- sum(Count)
    list(Count = Count, Total = total, p = Count / total, N0 = N0)
  },  by = j] -> y


  y
}

buildN <- function(x, N = 10, ...) {
  Ntab <- as.data.table(lapply((N-1):0, \(n) shift(x, n)))

  colnames(Ntab) <- c(paste0('Lag-', (N-1):1), 'Original')

  Ntab <- cbind(Ntab, as.data.table(list(...)))

  Ntab[ , I := 1:nrow(Ntab)]

  Ntab[]
}




model <- function(x, N = 5, ..., escape = 0, prior = NULL) {

  Nt <- buildN(x, N, ...)


  seqs <- colnames(Nt)[grepl('Lag-', colnames(Nt))]

  # enumerate all unique N grams
  for (n in 0:(N-1)) {
    cur <- tail(seqs, n)

    #
    counts <- Nt[ , list(I, Count = seq_along(I)  - 1 + escape), by = c(cur, 'Original')][ , c('I', "Count"), with = FALSE]
    colnames(counts) <- c('I', paste0('N', n))
    Nt <- counts[Nt, on = 'I']

    #
    # conditionCounts <- Nt[ , list(I, Condition = seq_along(I) -1 + escape), by = cur][ , c('I', "Condition"), with = FALSE]
    # colnames(conditionCounts) <- c('I', condition)
    # Nt <- conditionCounts[Nt, on = 'I']

    #
    if (!is.null(prior)) {
      jointPrior <- unique(prior, by = c(cur, 'Sequence'), fromLast = TRUE)
      Nt[jointPrior, on = c(cur, 'Sequence'), (joint) := get(joint) + get(paste0('i.', joint))]

      # conditionPrior <- unique(prior, by = cur, fromLast = TRUE)
      # Nt[conditionPrior, on = cur, (condition) := get(condition) + get(paste0('i.', condition))]

    }
  }

  setorder(Nt, I)

  Nt[]
}

pullFinal <- function(Nt) {
  seqs <- colnames(Nt)[grepl('Seq-', colnames(Nt))]

  Nt <- unique(Nt, by = seqs, fromLast = TRUE)

  setorderv(Nt, seqs)
  Nt


}



test <- c('I', 'IV', 'V', 'I', 'I', 'IV', 'I', 'I', 'IV', 'ii', 'V', 'I', 'ii','V','vi','IV','V','I')

test2 <- c('I', 'IV', 'V', 'I', 'I', 'IV', 'V', 'vi', 'I', 'ii', 'IV', 'V', 'I', 'I','IV', 'V', 'I')


# Joint / Condition or (Joint - 1) / (Condition - 1)

nth <- function(x) {
  na <- is.na(x)
  .x <- x[!na]
  i <- tapply(seq_along(.x), .x, force)
  n <- tapply(.x, .x, seq_along)
  .x <- unlist(n)[order(unlist(i))]

  output <- rep(NA_integer_, length(x))
  output[!na] <- .x
  output

}


ngrams <- function(..., maxN = 10) {

  vecs <- list(...)
  if (is.null(names(vecs)) || any(names(vecs) == '')) stop('Vecs must all be named')

  tables <-  Map(\(vec, name) {
                     table <- as.data.table(lapply(0:(maxN - 1), \(n) shift(vec, n)))
                     colnames(table) <- as.character(0:(maxN - 1))
                     table$I <- seq_len(nrow(table))
                     table
                     },
                 vecs, names(vecs))


  table <- do.call('cbind', tables)
  class(tables) <- c('ngram.table', class(tables))
  tables
}


sums <- function(tables, var, N = 2) {
  table <- do.call('cbind', tables[var])
  byjoint <- unlist(lapply(var, \(v) paste0(var, '.', 0:(N - 1))))
  bycond <- unlist(lapply(var, \(v) paste0(var, '.', 1:(N - 1))))

  table <- table[ , Joint := seq_along(.I), by = byjoint]
  table <- table[ , Condition := seq_along(.I), by = bycond]
  table[]
}
