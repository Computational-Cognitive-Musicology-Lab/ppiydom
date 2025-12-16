library(data.table)




lagMatrix <- function(x, N = 10, ...) {
  # x is a sequence of token
  # N is the maximum N-gram length
  # ... can be additional vectors, all the same length as x

  Ntab <- as.data.table(lapply((N - 1):0, \(n) shift(x, n)))

  colnames(Ntab) <- paste0('Lag', (N-1):0)

  Ntab <- cbind(Ntab, as.data.table(list(...)))

  Ntab[ , index := 1:nrow(Ntab)]

  Ntab[]
}


dynamicModel <- function(x, N = 5, escape = 0, prior = NULL) {
  # x is a sequence of token
  # N is the maximum N-gram length
  # ... can be additional vectors, all the same length as x
  # escape is a prior "escape" probability (a single integer) to be added to all gram counts
  # prior is Null, or a previous model created by dynamicModel

  countMatrix <- lagMatrix(x, N)

  lags <- setdiff(colnames(countMatrix), 'index')

  for (n in 1:N) {

    Ncolname <- paste0('N', n)

    groupby <- tail(lags, n)

    countMatrix <- countMatrix[ , (Ncolname) := seq_along(index) - 1 + escape, by = groupby]


    #
    if (!is.null(prior)) {
      priorFinal <- prior[, list(PriorCount = max(get(Ncolname)) + 1), by = groupby] # + 1 because prior countMatrix doesn't count the LAST time of each thing

      countMatrix <- merge(countMatrix, priorFinal, by = groupby, all.x = TRUE)
      countMatrix[ , (Ncolname) := get(Ncolname) + ifelse(is.na(PriorCount), 0, PriorCount)]
      countMatrix[, PriorCount := NULL]
      #
    }
  }

  setorder(countMatrix, index)

  countMatrix[]
}





test1 <- c('I', 'IV', 'V', 'I', 'I', 'IV', 'I', 'I', 'IV', 'ii', 'V', 'I', 'ii','V','vi','IV','V','I')
test2 <- c('I', 'IV', 'V', 'I', 'I', 'IV', 'V', 'vi', 'I', 'ii', 'IV', 'V', 'I', 'I','IV', 'V', 'I')
test3 <- sample(letters, 1e5, replace = TRUE)



mod1 <- dynamicModel(test1, N = 3)
mod2 <- dynamicModel(test2, N = 3, prior = mod1)
mod12 <- dynamicModel(c(test1, test2), N = 3)
