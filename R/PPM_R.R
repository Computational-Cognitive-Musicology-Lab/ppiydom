library(data.table)




lagMatrix <- function(x, N = 10, ...) {
  # x is a sequence of token
  # N is the maximum N-gram length
  # ... can be additional vectors, all the same length as x

  Ntab <- as.data.table(lapply((N - 1):0, \(n) shift(x, n)))

  setnames(Ntab, names(Ntab), paste0('Lag', (N-1):0))

  # Ntab <- cbind(Ntab, as.data.table(list(...)))

  Ntab[ , index := 1:nrow(Ntab)]

  Ntab[]
}


dynamicModel <- function(x, N = 5, escape = c('a', 'b', 'c', 'd'), prior = NULL) {
  # x is a sequence of token
  # N is the maximum N-gram length
  # ... can be additional vectors, all the same length as x
  # escape is a prior "escape" probability (a single integer) to be added to all gram counts
  # prior is Null, or a previous model created by dynamicModel

  escape <- match.arg(escape)

  countMatrix <- lagMatrix(x, N)

  lags <- setdiff(colnames(countMatrix), 'index')

  countMatrix[ , `C-1` := seq_along(index)  + if (!is.null(prior)) nrow(prior) else 0]

  for (n in 1:N) {

    Ccol <- paste0('C', n - 1)
    Tcol <- paste0('t', n - 1)
    groupby <- tail(lags, n)

    countMatrix <- countMatrix[ , (Ccol) := seq_along(index) - 1, by = groupby]
    if (n > 1)  {
      countMatrix[ , (Tcol) := c(0, head(cumsum(get(Ccol) == 0), -1)), by = eval(head(groupby, -1))]
    } else {
      countMatrix[ , (Tcol) := c(0, head(cumsum(get(Ccol) == 0), -1))]
    }
    countMatrix <- doescape(escape, countMatrix, n)

    #
    if (!is.null(prior)) {
      priorFinal <- prior[, list(PriorCount = max(get(Ncolname)) + 1), by = groupby] # + 1 because prior countMatrix doesn't count the LAST time of each thing


      countMatrix <- merge(countMatrix, priorFinal, by = groupby, all.x = TRUE)
      countMatrix[ , (Ncolname) := get(Ncolname) + ifelse(is.na(PriorCount), 0, PriorCount)]
      countMatrix[, PriorCount := NULL]
      #
    }
  }

  setnames(countMatrix, 'Lag0', 'Seq')

  setorder(countMatrix, index)
  countMatrix[ , (paste0('Lag', (N - 1):1)) := NULL]

  countMatrix[]
}


doescape <- function(escape = 'a', countMatrix, n) {
  C <- countMatrix[[paste0('C', n - 1)]]
  c <- countMatrix[[paste0('C', n)]]
  t <- countMatrix[[paste0('t', n)]]
  switch(escape,
         a = {
           c[c == 0] <- 1
           C <- C + 1
         },
         b = {
           c <- ifelse(c == 0, t, c - 1)
         },
         c = {
           c[c == 0] <- t[c == 0]
           C <- C + t
         },
         d = {
           c <- ifelse(c == 0, .5 * t, c - .5)
         })

  countMatrix[ , (paste0('C', n - 1)) := C]
  countMatrix[ , (paste0('C', n)) := c]
  countMatrix

}



test1 <- c('I', 'IV', 'V', 'I', 'I', 'IV', 'I', 'I', 'IV', 'ii', 'V', 'I', 'ii','V','vi','IV','V','I')
test2 <- c('I', 'IV', 'V', 'I', 'I', 'IV', 'V', 'vi', 'I', 'ii', 'IV', 'V', 'I', 'I','IV', 'V', 'I')
test3 <- sample(letters, 1e5, replace = TRUE)



# mod <- dynamicModel(test1, N=0)
# mod1 <- dynamicModel(test1, N = 3)
# mod2 <- dynamicModel(test2, N = 3, prior = mod1)
# mod12 <- dynamicModel(c(test1, test2), N = 3)
#
# ##
## escape probs----



## Cummulative entropy

cumH <- function(counts, base = 2) {

  l <- counts * log(counts, base)
  l_0 <- ifelse(counts <= 1L, 0, (counts - 1L) * log(counts - 1L, base))

  index <- seq_along(counts)

  -((cumsum(l - l_0) - index*log(index, base)) / index)
}
