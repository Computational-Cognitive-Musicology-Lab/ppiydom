library(data.table)




lagMatrix <- function(x, N = 10, ...) {
  # x is a sequence of token
  # N is the maximum N-gram length
  # ... can be additional vectors, all the same length as x

  Ntab <- as.data.table(lapply((N - 1):0, \(n) shift(x, n)))

  setnames(Ntab, names(Ntab), paste0('Lag', (N-1):0))

  # Ntab <- cbind(Ntab, as.data.table(list(...)))

  # Ntab[ , index := 1:nrow(Ntab)]

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

  lags <- rev(colnames(countMatrix))


  lapply(1:N,
         \(n) {

           countMatrix[ , (paste0('C', n)) := seq_along(.I) - 1, by = eval(lags[1:n])]
           if (n > 1) {
             countMatrix[ , (paste0('t', n)) := cumsum(c(0, !duplicated(Lag0)))[1:length(Lag0)]]
           } else {
             countMatrix[ , (paste0('t', n)) := cumsum(c(0, !duplicated(Lag0)))[1:length(Lag0)]]
           }


         })
  Cmatrix <- as.matrix(countMatrix[, grepl('^C', colnames(countMatrix)), with = FALSE])
  tmatrix <- as.matrix(countMatrix[, grepl('^t', colnames(countMatrix)), with = FALSE])
  rownames(Cmatrix) <- rownames(tmatrix) <- x

  denominator <- cbind(length(unique(x)),
                     rbind(0, Cmatrix[1:(nrow(Cmatrix) - 1), 1:(N - 1), drop = FALSE]))
  numerator <- Cmatrix

  switch(escape,
         a = {
           numerator <- pmax(numerator, 1)
           denominator <- denominator + 1
         },
         b = {
           numerator <- fifelse(numerator == 0, tmatrix, numerator - 1)
         })

  # denominator[Cmatrix == 0] <- numerator[Cmatrix == 0]
  # When back off is NOT happening, we can make the ratio 1 so the product has no effect.
  zeros <- which(Cmatrix != 0, arr.ind = TRUE)
  zeros[ , 'col'] <- zeros[ , 'col'] - 1
  zeros <- zeros[zeros[ , 'col'] > 0, , drop = FALSE]
  denominator[zeros] <- numerator[zeros]

  apply(numerator / denominator, 1, prod)
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


test0 <- c('I', 'IV', 'I', 'V','I','IV','I',"V",'I','IV','I','V','vi','IV','V','I')
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
