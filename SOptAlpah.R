SOptAlpah <- function(data, eta) {
  # Genes in row and cells are in column
  Vin <- c(0)
  #eta <- 0.01
  #n <- nrow(data) #      number of genes
  #  Select the range of alpha (0.1 < alpha < 0.5)
  # seq(0.49, 0.01, -0.001) for absE
  
  for (alpha in seq(0.2, 0.001, -0.001)) {
    # For complete data
    if (is.matrix(data)) {
      data <- apply(data, 2, sort, decreasing = F)
      n <- nrow(data) #      number of genes
      Y <-
        1 / (n - round(alpha * n)) * colSums(data[(round(alpha * n) + 1):(n - round(alpha * n)) ,])
      
      V <-
        1 / (1 - 2 * alpha) ^ 2 * (1 / n * rowSums(t(data[(round(alpha * n) + 1):(n - round(alpha * n)) ,]) - Y) ^
                                     2
                                   + alpha * (t(data[round(alpha * n) + 1, ]) - Y) ^ 2
                                   + alpha * (t(data[n - round(alpha * n), ]) - Y) ^ 2)
    } else {
      # for a single sample
      data <- sort(data)
      #  n <- 200 #      number of genes
      # n <- 20* round(sqrt(length(data))) for generate datasets 2 code
      n<- length(data)
      Y <-
        1 / (n - 2*round(alpha * n)) * sum(data[(round(alpha * n) + 1):(n - round(alpha * n))])
      
      V <-
        1 / (1 - 2 * alpha) ^ 2 * (1 / n * sum(t(data[(round(alpha * n) + 1):(n - round(alpha * n))]) - Y) ^
                                     2
                                   + alpha * (t(data[round(alpha * n) + 1]) - Y) ^ 2
                                   + alpha * (t(data[n - round(alpha * n)]) - Y) ^ 2)
    }
    
    if (all(abs(V - Vin) < eta)) {
      optalpah <- alpha
      break
    } 
    
    Vin <- V
    
  }
  if (exists("optalpah")) {
    return(optalpah)
  } else {
    return(0.05)
  }
}
