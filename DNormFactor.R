DgNormFactor <- function (object, lib.size = NULL, method = c("TMM", "TMMwsp", 
                                              "RLE", "upperquartile", "none"), refColumn = NULL, logratioTrim = 0.3, 
          sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10, p = 0.75, optimal = TRUE,
          ...) 
{
  x <- as.matrix(object)
  if (any(is.na(x))) 
    stop("NA counts not permitted")
  nsamples <- ncol(x)
  if (is.null(lib.size)) {
    lib.size <- colSums(x)
  }
  else {
    if (anyNA(lib.size)) 
      stop("NA lib.sizes not permitted")
    if (length(lib.size) != nsamples) {
      if (length(lib.size) > 1L) 
        warning("normLibSizes: length(lib.size) doesn't match number of samples", 
                call. = FALSE)
      lib.size <- rep_len(lib.size, nsamples)
    }
  }
  if (length(method) == 1L && method == "TMMwzp") {
    method <- "TMMwsp"
    message("TMMwzp has been renamed to TMMwsp")
  }
  method <- match.arg(method)
  allzero <- .rowSums(x > 0, nrow(x), nsamples) == 0L
  if (any(allzero)) 
    x <- x[!allzero, , drop = FALSE]
  if (nrow(x) == 0 || nsamples == 1) 
    method = "none"
  f <- switch(method, TMM = {
    if (is.null(refColumn)) {
      f75 <- suppressWarnings(QuantileFactor(data = x, 
                                                  lib.size = lib.size, p = 0.75))
      if (median(f75) < 1e-20) {
        refColumn <- which.max(colSums(sqrt(x)))
      } else {
        refColumn <- which.min(abs(f75 - mean(f75)))
      }
    }
    f <- rep_len(NA_real_, nsamples)
    for (i in 1:nsamples) f[i] <- DWNormFactor(obs = x[, 
                                                         i], ref = x[, refColumn], libsize.obs = lib.size[i], 
                                                 libsize.ref = lib.size[refColumn], logratioTrim = logratioTrim, 
                                                 sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff, optimal = optimal)
    f
  }, TMMwsp = {
    if (is.null(refColumn)) refColumn <- which.max(colSums(sqrt(x)))
    f <- rep_len(NA_real_, nsamples)
    for (i in 1:nsamples) f[i] <- .calcFactorTMMwsp(obs = x[, 
                                                            i], ref = x[, refColumn], libsize.obs = lib.size[i], 
                                                    libsize.ref = lib.size[refColumn], logratioTrim = logratioTrim, 
                                                    sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff)
    f
  }, RLE = .calcFactorRLE(x)/lib.size, upperquartile = .calcFactorQuantile(x, 
                                                                           lib.size, p = p), none = rep_len(1, nsamples))
  f <- f/exp(mean(log(f)))
  names(f) <- colnames(x)
  f
}
