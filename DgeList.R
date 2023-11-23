DgEList <- function (object, method = c("TMM", "TMMwsp", "RLE", "upperquartile", 
                             "none"), refColumn = NULL, logratioTrim = 0.3, sumTrim = 0.05, 
          doWeighting = TRUE, Acutoff = -1e+10, p = 0.75, optimal = TRUE, ...) 
{
  if (!is.null(object$offset)) 
    warning("object contains offsets, which take precedence over library\nsizes and norm factors (and which will not be recomputed).")
  
  object$samples$norm.factors <- mat.or.vec(ncol(object$counts), ncol(object$counts))
  if (is.null(refColumn)) {
    object$samples$norm.factors <- DgNormFactor(object = object$counts, 
                                                     lib.size = object$samples$lib.size, method = method, 
                                                     refColumn = refColumn, logratioTrim = logratioTrim, 
                                                     sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff, 
                                                     p = p, optimal = optimal)
  } else{
  for (r in 1: ncol(object$counts)) {
    
  object$samples$norm.factors[r, ] <- DgNormFactor(object = object$counts, 
                                              lib.size = object$samples$lib.size, method = method, 
                                              refColumn = r, logratioTrim = logratioTrim, 
                                              sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff, 
                                              p = p, optimal = optimal)
  }
  ifelse(object$samples$norm.factors==0,  object$samples$norm.factors <- colMeans(object$samples$norm.factors), object$samples$norm.factors <- geometricmeanCol(object$samples$norm.factors) )
  }
  object
}
