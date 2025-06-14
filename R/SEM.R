SEM = function(omega) {
  
  Family(
    ngradient = function(y, f, w = 1) {
      omega %*% (y - f)
    },
    loss = function(y, f) {
      sum(omega %*% (y - f)^2) 
    },
    offset = weighted.mean,
    check_y = function(y) {
      if (!is.numeric(y) || !is.null(dim(y)))
        stop("response is not a numeric vector but ", sQuote("family = SEM"))
      y
    },
    name = "Spatial (Durbin) Error Model",
    fW = function(f) rep(1, length(f)),
    response = function(f) f
  )
}