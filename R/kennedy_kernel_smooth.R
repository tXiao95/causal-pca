#' @description
#' univariate outcome on univariate treatment
#' 
kernel_smooth <- function(Y, Z, bw.range = c(.01, 50), tol = .01){
  n <- length(Y)
  
  kern <- function(x) stats::dnorm(x)
  
  ## hat-value function w(a; bw) for *all* a in Z
  w.fn <- function(bw) {
    w.vec <- numeric(0)
    for (a.val in Z) {
      z.std   <- (Z - a.val) / bw
      k.std   <- kern(z.std) / bw
      m1      <- mean(k.std)
      m2      <- mean(z.std * k.std)
      m3      <- mean(z.std^2 * k.std)
      w.a     <- m3 * kern(0) / bw /
        (m1 * m3 - m2^2)
      w.vec   <- c(w.vec, w.a / n)   # divide by n: hat-matrix diagonal element
    }
    w.vec
  }
  
  ## function returning diagonal of hat matrix at sample points
  hatvals <- function(bw){approx(Z, w.fn(bw), xout=Z)$y }
  
  ## local-linear smoother of pseudo.out at sample points
  cts.eff <- function(out, bw) {
    approx(KernSmooth::locpoly(x=Z, y=out, bandwidth = bw), xout=Z)$y
  }
  
  ## note: choice of bandwidth range depends on problem
  h.opt <- optimize( f = function(h){
    hats <- hatvals(h)
    mean( ((Y - cts.eff(out = Y, bw=h)) / (1-hats))^2)
  }, 
  interval = bw.range, 
  tol = tol)$minimum
  
  ## ------------------------------------------------------------------ ##
  ## 7.  Final dose–response estimate at original Z values
  ## ------------------------------------------------------------------ ##
  est <- approx(locpoly(Z, Y, bandwidth = h.opt), xout = Z)$y
  
  return(list(est=est, h.opt=h.opt))
}
