# Computes the projection matrix Pi(beta)
Pi <- function(beta) {
  beta %*% solve(t(beta) %*% beta) %*% t(beta)
}

Delta <- function(beta1, beta2, type = "F"){
  Pi1 <- Pi(beta1)
  Pi2 <- Pi(beta2)
  
  return( norm(Pi1 - Pi2, type = type) )
}
