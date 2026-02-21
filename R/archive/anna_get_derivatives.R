estimate_y_aTb <- function(Y,A,beta,wt,lib.y_aTb=c('SL.glm','SL.ranger'), SL=F){
  
  g <- as.matrix(A)%*%beta # A*b: the dimension-reduced version of A
  
  g <- as.data.frame(g, drop=F)
  
  if (SL){
    # fit model with SuperLearner (it takes too much time, I will use glm instead in selecting dimension)
    y_aTb_fit <- SuperLearner(Y=Y, X=g, family = quasipoisson(), SL.library = lib.y_aTb, obsWeights=wt)
    
    # make prediction: E[Y|A*beta]
    y_aTb <- predict(y_aTb_fit, type = "response")[[1]] %>% as.vector()
    
  }else{
    
    y_aTb_fit <- glm(formula=Y~.,data = data.frame(cbind(Y,g)), weights=wt, family = quasipoisson()) # fit a glm model to estimate E[Y|A*beta]
    
    y_aTb <- predict(y_aTb_fit, type = "response")
    
  }
  
  
  
  
  return(list(y_aTb_fit = y_aTb_fit, y_aTb = y_aTb))
}



get_r_beta <- function(beta, A_sample, p, d, A, C, Y, Y_ac, Yrep, 
                       Sigma, alpha, fA_fAC){
  
  n = length(Y) 
  beta = matrix(beta, nrow=p)
  d = ncol(beta)
  
  # Get U and E[U | A, C]
  y_aTb_estimate <- estimate_y_aTb(Y, A, beta, wt = fA_fAC, c('SL.glm','SL.ranger'))
  y_aTb = y_aTb_estimate$y_aTb
  y_aTb_fit = y_aTb_estimate$y_aTb_fit
  
  # i-th element of U and U_AC is the realization for the i-th row of the data
  
  U <- sweep(alpha, 1, Y - y_aTb, "*") # times each column of alpha with the difference between Y and E[Y|A*beta]
  # result in a matrix of n*(d x p)
  
  U_AC <- sweep(alpha, 1, Y_ac - y_aTb, "*") # times each column of alpha with the difference between Y_ac and E[Y|A*beta]
  # result in a matrix of n*(d x p)
  
  # Get Eq[ E[U|A,C] | C ] = \sum E[U | A_i, C] p(A_i) via summation over sampled A
  samp_size = nrow(A_sample)
  
  alpha_samp <- A_sample
  
  if (d>1){for (i in 2:d) {
    alpha_samp <- cbind(alpha_samp, A_sample^i)
  }}
  
  y_aTb_samp = predict(y_aTb_fit, newdata = as.data.frame(as.matrix(A_sample)%*%beta))[[1]] %>% as.vector()
  
  q_U_AC = matrix(NA, nrow = n, ncol = p*d)
  
  for (i in 1:n){
    y_ac_sample = Yrep[((i-1)*samp_size+1):(i*samp_size)] # extract E[Y|A_i,C], where A_i is the sampled A
    q_U_AC[i,] = t(y_ac_sample - y_aTb_samp)%*%as.matrix(alpha_samp)/samp_size # \sum E[U | A_i, C] p(A_i)
  }
  
  
  # Get r(beat)
  r_beta <- fA_fAC*(U - U_AC) + q_U_AC
  r_beta <- colMeans(r_beta) # sample mean of the IF: E[IF]
  
  return(r_beta)
}

get_first_derivatives <- function(beta, A_sample, p, d, A, C, Y, Y_ac, Yrep,
                                  Sigma, alpha, fA_fAC, delta){
  beta = as.vector(beta)
  K = length(beta)
  
  deriv = matrix(rep(0, K^2), nrow = K)
  
  for (k in 1:K){
    
    cat("Calculating the ", k, "-th/",K," first derivative of beta\n")
    
    e = matrix(rep(0, K), ncol = d) # the p*d matrix of beta
    e[k] = 1 #  indicator for updating the k-th beta
    
    bk_plus = beta + delta*e
    bk_minus = beta - delta*e
    
    r_bk_plus = get_r_beta(bk_plus, A_sample, p, d, A, C, Y, Y_ac, Yrep,
                           Sigma, alpha, fA_fAC)
    r_bk_minus = get_r_beta(bk_minus, A_sample, p, d, A, C, Y, Y_ac, Yrep,
                            Sigma, alpha, fA_fAC)
    first_d = (r_bk_plus - r_bk_minus)/(2*delta) # note that the r_beta is not squared here
    
    deriv[k,] = {r_bk_plus - r_bk_minus}/(2*delta)
  }
  
  return(deriv)
}



# get_first_derivatives <- function(beta, A_sample, p, d, A, C, Y, Y_ac, Yrep, 
#                                   Sigma, alpha, fA_fAC, delta){
#   beta = as.vector(beta)
#   K = length(beta)
#   
#   deriv = c()
#   
#   for (k in 1:K){
#     
#     cat("Calculating the ", k, "-th/",K," first derivative of beta\n")
#     
#     e = matrix(rep(0, K), ncol = d) # the p*d matrix of beta
#     e[k] = 1 #  indicator for updating the k-th beta
#     
#     bk_plus = beta + delta*e 
#     bk_minus = beta - delta*e 
#     
#     r_bk_plus = get_r_beta(bk_plus, A_sample, p, d, A, C, Y, Y_ac, Yrep, 
#                            Sigma, alpha, fA_fAC)
#     r_bk_minus = get_r_beta(bk_minus, A_sample, p, d, A, C, Y, Y_ac, Yrep, 
#                             Sigma, alpha, fA_fAC)
#     first_d = (sum((r_bk_plus)^2) - sum((r_bk_minus)^2))/(2*delta) # note that the r_beta should be squared here as suggested by Ma (2012) paper
#     deriv = c(deriv, first_d)
#   }
#   
#   return(deriv)
# }


# get_second_derivatives <- function(beta, A_sample, p, d, A, C, Y, Y_ac, Yrep, 
#                                    Sigma, alpha, fA_fAC, delta){
#   beta = as.vector(beta) 
#   K = length(beta)
#   
#   second_deriv = matrix(rep(0, K^2), nrow = K)
#   
#   for (k in 1:K){
#     
#     cat("Calculating the ", k, "-th/",K," second derivative of beta\n")
#     
#     e = matrix(rep(0, K), ncol = d)
#     e[k] = 1
#     
#     bk_plus = beta + delta*e
#     bk_minus = beta - delta*e 
#     
#     r_bk_plus = get_first_derivatives(bk_plus, A_sample, p, d, A, C, Y, Y_ac, Yrep, 
#                                       Sigma, alpha, fA_fAC, delta)
#     r_bk_minus = get_first_derivatives(bk_minus, A_sample, p, d, A, C, Y, Y_ac, Yrep, 
#                                        Sigma, alpha, fA_fAC, delta)
#     
#     second_deriv[k,] = {r_bk_plus - r_bk_minus}/(2*delta)
#   }
#   
#   return(second_deriv)
# }

# update beta via the Newton-Raphson method applied to r^2
# update_beta <- function(beta, first_deriv, second_deriv, lambda){
#   
#   vec_beta_new= as.vector(beta) - lambda*solve(second_deriv, tol = 1e-20)%*%first_deriv
#   beta_new = matrix(vec_beta_new, nrow = nrow(beta))
#   
#   return(beta_new)
# }

update_beta <- function(beta, first_deriv, lambda){
  
  r_bk = get_r_beta(as.vector(beta), A_sample, p, d, A, C, Y, Y_ac, Yrep, 
                    Sigma, alpha, fA_fAC)
  
  if(det(first_deriv) < 10^-76){
    
    vec_beta_new= as.vector(beta) - lambda*r_bk%*%ginv(first_deriv)
    
  }else{
    
    vec_beta_new= as.vector(beta) - lambda*r_bk%*%solve(first_deriv)
    
  }
  
  
  
  beta_new = matrix(vec_beta_new, nrow = nrow(beta))
  
  return(beta_new)
  
}


