get_fA_fAC <- function(C, A, Y, dat, method, samp_size){
  
  n = dim(C)[1]
  p = dim(A)[2]
  k = dim(C)[2]
  
  #------------------------------------------
  # Pick samples from p(A) used for calculate integration regarding p(A|C) later on
  #------------------------------------------
  idx = sample(1:n, samp_size, replace = TRUE)
  A_sample = A[idx, ]
  
  #------------------------------------------
  # calculate p(A)
  #------------------------------------------
  # p(A) can be any valid distribution of A
  # pA = rep(1/n, n)
  pA = rgamma(n, shape=1, rate=10^3)
  
  cat('summary of pA: ',summary(pA))
  
  #------------------------------------------
  # calculate p(A|C) via PCA
  #------------------------------------------
  
  # Perform PCA on A 
  pca_A <- prcomp(A, scale. = TRUE)
  summary(pca_A)
  
  # Retain the top k principal components
  k.A <- 3  # Choose a reasonable number of components
  A_reduced <- pca_A$x[, 1:k.A,drop=F]; colnames(A_reduced) <- paste0('A_PC.',1:k.A)
  
  # Use the reduced dimensions for KDE
  AC_reduced <- cbind(A_reduced, C)
  
  #------------------------------------------
  # calculate p(A|C) via fit_density
  #------------------------------------------
  
  if (method=='pca_condensier'){
    
    dens_fit <- fit_density(
      X = colnames(C),
      Y = colnames(A_reduced),
      input_data = as.data.frame(AC_reduced),
      nbins = 20,
      bin_method = "equal.len",
      bin_estimator = speedglmR6$new())
    
    pA_C <- predict_probability(dens_fit, newdata=as.data.frame(AC_reduced))
    
    cat('summary of pA_C: ',summary(pA_C))
    
  }
  
  #------------------------------------------
  # calculate p(A|C) via np
  #------------------------------------------
  
  if (method=='pca_np'){
    
    # Construct the formula using paste0
    formula_str <- paste0(paste(colnames(A_reduced), collapse = " + "), " ~ ", paste(colnames(C_reduced), collapse = " + "))
    
    # Create the formula
    formula <- as.formula(formula_str)
    
    cat("\n ( ****** fit p(A|C) ****** ) \n")
    start_time <- Sys.time()
    
    bw <- npcdensbw(formula=formula, data=as.data.frame(AC_reduced))
    pA_C_fit <- npcdens(bws=bw)
    
    end_time <- Sys.time()
    
    cat("\n ( ****** fit p(A|C) done in ", end_time - start_time, " ****** ) \n")
    
    # prediction
    pA_C <- predict(pA_C_fit)
    
  }
  
  
  
  #------------------------------------------
  # calculate p(A|C) via assuming normal distribution
  # #------------------------------------------
  
  if (method=='dmvnorm'){
    
    model_errors <- matrix(NA, nrow = n, ncol = p)
    predict.A_C <- matrix(NA, nrow = n, ncol = p)
    
    for (i in 1:p){ # loop over each variable in A
      
      cat("\n ( ****** predict at A #", i, "****** ) \n")
      
      model <- SuperLearner(Y=A[, i], X=as.data.frame(C), family = gaussian(), SL.library = c('SL.glm','SL.ranger'))
      predicted <- predict(model, type = "response")[[1]] %>% as.vector()
      model_errors[,i] <- A[,i] - predicted
      
      # model <- lm(A[, i] ~ . , data=as.data.frame(C)) # fit a linear model for each variable in A
      # model_errors[, i] <- A[,i] - predict(model) # errors of the model
      
      predict.A_C[,i] <- predict(model) # store the prediction result for E[A|C]
      
    } # end of for loop over variables
    
    # compute the variance-covariance matrix of the errors
    varcov <- cov(data.frame(model_errors))
    # varcov <- cov.shrink(data.frame(model_errors))
    
    # predict p(A|C) assuming is follows a normal distribution
    pA_C = rep(0, n)
    for (j in 1:n){ pA_C[j] <- dmvnorm(x = A[j, ], mean = predict.A_C[j,], sigma = varcov)}
    
  }
  
  #------------------------------------------
  # other method stop!
  # #------------------------------------------
  
  if (method %in% c('pca_condensier', 'pca_np', 'dmvnorm') == FALSE){
    
    stop("method for propensity not supported! select among 'pca_condensier', 'pca_np', 'dmvnorm'")
    
  }
  
  
  #------------------------------------------
  # calculate p(A)/p(A|C)
  #------------------------------------------
  pA_C[pA_C<=0.001] <- 0.001
  
  fA_fAC = pA/pA_C
  
  ## remove the outliers in this weight ratio by truncating applied by identifying lower bound = Q1 - 1.5 * IQR & upper bound = Q3 + 1.5 * IQR
  bp_stats = boxplot.stats(fA_fAC)
  
  # Identify lower and upper bounds for outliers
  lower_bound <- bp_stats$stats[1]
  upper_bound <- bp_stats$stats[5]
  
  fA_fAC[fA_fAC <= lower_bound ] <- lower_bound
  fA_fAC[fA_fAC >= upper_bound ] <- upper_bound
  
  cat('\n pA/pA_C: ', summary(round(fA_fAC, 3)), '\n')
  
  return(list(fA_fAC=fA_fAC, 
              A_sample=A_sample,
              pA_C=pA_C,
              pA=pA,
              idx=idx)) 
  
}

