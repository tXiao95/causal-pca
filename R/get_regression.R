
get_regression <- function(C, A, Y, dat, method ,samp_size, A_sample, verbose=T){
  
  n = dim(C)[1]
  p = dim(A)[2]
  
  Crep = C[rep(seq_len(n), each=samp_size), ]
  Arep = coredata(A_sample)[rep(seq(nrow(A_sample)), n), ]
  
  rownames(Crep) = 1:nrow(Crep)
  rownames(Arep) = 1:nrow(Arep)
  
  datY = data.frame(A, C, Y)
  
  if (method=='SL'){
    
    model <- SuperLearner(Y=datY[,'Y'], X=datY[,c(names(A), colnames(C))], family = quasipoisson(), SL.library = c('SL.glm','SL.ranger'))
    Y_ac <- predict(model, type = "response")[[1]] %>% as.vector()
    
    # Estimate E[Y | Asamp, C] 
    datACrep = data.frame(Arep, Crep); names(datACrep) = c(colnames(A),colnames(C))
    Yrep = predict(model, newdata=datACrep, type = "response")[[1]] %>% as.vector()
    
  }else if (method=='GLM'){
    
    model <- glm("Y~.", data=datY[,c("Y",names(A), colnames(C))], family = quasipoisson())
    Y_ac <- predict(model, type = "response")
    
    # Estimate E[Y | Asamp, C] 
    datACrep = data.frame(Arep, Crep); names(datACrep) = c(colnames(A),colnames(C))
    Yrep = predict(model, newdata=datACrep, type = "response")
    
  }
  
  
  if(verbose){
    cat("outcome regression fitted with GLM")
    cat("\n ( ****** E[Y | A, C] ****** ) \n")
    cat("summary of E[Y|A,C]: ", summary(Y_ac), "\n")
  }
  
  return(list(Y_ac = Y_ac, 
              Yrep = Yrep))
}

