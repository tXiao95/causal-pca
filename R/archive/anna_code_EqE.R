#------------------------------------------
# Pick samples from p(A) used for calculate integration regarding p(A|C) later on
#------------------------------------------
idx = sample(1:n, samp_size, replace = TRUE) # samp_size is the number of A samples draw to approximate integration
A_sample = A[idx, ]

#------------------------------------------
# Outcome regression
#------------------------------------------
model <- glm("Y~.", data=datY[,c("Y",names(A), colnames(C))], family = quasipoisson())
Y_ac <- predict(model, type = "response")

data_sample <- crossing(A_sample, dat[,colnames(C)])
Yrep = predict(model, newdata=data_sample, type = "response") # Outcome regression prediction for sampled A

#------------------------------------------
# E(Y|g(A*beta))
#------------------------------------------
g <- as.matrix(A)%*%beta # A*b: the dimension-reduced version of A  
g <- as.data.frame(g, drop=F)

y_aTb_fit <- glm(formula=Y~.,data = data.frame(cbind(Y,g)), weights=wt, family = quasipoisson())
y_aTb_samp = predict(y_aTb_fit, newdata = as.data.frame(as.matrix(A_sample)%*%beta))[[1]] %>% as.vector()

#------------------------------------------
# Eq(E(U|A,C)|C)
#------------------------------------------

# choose alpha(A) as (A,A^2,...)
# alpha_samp is alpha(A_sampled)
alpha_samp <- A_sample

if (d>1){for (i in 2:d) {
  alpha_samp <- cbind(alpha_samp, A_sample^i)
}}

q_U_AC = matrix(NA, nrow = n, ncol = p*d)

for (i in 1:n){
  y_ac_sample = Yrep[((i-1)*samp_size+1):(i*samp_size)] # extract E[Y|A_i,C], where A_i is the sampled A
  q_U_AC[i,] = t(y_ac_sample - y_aTb_samp)%*%as.matrix(alpha_samp)/samp_size # \sum E[U | A_i, C] p(A_i)
}