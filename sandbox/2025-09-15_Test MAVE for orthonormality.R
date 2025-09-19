
# We verify whether the MAVE::mave function enforces an orthonorma --------

library(MAVE)

n <- 400

# Create a X matrix with correlation and differing variances
Sigma <- matrix(c(3, 1, 1, 1), nrow = 2)
Sigma <- matrix(c(1, 0, 0, 1), nrow = 2)

mu <- c(1,2)
b1 <- c(2,1) / sqrt( sum(c(2,1)^2))

X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)

eps <- rnorm(n = n)

y <-  X %*% b1 + eps


dr.meanopg <- mave(y ~ X, method = 'meanOPG')
a <- dr.meanopg$dir[[2]]

t(a) %*% a
t(a) %*% solve(Sigma) %*% a
t(a) %*% (Sigma) %*% a


# Conclusion --------------------------------------------------------------

# The mave code is not "whitening" the variables. It does column-wise standardization using
# x_scaled <- scale(x). 

# After performing mave on "y ~ x_scaled", it transforms back into the original scale and then
# normalizes each column. This means there is no guarantee the returned beta is orthonormal. 

# CRAn code here

# mave.compute<-function(x, y, method='CSOPG', max.dim = 10, screen=nrow(x)/log(nrow(x))){
#   
#   method=toupper(method);
#   methodvec=c('CSOPG','CSMAVE','MEANOPG','MEANMAVE','KSIR')
#   
#   if(!(method %in% methodvec)){
#     stop('method should be one of CSMAVE, CSOPG, MEANOPG, MEANMAVE, KSIR')
#   }
#   y = as.matrix(y)
#   x = as.matrix(x)
#   if(nrow(x)!=nrow(y)){
#     stop('the row of x and the row of y is not compatible')
#   }
#   if(ncol(x)==1){
#     stop('x is one dimensional, no need do dimension reduction')
#   }
#   if(screen<ncol(x)){
#     cat("screening method is using to select import variables.")
#   }
#   
#   screen <- min(ncol(x),ceiling(screen))
#   max.dim <- min(max.dim,ncol(x),screen)
#   
#   x.scaled <- scale(x)  #scale x
#   dr <- MAVEfastCpp(x.scaled,y,method,max.dim,screen)
#   
#   dr$x <- x
#   dr$call <- match.call()
#   dr$method <- method
#   dr$dir <- M3d2list(dr$dir,colnames(x))
#   dr$dir <- dr$dir[1:max.dim]
#   len = apply(x,2,sd)
#   for(i in 1:max.dim){
#     dir <- dr$dir[[i]]
#     dir <- dir*matrix(1/len,ncol(x),i) #since x is scaled so we need to transform dir
#     lendir = sqrt(apply(dir^2,2,sum))
#     dr$dir[[i]] <- dir/matrix(lendir,ncol(x),i,byrow=T)
#   }
#   dr$max.dim = max.dim
#   dr$y=y;
#   class(dr) <- 'mave'
#   return(dr)
# }
