n <- 100
p <- 10
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
betas <- runif((p+1), -5,5)
beta0 <- betas[1]
betas_rest <- betas[2:(p+1)]
Y <- X%*%betas_rest
Y <- Y + beta0
eps <- rnorm(n,sd(Y)*0.1)
Y <- Y + eps

write.table(X, file="X.txt", row.names=FALSE, col.names=FALSE)
write.table(Y, file = "Y.txt", row.names=FALSE, col.names=FALSE)
print('done.')