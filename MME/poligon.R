# this script is there to check the dimensions

library(MASS)

mme = function(y, X, Z, A, sigma2_a, sigma2_e) {
    alpha = sigma2_e / sigma2_a
    invA = ginv(A)
    C = rbind(
        cbind(t(X)%*%X, t(X)%*%Z),
        cbind(t(Z)%*%X, t(Z)%*%Z+invA*c(alpha))
        )
    rhs = rbind(t(X)%*%y, t(Z)%*%y)
    invC = ginv(C)
    estimators = invC%*%rhs

    dimC11 <- dim(t(X)%*%X)
    dimC12 <- dim(t(X)%*%Z)
    dimC21 <- dim(t(Z)%*%X)
    dimC22 <- dim(t(Z)%*%Z+invA*c(alpha))
    dimC <- dim(C)

    list('C11' = dimC11,
        'C12' = dimC12, 'C21' = dimC21,
            'C22' = dimC22, 'C' = dimC)
    #list(C = C, est = estimators)
}


N_pheno <- 5
N <- 8
y <- c(4.5, 2.9, 3.9, 3.5, 5.0)
sex <- c(1, 0, 0, 1, 1)
X <- matrix(rep(0,N*N_pheno), nrow = N_pheno, ncol = 2)
X[,1] <- sex
X[,2] <- 1- sex
Z <- matrix(rep(0,N*N_pheno), nrow = N_pheno, ncol = N)
Z[,(N-N_pheno+1):N] <- diag(N_pheno)  # the error is here
A <- matrix(c(1.00,	0.00,	0.00,	0.500,	0.000,	0.50,	0.25,	0.250,
0.00,	1.00,	0.00,	0.000,	0.500,	0.50,	0.25,	0.250,
0.00,	0.00,	1.00,	0.000,	0.500,	0.00,	0.25,	0.500,
0.50,	0.00,	0.00,	1.000,	0.000,	0.25,	0.50,	0.125,
0.00,	0.50,	0.50,	0.000,	1.000,	0.25,	0.50,	0.375,
0.50,	0.50,	0.00,	0.250,	0.250,	1.00,	0.25,	0.500,
0.25,	0.25,	0.25,	0.500,	0.500,	0.25,	1.00,	0.250,
0.25,	0.25,	0.50,	0.125,	0.375,	0.50,	0.25,	1.000), nrow = 8, ncol = 8)

res <- mme(y, X, Z, A, 0.25, 0.25)

res