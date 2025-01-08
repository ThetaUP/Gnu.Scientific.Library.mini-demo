#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#define N 2
#define P 3

int main(){
    int i,j,status;
    gsl_vector *u = gsl_vector_alloc(P); // vector
    gsl_vector *res = gsl_vector_alloc(N); // result vector
    gsl_matrix *M = gsl_matrix_alloc(N,P); // matrix
    double alpha_param = 1.0;
    double beta_param = 0.0;
    
    for(i=0;i<P;i++){
        gsl_vector_set(u,i,i);
    }

    for(i=0;i<N;i++){
        for(j=0;j<P;j++){
            gsl_matrix_set(M,i,j,i);
        }
    }

    // perform matrix vector multiplication
    status = gsl_blas_dgemv(CblasNoTrans, alpha_param, M, u, beta_param, res);
    printf("Exitcode: %d\n", status);

    printf("Result: \n");
    for(i=0;i<N;i++){
        printf("%g ", gsl_vector_get(res,i));
    }
    printf("\n");

    return(0);
}



