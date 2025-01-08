#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

int main(){
    int i,j,N,P;
    double val;
    N = 100;   // number of samples
    P = 10;    // number of features
    FILE *input_X_file;
    FILE *input_Y_file;
    FILE *estimates_file;
    gsl_matrix *X_orig = gsl_matrix_alloc(N,P);
    
    input_X_file = fopen("X.txt", "r");
    gsl_matrix_fscanf(input_X_file, X_orig); // X without the intercept column
    P = P + 1;  // for the intercept
    gsl_matrix *X = gsl_matrix_alloc(N,P);  // X with the intercept column

    // add 1st col of 1s to X for the intercept
    for(i=0;i<N;i++){
        for(j=0;j<P;j++){
            if(j==0){
                gsl_matrix_set(X, i, j, 1.0); // for the intercept
            }
            else{
                val = gsl_matrix_get(X_orig, i, j-1);
                gsl_matrix_set(X, i, j, val);
            }
        }
    }

    // allocate space
    gsl_vector *Y = gsl_vector_alloc(N);
    gsl_matrix *X_T = gsl_matrix_alloc(P,N);      // t(X)
    gsl_matrix *mat_to_inv = gsl_matrix_alloc(P,P);  // t( X ) @ X
    gsl_matrix *inverted_mat = gsl_matrix_alloc(P,P); // inv(t( X ) @ X)
    gsl_matrix *inverted_mat_X = gsl_matrix_alloc(P,N); // inv(t( X ) @ X) @ X
    gsl_vector *results = gsl_vector_alloc(P); // inv(t( X ) @ X) @ X @ Y
    
    // open files
    input_Y_file = fopen("Y.txt", "r");
    estimates_file = fopen("OLS_estimates.txt", "w");

    // read data
    gsl_vector_fread(input_Y_file, Y);

    // perform matrix manipulations
    gsl_matrix_transpose_memcpy(X_T, X);  // t(X)
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X_T, X, 0.0, mat_to_inv); // t(X)@X
    gsl_matrix *LU_decomp_MAT = gsl_matrix_alloc(P,P); // for inv()
    gsl_matrix_memcpy(LU_decomp_MAT, mat_to_inv); // for inv()
    int sign; // for inv()
    gsl_permutation * perm_vec = gsl_permutation_alloc(P); // for inv()
    gsl_linalg_LU_decomp(LU_decomp_MAT, perm_vec, &sign);  // for inv()
    gsl_linalg_LU_invert(LU_decomp_MAT, perm_vec, inverted_mat); // inv(t( X ) @ X)
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, inverted_mat, X_T, 0.0, inverted_mat_X); // inv(t( X ) @ X) @ X
    gsl_blas_dgemv(CblasNoTrans, 1.0, inverted_mat_X, Y, 0.0, results);


    for(i=0;i<P;i++){
        val = gsl_vector_get(results, i);
        printf("%g\n", val);
    }


    return(0);
}