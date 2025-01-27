// this is unfinished

#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

int main(){
    
    // define some input parameters
    int N_total = 8; // total number of cows
    int N_pheno = 5; // number of phenotyped cows
    int N_fixed_effects = 2; // number of fixed effects

    // define input matrix/vector files
    FILE *y_f,*X_f,*A_f,*Z_f,*variance_components_f;
    y_f = fopen("y_vec.txt", "r");
    X_f = fopen("X_mat.txt", "r");
    A_f = fopen("A_mat.txt", "r");
    Z_f = fopen("Z_mat.txt", "r");
    variance_components_f = fopen("variance_components.txt", "r");

    // allocate space for input matrices and vectors
    gsl_vector *y = gsl_vector_alloc(N_pheno);
    gsl_matrix *X = gsl_matrix_alloc(N_pheno, N_fixed_effects);
    gsl_matrix *A = gsl_matrix_alloc(N_total, N_total);
    gsl_matrix *Z = gsl_matrix_alloc(N_pheno, N_total);
    
    
    // initialiaze vectors and matrices
    double sigma_a, sigma_e;
    fscanf(variance_components_f, "%lf %lf", &sigma_a, &sigma_e);
    gsl_matrix_fscanf(X_f, X);
    gsl_matrix_fscanf(A_f, A);
    gsl_matrix_fscanf(Z_f, Z);
    gsl_vector_fread(y_f, y);

    //allocate space for C matrix
    int dim_C = N_total + N_fixed_effects;
    gsl_matrix *C = gsl_matrix_alloc(dim_C, dim_C);

    // allocate space for the submratrices of C
    gsl_matrix *C11 = gsl_matrix_alloc(N_fixed_effects, N_fixed_effects);
    gsl_matrix *C12 = gsl_matrix_alloc(N_fixed_effects, N_total);
    gsl_matrix *C21 = gsl_matrix_alloc(N_total, N_fixed_effects);
    gsl_matrix *C22 = gsl_matrix_alloc(N_total, N_total);

    // perform the arithemtics
    /*
        t(X)%*%X -->  for C11   (DONE)
        t(X)%*%Z -->  for C22    (DONE)
        t(Z)%*%X --> for C21     (DONE)
        t(Z)%*%Z+invA*c(alpha) -->  for C22
    */
    gsl_matrix *X_T = gsl_matrix_alloc(N_fixed_effects, N_pheno);
    int status;
    status = gsl_matrix_transpose_memcpy(X_T, X);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X_T, X, 0.0, C11); //t(X)%*%X -->  for C11

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X_T, Z, 0.0, C12); //t(X)%*%Z -->  for C12

    gsl_matrix *Z_T = gsl_matrix_alloc(N_total, N_pheno);
    status = gsl_matrix_transpose_memcpy(Z_T, Z);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Z_T, X, 0.0, C21); //t(Z)%*%X -->  for C21
    

    // populate the C matrix
    /*
        Do it as shown in the pdf which is in the same directory.
    */

        

    return(0);
}