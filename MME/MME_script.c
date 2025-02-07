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
        t(Z)%*%Z+invA*c(alpha) -->  for C22 (DONE)
    */
    gsl_matrix *X_T = gsl_matrix_alloc(N_fixed_effects, N_pheno);
    int status;
    status = gsl_matrix_transpose_memcpy(X_T, X);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X_T, X, 0.0, C11); //t(X)%*%X -->  for C11

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X_T, Z, 0.0, C12); //t(X)%*%Z -->  for C12

    gsl_matrix *Z_T = gsl_matrix_alloc(N_total, N_pheno);
    status = gsl_matrix_transpose_memcpy(Z_T, Z);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Z_T, X, 0.0, C21); //t(Z)%*%X -->  for C21

    gsl_matrix *ZTZ = gsl_matrix_alloc(N_total, N_total);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Z_T, Z, 0.0, ZTZ);  //t(Z)%*%Z
    
    // invert the A matrix
    gsl_matrix *LU_decomp_MAT = gsl_matrix_alloc(N_total,N_total);
    gsl_matrix_memcpy(LU_decomp_MAT, A);
    int sign;
    gsl_permutation * perm_vec = gsl_permutation_alloc(N_total);
    gsl_linalg_LU_decomp(LU_decomp_MAT, perm_vec, &sign);
    gsl_matrix *A_inv = gsl_matrix_alloc(N_total, N_total);
    gsl_linalg_LU_invert(LU_decomp_MAT, perm_vec, A_inv); 

    // multiply A_inv by the scalar sigma_a
    int i,j;
    for(i=0;i<N_total;i++){
        for(j=0;j<N_total;j++){
            gsl_matrix_set(A_inv, i, j, gsl_matrix_get(A_inv, i, j)*sigma_a);
        }
    }

    // C22 = t(Z)%*%Z %*% A_inv*sigma_a
    gsl_matrix_add(ZTZ, A_inv); // the results of this matrix addition is stored in ZTZ, we will rename it for clarity
    for(i=0;i<N_total;i++){
        for(j=0;j<N_total;j++){
            gsl_matrix_set(C22,i,j,gsl_matrix_get(ZTZ,i,j));
        }
    }

    int k,l;
    // populate the C matrix
        // populate the C11 region
    k = 0;
    l = 0;
    for(i=0;i<N_fixed_effects;i++){
        for(j=0;j<N_fixed_effects;j++){
            gsl_matrix_set(C,i,j,gsl_matrix_get(C11,k,l));
            l = l + 1;
        }
        k = k + 1;
        l = 0;
    }

        // populate the C12 region
    k = 0;
    l = 0;
    for(i=0;i<N_fixed_effects;i++){
        for(j=N_fixed_effects;j<dim_C;j++){
            gsl_matrix_set(C,i,j,gsl_matrix_get(C12,k,l));
            l = l + 1;
        }
        k = k + 1;
        l = 0;
    }

        //populate the C21 region
    k = 0;
    l = 0;
    for(i=N_fixed_effects;i<dim_C;i++){
        for(j=0;j<N_fixed_effects;j++){
            gsl_matrix_set(C,i,j,gsl_matrix_get(C21,k,l));
            l = l + 1;
        }
        k = k + 1;
        l = 0;
    }

        //populate the C22 region
    k = 0;
    l = 0;
    for(i=N_fixed_effects;i<dim_C;i++){
        for(j=N_fixed_effects;j<dim_C;j++){
            gsl_matrix_set(C,i,j,gsl_matrix_get(C22,k,l));
            l = l + 1;
        }
        k = k + 1;
        l = 0;
    }

        

    return(0);
}