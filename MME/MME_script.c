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
    /*
    dim(C11) = N_fixed_effects x N_fixed_effects
    dim(C12) = N_fixed_effects x N_total
    dim(C21) = N_total x N_fixed_effects
    dim(C22) = N_total x N_total
    dim(C) = (N_total + N_fixed_efdfects) x (N_total + N_fixed_efdfects)
    */

    int dim_C = N_total + N_fixed_effects;
    gsl_matrix *C = gsl_matrix_alloc(dim_C, dim_C);
    
    // initilize the C matrix
    int i,j;

    for(i=0;i<dim_C;i++){
        for(j=0;j<dim_C;j++){

            if(i>=0 && i<N_fixed_effects && j>=0 && j<N_fixed_effects){ // here we are in the C11 region
                // here do t(X)%*%X and save it to the correct region of C
            }

            if(i>=0 && i<N_fixed_effects && j>=N_fixed_effects && j<dim_C){ // here we are in the C12 region
                // here do t(X)%*%Z and save it to the correct region of C
            }

            if(i>=N_fixed_effects && i<dim_C && j>=0 && j<=N_fixed_effects){ // here we are in the C21 region
                // here do t(Z)%*%X and save it to the correct region of C
            }

            if(i>=N_fixed_effects && i<dim_C && j>=N_fixed_effects && j<dim_C){ // here we are in the C22 region
                // here do t(Z)%*%Z+invA*c(alpha) and save it to the correct region of C
            }
        }
    }
    

    return(0);
}