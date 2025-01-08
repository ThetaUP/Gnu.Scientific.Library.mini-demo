#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#define DIM 3

int main(){
    double arr[DIM*DIM] = {4.0,7.0,2.0,3.0,5.0,1.0,1.0,2.0,3.0};
    double val;
    int i,j;
    gsl_matrix_view MAT_view;
    gsl_matrix *MAT;
    MAT_view = gsl_matrix_view_array(arr, DIM, DIM);
    MAT = &MAT_view.matrix;

    // LU decomposition specific
    gsl_matrix *LU_decomp_MAT = gsl_matrix_alloc(DIM,DIM); // matrix to hold the LU decomposition
    gsl_matrix_memcpy(LU_decomp_MAT, MAT); // copy the contents of MAT to LU_decomp_MAT
    int sign;
    gsl_permutation * perm_vec = gsl_permutation_alloc(DIM);
    gsl_linalg_LU_decomp(LU_decomp_MAT, perm_vec, &sign); // this function works inplace

    gsl_matrix *MAT_inv = gsl_matrix_alloc(DIM, DIM);
    gsl_linalg_LU_invert(LU_decomp_MAT, perm_vec, MAT_inv); // matrix inversion

    for(i=0;i<DIM;i++){
        for(j=0;j<DIM;j++){
            val = gsl_matrix_get(MAT_inv, i, j);
            printf("%g ", val);
        }
        printf("\n");
    }
    
    return(0);
}



