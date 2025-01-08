#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>  // BLAS support is needed
#define NROW 2
#define NCOL 3

int main(){
    int i,j;
    double val;
    double arr[6] = {1.0,2.0,3.0,4.0,5.0,6.0}; // our array
    gsl_matrix_view MAT_view; // define the matrix VIEW
    gsl_matrix *MAT; // define the MATRIX
    MAT_view = gsl_matrix_view_array(arr, NROW, NCOL); // create the view from the array
    MAT = &MAT_view.matrix; // convert the view to a standard gsl matrix

    // do something with the matrix
    for(i=0;i<NROW;i++){
        for(j=0;j<NCOL;j++){
            val = gsl_matrix_get(MAT,i,j);
            printf("%g ", val);
        }
        printf("\n");
    }
    return(0);
}