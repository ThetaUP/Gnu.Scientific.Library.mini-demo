#include <stdio.h>
#include <gsl/gsl_matrix.h>
#define N 10 // number of rows
#define P 5  // number of columns

int main(){
    int i,j;
    double val;
    gsl_matrix * M = gsl_matrix_alloc (N, P); // allocate space

    // 'populate' matrix (remember row-major ordering)
    for(i=0; i<N; i++){
        for(j=0; j<P; j++){
            gsl_matrix_set(M, i, j, (i+j));
            //gsl_matrix_set(matrix, row, col, value)    
        }
    }

    // get elements of matrix
    for(i=0;i<N;i++){
        for(j=0;j<P;j++){
            val = gsl_matrix_get(M, i, j); // returns a double
            printf("%g\n", val);
        }
    }

     gsl_matrix_free(M); // free allocated space

    return(0);
}




