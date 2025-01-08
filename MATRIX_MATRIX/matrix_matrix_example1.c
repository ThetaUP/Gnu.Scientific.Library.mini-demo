#include <stdio.h>
#include <gsl/gsl_matrix.h>
#define N 10 // rows
#define P 4 // columns

int main(){
    int i,j,status;
    double val;
    gsl_matrix * M = gsl_matrix_alloc(N, P);
    gsl_matrix * M_T = gsl_matrix_alloc(P,N);  // t(M)

    for(i=0;i<N;i++){
        for(j=0;j<P;j++){
            gsl_matrix_set(M,i,j,i+j);
        }
    }

    //transpose
    status = gsl_matrix_transpose_memcpy(M_T, M);
    printf("Exitcode: %d\n", status);
    printf("_____________________________________\n");

    printf("Transposed matrix:\n");
    for(i=0;i<P;i++){
        for(j=0;j<N;j++){
            val = gsl_matrix_get(M_T, i, j);
            printf("%g ", val);
        }
        printf("\n");
    }
    return(0);
}

