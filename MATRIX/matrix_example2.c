/*
How to read a matrix from a file.
*/

#include <stdio.h>
#include <gsl/gsl_matrix.h>

int main(){
    FILE *input_file;
    gsl_matrix *mat;
    int dim, status_read, i, j;
    dim = 3;
    input_file = fopen("input_file.txt", "r");
    mat = gsl_matrix_alloc(dim,dim);
    status_read = gsl_matrix_fscanf(input_file, mat);  // read from an ASCII file, THE FILE DIMENSIONS MUST MATCH EXACTLY
    printf("reading exitcode: %d\n", status_read);

    printf("File contents: \n");
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            printf("%g ", gsl_matrix_get(mat, i, j));
        }
        printf("\n");
    }
    return(0);
}


