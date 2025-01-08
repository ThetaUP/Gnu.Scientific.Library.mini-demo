/*
How to write a matrix to a file.
*/

#include <stdio.h>
#include <gsl/gsl_matrix.h>

int main(){
    /*
    To ensure that the dimensions of the output file match the dimensions
    of the matrix we must loop through the matrix, and write each element
    manuall to the file.
    */
    int dim,i,j,value;
    double element;
    FILE *output_file;
    gsl_matrix *mat;
    dim = 3;

    mat = gsl_matrix_alloc(dim,dim);
    value = 1;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            gsl_matrix_set(mat,i,j,value);
            value = value + 1;
        }
    }
    output_file = fopen("output_file.txt", "w");

    // loop through the matrix and write each single element
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            element = gsl_matrix_get(mat,i,j);
            fprintf(output_file, "%g ", element);
        }
        fprintf(output_file,"\n"); // add a newline after the line has ended

    }
    return(0);
}


