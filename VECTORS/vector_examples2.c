#include <stdio.h>
#include <gsl/gsl_blas.h>
#define SIZE 4

int main(){
    
    /*
    int gsl_blas_dsdot(const gsl_vector_float *x, const gsl_vector_float *y, double *result)

    input:
        - two vectors of type gsl_vector_float (struct)
        - a scalar value of type float (the result will be stored there)
    
    output:
        - an integer: 0 --> success, 1 --> failure
    */

    gsl_vector_float *u , *v;        // new type!
    int i, exit_code;
    float j;
    double result;

    u = gsl_vector_float_alloc(SIZE); // designated allocation function for this type
    v = gsl_vector_float_alloc(SIZE);

    j = 0.0;
    for(i=0;i<SIZE;i++){
        gsl_vector_float_set(u, i, j*10); // designated value set function for this type
        gsl_vector_float_set(v, i, j*100);
        j = j + 1.0;
    }

    exit_code = gsl_blas_dsdot(u, v, &result); // t(u)*v, '&result' reffers to the address of the variable 'result'

    printf("Exitcode: %d\n", exit_code);
    printf("t(u)*v = %g\n", result);

    return(0);
}


