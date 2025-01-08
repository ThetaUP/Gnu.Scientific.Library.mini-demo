#include <stdio.h>
#include <gsl/gsl_vector.h>

int main(){
    int i,a;
    a = 4;
    double max,min,argmax,argmin;

    // allocate space for vectors, these are structs and not standard C arrays
    gsl_vector *u = gsl_vector_alloc(5);

    for(i=0; i<5; i++){
        gsl_vector_set(u, i, i*10); // write the value i*10 to the i-th position of u
    }

    // calculate max,min,argmax,argmin, the return type of these functions is double
    max = gsl_vector_max(u);
    min = gsl_vector_min(u);
    argmax = gsl_vector_max_index(u);
    argmin = gsl_vector_min_index(u);

    printf("max: %g\n", max);
    printf("min: %g\n", min);
    printf("argmax: %g\n", argmax);
    printf("argmin: %g\n", argmin);

    return(0);
}

