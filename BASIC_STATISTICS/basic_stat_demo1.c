#include<stdio.h>
#include<gsl/gsl_statistics.h>   // load the 'gsl_statistics' header file

#define SIZE 4

int main(){

	double data[SIZE] = {11.00, 12.3, -9.10, 4.00};  // define data array of type double
	double mean, mean_stride2, variance, min, max;     // these variables will hold the output of the function, which return doubles
	int stride = 1;   // stride parameter
	int bigger_stride = 2;

	mean     = gsl_stats_mean(data, stride, SIZE);
	variance = gsl_stats_variance(data, stride, SIZE);
	max  = gsl_stats_max(data, stride, SIZE);
	min = gsl_stats_min(data, stride, SIZE);
	
	mean_stride2 = gsl_stats_mean(data, bigger_stride, SIZE);


	printf("mean with stride 1: %g\n", mean);
	printf("mean with stride 2: %g\n", mean_stride2);
    printf("variance: %g\n", variance);
    printf("max: %g\n", max);
    printf("min: %g\n", min);

	
	return(0);
}
