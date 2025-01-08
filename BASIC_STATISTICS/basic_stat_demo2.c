#include<stdio.h>
#include<gsl/gsl_statistics.h>

#define SIZE 20

int main(){

    double sample1[SIZE] = {4.236, 7.895, 2.467, 9.032, 5.364, 1.907, 8.123, 6.548, 0.483, 3.792, 7.024, 4.691, 9.865, 2.140, 5.378, 7.219, 0.837, 3.908, 6.255, 1.638};
    double sample2[SIZE] = {6.483, 2.917, 8.256, 4.732, 1.398, 7.031, 9.128, 0.573, 3.684, 5.439, 2.741, 8.395, 6.059, 4.820, 0.694, 7.567, 1.274, 9.342, 3.018, 5.768};
    double lag1_autocorr, correlation, covariance;
    int stride = 1;

    lag1_autocorr = gsl_stats_lag1_autocorrelation(sample1, stride, SIZE); // calculate lag 1 autocorrelation
    correlation = gsl_stats_correlation(sample1, stride, sample2, stride, SIZE); // calculate correlation
    covariance = gsl_stats_covariance(sample1, stride, sample2, stride, SIZE); // calculate covariance

    printf("lag 1 autocorrelation %g\n", lag1_autocorr);
    printf("correlation: %g\n", correlation);
    printf("covariance: %g\n", covariance);


    return(0);
}

