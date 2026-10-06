#include "gsl_support.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
int main(void) {
    double **a = dmatrix(0, 4, 0, 7);
    a[4][7] = 3.5;
    assert(a[4][7] == 3.5);
    free_dmatrix(a, 0, 4, 0, 7);
    int **b = imatrix(1, 4, 1, 7);
    b[1][1] = 17;
    b[4][7] = -9;
    assert(b[1][1] == 17 && b[4][7] == -9);
    free_imatrix(b, 1, 4, 1, 7);
    long seed = -494;
    double first = gsl_uniform_sample(&seed);
    seed = -494;
    assert(first == gsl_uniform_sample(&seed));
    seed = 494;
    gsl_random_cleanup();
    assert(first == gsl_uniform_sample(&seed));
    double sum = 0, squares = 0;
    const int count = 200000;
    for (int i = 0; i < count; i++) {
        double x = gsl_gaussian_sample(&seed);
        sum += x; squares += x*x;
    }
    double mean = sum / count;
    double variance = squares / count - mean*mean;
    assert(fabs(mean) < .01 && fabs(variance - 1) < .02);
    gsl_random_cleanup();
    printf("GSL bounds, repeatability and Gaussian moments passed: mean=%g variance=%g\n",mean,variance);
}
