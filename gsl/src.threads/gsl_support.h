#ifndef GSL_SUPPORT_H
#define GSL_SUPPORT_H
/* GSL storage with row pointers for the existing a[i][j] model equations.
   Bounds are inclusive and must be nonnegative. */
double **dmatrix(long row_first, long row_last, long col_first, long col_last);
void free_dmatrix(double **rows, long row_first, long row_last, long col_first, long col_last);
int **imatrix(long row_first, long row_last, long col_first, long col_last);
void free_imatrix(int **rows, long row_first, long row_last, long col_first, long col_last);
double gsl_uniform_sample(long *seed);
double gsl_gaussian_sample(long *seed);
void gsl_random_cleanup(void);
#endif
