/* GSL replacements written for these recovered research codes.
   No Numerical Recipes implementation is used here. */
#include "gsl_support.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_int.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef struct { gsl_matrix *values; double *rows[]; } DoubleMatrix;
typedef struct { gsl_matrix_int *values; int *rows[]; } IntegerMatrix;

static void check_bounds(long first_row, long last_row, long first_col, long last_col) {
    if (first_row < 0 || first_col < 0 || last_row < first_row || last_col < first_col) {
        fprintf(stderr, "Matrix bounds must be ordered and nonnegative.\n");
        exit(EXIT_FAILURE);
    }
}

double **dmatrix(long first_row, long last_row, long first_col, long last_col) {
    check_bounds(first_row, last_row, first_col, last_col);
    size_t row_count = (size_t)last_row + 1;
    size_t col_count = (size_t)last_col + 1;
    DoubleMatrix *matrix = calloc(1, sizeof(*matrix) + row_count * sizeof(double *));
    if (!matrix) { perror("Matrix row allocation"); exit(EXIT_FAILURE); }
    matrix->values = gsl_matrix_calloc(row_count, col_count);
    if (!matrix->values) { fprintf(stderr,"GSL matrix allocation failed.\n"); exit(EXIT_FAILURE); }
    for (size_t row = 0; row < row_count; row++)
        matrix->rows[row] = gsl_matrix_ptr(matrix->values, row, 0);
    return matrix->rows;
}

void free_dmatrix(double **rows, long first_row, long last_row, long first_col, long last_col) {
    (void)first_row; (void)last_row; (void)first_col; (void)last_col;
    if (!rows) return;
    DoubleMatrix *matrix = (DoubleMatrix *)((char *)rows - offsetof(DoubleMatrix, rows));
    gsl_matrix_free(matrix->values);
    free(matrix);
}

int **imatrix(long first_row, long last_row, long first_col, long last_col) {
    check_bounds(first_row, last_row, first_col, last_col);
    size_t row_count = (size_t)last_row + 1;
    size_t col_count = (size_t)last_col + 1;
    IntegerMatrix *matrix = calloc(1, sizeof(*matrix) + row_count * sizeof(int *));
    if (!matrix) { perror("Matrix row allocation"); exit(EXIT_FAILURE); }
    matrix->values = gsl_matrix_int_calloc(row_count, col_count);
    if (!matrix->values) { fprintf(stderr,"GSL integer matrix allocation failed.\n"); exit(EXIT_FAILURE); }
    for (size_t row = 0; row < row_count; row++)
        matrix->rows[row] = gsl_matrix_int_ptr(matrix->values, row, 0);
    return matrix->rows;
}

void free_imatrix(int **rows, long first_row, long last_row, long first_col, long last_col) {
    (void)first_row; (void)last_row; (void)first_col; (void)last_col;
    if (!rows) return;
    IntegerMatrix *matrix = (IntegerMatrix *)((char *)rows - offsetof(IntegerMatrix, rows));
    gsl_matrix_int_free(matrix->values);
    free(matrix);
}

static gsl_rng *generator = NULL;
static gsl_rng *random_generator(long *seed) {
    int needs_seed = !generator || *seed <= 0;
    if (!generator) generator = gsl_rng_alloc(gsl_rng_mt19937);
    if (needs_seed) {
        unsigned long magnitude = *seed < 0 ? (unsigned long)(-(*seed + 1)) + 1UL : (*seed > 0 ? (unsigned long)*seed : 1UL);
        gsl_rng_set(generator, magnitude);
        *seed = 1;
    }
    return generator;
}

double gsl_uniform_sample(long *seed) { return gsl_rng_uniform(random_generator(seed)); }
double gsl_gaussian_sample(long *seed) { return gsl_ran_gaussian(random_generator(seed), 1.0); }
void gsl_random_cleanup(void) {
    if (generator) gsl_rng_free(generator);
    generator = NULL;
}
