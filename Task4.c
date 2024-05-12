#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define INITIAL_SIZE 3
#define OUTPUT_FILE "out_linalsys.csv"
#define HEADER 7L

/* to be more concise, define a row type that aggregates all a*, b* and Q* values altogether */
typedef struct {
    double aStar;
    double bStar;
    double QStar;
} row_t;

FILE* safe_fopen(const char* path, const char* mode)
{
	FILE* fp = fopen(path, mode);
	if (fp == NULL) {
		perror("file open error.");
		exit(EXIT_FAILURE);
	}
	return fp;
}

void exit_if_null(void *ptr, char *msg) {
    if (!ptr) {
        printf("unexpected null pointer: %s\n", msg);
        exit(EXIT_FAILURE);
    }
}

/* this function takes a,b,c and Q vectors as arguement and perform Thomas algorithm, returns a vector "results" */
double *thomas(double *a, double *b, double *c, double *Q, int dataNum) {
    /* perform Guass-Elimination, calculate "a*"s and "Q*"s */
    double *results = (double *)malloc(sizeof(double) * dataNum);
    row_t *matrix = (row_t *)malloc(sizeof(row_t) * dataNum);

    /* special case: for i = 1 */
    matrix[0].aStar = a[0];
    matrix[0].QStar = Q[0];
    for (int i = 1; i<dataNum; i++) {
        matrix[i].aStar = a[i] - (c[i] * b[i-1] / matrix[i - 1].aStar); // calculat all a*'s
        matrix[i].QStar = Q[i] - (c[i] * matrix[i - 1].QStar / matrix[i - 1].aStar); //  calculate all Q*'s
    }

    /* now calculate the solution */
    /* special case: for i = N */
    results[dataNum - 1] = matrix[dataNum - 1].QStar / matrix[dataNum - 1].aStar;

    for (int i = dataNum - 2; i >= 0; i--) {
        results[i] = (matrix[i].QStar - b[i] * results[i + 1]) / matrix[i].aStar;
    }

    return results;
}

int main(int argc, char *argv[]) {
    int dataNum = 0;
    double *a, *b, *c, *Q, *results;
    char *inputFile = argv[1];
    size_t currentSize = INITIAL_SIZE;

    FILE *fp = safe_fopen(inputFile, "r");

    a = (double *)malloc(sizeof(double) * INITIAL_SIZE);
    b = (double *)malloc(sizeof(double) * INITIAL_SIZE);
    c = (double *)malloc(sizeof(double) * INITIAL_SIZE);
    Q = (double *)malloc(sizeof(double) * INITIAL_SIZE);

    /* load vectors into matrix */
    fseek(fp, HEADER, SEEK_SET);
    while (fscanf(fp, "%lf,%lf,%lf,%lf", &a[dataNum], &b[dataNum], &c[dataNum], &Q[dataNum]) == 4) {
        if (dataNum == currentSize-1) {
            currentSize *= 2;
            a = realloc(a, currentSize * sizeof(double));
            b = realloc(b, currentSize * sizeof(double));
            c = realloc(c, currentSize * sizeof(double));
            Q = realloc(Q, currentSize * sizeof(double));
            exit_if_null(a, "reallocation");
        }
        dataNum++;
    }

    results = thomas(a, b, c, Q, dataNum);

    FILE *fp2 = safe_fopen(OUTPUT_FILE, "w");

    for (int i = 0; i < dataNum; i++) {
        fprintf(fp2, "%lf\n", results[i]);
    }
}