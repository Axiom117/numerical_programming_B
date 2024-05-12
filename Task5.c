#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#define I_FIRST_LINE 6L
#define I_INITIAL_SIZE 10 //  the initial size of the parameters of interpolated functions

/* here are the parameters of one function in cubic spline */
typedef struct {
    double xo;
    double h;
    double a;
    double b;
    double c;
    double d;
} func_t;

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
    int dataNum = 0;    //  record how many data points will be readed
    double *C;  // record the "C" values in RHS in the linear system of equations
    double *alpha, *beta, *gamma;  // arguements for performing Thomas algorithm
    size_t currentSize = I_INITIAL_SIZE;
    char *inputFile = NULL;
    func_t buf; // file reading buffer
    func_t *interp;  // this list stores all the interpolated functions

    inputFile = argv[1];
    FILE *fp = safe_fopen(inputFile, "r");
    /* skip the header */
    fseek(fp, I_FIRST_LINE, SEEK_SET);

    interp = (func_t *)malloc(sizeof(func_t) * I_INITIAL_SIZE);

    /* step 1, start reading data into the list, set all a(i) = fx(i) */
    while (fscanf(fp, "%lf,%lf\n", &buf.xo, &buf.a) == 2) {
        if (dataNum == currentSize) {
            currentSize *= 2;
            interp = realloc(interp, currentSize * sizeof(func_t));
            exit_if_null(interp, "reallocation");
            assert(interp);
        }
        interp[dataNum] = buf;
        dataNum++;
    }

    C = (double *)malloc(sizeof(double) * dataNum);

    /* calculate the h's */
    for (int i = 0; i < dataNum - 2; i++) {
        interp[i].h = interp[i + 1].xo - interp[i].xo;
    }

    /* step 2, solve the maxtrix to obtain all values of c(i)'s by solving [A]{X] = {C}*/
    C[0] = C[dataNum] = 0.0;    // set the first and the last second derivative to 0
    for (int i = 1; i < dataNum - 2; i++) {
        C[i] = 3 * (interp[i + 1].a - interp[i].a) / interp[i].h + 3 * (interp[i].a - interp[i-1].a) / interp[i - 1].h;
    }
    /* now comes to the hardest part, I employed the Thomas algorithm where the [A] maxtrix
       is about to be eliminated since it's a tri-diagonal system, hence, things are getting much
       easier here */
    /* initialize a, b and c parameters */
    alpha = (double *)malloc(sizeof(double) * dataNum);
    beta = (double *)malloc(sizeof(double) * dataNum);
    gamma = (double *)malloc(sizeof(double) * dataNum);
    /* initialize the first and last row of [A] matrix */
    alpha[0] = 1;
    beta[0] = 0;
    gamma[0] = 0;
    alpha[dataNum - 1] = 0;
    beta[dataNum - 1] = 0;
    gamma[dataNum - 1] = 1;
    /* set the first row in [A] */
    alpha[0] = 1;
    beta[0] = 0;
    gamma[0] = 0;
    /* iterate to obtain the matrix */
    for (int i = 1; i <= dataNum - 2; ++i) {
        alpha[i] = 2 * (interp[i + 1].xo - interp[i - 1].xo) - interp[i - 1].h * beta[i - 1];
        beta[i] = interp[i].h / alpha[i];
        gamma[i] = (C[i] - interp[i - 1].h * gamma[i - 1]) / alpha[i];
    }
    /* set the last row in matrix */
    alpha[dataNum-1] = 1;
    alpha[dataNum-1] = 0;
    alpha[dataNum-1] = 0;

    /* now obtain c, b, d values */
    for (int j = dataNum - 2; j >= 0; --j) {
        interp[j].c = gamma[j] - beta[j] * interp[j + 1].c;
        interp[j].b = (interp[j + 1].a - interp[j].a) / interp[j].h - interp[j].h * (interp[j + 1].c + 2 * interp[j].c) / 3;
        interp[j].d = (interp[j + 1].c - interp[j].c) / (3 * interp[j].h);
    }

    /* now we have all the interpolated functions, it's time to find all fx values */
    double target = 4.0;
    FILE *fp2 = safe_fopen("out_interp.csv", "w");
    fprintf(fp2, "xo,f(xo)\n");
    for (int i = 0; i < dataNum - 2; i++) {
        /* catch the target within the interval */
        if (((target>=interp[i].xo) && (target<=(interp[i].xo+interp[i].h))) || ((target<=interp[i].xo) && (target>=(interp[i].xo+interp[i].h)))) {
            fprintf(fp2, "%lf,%lf\n", target, interp[i].a+interp[i].b*(target-interp[i].xo)+interp[i].c*(target-interp[i].xo)*(target-interp[i].xo) +
            interp[i].d*(target-interp[i].xo)*(target-interp[i].xo)* (target-interp[i].xo));
        }
    }
    free(C);
    free(alpha);
    free(gamma);
    free(beta);
    free(interp);
    return 0;
}

