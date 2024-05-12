/***************************************************************************
 *
 *   File        : tasks.c
 *   Student Id  : 926204
 *   Name        : HAORAN YAO
 *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "tasks.h"
#include <assert.h>

#define INITIAL_1 0.125  // the first initial interval boundary
#define INITIAL_2 0.375  // the second initial interval boundary
#define MAX_X 1.0 // the max value x can be
#define MAX_T 0.2   // the max time 
#define PI 3.1415926535 // pi
#define D_HEADER 17L
#define INITIAL_FILE "out_waveeqn_initial.csv"
#define D_OUTPUT1 "out_waveeqn_1U.csv"
#define D_OUTPUT2 "out_waveeqn_2C.csv"
#define S_OUTPUT_NAME "out_shock.csv"
#define S_FIRST_LINE 27L  // skip the header of the input file
#define E 0.00000001    // value that determines accuracy of the root
#define MAXITER 50  // maximum number of iterations
#define FACTOR 0.01745329252 // degree to radian conversion factor
#define IBETAU 90 // initial guess of the strong shock
#define INITIAL 1 // initial size of the array that stores the data for each Mach number
#define LBOUND 0.0  // to make sense of the solutions of beta, define the meaningful range of roots.
#define UBOUND 90.0 // same
#define LHEADER "x\n" // the header in x vector file
#define BUF_SIZE 10 //  the size of buffer in file reading
#define L_INITIAL_SIZE 3 // set the initial number of rows
#define L_OUTPUT_FILE "out_linalsys.csv"
#define HEADER 7L
#define I_FIRST_LINE 6L
#define I_INITIAL_SIZE 10 //  the initial size of the parameters of interpolated functions
#define I_OUTPUT_FILE "out_interp.csv"

/* to be more concise, define a row type that aggregates all a*, b* and Q* values altogether */
typedef struct {
    double aStar;
    double bStar;
    double QStar;
} row_t;

/* here are the parameters of one function in cubic spline */
typedef struct {
    double xo;
    double h;
    double a;
    double b;
    double c;
    double d;
} func_t;

/* some file manipulating functions */
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

/* customize trigonometric functions */
double cot(double x) {
    return 1.0 / tan(x);
}

double csc(double x) {
    return 1.0 / sin(x);
}

/* newton-raphson method function */
double newton(double iGuess, double theta, double m, double gamma) {
    double x = iGuess, f, df;
    for (int i = 0; i < MAXITER; i++) {
        f = 2*cot(x)*(m*m*sin(x)*sin(x)-1)/(m*m*(gamma+cos(2*x))+2)-tan(FACTOR*theta);
        df = (4.0 * m * m * cot(x) * (m * m * sin(x) * sin(x) - 1) * sin(2 * x)) / ((m * m * (cos(2 * x) + gamma) + 2.0) * (m * m * (cos(2 * x) + gamma) + 2.0))
            - 2.0 * csc(x) * csc(x) * (m * m * sin(x) * sin(x) - 1) / (m * m * (cos(2.0 * x) + gamma) + 2)
            + (4.0 * m * m * cos(x) * cot(x) * sin(x)) / (m * m * (cos(2.0 * x) + gamma) + 2.0);
        x -= f / df;
        if (fabs(f) < E) {
            break;
        }
    }
    return x;
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
    free(matrix);
    return results;
}

double rhs_u(double fa, double fb, double dX, double c) {
    return -c * (fa - fb) / dX;
}

double rhs_c(double fc, double fd, double dX, double c) {
    return -c * (fc - fd) / (dX * 2.0);
}

void array_cpy(double *a1, double *a2, size_t length) {
    for (int i = 0; i < length; i++) {
        a1[i] = a2[i];
    }
}

void shockwave(const char* q2_file)
{
    double ibeta_l, ibeta_u; // initial guesses of beta
    double *mach, beta_l, beta_u, theta, m, gamma;
    size_t currentSize = INITIAL;
    int validMach = 0;  // record the number of Mach numbers readed
    const char *inputFile = NULL;
    char buf[BUF_SIZE]; //  this array will be used as a buffer to skip a line in .csv file

    inputFile = q2_file; 

    FILE *fp = safe_fopen(inputFile, "r");

    /* task a) */
    fseek(fp, S_FIRST_LINE, SEEK_SET);
    fscanf(fp, "%lf,%lf,%lf,%lf,%lf\n", &m, &theta, &ibeta_l, &ibeta_u, &gamma);
    /* set the initial guess */
    ibeta_l = asin(1.0 / m);
    ibeta_u = IBETAU * FACTOR;
    /* apply newton-raphson's method */
    beta_l = (newton(ibeta_l, theta, m, gamma)) / FACTOR;
    beta_u = (newton(ibeta_u, theta, m, gamma)) / FACTOR;

    /*task b) */
    fgets(buf, BUF_SIZE, fp); // skip a line
    fscanf(fp, "%lf\n", &m);
    /* iterate to calculate roots from theta=0 to theta max */
    for (theta = 0.0; theta <= 90; theta += 1.0)
        {
            /* apply newton-raphson method, don't forget to convert back to degree */
            beta_l = (newton(ibeta_l, theta, m, gamma)) / FACTOR;
            beta_u = (newton(ibeta_u, theta, m, gamma)) / FACTOR;

            /* once found the roots, check if they make sense */
            if (beta_l < theta || beta_l > UBOUND || beta_u < theta || beta_u > UBOUND)
            {
                break;
            }
        }
    
    /* task c) */

    fgets(buf, BUF_SIZE, fp);

    /* dynamic allocate the space for Mach number */
    mach = (double *)malloc(sizeof(double) * currentSize);

    while (fscanf(fp, "%lf\n", &m) == 1) {
        if (validMach == currentSize) {
            currentSize *= 2;
            mach = realloc(mach, currentSize * sizeof(double));
            exit_if_null(mach, "reallocation");
        }
        mach[validMach] = m;
        validMach++;
    }

    FILE *fp2 = safe_fopen(S_OUTPUT_NAME, "w"); // this creates a new file where the data will be written into

    fprintf(fp2, "M,theta,beta_lower,beta_upper\n");

    for (int i = 0; i < validMach; i++) {
        /* set the initial guess for each corresponding Mach number */
        ibeta_l = asin(1.0 / mach[i]);
        ibeta_u = IBETAU * FACTOR;
        for (theta = 0.0; theta <= 90.0; theta += 1.0)
        {
            /* apply newton-raphson method, don't forget to convert back to degree */
            beta_l = (newton(ibeta_l, theta, mach[i], gamma)) / FACTOR;
            beta_u = (newton(ibeta_u, theta, mach[i], gamma)) / FACTOR;

            /* once found the roots, check if they make sense */
            if (beta_l < theta || beta_l > UBOUND || beta_u < theta || beta_u > UBOUND)
            {
                break;
            }
            fprintf(fp2, "%lf,%d,%lf,%lf\n", mach[i], (int)theta, beta_l, beta_u);
        }
    }
    free(mach);
    fclose(fp);
    fclose(fp2);
}

void linalgbsys(const char* q4_file)
{
    int dataNum = 0;
    double *a, *b, *c, *Q, *results;
    const char *inputFile = q4_file;
    size_t currentSize = L_INITIAL_SIZE;

    FILE *fp = safe_fopen(inputFile, "r");

    a = (double *)malloc(sizeof(double) * L_INITIAL_SIZE);
    b = (double *)malloc(sizeof(double) * L_INITIAL_SIZE);
    c = (double *)malloc(sizeof(double) * L_INITIAL_SIZE);
    Q = (double *)malloc(sizeof(double) * L_INITIAL_SIZE);

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

    FILE *fp2 = safe_fopen(L_OUTPUT_FILE, "w");

    fprintf(fp2, LHEADER);
    for (int i = 0; i < dataNum; i++)
    {
        fprintf(fp2, "%lf\n", results[i]);
    }
    free(a);
    free(b);
    free(c);
    free(Q);
    free(results);
    fclose(fp);
    fclose(fp2);
}

void interp(const char* q5_file, const double xo)
{
    int dataNum = 0;    //  record how many data points will be readed
    double *alpha, *beta, *gamma;  // arguements for performing Thomas algorithm
    size_t currentSize = I_INITIAL_SIZE;
    const char *inputFile = NULL;
    func_t buf; // file reading buffer
    func_t *interp;  // this list stores all the interpolated functions

    inputFile = q5_file;
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

    double C[dataNum];    // record the "C" values in RHS in the linear system of equations

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

    /** Numerical Analysis 9th ed - Burden, Faires (Ch. 3 Natural Cubic Spline, Pg. 149) */
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
    double target = xo;
    FILE *fp2 = safe_fopen("out_interp.csv", "w");
    fprintf(fp2, "xo,f(xo)\n");
    for (int i = 0; i < dataNum - 2; i++) {
        /* catch the target within the interval */
        if (((target>=interp[i].xo) && (target<=(interp[i].xo+interp[i].h))) || ((target<=interp[i].xo) && (target>=(interp[i].xo+interp[i].h)))) {
            fprintf(fp2, "%lf,%lf\n", target, interp[i].a+interp[i].b*(target-interp[i].xo)+interp[i].c*(target-interp[i].xo)*(target-interp[i].xo) +
            interp[i].d*(target-interp[i].xo)*(target-interp[i].xo)* (target-interp[i].xo));
        }
    }

    free(alpha);
    free(gamma);
    free(beta);
    free(interp);
    fclose(fp);
    fclose(fp2);
}

void waveeqn(const char* q6_file)
{
    double *f, *fNext, *fInter;   // f for fi with n, fNext means fi with (n+1), fInter stands for the intermediate level
    unsigned Nx, out_iter;
    double cfl, c, x = 0.0;

    FILE *fp0 = safe_fopen(q6_file, "r");
    fseek(fp0, D_HEADER, SEEK_SET);
    fscanf(fp0, "%lf,%u,%lf,%u\n", &c, &Nx, &cfl, &out_iter);

    double dX = MAX_X/(double)Nx;

    /* to satisfy CFL condition, dT must be calculated from the equation */
    double dT = dX*cfl/c;

    /* file operators */
    // FILE *fp1 = safe_fopen(INITIAL_FILE, "w");
    FILE *fp2 = safe_fopen(D_OUTPUT1, "w");
    FILE *fp3 = safe_fopen(D_OUTPUT2, "w");

    f = (double *)malloc(sizeof(double) * (Nx + 1));
    fNext = (double *)malloc(sizeof(double) * (Nx + 1));
    fInter = (double *)malloc(sizeof(double) * (Nx + 1));

    /* write the initial condition */
    for (int i=0; i <= Nx; i++) {
        if (x<INITIAL_1) {
            f[i] = 0.0;
        }
        else if (x <= INITIAL_2) {
            f[i] = 0.5 * (1.0 - cos(8.0 * PI * (x - INITIAL_1)));
        }
        else {
            f[i] = 0.0;
        }
        x += dX;
        // fprintf(fp1, "%lf,%lf\n", x, f[i]);
    }

    /* loop over the number of timesteps, upwind scheme */
    for (int n = 0; n < out_iter; n++) {
        /* first apply the boundary stencil, i=0 */
        fInter[0] = f[0] + dT * rhs_u(f[1], f[0], dX, c);
        for (int i = 1; i <= Nx; i++) {
            /* RK2 routines */
            fInter[i] = f[i] + dT * rhs_u(f[i], f[i-1], dX, c);
            fNext[i] = f[i] + 0.5 * dT * (rhs_u(f[i], f[i-1], dX, c) + rhs_u(fInter[i], fInter[i-1], dX, c));
        }
        /* update the f(), means increment by one timestep */
        fNext[0] = f[0] + 0.5 * dT * (rhs_u(f[1], f[0], dX, c) + rhs_u(fInter[1], fInter[0], dX, c));
        array_cpy(f, fNext, Nx + 1);
    }

    fprintf(fp2, "x,f(x)\n");
    for (int i = 0; i <= Nx; i++) {
        fprintf(fp2, "%lf,%lf\n", (double)i * dX, f[i]);
    }

    /* reset the array */
    x = 0.0;
    for (int i = 0; i <= Nx; i++)
        {
            if (x < INITIAL_1)
            {
                f[i] = 0.0;
            }
            else if (x <= INITIAL_2)
            {
                f[i] = 0.5 * (1.0 - cos(8.0 * PI * (x - INITIAL_1)));
            }
            else
            {
                f[i] = 0.0;
            }
            x += dX;   
    }

    /* loop over the number of timesteps central scheme */
    for (int n = 0; n < out_iter; n++) {
        /* first apply the boundary stencil, i=0 */
        fInter[0] = f[0] + dT * rhs_c(f[1], f[0], dX, c);
        for (int i = 1; i < Nx; i++) {
            /* RK2 routines part 1*/
            fInter[i] = f[i] + dT * rhs_c(f[i+1], f[i-1], dX, c);
        }
        /* continue with the stencil, i=0 */
        fNext[0] = f[0] + 0.5 * dT * (rhs_c(f[1], f[0], dX, c) + rhs_c(fInter[1], fInter[0], dX, c));
        
        for (int i = 1; i < Nx; i++) {
            /* RK2 routines part 2*/
            fNext[i] = f[i] + 0.5 * dT * (rhs_c(f[i+1], f[i-1], dX, c) + rhs_c(fInter[i+1], fInter[i-1], dX, c));
        }
        /* finish up the last term */
        fInter[Nx] = f[Nx] + dT * rhs_c(f[Nx], f[Nx - 1], dX, c);
        fNext[Nx] = f[Nx] + 0.5 * dT * (rhs_c(f[Nx], f[Nx-1], dX, c) + rhs_c(fInter[Nx], fInter[Nx-1], dX, c));
        /* update the f(), means increment by one timestep */        
        array_cpy(f, fNext, Nx + 1);
    }

    /* write data into a file */
    fprintf(fp3, "x,f(x)\n");
    for (int i = 0; i <= Nx; i++) {
        fprintf(fp3, "%lf,%lf\n", (double)i * dX, f[i]);
    }

    free(fInter);
    free(f);
    free(fNext);
    fclose(fp0);
    fclose(fp2);
    fclose(fp3);
}
