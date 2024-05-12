#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define FILE_NAME "out_shock.csv"
#define FIRST_LINE 27L  // skip the header of the input file
#define E 0.00000001    // value that determines accuracy of the root
#define MAXITER 50  // maximum number of iterations
#define FACTOR 0.01745329252 // degree to radian conversion factor
#define IBETAU 90 // initial guess of the strong shock
#define INITIAL 1 // initial size of the array that stores the data for each Mach number
#define LBOUND 0.0  // to make sense of the solutions of beta, define the meaningful range of roots.
#define UBOUND 90.0 // same
#define BUF_SIZE 10 //  the size of buffer in file reading

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

int main(int argc, char *argv[]) {
    double ibeta_l, ibeta_u; // initial guesses of beta
    double *mach, beta_l, beta_u, theta, m, gamma;
    size_t currentSize = INITIAL;
    int validMach = 0;  // record the number of Mach numbers readed
    char *inputFile = NULL;
    char buf[BUF_SIZE]; //  this array will be used as a buffer to skip a line in .csv file

    inputFile = argv[1]; 

    FILE *fp = safe_fopen(inputFile, "r");

    /* task a) */
    fseek(fp, FIRST_LINE, SEEK_SET);
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

    FILE *fp3 = safe_fopen("plot23.csv", "w");
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
            fprintf(fp3, "%lf,%lf,%lf\n", theta, beta_l, beta_u);
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

    FILE *fp2 = safe_fopen(FILE_NAME, "w"); // this creates a new file where the data will be written into

    fprintf(fp2, "M,theta,beta_lower,beta_upper\n");

    for (int i = 0; i < validMach; i++) {
        /* set the initial guess for each corresponding Mach number */
        ibeta_l = asin(1.0 / mach[i]);
        ibeta_u = IBETAU * FACTOR;
        for (theta = 0.0; theta <= 90; theta += 1.0)
        {
            /* apply newton-raphson method, don't forget to convert back to degree */
            beta_l = (newton(ibeta_l, theta, mach[i], gamma)) / FACTOR;
            beta_u = (newton(ibeta_u, theta, mach[i], gamma)) / FACTOR;

            /* once found the roots, check if they make sense */
            if (beta_l < theta || beta_l > UBOUND || beta_u < theta || beta_u > UBOUND)
            {
                break;
            }
            fprintf(fp2, "%lf,%lf,%lf,%lf\n", mach[i], theta, beta_l, beta_u);
        }
    }
    
    fclose(fp);
    fclose(fp2);
    fclose(fp3);

    return 0;
}
