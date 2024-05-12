#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#define INITIAL_1 0.125  // the first initial interval boundary
#define INITIAL_2 0.375  // the second initial interval boundary
#define MAX_X 1.0 // the max value x can be
#define MAX_T 0.2   // the max time 
#define PI 3.1415926535 // pi
#define INPUT "in_waveeqn.csv"
#define D_HEADER 17L
#define INITIAL_FILE "out_waveeqn_initial.csv"
#define D_OUTPUT1 "out_waveeqn_1U.csv"
#define D_OUTPUT2 "out_waveeqn_2C.csv"

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

int main(int argc, char *argv[]) {
    double *f, *fNext, *fInter;   // f for fi with n, fNext means fi with (n+1), fInter stands for the intermediate level
    unsigned Nx, out_iter;
    double cfl, c, x = 0.0;

    FILE *fp0 = safe_fopen(INPUT, "r");
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