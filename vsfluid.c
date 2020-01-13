#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define N 1000
#define SIZE 100
#define H 0.02
#define U 0.0
#define dt 0.01
#define nu 0.006
#define g 0.0

// Central difference derivatives
#define ddx(f) ((f[i][j+1] - f[i][j-1]) / (2*H))
#define ddy(f) ((f[i+1][j] - f[i-1][j]) / (2*H))
#define d2dx2(f) ((f[i][j+1] - 2*f[i][j] + f[i][j-1]) / (H*H))
#define d2dy2(f) ((f[i+1][j] - 2*f[i][j] + f[i-1][j]) / (H*H))

// First order upwind
#define fx(f) (sign(u) * (-f[i][j-sign(u)] + f[i][j]) / H)
#define fy(f) (sign(v) * (-f[i-sign(v)][j] + f[i][j]) / H)

// Second order upwind
#define fxx(f) ((f[i][j] + f[i][j-2*sign(u)] - 2 * f[i][j-sign(u)]) / (H*H))
#define fyy(f) ((f[i][j] + f[i-2*sign(u)][j] - 2 * f[i-sign(u)][j]) / (H*H))

#define sign(x) ((x > 0) - (x < 0))

#define EPSILON (1.5 * H)
double S(double f) {
    if (f < -EPSILON) {
        return -1;
    } else if (f > EPSILON) {
        return 1;
    } else {
        return f / EPSILON + (1 / M_PI) * sin(M_PI * f / EPSILON);
    }
}

typedef double field[SIZE][SIZE];

int main() {
    clock_t begin = clock();
    
    field omega = {0};
    field omega_new = {0};
    field psi = {0};
    field psi_new = {0};
    field phi = {0};
    field phi_new = {0};

    double *output = malloc(sizeof(field) * N);
    
    int n, i, j;

    // Initialise level set
    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++) {
            if (sqrtf((i-50)*(i-50)+(j-50)*(j-50))<10) {
                phi[i][j] = phi_new[i][j] = H/2;
            } else {
                phi[i][j] = phi_new[i][j] = -H/2;
            }
        }
    }
    
    for (n = 0; n < N; n++) {

        // Set boundaries
        for (i = 1; i < SIZE - 1; i++) {
            omega[0][i] = 2.0 * -psi[1][i] / (H*H);
            omega[SIZE-1][i] = 2.0 * -psi[SIZE-2][i] / (H*H);
            omega[i][0] = 2.0 * -psi[i][1] / (H*H);
            omega[i][SIZE-1] = 2.0 * -psi[i][SIZE-2] / (H*H);
        }

        // Time step
        for (i = 1; i < SIZE - 1; i++) {
            for (j = 1; j < SIZE - 1; j++) {
                double u = ddy(psi);
                double v = -ddx(psi);

                // Vorticity transport equation
                omega_new[i][j] = omega[i][j] + (- u * fx(omega) - v * fy(omega) + nu * (d2dx2(omega) + d2dy2(omega))) * dt;
                
                // Elliptic equation (one iteration)
                psi_new[i][j] = 0.25 * (psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] + H * H * omega_new[i][j]);
            }
        }

        memcpy(omega, omega_new, sizeof(field));
        memcpy(psi, psi_new, sizeof(field));
        
        for (i = 0; i < SIZE; i++) {
            for (j = 0; j < SIZE; j++) {
                double u = 0;
                double v = 0;
                if (i == 0 || j == 0 || i == SIZE - 1 || j == SIZE - 1) {
                    double u = 0;
                    double v = 0;
                } else {
                    double u = ddy(psi);
                    double v = -ddx(psi);
                }

                int sgn = sign(phi[i][j]);
                double dx_plus = (phi[i][j+(j<SIZE-1)] - phi[i][j]) / H;
                double dx_minus = (phi[i][j] - phi[i][j-(j>0)]) / H;
                double dy_plus = (phi[i+(i<SIZE-1)][j] - phi[i][j]) / H;
                double dy_minus = (phi[i][j] - phi[i-(i>0)][j]) / H;

                double dx = dx_minus * (dx_minus * sgn > 0 & dx_plus * sgn > -dx_minus * sgn) + dx_plus * (dx_plus * sgn < 0 | dx_minus * sgn < -dx_plus * sgn);
                double dy = dy_minus * (dy_minus * sgn > 0 & dy_plus * sgn > -dy_minus * sgn) + dy_plus * (dy_plus * sgn < 0 | dy_minus * sgn < -dy_plus * sgn);
                
                // Artificial timestep from CFL stability criteria
                double dtau = H / 2;
                double lambda = dtau / dt;

                // Convected level set
                phi_new[i][j] = phi[i][j] + ( -u * fx(phi) -v * fy(phi) ) * dt + sign(phi[i][j]) * (1 - sqrtf(dx * dx + dy * dy)) * dt * lambda;
            }
        }
        memcpy(phi, phi_new, sizeof(field));
    
        // Copy to correct offset in output buffer
        memcpy(output + n * SIZE * SIZE, phi, sizeof(field));
    }

    // Save to disk
    FILE *file = fopen("./output", "w");
    fwrite(output, sizeof(double), SIZE * SIZE * N, file);
    fclose(file);

    // Print elapsed time
    clock_t end = clock();
    double elapsed_time = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Done in %f s\n", elapsed_time);

    return 0;
}
