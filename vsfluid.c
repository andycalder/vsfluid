#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define N 300
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
#define fx(f) ((u > 0) * (f[i][j] - f[i][j-1]) / H + (u < 0) * (f[i][j+1] - f[i][j]) / H)
#define fy(f) ((v > 0) * (f[i][j] - f[i-1][j]) / H + (v < 0) * (f[i+1][j] - f[i][j]) / H)

typedef double field[SIZE][SIZE];

int main() {
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
            phi[i][j] = 15 * H - sqrtf((i-50)*H*(i-50)*H+(j-50)*H*(j-50)*H);
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
                if (i > 0 && j > 0 && i < SIZE - 1 && j < SIZE - 1) {
                    double u = ddy(psi);
                    double v = -ddx(psi);
                }

                double sgn = phi[i][j] / sqrtf(phi[i][j] * phi[i][j] + 0.00001);
                
                double dx_plus = (phi[i][j+(j<SIZE-1)] - phi[i][j]) / H;
                double dx_minus = (phi[i][j] - phi[i][j-(j>0)]) / H;
                double dy_plus = (phi[i+(i<SIZE-1)][j] - phi[i][j]) / H;
                double dy_minus = (phi[i][j] - phi[i-(i>0)][j]) / H;
                
                double dx = dx_plus * (dx_plus * sgn < 0 & (dx_plus + dx_minus) * sgn < 0) + dx_minus * (dx_minus * sgn > 0 & (dx_plus + dx_minus) * sgn > 0);
                double dy = dy_plus * (dy_plus * sgn < 0 & (dy_plus + dy_minus) * sgn < 0) + dy_minus * (dy_minus * sgn > 0 & (dy_plus + dy_minus) * sgn > 0);

                // Deal with local maxima/minima
                if (dx == 0 && dy == 0) {
                    dx = dx_plus;
                    dy = dy_plus;
                }
                
                // Artificial timestep from CFL stability criteria
                double dtau = H / 2;
                double lambda = dtau / dt;

                // Convect and redistance level set
                phi_new[i][j] = phi[i][j] + ( -u * fx(phi) -v * fy(phi) ) * dt + sgn * (1 - sqrtf(dx * dx + dy * dy)) * dt * lambda;
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

    return 0;
}
