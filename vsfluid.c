#include <stdlib.h>
#include <string.h>
#include <math.h>

// Central difference derivatives
#define ddx(f) ((f[i][j+1] - f[i][j-1]) / (2*h))
#define ddy(f) ((f[i+1][j] - f[i-1][j]) / (2*h))
#define d2dx2(f) ((f[i][j+1] - 2*f[i][j] + f[i][j-1]) / (h*h))
#define d2dy2(f) ((f[i+1][j] - 2*f[i][j] + f[i-1][j]) / (h*h))

typedef double field[100][100];

double* solve(field phi, int size, int n, double h, double dt, double nu, double U) {
    // Allocate working memory and output array
    field omega = {0};
    field omega_new = {0};
    field psi = {0};
    field psi_new = {0};
    field phi_new = {0};
    double* output = malloc(sizeof(field) * n);
    
    int i, j, k;
    for (k = 0; k < n; k++) {

        // Set boundaries
        for (i = 1; i < size - 1; i++) {
            omega[0][i] = 2.0 * -psi[1][i] / (h*h) + 2 * U / h;
            omega[size-1][i] = 2.0 * -psi[size-2][i] / (h*h);
            omega[i][0] = 2.0 * -psi[i][1] / (h*h);
            omega[i][size-1] = 2.0 * -psi[i][size-2] / (h*h);
        }

        // Time step
        for (i = 1; i < size - 1; i++) {
            for (j = 1; j < size - 1; j++) {
                // Vorticity transport equation
                omega_new[i][j] = omega[i][j] + (- ddy(psi) * ddx(omega) + ddx(psi) * ddy(omega) + nu * (d2dx2(omega) + d2dy2(omega))) * dt;
                
                // Elliptic equation (one iteration)
                psi_new[i][j] = 0.25 * (psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] + h * h * omega_new[i][j]);
            }
        }

        memcpy(omega, omega_new, sizeof(field));
        memcpy(psi, psi_new, sizeof(field));
        
        // Convect and redistance level set
        for (i = 0; i < size; i++) {
            for (j = 0; j < size; j++) {
                double u = 0;
                double v = 0;
                if (i > 0 && j > 0 && i < size - 1 && j < size - 1) {
                    double u = ddy(psi);
                    double v = -ddx(psi);
                }
                if (i == 0) {
                    double u = U;
                }

                // Smoothed sign function
                double sgn = phi[i][j] / sqrtf(phi[i][j] * phi[i][j] + 0.00001);
                
                double dx_plus = (phi[i][j+(j<size-1)] - phi[i][j]) / h;
                double dx_minus = (phi[i][j] - phi[i][j-(j>0)]) / h;
                double dy_plus = (phi[i+(i<size-1)][j] - phi[i][j]) / h;
                double dy_minus = (phi[i][j] - phi[i-(i>0)][j]) / h;
                
                double dx = dx_plus * (dx_plus * sgn < 0 & (dx_plus + dx_minus) * sgn < 0) + dx_minus * (dx_minus * sgn > 0 & (dx_plus + dx_minus) * sgn > 0);
                double dy = dy_plus * (dy_plus * sgn < 0 & (dy_plus + dy_minus) * sgn < 0) + dy_minus * (dy_minus * sgn > 0 & (dy_plus + dy_minus) * sgn > 0);

                // Deal with local maxima/minima
                if (dx == 0 && dy == 0) {
                    dx = dx_plus;
                    dy = dy_plus;
                }
                
                // Artificial timestep from CFL stability criteria
                double dtau = h / 2;
                double lambda = dtau / dt;

                // Convect and redistance level set
                phi_new[i][j] = phi[i][j] + ( -u * ddx(phi) -v * ddy(phi) ) * dt;// + sgn * (1 - sqrtf(dx * dx + dy * dy)) * dt * lambda;
            }
        }
        memcpy(phi, phi_new, sizeof(field));
    
        // Copy to correct offset in output buffer
        memcpy(output + k * size * size, psi, sizeof(field));
    }

    return output;
}
