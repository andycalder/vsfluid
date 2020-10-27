#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fluid.h"

double* solve(struct Params params) {
    int nx = params.nx;
    int ny = params.ny;
    int nt = params.nt;
    double h = params.h;
    double dt = params.dt;
    double nu = params.nu;
    double epsilon = params.epsilon;
    double sigma = params.sigma;
    double g = params.g;
    double* phi0 = params.phi0;

    double level = 2 * epsilon / M_PI;

    // Allocate working memory and output buffer
    typedef double field[nx][ny];
    field* omega = malloc(sizeof(field) * 2);
    field* psi = malloc(sizeof(field));
    field* phi = malloc(sizeof(field) * 2);
    double* output = malloc(nx * ny * nt * sizeof(double));
    
    if (!omega || !psi || !phi || !output) {
        printf("Memory allocation failed\n");
        return NULL;
    }

    memcpy(phi[0], phi0, sizeof(field));
    memcpy(phi[1], phi0, sizeof(field));

    // Lid velocity
    double U = 1.0;
    
    int n = 0;

    // Solve for k timesteps
    for (int k = 0; k < nt; k++) {

        // Set vorticity on boundaries (from Taylor series expansion)
        for (int i = 1; i < nx - 1; i++) {
            omega[n][0][i] = 2 * -psi[n][1][i] / (h*h) + 2 * U / h;
            omega[n][nx-1][i] = 2 * -psi[n][nx-2][i] / (h*h);
            omega[n][i][0] = 2 * -psi[n][i][1] / (h*h);
            omega[n][i][nx-1] = 2 * -psi[n][i][nx-2] / (h*h);
        }
        
        // Solve for psi iteratively (Gauss-Seidel method)
        for (int iter = 0; iter < 10; iter++) {
            for (int i = 1; i < ny - 1; i++) {
                for (int j = 1; j < nx - 1; j++) {
                    // Use elliptic equation for vorticity
                    psi[n][i][j] = 0.25 * (psi[n][i+1][j] + psi[n][i-1][j] + psi[n][i][j+1] + psi[n][i][j-1] + h * h * omega[n][i][j]);
                }
            }
        }
        
        // Update vorticity
        for (int i = 1; i < ny - 1; i++) {
            for (int j = 1; j < nx - 1; j++) {
                double u = ddy(psi[n]);
                double v = -ddx(psi[n]);

                // Vorticity transport equation
                omega[n+1][i][j] = omega[n][i][j] + (-u * ddx(omega[n]) - v * ddy(omega[n]) + nu * (d2dx2(omega[n]) + d2dy2(omega[n]))) * dt;
            }
        }
        
        memcpy(omega[0], omega[1], sizeof(field));
        
        // Advect and reinitialise level set
        for (int i = 0; i < ny; i++) {
            for (int j = 0; j < nx; j++) {
                // Use one-sided difference on boundaries, central difference in interior
                double dx = (j == 0) ? ddx_plus(phi[n]) : (j == nx-1) ? ddx_minus(phi[n]) : ddx(phi[n]);
                double dy = (i == 0) ? ddy_plus(phi[n]) : (i == ny-1) ? ddy_minus(phi[n]) : ddy(phi[n]);

                if (dx == 0 && dy == 0) {
                    continue;
                }

                double u = (i == 0 && j > 0 && j < nx-1) ? U : 0;
                double v = 0;
                
                if (!isBoundary) {
                    u = ddy(psi[n]);
                    v = -ddx(psi[n]);
                }

                // Advection component
                double phi_new = phi[n][i][j] + (-u * ddx_upwind(phi[n]) -v * ddy_upwind(phi[n])) * dt;

                // Smoothed sign function
                double sign = phi[n][i][j] / sqrtf(phi[n][i][j] * phi[n][i][j] + h * h);

                // Artificial timestep from CFL stability criteria
                double dtau = h / 2;
                double lambda = dtau / dt;

                double target = sqrtf(1 - pow(phi[n][i][j] / level, 2));
                double actual = sqrtf(dx * dx + dy * dy);

                // Reinitialisation component
                phi_new += sign * (target - actual) * dt * lambda;
                
                phi[n+1][i][j] = clamp(phi_new, -level, level);
            }
        }

        memcpy(phi[0], phi[1], sizeof(field));

        // Copy to correct offset in output buffer
        memcpy(output + k * nx * ny, phi[0], sizeof(field));
    }
    
    free(omega);
    free(psi);
    free(phi);
    
    return output;
}
