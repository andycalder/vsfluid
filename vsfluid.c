#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N 1000
#define SIZE 100
#define H 0.025
#define U 1.0
#define dt (1.0 / 60)
#define nu 0.005

// Central difference derivatives
#define ddx(f) ((f[i][j+1] - f[i][j-1]) / (2*H))
#define ddy(f) ((f[i+1][j] - f[i-1][j]) / (2*H))
#define d2dx2(f) ((f[i][j+1] - 2*f[i][j] + f[i][j-1]) / (H*H))
#define d2dy2(f) ((f[i+1][j] - 2*f[i][j] + f[i-1][j]) / (H*H))

// One sided differences
#define ddx_fwd(f) ((f[i][j+(j<SIZE-1)] - f[i][j]) / H)
#define ddx_back(f) ((f[i][j] - f[i][j-(j>0)]) / H)
#define ddy_fwd(f) ((f[i+(i<SIZE-1)][j] - f[i][j]) / H)
#define ddy_back(f) ((f[i][j] - f[i-(i>0)][j]) / H)

// First order upwind
#define fx(f) (sign(u) * (-phi[i][j-sign(u)] + phi[i][j]) / H)
#define fy(f) (sign(v) * (-phi[i-sign(v)][j] + phi[i][j]) / H)

#define sign(x) ((x > 0) - (x < 0))

typedef float field[SIZE][SIZE];

int main() {
    field omega = {0};
    field omega_new = {0};
    field psi = {0};
    field psi_new = {0};
    field phi = {0};
    field phi_new = {0};

    float *output = malloc(sizeof(field) * N);
    
    int n, i, j;

    // Initialise level set
    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++) {
            phi[i][j] = -1.0 * (j <= 50) + 1.0 * (j > 50);
            phi_new[i][j] = -1.0 * (j <= 50) + 1.0 * (j > 50);
        }
    }
    
    for (n = 0; n < N; n++) {

        // Set vorticity on boundaries
        for (i = 1; i < SIZE - 1; i++) {
            omega[0][i] = 2.0 * -psi[1][i] / (H*H) + 2 * U / H;
            omega[SIZE-1][i] = 2.0 * -psi[SIZE-2][i] / (H*H);
            omega[i][0] = 2.0 * -psi[i][1] / (H*H);
            omega[i][SIZE-1] = 2.0 * -psi[i][SIZE-2] / (H*H);
        }

        // Time step
        for (i = 1; i < SIZE - 1; i++) {
            for (j = 1; j < SIZE - 1; j++) {
                float u = ddy(psi);
                float v = -ddx(psi);
                
                // Vorticity transport equation
                omega_new[i][j] = omega[i][j] + ( -u * ddx(omega) - v * ddy(omega) + nu * (d2dx2(omega) + d2dy2(omega)) ) * dt;
                
                // Elliptic equation (one iteration)
                psi_new[i][j] = 0.25 * (psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] + H * H * omega_new[i][j]);
            }
        }

        memcpy(omega, omega_new, sizeof(field));
        memcpy(psi, psi_new, sizeof(field));

        // Advect level set on moving boundary
        for (j = 1; j < SIZE-1; j++) {
            i = 0;
            float u = U;
            phi_new[0][j] = phi[0][j] - u * fx(phi) * dt;
        }

        // Advect level set in interior
        for (i = 1; i < SIZE - 1; i++) {
            for (j = 1; j < SIZE - 1; j++) {
                float u = ddy(psi);
                float v = -ddx(psi);
                
                phi_new[i][j] = phi[i][j] + ( -u * fx(phi) -v * fy(phi) ) * dt;
            }
        }

        memcpy(phi, phi_new, sizeof(field));

        // Redistancing of level set function
        for (i = 0; i < SIZE; i++) {
            for (j = 0; j < SIZE; j++) {

                float dx = ddx_back(phi) * ((ddx_back(phi) * sign(phi[i][j])) > 0) + ddx_fwd(phi) * ((ddx_fwd(phi) * sign(phi[i][j])) < 0);
                float dy = ddy_back(phi) * ((ddy_back(phi) * sign(phi[i][j])) > 0) + ddy_fwd(phi) * ((ddy_fwd(phi) * sign(phi[i][j])) < 0);
                
                // Artificial timestep from CFL stability criteria
                float dtau = H / 2;

                // One step of redistancing equation
                phi_new[i][j] = phi[i][j] + sign(phi[i][j]) * (1 - sqrtf(dx * dx + dy * dy)) * dtau;
            }
        }

        memcpy(phi, phi_new, sizeof(field));
        
        // Copy to correct location in output buffer
        memcpy(output + n * SIZE * SIZE, phi, sizeof(field));
    }

    // Save to disk
    FILE *file = fopen("./output", "w");
    fwrite(output, sizeof(float), SIZE * SIZE * N, file);
    fclose(file);

    return 0;
}
