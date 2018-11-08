#include <stdio.h>
#include <string.h>

#define N 1000
#define SIZE 40
#define H 0.025
#define U 1.0
#define DT 0.00015
#define nu 0.002

// Finite difference derivatives
#define ddx(f) ((f[i][j+1] - f[i][j-1]) / (2*H))
#define ddy(f) ((f[i+1][j] - f[i-1][j]) / (2*H))
#define d2dx2(f) ((f[i][j+1] - 2*f[i][j] + f[i][j-1]) / (H*H))
#define d2dy2(f) ((f[i+1][j] - 2*f[i][j] + f[i-1][j]) / (H*H))

int main() {
    // Declare field variables
    float omega[SIZE][SIZE] = {0};
    float domega[SIZE][SIZE] = {0};
    float psi[SIZE][SIZE] = {0};
    float dpsi[SIZE][SIZE] = {0};

    float output[SIZE * SIZE * N];
    
    int n, i, j;
    
    for (n = 0; n < N; n++) {

        // Set boundary conditions
        for (i = 1; i < SIZE - 1; i++) {
            // Top
            omega[0][i] = 2.0 * -psi[1][i] / (H*H) + 2 * U / H;
            // Bottom
            omega[SIZE-1][i] = 2.0 * -psi[SIZE-2][i] / (H*H);
            // Left
            omega[i][0] = 2.0 * -psi[i][1] / (H*H);
            // Right
            omega[i][SIZE-1] = 2.0 * -psi[i][SIZE-2] / (H*H);
        }

        // Update
        for (i = 1; i < SIZE - 1; i++) {
            for (j = 1; j < SIZE - 1; j++) {
                domega[i][j] = -ddy(psi) * ddx(omega) + ddx(psi) * ddy(omega) + nu * (d2dx2(omega) + d2dy2(omega));
                dpsi[i][j] = 0.25 * (psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] + H * H * (omega[i][j] + domega[i][j] * DT));
            }
        }
        for (i = 1; i < SIZE - 1; i++) {
            for (j = 1; j < SIZE - 1; j++) {
                omega[i][j] += domega[i][j] * DT;
                psi[i][j] = dpsi[i][j];
            }
        }

        // Copy to correct location in output buffer
        memcpy(output + n * SIZE * SIZE, psi, SIZE * SIZE * sizeof(float));
    }

    // Save to disk
    printf("Writing file...");
    FILE *file = fopen("./output", "w");
    fwrite(output, sizeof(float), SIZE * SIZE * N, file);
    fclose(file);
    printf("done\n");

    return 0;
}
