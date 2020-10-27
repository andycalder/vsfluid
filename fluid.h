// Central differences
#define ddx(f) ((f[i][j+1] - f[i][j-1]) / (2*h))
#define ddy(f) ((f[i+1][j] - f[i-1][j]) / (2*h))
#define d2dx2(f) ((f[i][j+1] - 2*f[i][j] + f[i][j-1]) / (h*h))
#define d2dy2(f) ((f[i+1][j] - 2*f[i][j] + f[i-1][j]) / (h*h))

// First order one-sided differences
#define ddx_plus(f) ((f[i][j+1] - f[i][j]) / h)
#define ddx_minus(f) ((f[i][j] - f[i][j-1]) / h)
#define ddy_plus(f) ((f[i+1][j] - f[i][j]) / h)
#define ddy_minus(f) ((f[i][j] - f[i-1][j]) / h)

// Second order one-sided differences
#define d2dx2_plus(f) (f[i][j+2] - 2 * f[i][j+1] + f[i][j] / (h*h))
#define d2dx2_minus(f) (f[i][j] - 2 * f[i][j-1] + f[i][j-2] / (h*h))
#define d2dy2_plus(f) (f[i+2][j] - 2 * f[i+1][j] + f[i][j] / (h*h))
#define d2dy2_minus(f) (f[i][j] - 2 * f[i-1][j] + f[i-2][j] / (h*h))

// First order upwind differences
#define ddx_upwind(f) ((u > 0) * ddx_minus(f) + (u < 0) * ddx_plus(f))
#define ddy_upwind(f) ((v > 0) * ddy_minus(f) + (v < 0) * ddy_plus(f))

#define isBoundary (j == 0 || j == nx - 1 || i == 0 || i == ny - 1)

#define min(a, b) (a < b ? a : b)
#define max(a, b) (a > b ? a : b)
#define clamp(x, lower, upper) (min(upper, max(x, lower)))

struct Params {
    int nx;
    int ny;
    int nt;
    double h;
    double dt;
    double nu;
    double epsilon;
    double sigma;
    double g;
    double* phi0;
};
