
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#define VEC_SIZE 2

typedef double Vec[VEC_SIZE];

static void rhs(double t, const double *y, double *f) {
    (void)t;
    f[0] = 0.1 * y[0] - 0.3 * y[0] * y[1];
    f[1] = 0.5 * (y[0] - 1.0) * y[1];
}

static void explicit_euler_step(const double *y, double dt, double *ynext) {
    Vec k;
    rhs(0.0, y, k);
    for (int i = 0; i < VEC_SIZE; i++) ynext[i] = y[i] + dt * k[i];
}
static int implicit_euler_step(const double *y, double dt, const double *guess, double *ynext) {
    double x = y[0], z = y[1];
    double A = 1.0 - 0.1 * dt, D = 1.0 + 0.5 * dt;
    double a = -0.5 * dt * A;
    double b = A * D + 0.5 * dt * x + 0.3 * dt * z;
    double c = -x * D;
    double roots[2];
    int nroots = 0;
    if (fabs(dt) < 1.0e-15) { ynext[0] = x; ynext[1] = z; return 0; }
    if (fabs(a) < 1.0e-15) {
        if (fabs(b) < 1.0e-15) return 1;
        roots[nroots++] = -c / b;
    } else {
        double discriminant = b*b - 4.0*a*c;
        if (discriminant < 0.0) {
            if (discriminant > -1.0e-14) discriminant = 0.0;
            else return 1;
        }
	double sqrt_discriminant = sqrt(discriminant);

        if (sqrt_discriminant == 0.0)
        {
           roots[nroots++] = -b / (2.0 * a);
        }
        else
        {
           double q = -0.5 * (b + copysign(sqrt_discriminant, b));

        if (fabs(q) < 1.0e-300)
        {
           return 1;
        }

        roots[nroots++] = q / a;
        roots[nroots++] = c / q;
	}
    }
    double reference = guess != NULL ? guess[0] : x + dt*(0.1*x - 0.3*x*z);
    double best = roots[0], best_distance = fabs(roots[0] - reference);
    for (int i=1; i<nroots; i++) {
        double denominator = D - 0.5*dt*roots[i];
        if (fabs(denominator) < 1.0e-14) continue;
        double distance = fabs(roots[i] - reference);
        if (distance < best_distance) { best = roots[i]; best_distance = distance; }
    }
    double denominator = D - 0.5*dt*best;
    if (fabs(denominator) < 1.0e-14) return 1;
    ynext[0] = best;
    ynext[1] = z / denominator;
    return 0;
}
static void richardson_step(const double *y, double dt, double *ynext, double *error) {
    Vec one_step, half_step, two_half_steps, guess_one, guess_half;
    double error_squared = 0.0;
    explicit_euler_step(y, dt, guess_one);
    explicit_euler_step(y, 0.5*dt, guess_half);
    if (implicit_euler_step(y, dt, guess_one, one_step) != 0 ||
        implicit_euler_step(y, 0.5*dt, guess_half, half_step) != 0 ||
        implicit_euler_step(half_step, 0.5*dt, one_step, two_half_steps) != 0) {
        *error = INFINITY;
        return;
    }
    for (int k=0; k<VEC_SIZE; k++) {
        double difference = 2.0*(two_half_steps[k] - one_step[k]);
        ynext[k] = one_step[k];
        error_squared += difference*difference;
    }
    *error = sqrt(error_squared);
}
static double determine_initial_step(const double *y, double tolerance) {
    Vec f, y1, f1;
    double norm_y = 0.0, norm_f = 0.0, norm_difference = 0.0;
    rhs(0.0, y, f);
    for (int k = 0; k < VEC_SIZE; k++) {
        norm_y += y[k] * y[k];
        norm_f += f[k] * f[k];
    }
    double d0 = sqrt(norm_y) / (sqrt((double)VEC_SIZE) * tolerance);
    double d1 = sqrt(norm_f) / (sqrt((double)VEC_SIZE) * tolerance);
    double h1 = (d0 < 1.0e-5 || d1 < 1.0e-5) ? 1.0e-6 : 0.01 * d0 / d1;
    for (int k = 0; k < VEC_SIZE; k++) y1[k] = y[k] + h1 * f[k];
    rhs(h1, y1, f1);
    for (int k = 0; k < VEC_SIZE; k++) {
        double difference = f1[k] - f[k];
        norm_difference += difference * difference;
    }
    double d2 = sqrt(norm_difference) / (h1 * sqrt((double)VEC_SIZE) * tolerance);
    double dmax = fmax(d1, d2);
    double h2 = dmax <= 1.0e-15 ? fmax(1.0e-6, 0.001 * h1) : sqrt(0.01 / dmax);
    return fmin(100.0 * h1, h2);
}

static int generate_solution(double tstart, double tstop, double tolerance, const char *filename, Vec solution, int *intervals_ptr) {
    double t = tstart;
    double dt, error;
    int intervals = 0;
    Vec y = {1.0, 1.0}, ynext;
    FILE *file = fopen(filename, "w");
    if (file == NULL) return 1;
    fprintf(file, "%.14e %.14e %.14e %d\n", y[0], y[1], t, intervals);
    dt = determine_initial_step(y, tolerance);
    richardson_step(y, dt, ynext, &error);
    dt *= 0.9 * sqrt(tolerance / (error + 1.0e-300));
    while (t < tstop) {
        if (t + dt > tstop) dt = tstop - t;
        richardson_step(y, dt, ynext, &error);
        if (error <= tolerance || dt <= 1.0e-14 * (tstop - tstart)) {
            t += dt;
            memcpy(y, ynext, sizeof(Vec));
            intervals++;
            fprintf(file, "%.14e %.14e %.14e %d\n", y[0], y[1], t, intervals);
            if (t >= tstop) break;
            double factor = 0.9 * sqrt(tolerance / (error + 1.0e-300));
            dt *= fmin(4.0, fmax(0.25, factor));
        } else {
            double factor = 0.9 * sqrt(tolerance / (error + 1.0e-300));
            dt *= fmax(0.1, factor);
        }
    }
    memcpy(solution, y, sizeof(Vec));
    fclose(file);
    *intervals_ptr = intervals;
    return 0;
}
static int ensure_results_directory(void) {
    if (mkdir("results", 0777) == 0 || errno == EEXIST) return 0;
    fprintf(stderr, "Cannot create results directory: %s\n", strerror(errno));
    return 1;
}

int main(int argc, char **argv) {
    double tstart = 0.0, tstop = 62.0, tolerance = 1.0e-3;
    int intervals;
    Vec solution;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-tol") && i + 1 < argc) tolerance = atof(argv[++i]);
        else if (!strcmp(argv[i], "-tstop") && i + 1 < argc) tstop = atof(argv[++i]);
        else {
            fprintf(stderr, "Unknown or incomplete option: %s\n", argv[i]);
            return 1;
        }
    }
    if (tstop <= tstart || tolerance <= 0.0 || ensure_results_directory() != 0) return 1;
    char filename[512];
    snprintf(filename, sizeof(filename), "results/LV-sequential_IE-%e.txt", tolerance);
    clock_t start = clock();
    if (generate_solution(tstart, tstop, tolerance, filename, solution, &intervals) != 0) return 1;
    printf("Sequential Lotka-Volterra implicit Euler solve finished.\n");
    printf("Tolerance:              %.6e\n", tolerance);
    printf("Number of intervals:    %d\n", intervals);
    printf("Output file:            %s\n", filename);
    printf("CPU time:               %.6f seconds\n", (double)(clock() - start) / CLOCKS_PER_SEC);
    return 0;
}

