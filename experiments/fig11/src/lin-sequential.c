
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
    f[0] = -0.2 * y[0];
    f[1] = -0.4 * (y[1] - y[0]);
}

static void implicit_euler_step(const double *y, double dt, double *ynext) {
    double a11 = 1.0 + 0.2 * dt;
    double a22 = 1.0 + 0.4 * dt;
    double a21 = -0.4 * dt;
    ynext[0] = y[0] / a11;
    ynext[1] = (y[1] - a21 * ynext[0]) / a22;
}

static void richardson_step(const double *y, double dt, double *ynext, double *error) {
    Vec one_step, half_step, two_half_steps;
    double error_squared = 0.0;
    implicit_euler_step(y, dt, one_step);
    implicit_euler_step(y, 0.5 * dt, half_step);
    implicit_euler_step(half_step, 0.5 * dt, two_half_steps);
    for (int k = 0; k < VEC_SIZE; k++) {
        double difference = 2.0*(two_half_steps[k] - one_step[k]);
        ynext[k] = one_step[k];
        error_squared += difference * difference;
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

static int generate_solution(double tstart, double tstop, double tolerance, int *points_ptr, double **times_ptr, Vec **solution_ptr) {
    int capacity = 1024;
    int points = 1;
    double *times = malloc((size_t)capacity * sizeof(double));
    Vec *solution = malloc((size_t)capacity * sizeof(Vec));
    double t = tstart;
    double dt, error;
    Vec y = {0.3, 0.0}, ynext;
    if (times == NULL || solution == NULL) return 1;
    times[0] = t;
    memcpy(solution[0], y, sizeof(Vec));
    dt = determine_initial_step(y, tolerance);
    richardson_step(y, dt, ynext, &error);
    dt *= 0.9 * sqrt(tolerance / (error + 1.0e-300));
    while (t < tstop) {
        if (t + dt > tstop) dt = tstop - t;
        richardson_step(y, dt, ynext, &error);
        if (error <= tolerance || dt <= 1.0e-14 * (tstop - tstart)) {
            t += dt;
            memcpy(y, ynext, sizeof(Vec));
            if (points >= capacity) {
                capacity *= 2;
                double *new_times = realloc(times, (size_t)capacity * sizeof(double));
                Vec *new_solution = realloc(solution, (size_t)capacity * sizeof(Vec));
                if (new_times == NULL || new_solution == NULL) return 1;
                times = new_times;
                solution = new_solution;
            }
            times[points] = t;
            memcpy(solution[points], y, sizeof(Vec));
            points++;
            if (t >= tstop) break;
            double factor = 0.9 * sqrt(tolerance / (error + 1.0e-300));
            dt *= fmin(4.0, fmax(0.25, factor));
        } else {
            double factor = 0.9 * sqrt(tolerance / (error + 1.0e-300));
            dt *= fmax(0.1, factor);
        }
    }
    *points_ptr = points;
    *times_ptr = times;
    *solution_ptr = solution;
    return 0;
}

static int ensure_results_directory(void) {
    if (mkdir("results", 0777) == 0 || errno == EEXIST) return 0;
    fprintf(stderr, "Cannot create results directory: %s\n", strerror(errno));
    return 1;
}

static int write_solution(const char *filename, const double *times, Vec *solution, int points) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) return 1;
    for (int i = 0; i < points; i++) fprintf(file, "%.14e %.14e %.14e %d\n", solution[i][0], solution[i][1], times[i], i);
    fclose(file);
    return 0;
}

int main(int argc, char **argv) {
    double tstart = 0.0, tstop = 20.0, tolerance = 1.0e-3;
    int points;
    double *times = NULL;
    Vec *solution = NULL;
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
    snprintf(filename, sizeof(filename), "results/lin-sequential-%e.txt", tolerance);
    clock_t start = clock();
    if (generate_solution(tstart, tstop, tolerance, &points, &times, &solution) != 0) return 1;
    if (write_solution(filename, times, solution, points) != 0) return 1;
    printf("Sequential implicit Euler solve finished.\n");
    printf("Tolerance:              %.6e\n", tolerance);
    printf("Number of intervals:    %d\n", points - 1);
    printf("Output file:            %s\n", filename);
    printf("CPU time:               %.6f seconds\n", (double)(clock() - start) / CLOCKS_PER_SEC);
    free(times);
    free(solution);
    return 0;
}

