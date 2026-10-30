#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include "braid.h"

#define VEC_SIZE 2

typedef double Vec[VEC_SIZE];

typedef struct _braid_App_struct {
    int rank;
    double tstart;
    double eps;
    double eps0;
    int ntgrid;
    double *tgrid;
    Vec *ugrid;
} my_App;

typedef struct _braid_Vector_struct {
    Vec values;
} my_Vector;

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

static double determine_initial_step(const double *y, double eps0) {
    Vec f, y1, f1;
    double norm_y = 0.0, norm_f = 0.0, norm_difference = 0.0;
    rhs(0.0, y, f);
    for (int k = 0; k < VEC_SIZE; k++) {
        norm_y += y[k] * y[k];
        norm_f += f[k] * f[k];
    }
    double d0 = sqrt(norm_y) / (sqrt((double)VEC_SIZE) * eps0);
    double d1 = sqrt(norm_f) / (sqrt((double)VEC_SIZE) * eps0);
    double h1 = (d0 < 1.0e-5 || d1 < 1.0e-5) ? 1.0e-6 : 0.01 * d0 / d1;
    for (int k = 0; k < VEC_SIZE; k++) y1[k] = y[k] + h1 * f[k];
    rhs(h1, y1, f1);
    for (int k = 0; k < VEC_SIZE; k++) {
        double difference = f1[k] - f[k];
        norm_difference += difference * difference;
    }
    double d2 = sqrt(norm_difference) / (h1 * sqrt((double)VEC_SIZE) * eps0);
    double dmax = fmax(d1, d2);
    double h2 = dmax <= 1.0e-15 ? fmax(1.0e-6, 0.001 * h1) : sqrt(0.01 / dmax);
    return fmin(100.0 * h1, h2);
}

static void generate_adaptive_grid(double tstart, double tstop, double eps0, int *ntgrid_ptr, double **tgrid_ptr, Vec **ugrid_ptr) {
    int capacity = 1024;
    int ntgrid = 1;
    double *tgrid = malloc((size_t)capacity * sizeof(double));
    Vec *ugrid = malloc((size_t)capacity * sizeof(Vec));
    double t = tstart;
    double dt, error;
    Vec y = {0.3, 0.0}, ynext;
    if (tgrid == NULL || ugrid == NULL) {
        fprintf(stderr, "Allocation failure in generate_adaptive_grid.\n");
        exit(EXIT_FAILURE);
    }
    tgrid[0] = t;
    memcpy(ugrid[0], y, sizeof(Vec));
    dt = determine_initial_step(y, eps0);
    richardson_step(y, dt, ynext, &error);
    dt *= 0.9 * sqrt(eps0 / (error + 1.0e-300));
    while (t < tstop) {
        richardson_step(y, dt, ynext, &error);
        if (error <= eps0 || dt <= 1.0e-14 * (tstop - tstart)) {
            t += dt;
            memcpy(y, ynext, sizeof(Vec));
            if (ntgrid >= capacity) {
                capacity *= 2;
                double *new_tgrid = realloc(tgrid, (size_t)capacity * sizeof(double));
                Vec *new_ugrid = realloc(ugrid, (size_t)capacity * sizeof(Vec));
                if (new_tgrid == NULL || new_ugrid == NULL) {
                    fprintf(stderr, "Reallocation failure in generate_adaptive_grid.\n");
                    exit(EXIT_FAILURE);
                }
                tgrid = new_tgrid;
                ugrid = new_ugrid;
            }
            tgrid[ntgrid] = t;
            memcpy(ugrid[ntgrid], y, sizeof(Vec));
            ntgrid++;
            if (t >= tstop) break;
            double factor = 0.9 * sqrt(eps0 / (error + 1.0e-300));
            factor = fmin(4.0, fmax(0.25, factor));
            dt *= factor;
        } else {
            double factor = 0.9 * sqrt(eps0 / (error + 1.0e-300));
            dt *= fmax(0.1, factor);
        }
    }
    *ntgrid_ptr = ntgrid;
    *tgrid_ptr = tgrid;
    *ugrid_ptr = ugrid;
}

static double refined_time(const double *tgrid, int ntgrid, int i, double theta) {
    double h = tgrid[i + 1] - tgrid[i];
    double hm = i > 0 ? tgrid[i] - tgrid[i - 1] : h;
    double hp = i < ntgrid - 2 ? tgrid[i + 2] - tgrid[i + 1] : h;
    double d0 = 2.0 * hm * h / (hm + h);
    double d1 = 2.0 * h * hp / (h + hp);
    double a = d0;
    double b = 3.0 * h - 2.0 * d0 - d1;
    double c = d0 + d1 - 2.0 * h;
    return tgrid[i] + a * theta + b * theta * theta + c * theta * theta * theta;
}

static int refinement_factor(double eps0, double eps) {
    if (eps <= 0.0 || eps0 <= eps) return 1;
    int factor = (int)ceil(sqrt(eps0 / eps));
    return factor < 1 ? 1 : factor;
}

static void refine_grid(int ntcoarse, const double *tcoarse, const Vec *ucoarse, int factor, int *ntfine_ptr, double **tfine_ptr, Vec **ufine_ptr) {
    int ntfine = (ntcoarse - 1) * factor + 1;
    double *tfine = malloc((size_t)ntfine * sizeof(double));
    Vec *ufine = malloc((size_t)ntfine * sizeof(Vec));
    if (tfine == NULL || ufine == NULL) {
        fprintf(stderr, "Allocation failure in refine_grid.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < ntcoarse - 1; i++) {
        double h = tcoarse[i + 1] - tcoarse[i];
        for (int j = 0; j < factor; j++) {
            int index = i * factor + j;
            double theta = (double)j / (double)factor;
            tfine[index] = refined_time(tcoarse, ntcoarse, i, theta);
            double beta = fabs(h) < 1.0e-300 ? 0.0 : (tfine[index] - tcoarse[i]) / h;
            for (int k = 0; k < VEC_SIZE; k++) {
                ufine[index][k] = (1.0 - beta) * ucoarse[i][k] + beta * ucoarse[i + 1][k];
            }
        }
    }
    tfine[ntfine - 1] = tcoarse[ntcoarse - 1];
    memcpy(ufine[ntfine - 1], ucoarse[ntcoarse - 1], sizeof(Vec));
    *ntfine_ptr = ntfine;
    *tfine_ptr = tfine;
    *ufine_ptr = ufine;
}

static void truncate_grid(double tstop, int *ntgrid_ptr, double **tgrid_ptr, Vec **ugrid_ptr) {
    int ntgrid = *ntgrid_ptr;
    double *tgrid = *tgrid_ptr;
    Vec *ugrid = *ugrid_ptr;
    int right = 0;
    while (right < ntgrid && tgrid[right] < tstop) right++;
    if (right >= ntgrid) {
        fprintf(stderr, "The prerefined grid does not reach tstop.\n");
        exit(EXIT_FAILURE);
    }
    if (fabs(tgrid[right] - tstop) > 1.0e-14 * fmax(1.0, fabs(tstop))) {
        int left = right - 1;
        double beta = (tstop - tgrid[left]) / (tgrid[right] - tgrid[left]);
        for (int k = 0; k < VEC_SIZE; k++) {
            ugrid[right][k] = (1.0 - beta) * ugrid[left][k] + beta * ugrid[right][k];
        }
    }
    tgrid[right] = tstop;
    *ntgrid_ptr = right + 1;
    *tgrid_ptr = tgrid;
    *ugrid_ptr = ugrid;
}

static int find_time_index(const my_App *app, double t) {
    int left = 0;
    int right = app->ntgrid - 1;
    if (t <= app->tgrid[0]) return 0;
    if (t >= app->tgrid[right]) return right;
    while (left <= right) {
        int middle = left + (right - left) / 2;
        if (app->tgrid[middle] < t) left = middle + 1;
        else if (app->tgrid[middle] > t) right = middle - 1;
        else return middle;
    }
    return fabs(app->tgrid[left] - t) < fabs(t - app->tgrid[left - 1]) ? left : left - 1;
}

static int my_TimeGrid(braid_App app, braid_Real *times, braid_Int *lower_ptr, braid_Int *upper_ptr) {
    for (braid_Int i = *lower_ptr; i <= *upper_ptr; i++) times[i - *lower_ptr] = app->tgrid[i];
    return 0;
}

static int my_Step(braid_App app, braid_Vector ustop, braid_Vector fstop, braid_Vector u, braid_StepStatus status) {
    double tstart, tstop;
    Vec unew;
    braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
    implicit_euler_step(u->values, tstop - tstart, unew);
    memcpy(u->values, unew, sizeof(Vec));
    (void)app;
    (void)ustop;
    (void)fstop;
    return 0;
}

static int my_Init(braid_App app, double t, braid_Vector *u_ptr) {
    int index = find_time_index(app, t);
    my_Vector *u = malloc(sizeof(my_Vector));
    if (u == NULL) return 1;
    memcpy(u->values, app->ugrid[index], sizeof(Vec));
    *u_ptr = u;
    return 0;
}

static int my_Clone(braid_App app, braid_Vector u, braid_Vector *v_ptr) {
    my_Vector *v = malloc(sizeof(my_Vector));
    if (v == NULL) return 1;
    memcpy(v->values, u->values, sizeof(Vec));
    *v_ptr = v;
    (void)app;
    return 0;
}

static int my_Free(braid_App app, braid_Vector u) {
    free(u);
    (void)app;
    return 0;
}

static int my_Sum(braid_App app, double alpha, braid_Vector x, double beta, braid_Vector y) {
    for (int i = 0; i < VEC_SIZE; i++) y->values[i] = alpha * x->values[i] + beta * y->values[i];
    (void)app;
    return 0;
}

static int my_SpatialNorm(braid_App app, braid_Vector u, double *norm_ptr) {
    double sum = 0.0;
    for (int i = 0; i < VEC_SIZE; i++) sum += u->values[i] * u->values[i];
    *norm_ptr = sqrt(sum);
    (void)app;
    return 0;
}

static int ensure_results_directory(void) {
    if (mkdir("results", 0777) == 0 || errno == EEXIST) return 0;
    fprintf(stderr, "Cannot create results directory: %s\n", strerror(errno));
    return 1;
}

static int my_Access(braid_App app, braid_Vector u, braid_AccessStatus status) {
    int index, done;
    char filename[512];
    braid_AccessStatusGetTIndex(status, &index);
    braid_AccessStatusGetDone(status, &done);
    if (!done || index != app->ntgrid - 1) return 0;
    snprintf(filename, sizeof(filename), "results/lin-hermite_harmonic-%e-%e.txt", app->eps0, app->eps);
    FILE *file = fopen(filename, "w");
    if (file == NULL) return 1;
    fprintf(file, "%.14e %.14e\n", u->values[0], u->values[1]);
    fclose(file);
    return 0;
}
static int my_BufSize(braid_App app, int *size_ptr, braid_BufferStatus status) {
    *size_ptr = VEC_SIZE * (int)sizeof(double);
    (void)app;
    (void)status;
    return 0;
}

static int my_BufPack(braid_App app, braid_Vector u, void *buffer, braid_BufferStatus status) {
    memcpy(buffer, u->values, sizeof(Vec));
    braid_BufferStatusSetSize(status, VEC_SIZE * (int)sizeof(double));
    (void)app;
    return 0;
}

static int my_BufUnpack(braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus status) {
    my_Vector *u = malloc(sizeof(my_Vector));
    if (u == NULL) return 1;
    memcpy(u->values, buffer, sizeof(Vec));
    *u_ptr = u;
    (void)app;
    (void)status;
    return 0;
}

int main(int argc, char **argv) {
    MPI_Comm comm = MPI_COMM_WORLD;
    braid_Core core;
    int rank, ntcoarse, ntgrid, factor;
    double tstart = 0.0, tstop = 20.0;
    double eps0 = 1.0e-3, eps = 1.0e-5, tol = 1.0e-10;
    int max_levels = 15, min_coarse = 64, nrelax = 1, nrelax0 = -1, nrelaxc = 1;
    int cfactor = 2, max_iter = 100, skip = 0, storage = -1;
    double *tcoarse = NULL, *tgrid = NULL;
    Vec *ucoarse = NULL, *ugrid = NULL;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-eps0") && i + 1 < argc) eps0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "-eps") && i + 1 < argc) eps = atof(argv[++i]);
        else if (!strcmp(argv[i], "-ml") && i + 1 < argc) max_levels = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-mc") && i + 1 < argc) min_coarse = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-nu") && i + 1 < argc) nrelax = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-nu0") && i + 1 < argc) nrelax0 = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-nuc") && i + 1 < argc) nrelaxc = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-tol") && i + 1 < argc) tol = atof(argv[++i]);
        else if (!strcmp(argv[i], "-cf") && i + 1 < argc) cfactor = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-mi") && i + 1 < argc) max_iter = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-skip") && i + 1 < argc) skip = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-storage") && i + 1 < argc) storage = atoi(argv[++i]);
        else {
            if (rank == 0) fprintf(stderr, "Unknown or incomplete option: %s\n", argv[i]);
            MPI_Finalize();
            return 1;
        }
    }
    if (ensure_results_directory() != 0) {
        MPI_Finalize();
        return 1;
    }
    if (rank == 0) generate_adaptive_grid(tstart, tstop, eps0, &ntcoarse, &tcoarse, &ucoarse);
    MPI_Bcast(&ntcoarse, 1, MPI_INT, 0, comm);
    if (rank != 0) {
        tcoarse = malloc((size_t)ntcoarse * sizeof(double));
        ucoarse = malloc((size_t)ntcoarse * sizeof(Vec));
    }
    MPI_Bcast(tcoarse, ntcoarse, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&ucoarse[0][0], ntcoarse * VEC_SIZE, MPI_DOUBLE, 0, comm);
    factor = refinement_factor(eps0, eps);
    refine_grid(ntcoarse, tcoarse, ucoarse, factor, &ntgrid, &tgrid, &ugrid);
    truncate_grid(tstop, &ntgrid, &tgrid, &ugrid);
    free(tcoarse);
    free(ucoarse);
    my_App *app = malloc(sizeof(my_App));
    app->rank = rank;
    app->tstart = tstart;
    app->eps0 = eps0;
    app->eps = eps;
    app->ntgrid = ntgrid;
    app->tgrid = tgrid;
    app->ugrid = ugrid;
    if (rank == 0) {
        char filename[512];
        snprintf(filename, sizeof(filename), "results/lin-hermite_harmonic-%e-%e.txt", eps0, eps);
        remove(filename);
        printf("Initial adaptive intervals: %d\n", ntcoarse - 1);
        printf("Refinement factor:          %d\n", factor);
        printf("Refined intervals:          %d\n", ntgrid - 1);
        printf("Output file:                %s\n", filename);
    }
    braid_Init(comm, comm, tgrid[0], tgrid[ntgrid - 1], ntgrid - 1, app, my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
    braid_SetTimeGrid(core, my_TimeGrid);
    braid_SetPrintLevel(core, 2);
    braid_SetMaxLevels(core, max_levels);
    braid_SetNRelax(core, -1, nrelax);
    if (nrelax0 >= 0) braid_SetNRelax(core, 0, nrelax0);
    braid_SetNRelax(core, max_levels - 1, nrelaxc);
    braid_SetMinCoarse(core, min_coarse);
    braid_SetAbsTol(core, tol);
    braid_SetCFactor(core, -1, cfactor);
    braid_SetMaxIter(core, max_iter);
    braid_SetSkip(core, skip);
    if (storage >= -2) braid_SetStorage(core, storage);
    clock_t start = clock();
    braid_Drive(core);
    double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
    if (rank == 0) printf("Braid CPU time:             %.6f seconds\n", elapsed);
    braid_Destroy(core);
    free(app->tgrid);
    free(app->ugrid);
    free(app);
    MPI_Finalize();
    return 0;
}

