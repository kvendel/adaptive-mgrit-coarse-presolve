#define _POSIX_C_SOURCE 200809L
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static const double alpha = 0.1;
static const double safety = 0.90;
static const double tiny = 1.0e-300;

typedef struct
{
   double tstart, tstop;
   double xstart, xstop;
   double xleft, xright;
   int nspace;
   double tolerance;
   double *tri_c;
   double *rich_one;
   double *rich_mid;
   double *rich_two;
   double *dt_F;
   double *dt_y1;
   double *dt_F1;
} Problem;

static void *checked_malloc(size_t bytes, const char *where)
{
   void *ptr = malloc(bytes);
   if (ptr == NULL)
   {
      fprintf(stderr, "Allocation failure in %s (%zu bytes).\n", where, bytes);
      exit(EXIT_FAILURE);
   }
   return ptr;
}

static void *checked_realloc(void *ptr, size_t bytes, const char *where)
{
   void *new_ptr = realloc(ptr, bytes);
   if (new_ptr == NULL)
   {
      fprintf(stderr, "Reallocation failure in %s (%zu bytes).\n", where, bytes);
      exit(EXIT_FAILURE);
   }
   return new_ptr;
}

static void allocate_workspace(Problem *problem)
{
   size_t bytes = (size_t)problem->nspace * sizeof(double);
   problem->tri_c = checked_malloc(bytes, "workspace/tri_c");
   problem->rich_one = checked_malloc(bytes, "workspace/rich_one");
   problem->rich_mid = checked_malloc(bytes, "workspace/rich_mid");
   problem->rich_two = checked_malloc(bytes, "workspace/rich_two");
   problem->dt_F = checked_malloc(bytes, "workspace/dt_F");
   problem->dt_y1 = checked_malloc(bytes, "workspace/dt_y1");
   problem->dt_F1 = checked_malloc(bytes, "workspace/dt_F1");
}

static void free_workspace(Problem *problem)
{
   free(problem->tri_c);
   free(problem->rich_one);
   free(problem->rich_mid);
   free(problem->rich_two);
   free(problem->dt_F);
   free(problem->dt_y1);
   free(problem->dt_F1);
}

static double wall_time(void)
{
   struct timespec ts;
   if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0)
   {
      perror("clock_gettime");
      exit(EXIT_FAILURE);
   }
   return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

static void initial_condition(const Problem *problem, double *u)
{
   int n = problem->nspace;
   u[0] = problem->xleft;
   u[n - 1] = problem->xright;
   for (int i = 1; i < n - 1; i++)
   {
      double x = problem->xstart + (problem->xstop - problem->xstart) *
                 (double)i / (double)(n - 1);
      double z = x - 0.435576;
      u[i] = 10.0 * exp(-z * z) - 7.27185;
   }
}

/* Identical implicit Euler propagator to the two XBraid versions. */
static void implicit_euler_step(Problem *problem, const double *y,
                                double dt, double *ynext)
{
   const int n = problem->nspace;
   const double dx = (problem->xstop - problem->xstart) / (double)(n - 1);
   const double r = alpha * dt / (dx * dx);
   const double lower = -r;
   const double diagonal = 1.0 + 2.0 * r;
   double *cp = problem->tri_c;

   ynext[0] = problem->xleft;
   cp[0] = 0.0;
   for (int i = 1; i < n - 1; i++)
   {
      double denominator = diagonal - lower * cp[i - 1];
      cp[i] = lower / denominator;
      ynext[i] = (y[i] - lower * ynext[i - 1]) / denominator;
   }
   ynext[n - 1] = problem->xright;
   for (int i = n - 2; i >= 1; i--) ynext[i] -= cp[i] * ynext[i + 1];
}

/* The accepted solution is the one-full-step implicit Euler value, while the
   two-half-step difference is used only as the Richardson error estimate. */
static void richardson_step(Problem *problem, const double *y, double dt,
                            double *ynext, double *error)
{
   int n = problem->nspace;
   double error_squared = 0.0;
   implicit_euler_step(problem, y, dt, problem->rich_one);
   implicit_euler_step(problem, y, 0.5 * dt, problem->rich_mid);
   implicit_euler_step(problem, problem->rich_mid, 0.5 * dt, problem->rich_two);
   for (int i = 0; i < n; i++)
   {
      double difference = 2.0*(problem->rich_two[i] - problem->rich_one[i]);
      ynext[i] = problem->rich_one[i];
      error_squared += difference * difference;
   }
   *error = sqrt(error_squared / (double)n);
}

static double det_dt0(const double *y, double eps0, Problem *problem)
{
   const int n = problem->nspace;
   const double dx = (problem->xstop - problem->xstart) / (double)(n - 1);
   const double scale = alpha / (dx * dx);
   double *F = problem->dt_F;
   double *y1 = problem->dt_y1;
   double *F1 = problem->dt_F1;
   double norm_y2 = 0.0, norm_F2 = 0.0, norm_diff2 = 0.0;

   F[0] = 0.0;
   F[n - 1] = 0.0;
   for (int i = 1; i < n - 1; i++)
      F[i] = scale * (y[i - 1] - 2.0 * y[i] + y[i + 1]);
   for (int i = 0; i < n; i++)
   {
      norm_y2 += y[i] * y[i];
      norm_F2 += F[i] * F[i];
   }
   double root_n_eps = sqrt((double)n) * eps0;
   double d0 = sqrt(norm_y2) / root_n_eps;
   double d1 = sqrt(norm_F2) / root_n_eps;
   double h1 = (d0 < 1.0e-5 || d1 < 1.0e-5) ? 1.0e-6 : 0.01 * d0 / d1;

   for (int i = 0; i < n; i++) y1[i] = y[i] + h1 * F[i];
   F1[0] = 0.0;
   F1[n - 1] = 0.0;
   for (int i = 1; i < n - 1; i++)
      F1[i] = scale * (y1[i - 1] - 2.0 * y1[i] + y1[i + 1]);
   for (int i = 0; i < n; i++)
   {
      double d = F1[i] - F[i];
      norm_diff2 += d * d;
   }
   double d2 = sqrt(norm_diff2) / (h1 * root_n_eps);
   double dmax = d1 > d2 ? d1 : d2;
   double h2 = dmax <= 1.0e-15 ? fmax(1.0e-6, h1 * 0.001)
                                : sqrt(0.01 / dmax);
   return fmin(100.0 * h1, h2);
}

static int generate_solution(Problem *problem, int *points_ptr,
                             double **times_ptr, double **solution_ptr)
{
   int capacity = 1024;
   int points = 1;
   int n = problem->nspace;
   double time = problem->tstart;
   double *times = checked_malloc((size_t)capacity * sizeof(double), "solution/times");
   double *solution = checked_malloc((size_t)capacity * n * sizeof(double), "solution/values");
   double *y = checked_malloc((size_t)n * sizeof(double), "solution/y");
   double *ynext = checked_malloc((size_t)n * sizeof(double), "solution/ynext");
   double error;

   initial_condition(problem, y);
   times[0] = time;
   memcpy(solution, y, (size_t)n * sizeof(double));
   double step = det_dt0(y, problem->tolerance, problem);
   richardson_step(problem, y, step, ynext, &error);
   step *= safety * sqrt(problem->tolerance / (error + tiny));

   while (time < problem->tstop)
   {
      if (time + step > problem->tstop) step = problem->tstop - time;
      richardson_step(problem, y, step, ynext, &error);
      double factor = safety * sqrt(problem->tolerance / (error + tiny));
      if (error <= problem->tolerance ||
          step <= 1.0e-14 * (problem->tstop - problem->tstart))
      {
         time += step;
         memcpy(y, ynext, (size_t)n * sizeof(double));
         if (points >= capacity)
         {
            capacity *= 2;
            times = checked_realloc(times, (size_t)capacity * sizeof(double), "solution/times");
            solution = checked_realloc(solution, (size_t)capacity * n * sizeof(double), "solution/values");
         }
         times[points] = time;
         memcpy(solution + (size_t)points * n, y, (size_t)n * sizeof(double));
         points++;
         if (factor > 4.0) factor = 4.0;
         if (factor < 0.25) factor = 0.25;
      }
      else if (factor < 0.1)
      {
         factor = 0.1;
      }
      step *= factor;
   }

   free(y);
   free(ynext);
   *points_ptr = points;
   *times_ptr = times;
   *solution_ptr = solution;
   return 0;
}

int main(int argc, char **argv)
{
   Problem problem = {0.0, 5.0, 0.0, 1.0, 1.0, 0.0, 17, 1.0e-5,
                      NULL, NULL, NULL, NULL, NULL, NULL, NULL};
   int points;
   double *times = NULL;
   double *solution = NULL;

   for (int i = 1; i < argc; i++)
   {
      if (!strcmp(argv[i], "-tol") && i + 1 < argc) problem.tolerance = atof(argv[++i]);
      else if (!strcmp(argv[i], "-nx") && i + 1 < argc) problem.nspace = atoi(argv[++i]);
      else
      {
         fprintf(stderr, "Unknown or incomplete option: %s\n", argv[i]);
         return EXIT_FAILURE;
      }
   }
   if (problem.nspace < 3 || problem.tolerance <= 0.0)
   {
      fprintf(stderr, "Invalid nx or tolerance.\n");
      return EXIT_FAILURE;
   }
   allocate_workspace(&problem);

   double total_start = wall_time();
   double solve_start = wall_time();
   generate_solution(&problem, &points, &times, &solution);
   double solve_time = wall_time() - solve_start;
   double total_time = wall_time() - total_start;

   printf("Sequential heat equation solve finished.\n");
   printf("Tolerance:             %.6e\n", problem.tolerance);
   printf("Number of time points: %d\n", points);
   printf("Number of intervals:   %d\n", points - 1);
   printf("Spatial points:        %d\n", problem.nspace);
   printf("Solve wall time:       %.8e seconds\n", solve_time);
   printf("Total wall time:       %.8e seconds\n", total_time);

   free(times);
   free(solution);
   free_workspace(&problem);
   return EXIT_SUCCESS;
}
