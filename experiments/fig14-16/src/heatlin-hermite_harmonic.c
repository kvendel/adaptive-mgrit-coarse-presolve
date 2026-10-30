#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "braid.h"

static const double alpha = 0.1;
static const double safety = 0.90;
static const double tiny = 1.0e-300;

typedef struct _braid_App_struct
{
   MPI_Comm comm;
   int rank;
   int nprocs;
   double tstart, tstop;
   double xstart, xstop, xLeft, xRight;
   int nspace;
   double eps, eps0;
   int ntgrid;
   double *tgrid;
   int ntcoarse;
   double *tcoarse;
   double *ucoarse;
   int mgrit_iterations;
   double *tri_c;
   double *rich_one;
   double *rich_mid;
   double *rich_two;
   double *step_work;
   double *dt_F;
   double *dt_y1;
   double *dt_F1;
} my_App;

typedef struct _braid_Vector_struct
{
   int size;
   double *values;
} my_Vector;

static void die(const char *message)
{
   fprintf(stderr, "%s\n", message);
   MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

static void *checked_malloc(size_t bytes, const char *where)
{
   void *ptr = malloc(bytes);
   if (ptr == NULL)
   {
      fprintf(stderr, "Allocation failure in %s (%zu bytes).\n", where, bytes);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
   }
   return ptr;
}

static void *checked_realloc(void *ptr, size_t bytes, const char *where)
{
   void *new_ptr = realloc(ptr, bytes);
   if (new_ptr == NULL)
   {
      fprintf(stderr, "Reallocation failure in %s (%zu bytes).\n", where, bytes);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
   }
   return new_ptr;
}

static void allocate_workspace(my_App *app)
{
   size_t bytes = (size_t)app->nspace * sizeof(double);
   app->tri_c = checked_malloc(bytes, "allocate_workspace/tri_c");
   app->rich_one = checked_malloc(bytes, "allocate_workspace/rich_one");
   app->rich_mid = checked_malloc(bytes, "allocate_workspace/rich_mid");
   app->rich_two = checked_malloc(bytes, "allocate_workspace/rich_two");
   app->step_work = checked_malloc(bytes, "allocate_workspace/step_work");
   app->dt_F = checked_malloc(bytes, "allocate_workspace/dt_F");
   app->dt_y1 = checked_malloc(bytes, "allocate_workspace/dt_y1");
   app->dt_F1 = checked_malloc(bytes, "allocate_workspace/dt_F1");
}

static void free_workspace(my_App *app)
{
   free(app->tri_c);
   free(app->rich_one);
   free(app->rich_mid);
   free(app->rich_two);
   free(app->step_work);
   free(app->dt_F);
   free(app->dt_y1);
   free(app->dt_F1);
}

static void initial_condition(double *u, int nspace, double xLeft, double xRight)
{
   u[0] = xLeft;
   u[nspace - 1] = xRight;
   for (int i = 1; i < nspace - 1; i++)
   {
      double x = (double)i / (double)(nspace - 1);
      double z = x - 0.435576;
      u[i] = 10.0 * exp(-z * z) - 7.27185;
   }
}

/* Optimized Thomas solve for the heat-equation matrix.  Only one reusable
   c-prime work array is needed; no allocation occurs inside a time step. */
static void ie_step(my_App *app, const double *y, double dt, double *ynext)
{
   const int n = app->nspace;
   const double dx = (app->xstop - app->xstart) / (double)(n - 1);
   const double r = alpha * dt / (dx * dx);
   const double lower = -r;
   const double diagonal = 1.0 + 2.0 * r;
   double *cp = app->tri_c;

   ynext[0] = app->xLeft;
   cp[0] = 0.0;
   for (int i = 1; i < n - 1; i++)
   {
      double denominator = diagonal - lower * cp[i - 1];
      cp[i] = lower / denominator;
      ynext[i] = (y[i] - lower * ynext[i - 1]) / denominator;
   }
   ynext[n - 1] = app->xRight;
   for (int i = n - 2; i >= 1; i--)
      ynext[i] -= cp[i] * ynext[i + 1];
}

static void ie_richardson_step(my_App *app, const double *y, double dt,
                               double *ynext, double *err_ptr)
{
   const int n = app->nspace;
   double err2 = 0.0;
   ie_step(app, y, dt, app->rich_one);
   ie_step(app, y, 0.5 * dt, app->rich_mid);
   ie_step(app, app->rich_mid, 0.5 * dt, app->rich_two);
   for (int i = 0; i < n; i++)
   {
      double e = 2.0*(app->rich_two[i] - app->rich_one[i]);
      ynext[i] = app->rich_one[i];
      err2 += e * e;
   }
   *err_ptr = sqrt(err2 / (double)n);
}

static double det_dt0(const double *y, double eps0, my_App *app)
{
   const int n = app->nspace;
   const double dx = (app->xstop - app->xstart) / (double)(n - 1);
   const double scale = alpha / (dx * dx);
   double *F = app->dt_F;
   double *y1 = app->dt_y1;
   double *F1 = app->dt_F1;
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

/* Deliberately generates the pilot grid beyond tstop.  Truncation is applied
   only after prerefinement, so only the final fine interval is shortened. */
static void generate_adaptive_grid(my_App *app, int *ntgrid_ptr,
                                   double **tgrid_ptr, double **ugrid_ptr)
{
   int capacity = 1024;
   int ntgrid = 1;
   int n = app->nspace;
   double *tgrid = checked_malloc((size_t)capacity * sizeof(double), "pilot/tgrid");
   double *ugrid = checked_malloc((size_t)capacity * n * sizeof(double), "pilot/ugrid");
   double *y = checked_malloc((size_t)n * sizeof(double), "pilot/y");
   double *ynext = checked_malloc((size_t)n * sizeof(double), "pilot/ynext");
   double t = app->tstart;
   double err;

   initial_condition(y, n, app->xLeft, app->xRight);
   tgrid[0] = t;
   memcpy(ugrid, y, (size_t)n * sizeof(double));
   double dt = det_dt0(y, app->eps0, app);
   ie_richardson_step(app, y, dt, ynext, &err);
   dt *= safety * sqrt(app->eps0 / (err + tiny));

   while (t < app->tstop)
   {
      ie_richardson_step(app, y, dt, ynext, &err);
      double fac = safety * sqrt(app->eps0 / (err + tiny));
      if (err <= app->eps0 || dt <= 1.0e-14 * (app->tstop - app->tstart))
      {
         t += dt;
         memcpy(y, ynext, (size_t)n * sizeof(double));
         if (ntgrid >= capacity)
         {
            capacity *= 2;
            tgrid = checked_realloc(tgrid, (size_t)capacity * sizeof(double), "pilot/tgrid");
            ugrid = checked_realloc(ugrid, (size_t)capacity * n * sizeof(double), "pilot/ugrid");
         }
         tgrid[ntgrid] = t;
         memcpy(ugrid + (size_t)ntgrid * n, y, (size_t)n * sizeof(double));
         ntgrid++;
         if (t >= app->tstop) break;
         if (fac > 4.0) fac = 4.0;
         if (fac < 0.25) fac = 0.25;
      }
      else if (fac < 0.1)
      {
         fac = 0.1;
      }
      dt *= fac;
   }
   free(y);
   free(ynext);
   *ntgrid_ptr = ntgrid;
   *tgrid_ptr = tgrid;
   *ugrid_ptr = ugrid;
}

static void truncate_grid_at_tstop(double tstop, int *ntgrid_ptr, double **tgrid_ptr)
{
   int ntgrid = *ntgrid_ptr;
   double *tgrid = *tgrid_ptr;
   int right = 0;
   while (right < ntgrid && tgrid[right] < tstop) right++;
   if (right >= ntgrid) die("The prerefined grid does not reach tstop.");
   if (fabs(tgrid[right] - tstop) > 1.0e-14 * fmax(1.0, fabs(tstop)))
      tgrid[right] = tstop;
   *ntgrid_ptr = right + 1;
   *tgrid_ptr = checked_realloc(tgrid, (size_t)(right + 1) * sizeof(double),
                                "truncate/tgrid");
}
static int compute_ref_factor(double eps0, double eps)
{
   if (eps <= 0.0 || eps0 <= eps) return 1;
   int rf = (int)ceil(sqrt(eps0 / eps));
   return rf < 1 ? 1 : rf;
}

static double hermite_time(const double *tgrid, int ntgrid, int i, double theta)
{
   double h = tgrid[i + 1] - tgrid[i];
   double hm = i > 0 ? tgrid[i] - tgrid[i - 1] : h;
   double hp = i < ntgrid - 2 ? tgrid[i + 2] - tgrid[i + 1] : h;
   double d0 = 2.0 * hm * h / (hm + h);
   double d1 = 2.0 * h * hp / (h + hp);
   double b = 3.0 * h - 2.0 * d0 - d1;
   double c = d0 + d1 - 2.0 * h;
   return tgrid[i] + d0 * theta + b * theta * theta + c * theta * theta * theta;
}

static void fill_refined_interval(const double *tcoarse, int ntcoarse, int rf, int i, double *tfine, int base)
{
   for (int j = 0; j < rf; j++)
   {
      int id = base + j;
      double theta = (double)j / (double)rf;
      tfine[id] = hermite_time(tcoarse, ntcoarse, i, theta);
   }
}
static void refine_grid_distributed(MPI_Comm comm, int rank, int nprocs, int ntcoarse,
                                    const double *tcoarse, int rf, int *ntfine_ptr,
                                    double **tfine_ptr)
{
   int intervals = ntcoarse - 1;
   int ibegin = rank * intervals / nprocs;
   int iend = (rank + 1) * intervals / nprocs;
   int local_points = (iend - ibegin) * rf;
   int ntfine = intervals * rf + 1;
   double *local_t = checked_malloc((size_t)(local_points > 0 ? local_points : 1) * sizeof(double), "dist/local_t");
   for (int i = ibegin; i < iend; i++)
      fill_refined_interval(tcoarse, ntcoarse, rf, i, local_t, (i - ibegin) * rf);
   int *counts = checked_malloc((size_t)nprocs * sizeof(int), "dist/counts_t");
   int *displs = checked_malloc((size_t)nprocs * sizeof(int), "dist/displs_t");
   for (int q = 0; q < nprocs; q++)
   {
      int qb = q * intervals / nprocs;
      int qe = (q + 1) * intervals / nprocs;
      counts[q] = (qe - qb) * rf;
      displs[q] = qb * rf;
   }
   double *tfine = checked_malloc((size_t)ntfine * sizeof(double), "dist/tfine");
   MPI_Allgatherv(local_t, local_points, MPI_DOUBLE, tfine, counts, displs, MPI_DOUBLE, comm);
   tfine[ntfine - 1] = tcoarse[ntcoarse - 1];
   free(local_t); free(counts); free(displs);
   *ntfine_ptr = ntfine;
   *tfine_ptr = tfine;
}
static int find_coarse_interval(const my_App *app, double t)
{
   if (t <= app->tcoarse[0]) return 0;
   if (t >= app->tcoarse[app->ntcoarse - 1]) return app->ntcoarse - 2;
   int lo = 0, hi = app->ntcoarse - 1;
   while (hi - lo > 1)
   {
      int mid = lo + (hi - lo) / 2;
      if (app->tcoarse[mid] <= t) lo = mid; else hi = mid;
   }
   return lo;
}

static int my_TimeGrid(braid_App app, braid_Real *ta,
                       braid_Int *ilower_ptr, braid_Int *iupper_ptr)
{
   braid_Int ilower = *ilower_ptr, iupper = *iupper_ptr;
   for (braid_Int i = ilower; i <= iupper; i++) ta[i - ilower] = app->tgrid[i];
   return 0;
}

static int my_Step(braid_App app, braid_Vector ustop, braid_Vector fstop,
                   braid_Vector u, braid_StepStatus status)
{
   double tstart, tstop;
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   ie_step(app, u->values, tstop - tstart, app->step_work);
   memcpy(u->values, app->step_work, (size_t)u->size * sizeof(double));
   (void)ustop;
   (void)fstop;
   return 0;
}

static int my_Init(braid_App app, double t, braid_Vector *u_ptr)
{
   int i = find_coarse_interval(app, t);
   double beta = (t - app->tcoarse[i]) / (app->tcoarse[i + 1] - app->tcoarse[i]);
   const double *ul = app->ucoarse + (size_t)i * app->nspace;
   const double *ur = app->ucoarse + (size_t)(i + 1) * app->nspace;
   my_Vector *u = checked_malloc(sizeof(my_Vector), "my_Init/vector");
   u->size = app->nspace;
   u->values = checked_malloc((size_t)u->size * sizeof(double), "my_Init/values");
   for (int k = 0; k < u->size; k++)
      u->values[k] = (1.0 - beta) * ul[k] + beta * ur[k];
   *u_ptr = u;
   return 0;
}

static int my_Clone(braid_App app, braid_Vector u, braid_Vector *v_ptr)
{
   my_Vector *v = checked_malloc(sizeof(my_Vector), "my_Clone/vector");
   v->size = u->size;
   v->values = checked_malloc((size_t)v->size * sizeof(double), "my_Clone/values");
   memcpy(v->values, u->values, (size_t)v->size * sizeof(double));
   *v_ptr = v;
   (void)app;
   return 0;
}

static int my_Free(braid_App app, braid_Vector u)
{
   free(u->values);
   free(u);
   (void)app;
   return 0;
}

static int my_Sum(braid_App app, double a, braid_Vector x, double b, braid_Vector y)
{
   for (int i = 0; i < y->size; i++) y->values[i] = a * x->values[i] + b * y->values[i];
   (void)app;
   return 0;
}

static int my_SpatialNorm(braid_App app, braid_Vector u, double *norm_ptr)
{
   double sum = 0.0;
   for (int i = 0; i < u->size; i++) sum += u->values[i] * u->values[i];
   *norm_ptr = sqrt(sum);
   (void)app;
   return 0;
}

static int my_Access(braid_App app, braid_Vector u, braid_AccessStatus status)
{
   braid_Int iter;
   braid_AccessStatusGetIter(status, &iter);
   if ((int)iter > app->mgrit_iterations) app->mgrit_iterations = (int)iter;
   (void)u;
   return 0;
}

static int my_BufSize(braid_App app, int *size_ptr, braid_BufferStatus status)
{
   *size_ptr = (app->nspace + 1) * (int)sizeof(double);
   (void)status;
   return 0;
}

static int my_BufPack(braid_App app, braid_Vector u, void *buffer, braid_BufferStatus status)
{
   double *d = buffer;
   d[0] = (double)u->size;
   memcpy(d + 1, u->values, (size_t)u->size * sizeof(double));
   braid_BufferStatusSetSize(status, (u->size + 1) * (int)sizeof(double));
   (void)app;
   return 0;
}

static int my_BufUnpack(braid_App app, void *buffer, braid_Vector *u_ptr,
                        braid_BufferStatus status)
{
   double *d = buffer;
   my_Vector *u = checked_malloc(sizeof(my_Vector), "my_BufUnpack/vector");
   u->size = (int)d[0];
   u->values = checked_malloc((size_t)u->size * sizeof(double), "my_BufUnpack/values");
   memcpy(u->values, d + 1, (size_t)u->size * sizeof(double));
   *u_ptr = u;
   (void)app;
   (void)status;
   return 0;
}

static double max_elapsed(MPI_Comm comm, double start)
{
   double local = MPI_Wtime() - start;
   double maximum = 0.0;
   MPI_Reduce(&local, &maximum, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
   return maximum;
}

int main(int argc, char *argv[])
{
   MPI_Comm comm = MPI_COMM_WORLD;
   braid_Core core;
   int rank, nprocs, ntcoarse = 0, ntgrid = 0, ntime, rf;
   double *tcoarse = NULL, *ucoarse = NULL, *tgrid = NULL;
   int max_levels = 15, min_coarse = 64, nrelax = 1, nrelax0 = -1, nrelaxc = 1;
   int cfactor = 2, max_iter = 100, skip = 0, storage = -1;
   double tol = 1.0e-9;

   MPI_Init(&argc, &argv);
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nprocs);

   my_App *app = calloc(1, sizeof(my_App));
   if (app == NULL) die("Allocation failure for app.");
   app->comm = comm;
   app->rank = rank;
   app->nprocs = nprocs;
   app->tstart = 0.0;
   app->tstop = 5.0;
   app->xstart = 0.0;
   app->xstop = 1.0;
   app->xLeft = 1.0;
   app->xRight = 0.0;
   app->nspace = 17;
   app->eps0 = 1.0e-5;
   app->eps = 1.0e-7;

   for (int i = 1; i < argc; i++)
   {
      if (!strcmp(argv[i], "-eps0") && i + 1 < argc) app->eps0 = atof(argv[++i]);
      else if (!strcmp(argv[i], "-eps") && i + 1 < argc) app->eps = atof(argv[++i]);
      else if (!strcmp(argv[i], "-nx") && i + 1 < argc) app->nspace = atoi(argv[++i]);
      else if (!strcmp(argv[i], "-xL") && i + 1 < argc) app->xLeft = atof(argv[++i]);
      else if (!strcmp(argv[i], "-xR") && i + 1 < argc) app->xRight = atof(argv[++i]);
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
      else
      {
         if (rank == 0) fprintf(stderr, "Unknown or incomplete option: %s\n", argv[i]);
         MPI_Abort(comm, EXIT_FAILURE);
      }
   }
   if (app->nspace < 3 || app->eps <= 0.0 || app->eps0 <= 0.0) die("Invalid nx, eps, or eps0.");
   allocate_workspace(app);

   MPI_Barrier(comm);
   double total_start = MPI_Wtime();

   double pilot_time = 0.0;
   if (rank == 0)
   {
      double start = MPI_Wtime();
      generate_adaptive_grid(app, &ntcoarse, &tcoarse, &ucoarse);
      pilot_time = MPI_Wtime() - start;
   }

   MPI_Barrier(comm);
   double prerefine_start = MPI_Wtime();
   MPI_Bcast(&ntcoarse, 1, MPI_INT, 0, comm);
   if (rank != 0)
   {
      tcoarse = checked_malloc((size_t)ntcoarse * sizeof(double), "pilot_bcast/tcoarse");
      ucoarse = checked_malloc((size_t)ntcoarse * app->nspace * sizeof(double), "pilot_bcast/ucoarse");
   }
   MPI_Bcast(tcoarse, ntcoarse, MPI_DOUBLE, 0, comm);
   MPI_Bcast(ucoarse, ntcoarse * app->nspace, MPI_DOUBLE, 0, comm);
   rf = compute_ref_factor(app->eps0, app->eps);
   refine_grid_distributed(comm, rank, nprocs, ntcoarse, tcoarse, rf, &ntgrid, &tgrid);
   truncate_grid_at_tstop(app->tstop, &ntgrid, &tgrid);
   double prerefine_time = max_elapsed(comm, prerefine_start);

   ntime = ntgrid - 1;
   app->ntgrid = ntgrid;
   app->tgrid = tgrid;
   app->ntcoarse = ntcoarse;
   app->tcoarse = tcoarse;
   app->ucoarse = ucoarse;

   MPI_Barrier(comm);
   double setup_start = MPI_Wtime();
   braid_Init(comm, comm, tgrid[0], tgrid[ntgrid - 1], ntime, app,
              my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm,
              my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
   braid_SetTimeGrid(core, my_TimeGrid);
   braid_SetPrintLevel(core, 0);
   braid_SetMaxLevels(core, max_levels);
   braid_SetNRelax(core, -1, nrelax);
   if (nrelax0 > -1) braid_SetNRelax(core, 0, nrelax0);
   braid_SetNRelax(core, max_levels - 1, nrelaxc);
   braid_SetMinCoarse(core, min_coarse);
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, cfactor);
   braid_SetMaxIter(core, max_iter);
   braid_SetSkip(core, skip);
   if (storage >= -2) braid_SetStorage(core, storage);
   double setup_time = max_elapsed(comm, setup_start);
   (void)setup_time;

   MPI_Barrier(comm);
   double drive_start = MPI_Wtime();
   braid_Drive(core);
   double drive_time = max_elapsed(comm, drive_start);
   int mgrit_iterations = 0;
   MPI_Reduce(&app->mgrit_iterations, &mgrit_iterations, 1, MPI_INT, MPI_MAX, 0, comm);

   braid_Destroy(core);
   MPI_Barrier(comm);
   double total_time = max_elapsed(comm, total_start);

   if (rank == 0)
   {
      printf("%d %d %.14e %.14e %.14e %.14e %d\n", rf, nprocs, pilot_time, prerefine_time, drive_time, total_time, mgrit_iterations);
      fflush(stdout);
   }

   free(app->tgrid);
   free(app->tcoarse);
   free(app->ucoarse);
   free_workspace(app);
   free(app);
   MPI_Finalize();
   return 0;
}
