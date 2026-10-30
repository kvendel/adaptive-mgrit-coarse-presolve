
#include <errno.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "braid.h"

#define VecSize 2
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define VDP_MU 2.0

typedef double VEC[VecSize];

typedef struct _braid_App_struct
{
   int       rank;
   double    tstart;
   double    eps;
   double    eps0;
   int       ntgrid;
   double   *tgrid;
   VEC      *ugrid;
} my_App;

typedef struct _braid_Vector_struct
{
   VEC values;
} my_Vector;

void f(double t, double *y, double *k)
{
   (void)t;
   k[0] = y[1];
   k[1] = VDP_MU*(1.0-y[0]*y[0])*y[1] - y[0];
}

double hermite_time(double *tgrid, int ntgrid, int i, double theta)
{
   double h, hm, hp, d0, d1, a, b, c;

   h = tgrid[i+1] - tgrid[i];

   if (i > 0) hm = tgrid[i] - tgrid[i-1];
   else       hm = h;

   if (i < ntgrid-2) hp = tgrid[i+2] - tgrid[i+1];
   else              hp = h;

   d0 = 2*hm*h/(hm+h);
   d1 = 2*h*hp/(h+hp);

   a = d0;
   b = 3*h-2*d0-d1;
   c = d0+d1-2*h;

   return tgrid[i] + a*theta + b*theta*theta + c*theta*theta*theta;
}

int find_time_index(const my_App *app, double t)
{
   int left = 0;
   int right = app->ntgrid - 1;

   if (t <= app->tgrid[0]) return 0;
   if (t >= app->tgrid[right]) return right;

   while (left <= right)
   {
      int mid = left + (right-left)/2;
      double tm = app->tgrid[mid];
      if (tm < t) left = mid + 1;
      else if (tm > t) right = mid - 1;
      else return mid;
   }

   if (fabs(app->tgrid[left] - t) < fabs(t - app->tgrid[left-1])) return left;
   return left-1;
}

int my_TimeGrid(braid_App app, braid_Real *ta,
                braid_Int *ilower_ptr, braid_Int *iupper_ptr)
{
   braid_Int ilower = *ilower_ptr;
   braid_Int iupper = *iupper_ptr;

   for (braid_Int i=ilower; i<=iupper; i++)
      ta[i-ilower] = app->tgrid[i];

   return 0;
}

void bogacki_shampine_step(double t, double *y, double dt, double *ynext, double *err_ptr)
{
   double ystage[VecSize];
   double k1[VecSize];
   double k2[VecSize];
   double k3[VecSize];
   double k4[VecSize];
   double err = 0.0;
   f(t, y, k1);
   for (int k=0; k<VecSize; k++)
   {
      ystage[k] = y[k] + 0.5*dt*k1[k];
   }
   f(t+0.5*dt, ystage, k2);
   for (int k=0; k<VecSize; k++)
   {
      ystage[k] = y[k] + 0.75*dt*k2[k];
   }
   f(t+0.75*dt, ystage, k3);
   for (int k=0; k<VecSize; k++)
   {
      ynext[k] = y[k] + dt*(2.0/9.0*k1[k] + 1.0/3.0*k2[k] + 4.0/9.0*k3[k]);
   }
   f(t+dt, ynext, k4);
   for (int k=0; k<VecSize; k++)
   {
      double difference = dt*(-5.0/72.0*k1[k] + 1.0/12.0*k2[k] + 1.0/9.0*k3[k] - 1.0/8.0*k4[k]);
      err += difference*difference;
   }
   *err_ptr = sqrt(err);
}

double det_dt0(double *y, double eps0)
{
   double F[VecSize];
   double y1[VecSize];
   double F1[VecSize];

   double h1;
   double d0, d1, d2;
   double normF = 0.0;
   double normDiff = 0.0;

   double p = 2.0;

   f(0.0, y, F);

   d0 = 0.0;
   for (int k=0; k<VecSize; k++)
      d0 += y[k]*y[k];
   d0 = sqrt(d0)/(sqrt((double)VecSize)*eps0);

   for (int k=0; k<VecSize; k++)
      normF += F[k]*F[k];
   d1 = sqrt(normF)/(sqrt((double)VecSize)*eps0);

   if (d0 < 1.0e-5 || d1 < 1.0e-5)
      h1 = 1.0e-6;
   else
      h1 = 0.01*d0/d1;

   for (int k=0; k<VecSize; k++)
      y1[k] = y[k] + h1*F[k];

   f(h1, y1, F1);

   for (int k=0; k<VecSize; k++)
      normDiff += (F1[k]-F[k])*(F1[k]-F[k]);

   d2 = sqrt(normDiff)/(h1*sqrt((double)VecSize)*eps0);

   double dmax = (d1 > d2) ? d1 : d2;

   double h2;
   if (dmax <= 1.0e-15)
      h2 = (1.0e-6 > h1*0.001) ? 1.0e-6 : h1*0.001;
   else
      h2 = pow(0.01/dmax, 1.0/(p+1.0));

   return (100.0*h1 < h2) ? 100.0*h1 : h2;
}

void generate_adaptive_grid(double tstart, double tstop, double eps0, double dt0, int *ntgrid_ptr, double **tgrid_ptr, VEC **ugrid_ptr)
{
   int capacity = 1024;
   int ntgrid = 1;
   int first_accepted = 1;

   double *tgrid = (double *)malloc(capacity * sizeof(double));
   VEC *ugrid = (VEC *)malloc(capacity * sizeof(VEC));

   double t = tstart;

   double y[VecSize];
   double ynext[VecSize];

   double q = 3.0;

   double kI = 0.70 / q;
   double kP = -0.40 / q;

   double safety = 0.90;
   double err_old = eps0;
   double tiny = 1.0e-300;

   double dt_min = 1.0e-14 * (tstop - tstart);

   if (tgrid == NULL || ugrid == NULL)
   {
      free(tgrid);
      free(ugrid);

      *ntgrid_ptr = 0;
      *tgrid_ptr = NULL;
      *ugrid_ptr = NULL;

      return;
   }

   y[0] = 2.0;
   y[1] = 0.0;

   tgrid[0] = tstart;

   for (int k = 0; k < VecSize; k++)
   {
      ugrid[0][k] = y[k];
   }

   double dt = det_dt0(y,eps0);
   double err;
   bogacki_shampine_step(t, y, dt, ynext, &err);
   dt = safety * pow(eps0 / (err + tiny), 1.0 / q) * dt;
   while (t < tstop)
   {
      double err;
      double fac;
      double dt_trial = dt;
      int final_step = 0;

      bogacki_shampine_step(t, y, dt_trial, ynext, &err);

      if (!isfinite(err))
      {
         dt = 0.2 * dt_trial;
         continue;
      }

      if (err <= eps0 || dt_trial <= dt_min)
      {
         t += dt_trial;

         for (int k = 0; k < VecSize; k++)
         {
            y[k] = ynext[k];
         }

         if (ntgrid >= capacity)
         {
            double *new_tgrid;
            VEC *new_ugrid;

            capacity *= 2;

            new_tgrid = (double *)realloc(tgrid,capacity * sizeof(double));
            new_ugrid = (VEC *)realloc(ugrid,capacity * sizeof(VEC));

            if (new_tgrid == NULL || new_ugrid == NULL)
            {
               fprintf(stderr, "Memory allocation failed in " "generate_adaptive_grid().\n");
               exit(EXIT_FAILURE);
            }

            tgrid = new_tgrid;
            ugrid = new_ugrid;
         }

         tgrid[ntgrid] = t;

         for (int k = 0; k < VecSize; k++)
         {
            ugrid[ntgrid][k] = y[k];
         }

         ntgrid++;

	 if (t >= tstop)
	 {
            break;
	 }

         if (final_step)
         {
            err_old = err;
            break;
         }

         if (first_accepted)
         {
            
            fac = safety * pow(eps0 / (err + tiny), 1.0 / q);
            first_accepted = 0;
         }
         else
         {
            fac = safety * pow(eps0 / (err + tiny), kI) * pow(eps0  / (err_old + tiny), kP);
	    fac = safety * pow(eps0 / (err + tiny), 1.0 / q);
         }

         if (!isfinite(fac) || fac <= 0.0)
         {
            fac = 0.5;
         }

         if (fac > 4.0)
         {
            fac = 4.0;
         }

         if (fac < 0.25)
         {
            fac = 0.25;
         }

         err_old = fmax(err, tiny);

         dt = dt_trial * fac;
      }
      else
      {
         
         fac = safety * pow(eps0 / (err + tiny), 1.0 / q);

         if (!isfinite(fac) || fac <= 0.0)
         {
            fac = 0.1;
         }

         if (fac < 0.1)
         {
            fac = 0.1;
         }

         dt = dt_trial * fac;
      }
   }

   *ntgrid_ptr = ntgrid;
   *tgrid_ptr = tgrid;
   *ugrid_ptr = ugrid;
}

int compute_ref_factor(double eps0, double eps)
{
   int rf;
   double exponent = 1.0/3.0;

   if (eps <= 0.0 || eps0 <= eps)
   {
      return 1;
   }
   rf = (int)ceil(pow(eps0/eps, exponent));
   if (rf < 1)
   {
      rf = 1;
   }
   return rf;
}

void refine_grid_hermite(int ntcoarse, double *tcoarse, VEC *ucoarse, int rf,
                         int *ntfine_ptr, double **tfine_ptr, VEC **ufine_ptr)
{
   int ntfine = (ntcoarse-1)*rf + 1;
   double *tfine = (double *) malloc(ntfine*sizeof(double));
   VEC *ufine = (VEC *) malloc(ntfine*sizeof(VEC));

   #pragma omp parallel for if(ntcoarse > 64)
   for (int i=0; i<ntcoarse-1; i++)
   {
      for (int j=0; j<rf; j++)
      {
         int id = i*rf + j;
         double theta = ((double)j)/((double)rf);
         double h = tcoarse[i+1] - tcoarse[i];
         double alpha;
         tfine[id] = hermite_time(tcoarse, ntcoarse, i, theta);
         alpha = (fabs(h) < 1.0e-300) ? 0.0 : (tfine[id]-tcoarse[i])/h;
         for (int k=0; k<VecSize; k++)
         {
            ufine[id][k] = (1.0-alpha)*ucoarse[i][k] + alpha*ucoarse[i+1][k];
         }
      }
   }

   tfine[ntfine-1] = tcoarse[ntcoarse-1];
   for (int k=0; k<VecSize; k++)
   {
      ufine[ntfine-1][k] = ucoarse[ntcoarse-1][k];
   }

   *ntfine_ptr = ntfine;
   *tfine_ptr = tfine;
   *ufine_ptr = ufine;
}

void truncate_grid_at_tstop(double tstop, int *ntgrid_ptr, double **tgrid_ptr, VEC **ugrid_ptr)
{
   int ntgrid = *ntgrid_ptr;
   double *tgrid = *tgrid_ptr;
   VEC *ugrid = *ugrid_ptr;
   int right = 0;
   double tolerance =1.0e-14 * fmax(1.0, fabs(tstop));

   while (right < ntgrid && tgrid[right] < tstop)
   {
      right++;
   }
   if (right >= ntgrid)
   {
      fprintf(stderr,
              "The prerefined grid does not reach tstop.\n");

      exit(EXIT_FAILURE);
   }

   if (fabs(tgrid[right] - tstop) <= tolerance)
   {
      double *new_tgrid;
      VEC *new_ugrid;

      tgrid[right] = tstop;
      ntgrid = right + 1;

      new_tgrid = (double *)realloc(tgrid, ntgrid * sizeof(double));

      new_ugrid = (VEC *)realloc(ugrid, ntgrid * sizeof(VEC));

      if (new_tgrid != NULL)
      {
         tgrid = new_tgrid;
      }

      if (new_ugrid != NULL)
      {
         ugrid = new_ugrid;
      }

      *ntgrid_ptr = ntgrid;
      *tgrid_ptr = tgrid;
      *ugrid_ptr = ugrid;

      return;
   }

   if (right == 0)
   {
      fprintf(stderr, "tstop lies before the first grid point.\n");

      exit(EXIT_FAILURE);
   }

   {
      int left = right - 1;
      double interval_length = tgrid[right] - tgrid[left];
      double beta = (tstop - tgrid[left]) / interval_length;

      for (int k = 0; k < VecSize; k++)
      {
         ugrid[right][k] = (1.0 - beta) * ugrid[left][k] + beta * ugrid[right][k];
      }
      tgrid[right] = tstop;
   }
   ntgrid = right + 1;

   {
      double *new_tgrid;
      VEC *new_ugrid;

      new_tgrid = (double *)realloc(tgrid,ntgrid * sizeof(double));

      new_ugrid = (VEC *)realloc(ugrid,ntgrid * sizeof(VEC));

      if (new_tgrid != NULL)
      {
         tgrid = new_tgrid;
      }

      if (new_ugrid != NULL)
      {
         ugrid = new_ugrid;
      }
   }

   *ntgrid_ptr = ntgrid;
   *tgrid_ptr = tgrid;
   *ugrid_ptr = ugrid;
}

int my_Step(braid_App app, braid_Vector ustop, braid_Vector fstop,
            braid_Vector u, braid_StepStatus status)
{
   double tstart;
   double tstop;
   int istart;
   int iter;
   int level;
   int nrefine;

   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   braid_StepStatusGetTIndex(status, &istart);
   braid_StepStatusGetIter(status, &iter);
   braid_StepStatusGetLevel(status, &level);
   braid_StepStatusGetNRefine(status, &nrefine);

   double deltaT = tstop - tstart;

   double u_new[VecSize];
   double error;
   bogacki_shampine_step(tstart, u->values, deltaT, u_new, &error);

   for (int k=0; k<VecSize; k++)
   {
      u->values[k] = u_new[k];
   }

   (void)fstop;
   (void)istart;
   (void)iter;
   (void)level;
   (void)nrefine;
   return 0;
}

int my_Init(braid_App app, double t, braid_Vector *u_ptr)
{
   int index = find_time_index(app, t);
   my_Vector *u = (my_Vector *) malloc(sizeof(my_Vector));
   if (u == NULL) return 1;

   for (int k=0; k<VecSize; k++) u->values[k] = app->ugrid[index][k];
   *u_ptr = u;
   return 0;
}

int my_Clone(braid_App app, braid_Vector u, braid_Vector *v_ptr)
{
   my_Vector *v = (my_Vector *) malloc(sizeof(my_Vector));
   for (int i=0; i<VecSize; i++)
   {
      v->values[i] = u->values[i];
   }
   *v_ptr = v;
   return 0;
}

int my_Free(braid_App app, braid_Vector u)
{
   free(u);
   return 0;
}

int my_Sum(braid_App app, double alpha, braid_Vector x, double beta, braid_Vector y)
{
   for (int i=0; i<VecSize; i++)
   {
      y->values[i] = alpha*x->values[i] + beta*y->values[i];
   }
   return 0;
}

int my_InnerProd(braid_App app, braid_Vector u, braid_Vector v, double *prod_ptr)
{
   double dot = 0.0;
   for (int i=0; i<VecSize; i++)
   {
      dot += u->values[i]*v->values[i];
   }
   *prod_ptr = dot;
   return 0;
}

int my_SpatialNorm(braid_App app, braid_Vector u, double *norm_ptr)
{
   my_InnerProd(app, u, u, norm_ptr);
   *norm_ptr = sqrt(*norm_ptr);
   return 0;
}

int my_Access(braid_App app, braid_Vector u, braid_AccessStatus status) {
    int index, done;
    char filename[512];
    braid_AccessStatusGetTIndex(status, &index);
    braid_AccessStatusGetDone(status, &done);
    if (!done || index != app->ntgrid - 1) return 0;
    snprintf(filename, sizeof(filename), "results/vdp-BS-%e-%e.txt", app->eps0, app->eps);
    FILE *file = fopen(filename, "w");
    if (file == NULL) return 1;
    fprintf(file, "%.14e %.14e\n", u->values[0], u->values[1]);
    fclose(file);
    return 0;
}
int my_BufSize(braid_App app, int *size_ptr, braid_BufferStatus bstatus)
{
   *size_ptr = VecSize*sizeof(double);
   return 0;
}

int my_BufPack(braid_App app, braid_Vector u, void *buffer, braid_BufferStatus bstatus)
{
   double *dbuffer = buffer;
   for (int i=0; i<VecSize; i++)
   {
      dbuffer[i] = u->values[i];
   }
   braid_BufferStatusSetSize(bstatus, VecSize*sizeof(double));
   return 0;
}

int my_BufUnpack(braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus bstatus)
{
   double *dbuffer = buffer;
   my_Vector *u = (my_Vector *) malloc(sizeof(my_Vector));
   for (int i=0; i<VecSize; i++)
   {
      u->values[i] = dbuffer[i];
   }
   *u_ptr = u;
   return 0;
}

int my_BufAlloc(braid_App app, void **buffer, braid_Int nbytes, braid_BufferStatus bstatus)
{
   *buffer = malloc(nbytes);
   return 0;
}

int my_BufFree(braid_App app, void **buffer)
{
   free((char *) *buffer);
   *buffer = NULL;
   return 0;
}

static int ensure_results_directory(void) {
   if (mkdir("results", 0777) == 0 || errno == EEXIST) return 0;
   return 1;
}
int main(int argc, char *argv[])
{
   braid_Core core;
   MPI_Comm comm;
   my_App *app;
   double tstart = 0.0;
   double tstop = 10.0;
   double dt0;
   int ntime;

   int max_levels = 15;
   int min_coarse = 64;
   int nrelax = 1;
   int nrelax0 = -1;
   int nrelaxc = 1;
   int cfactor = 2;
   int max_iter = 100;
   int skip = 0;
   int bufalloc = 0;
   int storage = -1;
   
   double tol = 1.0e-09;
   double eps = 0.01;
   double eps0 = 1.0e-3;

   int rank;
   int arg_index;
   int ntcoarse;
   int ntgrid;
   int rfactor;
   double *tcoarse = NULL;
   double *tgrid = NULL;
   VEC *ucoarse = NULL;
   VEC *ugrid = NULL;

   comm = MPI_COMM_WORLD;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(comm, &rank);
   if (ensure_results_directory() != 0) { MPI_Finalize(); return 1; }
   dt0 = (tstop-tstart)/32.0;

   arg_index = 1;
   while (arg_index < argc)
   {
      if (strcmp(argv[arg_index], "-help") == 0)
      {
         if (rank == 0)
         {
            printf("\nVDP prerefined solver: VDP mu=2, Hermite prerefined physical grid, PI initial controller\n\n");
            printf("  -dt0, -eps0, -eps, -ml, -mc, -nu, -nu0, -nuc, -tol, -cf, -mi, -skip, -storage\n");
         }
         exit(1);
      }
      else if (strcmp(argv[arg_index], "-dt0") == 0)
      {
         arg_index++;
         dt0 = atof(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-eps0") == 0)
      {
         arg_index++;
         eps0 = atof(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-eps") == 0)
      {
         arg_index++;
         eps = atof(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-ml") == 0)
      {
         arg_index++;
         max_levels = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-storage") == 0)
      {
         arg_index++;
         storage = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-mc") == 0)
      {
         arg_index++;
         min_coarse = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-nu") == 0)
      {
         arg_index++;
         nrelax = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-nu0") == 0)
      {
         arg_index++;
         nrelax0 = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-nuc") == 0)
      {
         arg_index++;
         nrelaxc = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-tol") == 0)
      {
         arg_index++;
         tol = atof(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-cf") == 0)
      {
         arg_index++;
         cfactor = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-mi") == 0)
      {
         arg_index++;
         max_iter = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-bufalloc") == 0)
      {
         arg_index++;
         bufalloc = 1;
      }
      else if (strcmp(argv[arg_index], "-skip") == 0)
      {
         arg_index++;
         skip = atoi(argv[arg_index++]);
      }
      else
      {
         arg_index++;
      }
   }

   if (rank == 0)
   {
      generate_adaptive_grid(tstart, tstop, eps0, dt0, &ntcoarse, &tcoarse, &ucoarse);
      printf("Initial adaptive VDP grid: eps0=%e, ntgrid=%d, intervals=%d\n",
             eps0, ntcoarse, ntcoarse-1);
   }

   MPI_Bcast(&ntcoarse, 1, MPI_INT, 0, comm);
   if (rank != 0)
   {
      tcoarse = (double *) malloc(ntcoarse*sizeof(double));
      ucoarse = (VEC *) malloc(ntcoarse*sizeof(VEC));
   }
   MPI_Bcast(tcoarse, ntcoarse, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&ucoarse[0][0], ntcoarse*VecSize, MPI_DOUBLE, 0, comm);

   rfactor = compute_ref_factor(eps0, eps);
   
   refine_grid_hermite(ntcoarse, tcoarse, ucoarse, rfactor, &ntgrid, &tgrid, &ugrid);
   truncate_grid_at_tstop(tstop,&ntgrid,&tgrid,&ugrid);
   ntime = ntgrid - 1;

   if (rank == 0)
   {
      char filename[400];
      printf("Hermite pre-refinement: eps0=%e, eps=%e, rfactor=%d, ntgrid=%d, intervals=%d\n",
             eps0, eps, rfactor, ntgrid, ntgrid-1);
      snprintf(filename, sizeof(filename), "results/vdp-BS-%e-%e.txt", eps0, eps);
      remove(filename);
   }

   free(tcoarse);
   free(ucoarse);

   app = (my_App *) malloc(sizeof(my_App));
   app->rank = rank;
   app->tstart = tstart;
   app->eps = eps;
   app->eps0 = eps0;
   app->ntgrid = ntgrid;
   app->tgrid = tgrid;
   app->ugrid = ugrid;

   braid_Init(comm, comm, tgrid[0], tgrid[ntgrid-1], ntime, app,
              my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm,
              my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
   braid_SetTimeGrid(core, my_TimeGrid);

   braid_SetPrintLevel(core, 2);
   braid_SetMaxLevels(core, max_levels);
   braid_SetNRelax(core, -1, nrelax);
   if (nrelax0 > -1)
   {
      braid_SetNRelax(core, 0, nrelax0);
   }
   braid_SetNRelax(core, max_levels-1, nrelaxc);
   braid_SetMinCoarse(core, min_coarse);
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, cfactor);
   braid_SetMaxIter(core, max_iter);
   braid_SetSkip(core, skip);

   if (bufalloc)
   {
      braid_SetBufAllocFree(core, my_BufAlloc, my_BufFree);
   }
   if (storage >= -2)
   {
      braid_SetStorage(core, storage);
   }

   clock_t run_time = clock();
   braid_Drive(core);
   run_time = clock() - run_time;
   if (rank == 0)
   {
      printf("Braid took %f seconds to execute \n", ((double)run_time)/CLOCKS_PER_SEC);
   }

   braid_Destroy(core);
   free(app->tgrid);
   free(app->ugrid);
   free(app);
   MPI_Finalize();
   return 0;
}


