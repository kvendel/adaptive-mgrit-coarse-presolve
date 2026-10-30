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

static double cbrt_real(double x)
{
   if (x >= 0.0)
   {
      return pow(x, 1.0/3.0);
   }
   return -pow(-x, 1.0/3.0);
}

static int cubic_real_roots(double a, double b, double c, double d, double roots[3])
{
   double tol = 1.0e-14;

   if (fabs(a) < tol)
   {
      if (fabs(b) < tol)
      {
         if (fabs(c) < tol)
         {
            return 0;
         }
         roots[0] = -d/c;
         return 1;
      }

      double disc = c*c - 4.0*b*d;
      if (disc < -tol)
      {
         return 0;
      }
      if (disc < 0.0)
      {
         disc = 0.0;
      }
      roots[0] = (-c + sqrt(disc))/(2.0*b);
      roots[1] = (-c - sqrt(disc))/(2.0*b);
      if (disc == 0.0)
      {
         return 1;
      }
      return 2;
   }

   double A = b/a;
   double B = c/a;
   double C = d/a;
   double p = B - A*A/3.0;
   double q = 2.0*A*A*A/27.0 - A*B/3.0 + C;
   double shift = -A/3.0;
   double Delta = 0.25*q*q + p*p*p/27.0;

   if (Delta > tol)
   {
      double sd = sqrt(Delta);
      roots[0] = shift + cbrt_real(-0.5*q + sd) + cbrt_real(-0.5*q - sd);
      return 1;
   }
   else if (fabs(Delta) <= tol)
   {
      double z = cbrt_real(-0.5*q);
      roots[0] = shift + 2.0*z;
      roots[1] = shift - z;
      if (fabs(roots[0]-roots[1]) < 1.0e-12)
      {
         return 1;
      }
      return 2;
   }
   else
   {
      double arg;
      double phi;
      double r;

      if (fabs(p) < tol)
      {
         roots[0] = shift;
         return 1;
      }

      arg = (3.0*q/(2.0*p))*sqrt(-3.0/p);
      if (arg > 1.0) arg = 1.0;
      if (arg < -1.0) arg = -1.0;
      phi = acos(arg)/3.0;
      r = 2.0*sqrt(-p/3.0);

      roots[0] = shift + r*cos(phi);
      roots[1] = shift + r*cos(phi - 2.0*M_PI/3.0);
      roots[2] = shift + r*cos(phi - 4.0*M_PI/3.0);
      return 3;
   }
}

void explicit_euler_step(double t, double *y, double dt, double *ynext)
{
   double k[VecSize];
   f(t, y, k);
   for (int i=0; i<VecSize; i++)
   {
      ynext[i] = y[i] + dt*k[i];
   }
}

int ie_step(double t, double *u_old, double dt, double *guess, double *u)
{
   (void)t;
   double xn = u_old[0];
   double yn = u_old[1];
   double a = dt*VDP_MU;
   double b = -dt*VDP_MU*xn;
   double c = 1.0 - dt*VDP_MU + dt*dt;
   double d = -(1.0 - dt*VDP_MU)*xn - dt*yn;
   double roots[3];
   int nroots;
   double xref;
   double best;
   double best_dist;

   if (fabs(dt) < 1.0e-15)
   {
      u[0] = u_old[0];
      u[1] = u_old[1];
      return 0;
   }

   nroots = cubic_real_roots(a, b, c, d, roots);
   if (nroots <= 0)
   {
      return 1;
   }

   if (guess != NULL)
   {
      xref = guess[0];
   }
   else
   {
      xref = xn + dt*yn;
   }

   best = roots[0];
   best_dist = fabs(roots[0] - xref);
   for (int i=1; i<nroots; i++)
   {
      double dist = fabs(roots[i] - xref);
      if (dist < best_dist)
      {
         best = roots[i];
         best_dist = dist;
      }
   }

   u[0] = best;
   u[1] = (yn - dt * best) / (1.0 - dt * VDP_MU * (1.0 - best * best));
   return 0;
}

void ie_richardson_step_with_error(double t, double *y, double dt, double *ynext, double *err_ptr)
{
   double y_one[VecSize];
   double y_mid[VecSize];
   double y_two[VecSize];
   double guess_one[VecSize];
   double guess_mid[VecSize];
   double err = 0.0;

   explicit_euler_step(t, y, dt, guess_one);
   explicit_euler_step(t, y, 0.5*dt, guess_mid);

   ie_step(t, y, dt, guess_one, y_one);
   ie_step(t, y, 0.5*dt, guess_mid, y_mid);
   ie_step(t+0.5*dt, y_mid, 0.5*dt, y_one, y_two);

   for (int k=0; k<VecSize; k++)
   {
      ynext[k] = y_one[k];
      double difference = y_two[k] - y_one[k];
      err += 4.0 * difference * difference;
   }
   *err_ptr = sqrt(err);
}

void implicit_euler_step_with_error(double t, double *y, double dt, double *ynext, double *err_ptr) {
   ie_richardson_step_with_error(t, y, dt, ynext, err_ptr);
}

void implicit_euler_step(double t, double *y, double dt, double *guess, double *ynext) {
   double fe_guess[VecSize];
   if (guess == NULL) {
      explicit_euler_step(t, y, dt, fe_guess);
      guess = fe_guess;
   }
   ie_step(t, y, dt, guess, ynext);
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

   double p = 1.0;

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
   (void)dt0;
   int capacity = 1024;
   int ntgrid = 1;

   double *tgrid = (double *)malloc(capacity * sizeof(double));
   VEC *ugrid = (VEC *)malloc(capacity * sizeof(VEC));

   double t = tstart;

   double y[VecSize];
   double ynext[VecSize];

   double q = 2.0;

   double kI = 0.70 / q;
   double kP = -0.40 / q;
   (void)kI;
   (void)kP;

   double safety = 0.90;
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
   implicit_euler_step_with_error(t, y, dt, ynext, &err);
   dt = safety * pow(eps0 / (err + tiny), 1.0 / q) * dt;
   while (t < tstop)
   {
      double err;
      double fac;
      double dt_trial = dt;
      int final_step = 0;

      implicit_euler_step_with_error(t, y, dt_trial, ynext, &err);

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
            break;
         }
         fac = safety * pow(eps0 / (err + tiny), 1.0 / q);

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
   double exponent = 0.5;

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
   double *guess = NULL;
   if (ustop != NULL)
   {
      guess = ustop->values;
   }
   implicit_euler_step(tstart, u->values, deltaT, guess, u_new);

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

int my_Access(braid_App app, braid_Vector u, braid_AccessStatus astatus)
{
   int index;
   double t;
   FILE *file;

   braid_AccessStatusGetTIndex(astatus, &index);
   braid_AccessStatusGetT(astatus, &t);

   char filename[400];
   snprintf(filename, sizeof(filename), "results/vdp-hermite_harmonic_I-%e-%e.txt", app->eps0, app->eps);
   file = fopen(filename, "a");
   if (!file)
   {
      return 1;
   }
   for (int i=0; i<VecSize; i++)
   {
      fprintf(file, " %.14e", u->values[i]);
   }
   fprintf(file, " %.14e %d\n", t, index);
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
      snprintf(filename, sizeof(filename), "results/vdp-hermite_harmonic_I-%e-%e.txt", eps0, eps);
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


