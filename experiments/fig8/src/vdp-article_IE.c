
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "braid.h"

#define VecSize 2
#define VDP_MU 2.0
typedef double VEC[VecSize];

/*--------------------------------------------------------------------------
 * User-defined routines and structures
 *--------------------------------------------------------------------------*/

/* App structure can contain anything, and be named anything as well */
typedef struct _braid_App_struct
{
   int       rank;
   double    tstart;
   double    eps;

} my_App;

/* Vector structure can contain anything, and be name anything as well */
typedef struct _braid_Vector_struct
{
   VEC values;

} my_Vector;

void f(double t, double* y, double* k)
{
   k[0] = y[1];
   k[1] = 2*(1-y[0]*y[0])*y[1]-y[0];
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


int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   double tstart;
   double tstop;
   int    istart;
   int iter;
   int level;
   int nrefine;
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   braid_StepStatusGetTIndex(status, &istart);
   braid_StepStatusGetIter (status, &iter);
   braid_StepStatusGetLevel (status, &level);
   braid_StepStatusGetNRefine(status, &nrefine);

   double deltaT = tstop - tstart;
   double u_new[VecSize];
   double loc_err;

   implicit_euler_step_with_error(tstart, u->values, deltaT, u_new, &loc_err);

   if (level == 0 && iter>3 && (nrefine+iter) % 4 == 0)
   {
      int rf = ceil(pow(loc_err/app->eps, 1./2.)-0.4);
      braid_StepStatusSetRFactor(status, rf);
   }

   for (int k=0; k<VecSize; k++)
   {
      u->values[k] = u_new[k];
   }

   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   u->values[0] = 2.0;
   u->values[1] = 0.0;
   *u_ptr = u;

   return 0;
}

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   my_Vector *v;

   v = (my_Vector *) malloc(sizeof(my_Vector));
   for (int i=0;i<VecSize;i++)
   {
   (v->values[i]) = (u->values[i]);
   }
   *v_ptr = v;

   return 0;
}

int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u);

   return 0;
}

int
my_Sum(braid_App app,
           double alpha,
           braid_Vector x,
           double beta,
           braid_Vector y)
{
   for (int i = 0; i < VecSize; i++)
   {
      y->values[i] = alpha * x->values[i] + beta * y->values[i];
   }

   return 0;
}


int my_InnerProd(braid_App app, braid_Vector u, braid_Vector v, double *prod_ptr)
{
   double dot = 0.;

   for (int i = 0; i < VecSize; i++)
   {
      dot += (u->values[i]) * (v->values[i]);
   }
   *prod_ptr = dot;
   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{

   my_InnerProd(app, u, u, norm_ptr);
   *norm_ptr = sqrt(*norm_ptr);

   return 0;
}

int my_Access(braid_App app, braid_Vector u, braid_AccessStatus status) {
    int done;
    double t;
    char filename[400];
    braid_AccessStatusGetDone(status, &done);
    braid_AccessStatusGetT(status, &t);
    if (!done || fabs(t - 10.0) > 1.0e-12) return 0;
    snprintf(filename, sizeof(filename), "results/vdp-article_IE-%e.txt", app->eps);
    FILE *file = fopen(filename, "w");
    if (file == NULL) return 1;
    fprintf(file, "%.14e %.14e\n", u->values[0], u->values[1]);
    fclose(file);
    return 0;
}
int
my_BufSize(braid_App          app,
           int                *size_ptr,
           braid_BufferStatus bstatus)
{
   *size_ptr = VecSize * sizeof(double);
   return 0;
}

int
my_BufPack(braid_App          app,
           braid_Vector       u,
           void               *buffer,
           braid_BufferStatus bstatus)
{
   double *dbuffer = buffer;

   for (int i = 0; i < VecSize; i++)
   {
      dbuffer[i] = (u->values[i]);
   }
   braid_BufferStatusSetSize(bstatus, VecSize * sizeof(double));

   return 0;
}

int
my_BufUnpack(braid_App          app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus bstatus)
{
   double    *dbuffer = buffer;
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   for (int i = 0; i < VecSize; i++)
   {
      (u->values[i]) = dbuffer[i];
   }
   *u_ptr = u;

   return 0;
}

int
my_BufAlloc(braid_App          app,
            void               **buffer,
            braid_Int          nbytes,
            braid_BufferStatus bstatus)
{
   *buffer = malloc(nbytes);
   return 0;
}

int
my_BufFree(braid_App          app,
           void               **buffer)
{
   free((char *) *buffer);
   *buffer = NULL;
   return 0;
}

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
   braid_Core    core;
   MPI_Comm      comm;
   my_App       *app;
   double        tstart, tstop;
   int           ntime;

   int           max_levels    = 15;
   int           min_coarse    = 128;
   int           nrelax        = 1;
   int           nrelax0       = -1;
   int           nrelaxc       = 1;
   double        tol           = 1.0e-11;
   int           cfactor       = 2;
   int           max_iter      = 100;
   int           skip          = 0;
   int           bufalloc      = 0;
   int           storage       = -1;
   double        eps           = 0.01;

   int           arg_index;
   int           rank;

   /* Initialize MPI */
   comm   = MPI_COMM_WORLD;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(comm, &rank);

   /* Define time domain: ntime intervals */
   ntime  = 32;
   tstart = 0.0;
   tstop  = 10;
   
   /* Parse command line */
   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         if ( rank == 0 )
         {
            printf("\nExample 1: Solve a scalar ODE \n\n");
            printf("  -ntime <ntime>    : set num time points\n");
            printf("  -ml  <max_levels> : set max levels\n");
            printf("  -nu  <nrelax>     : set num F-C relaxations\n");
            printf("  -nu0 <nrelax>     : set num F-C relaxations on level 0\n");
            printf("  -nuc <nrelax>     : set num F-C relaxations on coarsest grid\n");
            printf("  -tol <tol>        : set stopping tolerance\n");
            printf("  -cf  <cfactor>    : set coarsening factor\n");
            printf("  -mi  <max_iter>   : set max iterations\n");
            printf("  -skip <set_skip>  : set skip relaxations on first down-cycle; 0: no skip;  1: skip\n");
            printf("  -fmg              : use FMG cycling\n");
            printf("  -bufalloc         : user-defined MPI buffer allocation\n");
            printf("  -st <stepper>     : set which stepper in use\n");
         }
         exit(1);
      } 
      else if ( strcmp(argv[arg_index], "-ntime") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ml") == 0 )
      {
         arg_index++;
         max_levels = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-storage") == 0 ){
         arg_index++;
         storage = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-mc") == 0 ){
          arg_index++;
          min_coarse = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu") == 0 )
      {
         arg_index++;
         nrelax = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu0") == 0 )
      {
         arg_index++;
         nrelax0 = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nuc") == 0 )
      {
         arg_index++;
         nrelaxc = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tol") == 0 ) 
      {
         arg_index++;
         tol = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-cf") == 0 )
      {
         arg_index++;
         cfactor = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mi") == 0 )
      {
         arg_index++;
         max_iter = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-bufalloc") == 0 )
      {
         arg_index++;
         bufalloc = 1;
      }
      else if ( strcmp(argv[arg_index], "-skip") == 0 )
      {
         arg_index++;
         skip = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-eps") == 0 )
      {
         arg_index++;
         eps = atof(argv[arg_index++]);
      }
      else
      {
         arg_index++;
      }
   }

   if (rank == 0)
   {
      char filename[400];
      snprintf(filename, sizeof(filename), "results/vdp-article_IE-%e.txt", eps);
      remove(filename);
   }

   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->rank)   = rank;
   (app->tstart) = tstart;
   (app->eps)    = eps;

   /* initialize XBraid and set options */
   braid_Init(comm, comm, tstart, tstop, ntime, app,
          my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
          my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

   braid_SetPrintLevel( core, 2);
   braid_SetMaxLevels(core, max_levels);
   braid_SetNRelax(core, -1, nrelax);
   if (nrelax0 > -1)
   {
      braid_SetNRelax(core,  0, nrelax0);
   }
   braid_SetNRelax(core, max_levels-1, nrelaxc);
   braid_SetMinCoarse(core, min_coarse);
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, cfactor);
   braid_SetMaxIter(core, max_iter);
   braid_SetSkip(core, skip);
   braid_SetRefine(core, 1);
   if (bufalloc)
   {
      braid_SetBufAllocFree(core, my_BufAlloc, my_BufFree);
   }
   if (storage >= -2) {
         braid_SetStorage(core, storage);
      }

   /* Run simulation, and then clean up */
   clock_t t; 
   t = clock(); 
   braid_Drive(core);
   t = clock() - t; 
   double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 

   if ( rank == 0)
   {
      printf("Braid took %f seconds to execute \n", time_taken);
   }

   braid_Destroy(core);
   free(app);
   MPI_Finalize();

   return (0);
}

