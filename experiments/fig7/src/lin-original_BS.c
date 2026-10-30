
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "braid.h"

#define VecSize 2
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
   k[0] = -0.2*y[0];
   k[1] = -0.4*(y[1]-y[0]);
}

int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   int    istart;   /* time point index value corresponding to tstart on (global) fine grid */
   int iter;
   int level;
   int nrefine;
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   braid_StepStatusGetTIndex(status, &istart);
   braid_StepStatusGetIter (status, &iter);
   braid_StepStatusGetLevel (status, &level);
   braid_StepStatusGetNRefine(status, &nrefine);

   /* Use backward Euler to propagate solution */
   double deltaT = tstop - tstart;
   double *u_new;
   double *k1;
   double *k2;
   double *k3;
   double *k4;

   u_new=(double *) malloc(sizeof(double)*VecSize);
   k1=(double *) malloc(sizeof(double)*VecSize);
   k2=(double *) malloc(sizeof(double)*VecSize);
   k3=(double *) malloc(sizeof(double)*VecSize);
   k4=(double *) malloc(sizeof(double)*VecSize);

   for (int k=0; k<VecSize; k++)
   {
      u_new[k]=u->values[k];
   }
   f(tstart, u_new, k1);

   for (int k=0; k<VecSize; k++)
   {
      u_new[k]=u->values[k]+1./2.*deltaT*k1[k];
   }
   f(tstart+deltaT/2, u_new, k2);

   for (int k=0; k<VecSize; k++)
   {
      u_new[k]=u->values[k]+3./4.*deltaT*k2[k];
   }
   f(tstart+3./4.*deltaT, u_new, k3);

   for (int k=0; k<VecSize; k++)
   {
      u_new[k]=u->values[k]+deltaT*(2./9.*k1[k]+1./3.*k2[k]+4./9.*k3[k]);
   }
   f(tstart+deltaT, u_new, k4);

   for (int k=0; k<VecSize; k++)
   {
      u->values[k]=u->values[k]+deltaT*(2./9.*k1[k]+1./3.*k2[k]+4./9.*k3[k]);
   }


   if (level == 0 && iter>3 && (iter+nrefine) % 3 == 1)
   {
      double loc_err = 0;
      for (int k=0;k<VecSize;k++)
      {
         loc_err += pow(deltaT*(-5./72.*k1[k]+1./12.*k2[k]+1./9.*k3[k]-1./8.*k4[k]),2);
      }
      loc_err = sqrt(loc_err);
      int rf = ceil(pow(loc_err/app->eps,1./3.));
      braid_StepStatusSetRFactor(status, rf);
   }

   free(u_new);
   free(k1);
   free(k2);
   free(k3);
   free(k4);

   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   u->values[0] = 0.3;
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
    int index;
    double t;
    char filename[512];
    braid_AccessStatusGetTIndex(status, &index);
    braid_AccessStatusGetT(status, &t);
    snprintf(filename, sizeof(filename), "results/lin-original_BS-%e.txt", app->eps);
    FILE *file = fopen(filename, "a");
    if (file == NULL) return 1;
    for (int i = 0; i < VecSize; i++) fprintf(file, " %.14e", u->values[i]);
    fprintf(file, " %.14e %d\n", t, index);
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
   int           min_coarse    = 64;
   int           nrelax        = 1;
   int           nrelax0       = -1;
   int           nrelaxc       = 1;
   double        tol           = 1.0e-09;
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
   tstop  = 20;
   
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
      snprintf(filename, sizeof(filename), "results/lin-original_BS-%e.txt", eps);
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

