#include <errno.h>
#include <sys/stat.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define VecSize 2
#define VDP_MU 2.0

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef double VEC[VecSize];

void f(double t, const double *y, double *k)
{
   (void)t;

   k[0] = y[1];
   k[1] = VDP_MU * (1.0 - y[0] * y[0]) * y[1] - y[0];
}

void bogacki_shampine_step(double t, const double *y, double dt, double *ynext, double *error_ptr)
{
   VEC ystage, k1, k2, k3, k4;
   double error = 0.0;
   f(t, y, k1);
   for (int k=0; k<VecSize; k++) ystage[k] = y[k] + 0.5*dt*k1[k];
   f(t+0.5*dt, ystage, k2);
   for (int k=0; k<VecSize; k++) ystage[k] = y[k] + 0.75*dt*k2[k];
   f(t+0.75*dt, ystage, k3);
   for (int k=0; k<VecSize; k++) ynext[k] = y[k] + dt*(2.0/9.0*k1[k] + 1.0/3.0*k2[k] + 4.0/9.0*k3[k]);
   f(t+dt, ynext, k4);
   for (int k=0; k<VecSize; k++)
   {
      double difference = dt*(-5.0/72.0*k1[k] + 1.0/12.0*k2[k] + 1.0/9.0*k3[k] - 1.0/8.0*k4[k]);
      error += difference*difference;
   }
   *error_ptr = sqrt(error);
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

int generate_adaptive_solution(double tstart, double tstop, double tolerance, double dt0, const char *filename,
                               VEC solution, int *number_of_intervals_ptr)
{
   int number_of_intervals = 0;
   int first_accepted = 1;
   double t = tstart;
   VEC y = {2.0, 0.0}, ynext;
   FILE *file = fopen(filename, "w");
   if (file == NULL) return 1;
   fprintf(file, "%.14e %.14e %.14e %d\n", y[0], y[1], t, number_of_intervals);
   const double q = 3.0;
   const double kI = 0.70 / q;
   const double kP = -0.40 / q;
   const double safety = 0.90;
   const double tiny = 1.0e-300;
   const double dt_min = 1.0e-14 * (tstop - tstart);
   double error_old = tolerance;
   double dt = det_dt0(y, tolerance);
   double err;
   (void)dt0;
   bogacki_shampine_step(t, y, dt, ynext, &err);
   dt = 0.9 * pow(tolerance / (err + tiny), 1.0 / q) * dt;
   while (t < tstop)
   {
      double error;
      double factor;
      double dt_trial = dt;
      if (t + dt_trial > tstop) dt_trial = tstop - t;
      bogacki_shampine_step(t, y, dt_trial, ynext, &error);
      if (!isfinite(error))
      {
         dt = 0.2 * dt_trial;
         continue;
      }
      if (error <= tolerance || dt_trial <= dt_min)
      {
         t += dt_trial;
         y[0] = ynext[0];
         y[1] = ynext[1];
         number_of_intervals++;
         fprintf(file, "%.14e %.14e %.14e %d\n", y[0], y[1], t, number_of_intervals);
         if (t >= tstop) break;
         if (first_accepted)
         {
            factor = safety * pow(tolerance / (error + tiny), 1.0 / q);
            first_accepted = 0;
         }
         else
         {
            factor = safety * pow(tolerance / (error + tiny), kI)
                            * pow(tolerance / (error_old + tiny), kP);
	    factor = safety * pow(tolerance / (error + tiny), 1.0 / q);
         }
         if (!isfinite(factor) || factor <= 0.0) factor = 0.5;
         if (factor > 4.0) factor = 4.0;
         if (factor < 0.25) factor = 0.25;
         error_old = fmax(error, tiny);
         dt = dt_trial * factor;
      }
      else
      {
         factor = safety * pow(tolerance / (error + tiny), 1.0 / q);
         if (!isfinite(factor) || factor <= 0.0) factor = 0.1;
         if (factor < 0.1) factor = 0.1;
         dt = dt_trial * factor;
      }
   }
   solution[0] = y[0];
   solution[1] = y[1];
   *number_of_intervals_ptr = number_of_intervals;
   return 0;
}


int main(int argc, char *argv[])
{
   double tstart = 0.0;
   double tstop = 10.0;
   double tolerance = 1.0e-3;
   double dt0 = -1.0;
   char output_filename[400];
   snprintf(output_filename, sizeof(output_filename), "results/vdp-sequential-%e.txt", tolerance);
   int number_of_intervals;
   VEC solution;
   clock_t start_time = clock();

   for (int i = 1; i < argc; i++)
   {
      if (strcmp(argv[i], "-tol") == 0 && i + 1 < argc)
      {
         tolerance = atof(argv[++i]);
      }
      else if (strcmp(argv[i], "-dt0") == 0 && i + 1 < argc)
      {
         dt0 = atof(argv[++i]);
      }
      else if (strcmp(argv[i], "-tstop") == 0 && i + 1 < argc)
      {
         tstop = atof(argv[++i]);
      }
      else
      {
         fprintf(stderr, "Unknown or incomplete option: %s\n", argv[i]);
         return 1;
      }
   }

   if (mkdir("results", 0777) != 0 && errno != EEXIST) return 1;
   snprintf(output_filename, sizeof(output_filename), "results/vdp-sequential_BS-%e.txt", tolerance);
   if (dt0 <= 0.0)
   {
      dt0 = (tstop - tstart) / 32.0;
   }

   if (generate_adaptive_solution(tstart, tstop, tolerance, dt0, output_filename, solution, &number_of_intervals) != 0)
   {
      return 1;
   }


   printf("Sequential VDP Bogacki-Shampine solve finished.\n");
   printf("VDP mu:                %.1f\n", VDP_MU);
   printf("Tolerance:             %.6e\n", tolerance);
   printf("Initial trial step:    %.14e\n", dt0);
   printf("Number of intervals:   %d\n", number_of_intervals);
   printf("Output file:           %s\n", output_filename);
   printf("CPU time:              %.6f seconds\n",
          (double)(clock() - start_time) / CLOCKS_PER_SEC);


   return 0;
}


