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

static double cbrt_real(double x)
{
   if (x >= 0.0)
   {
      return pow(x, 1.0 / 3.0);
   }

   return -pow(-x, 1.0 / 3.0);
}

static int cubic_real_roots(double a,
                            double b,
                            double c,
                            double d,
                            double roots[3])
{
   const double tolerance = 1.0e-14;

   if (fabs(a) < tolerance)
   {
      if (fabs(b) < tolerance)
      {
         if (fabs(c) < tolerance)
         {
            return 0;
         }

         roots[0] = -d / c;
         return 1;
      }

      double discriminant = c * c - 4.0 * b * d;

      if (discriminant < -tolerance)
      {
         return 0;
      }

      if (discriminant < 0.0)
      {
         discriminant = 0.0;
      }

      roots[0] = (-c + sqrt(discriminant)) / (2.0 * b);
      roots[1] = (-c - sqrt(discriminant)) / (2.0 * b);

      if (discriminant == 0.0)
      {
         return 1;
      }

      return 2;
   }

   double A = b / a;
   double B = c / a;
   double C = d / a;
   double p = B - A * A / 3.0;
   double q = 2.0 * A * A * A / 27.0 - A * B / 3.0 + C;
   double shift = -A / 3.0;
   double discriminant = 0.25 * q * q + p * p * p / 27.0;

   if (discriminant > tolerance)
   {
      double square_root = sqrt(discriminant);

      roots[0] = shift
               + cbrt_real(-0.5 * q + square_root)
               + cbrt_real(-0.5 * q - square_root);

      return 1;
   }
   else if (fabs(discriminant) <= tolerance)
   {
      double z = cbrt_real(-0.5 * q);

      roots[0] = shift + 2.0 * z;
      roots[1] = shift - z;

      if (fabs(roots[0] - roots[1]) < 1.0e-12)
      {
         return 1;
      }

      return 2;
   }
   else
   {
      double argument;
      double phi;
      double radius;

      if (fabs(p) < tolerance)
      {
         roots[0] = shift;
         return 1;
      }

      argument = (3.0 * q / (2.0 * p)) * sqrt(-3.0 / p);

      if (argument > 1.0)
      {
         argument = 1.0;
      }

      if (argument < -1.0)
      {
         argument = -1.0;
      }

      phi = acos(argument) / 3.0;
      radius = 2.0 * sqrt(-p / 3.0);

      roots[0] = shift + radius * cos(phi);
      roots[1] = shift + radius * cos(phi - 2.0 * M_PI / 3.0);
      roots[2] = shift + radius * cos(phi - 4.0 * M_PI / 3.0);

      return 3;
   }
}

void explicit_euler_step(double t,
                         const double *y,
                         double dt,
                         double *ynext)
{
   double derivative[VecSize];

   f(t, y, derivative);

   for (int i = 0; i < VecSize; i++)
   {
      ynext[i] = y[i] + dt * derivative[i];
   }
}

int ie_step(double t,
            const double *u_old,
            double dt,
            const double *guess,
            double *u)
{
   (void)t;

   double xn = u_old[0];
   double yn = u_old[1];
   double a = dt * VDP_MU;
   double b = -dt * VDP_MU * xn;
   double c = 1.0 - dt * VDP_MU + dt * dt;
   double d = -(1.0 - dt * VDP_MU) * xn - dt * yn;
   double roots[3];
   int number_of_roots;
   double reference;
   double best;
   double best_distance;

   if (fabs(dt) < 1.0e-15)
   {
      u[0] = u_old[0];
      u[1] = u_old[1];
      return 0;
   }

   number_of_roots = cubic_real_roots(a, b, c, d, roots);

   if (number_of_roots <= 0)
   {
      return 1;
   }

   if (guess != NULL)
   {
      reference = guess[0];
   }
   else
   {
      reference = xn + dt * yn;
   }

   best = roots[0];
   best_distance = fabs(roots[0] - reference);

   for (int i = 1; i < number_of_roots; i++)
   {
      double distance = fabs(roots[i] - reference);

      if (distance < best_distance)
      {
         best = roots[i];
         best_distance = distance;
      }
   }

   u[0] = best;
   u[1] = (yn - dt * best) / (1.0 - dt * VDP_MU * (1.0 - best * best));

   return 0;
}

void ie_richardson_step_with_error(double t,
                                   const double *y,
                                   double dt,
                                   double *ynext,
                                   double *error_ptr)
{
   VEC y_one;
   VEC y_mid;
   VEC y_two;
   VEC guess_one;
   VEC guess_mid;
   double error = 0.0;

   explicit_euler_step(t, y, dt, guess_one);
   explicit_euler_step(t, y, 0.5 * dt, guess_mid);

   if (ie_step(t, y, dt, guess_one, y_one) != 0 ||
       ie_step(t, y, 0.5 * dt, guess_mid, y_mid) != 0 ||
       ie_step(t + 0.5 * dt, y_mid, 0.5 * dt, y_one, y_two) != 0)
   {
      *error_ptr = INFINITY;
      return;
   }

   for (int k = 0; k < VecSize; k++)
   {
      double difference = 2.0*(y_two[k] - y_one[k]);

      ynext[k] = y_one[k];

      error += difference * difference;
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

int generate_adaptive_solution(double tstart, double tstop, double tolerance, double dt0,
                               VEC solution, int *number_of_intervals_ptr)
{
   int number_of_intervals = 0, first_accepted = 1;
   double t = tstart; VEC y = {2.0, 0.0}, ynext;
   const double q = 2.0, kI = 0.70 / q, kP = -0.40 / q, safety = 0.90, tiny = 1.0e-300;
   const double dt_min = 1.0e-14 * (tstop - tstart);
   double error_old = tolerance, dt = det_dt0(y, tolerance), err; (void)dt0;
   ie_richardson_step_with_error(t, y, dt, ynext, &err);
   dt = 0.9 * pow(tolerance / (err + tiny), 1.0 / q) * dt;
   while (t < tstop) {
      double error, factor, dt_trial = dt;
      if (t + dt_trial > tstop) dt_trial = tstop - t;
      ie_richardson_step_with_error(t, y, dt_trial, ynext, &error);
      if (!isfinite(error)) { dt = 0.2 * dt_trial; continue; }
      if (error <= tolerance || dt_trial <= dt_min) {
         t += dt_trial; y[0] = ynext[0]; y[1] = ynext[1]; number_of_intervals++;
         if (t >= tstop) break;
         if (first_accepted) { factor = safety * pow(tolerance/(error+tiny),1.0/q); first_accepted=0; }
         else factor = safety * pow(tolerance/(error+tiny),1.0/q);
         if (!isfinite(factor) || factor <= 0.0) factor=0.5;
         if (factor > 4.0) factor=4.0;
         if (factor < 0.25) factor=0.25;
         error_old=fmax(error,tiny); dt=dt_trial*factor;
      } else {
         factor=safety*pow(tolerance/(error+tiny),1.0/q);
         if (!isfinite(factor) || factor<=0.0) factor=0.1;
         if (factor < 0.1) factor = 0.1;
         dt = dt_trial * factor;
      }
   }
   solution[0]=y[0]; solution[1]=y[1]; *number_of_intervals_ptr=number_of_intervals; return 0;
}

int write_solution(const char *filename, const VEC solution) {
   FILE *file=fopen(filename,"w"); if(file==NULL) return 1;
   fprintf(file,"%.14e %.14e\n",solution[0],solution[1]); fclose(file); return 0;
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
   snprintf(output_filename, sizeof(output_filename), "results/vdp-sequential_IE-%e.txt", tolerance);
   if (dt0 <= 0.0)
   {
      dt0 = (tstop - tstart) / 32.0;
   }

   if (generate_adaptive_solution(tstart, tstop, tolerance, dt0, solution, &number_of_intervals) != 0)
   {
      return 1;
   }

   if (write_solution(output_filename, solution) != 0) return 1;

   printf("Sequential VDP solve finished.\n");
   printf("VDP mu:                %.1f\n", VDP_MU);
   printf("Tolerance:             %.6e\n", tolerance);
   printf("Initial trial step:    %.14e\n", dt0);
   printf("Number of intervals:   %d\n", number_of_intervals);
   printf("Output file:           %s\n", output_filename);
   printf("CPU time:              %.6f seconds\n",
          (double)(clock() - start_time) / CLOCKS_PER_SEC);


   return 0;
}


