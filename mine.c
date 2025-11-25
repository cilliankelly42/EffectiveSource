#include <stdlib.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <gsl/gsl_errno.h>
#include "effsource.h"
#include "./test/decompose.h"
#include <sys/stat.h>

/* Return the value of the singular field at the point x */
double PhiS_calc(struct coordinate * x)
{
  double PhiS;
  effsource_PhiS(x, &PhiS);
  return PhiS;
}

/* Return the value of the effective source at the point x */
double src_calc(struct coordinate * x)
{
  double PhiS, dPhiS_dx[4], d2PhiS_dx2[10], src;
  effsource_calc(x, &PhiS, dPhiS_dx, d2PhiS_dx2, &src);
  return src;
}

int main()
{
  FILE *parameter_file, *input_file, *output_file;

  // Load basic parameters energy, angular_momentum, array_length
  char buffer[256];
  double e, p, E, L;
  int arrData_length;

  parameter_file = fopen("/home/cillian/.eff_source_J/notebooks/parameters.dat", "r");
  if(parameter_file == NULL){
    printf("Error finding the parameter file");
  }

  fgets(buffer, sizeof(buffer), parameter_file);
  sscanf(buffer, "%lf %lf %lf %lf %d", &e, &p, &E, &L, &arrData_length);
  fclose(parameter_file);

  /* Mass and spin of the central black hole */
  const double a = 0.5;
  const double M = 1.0;
  effsource_init(M, a);

  double t_p, r_p, phi_p, ur; 
  struct coordinate xp;

  /* The point where we measure the singular field/effective source (ordered like {r,theta,phi,t})*/
  struct coordinate x = {20,  M_PI_2 + 0.00001, 0.0, 0.0}; // It seems this definition is unnecessary as the r and theta values are overwritten in the loop below and the phi & t values are not used
  
  int m = 4;
  char output_path[256];
  sprintf(output_path, "./output_data/eccentric_%.2f", e);
  mkdir(output_path, 0777);

  /* Load the orbit data */
  input_file = fopen("/home/cillian/.eff_source_J/notebooks/trajectory_data_t.dat", "r");
  if(input_file == NULL){
    printf("Couldn't open input_file");
  }

  char buf[0x100];
  sprintf(buf, "./output_data/eccentric_%.2f/output_t_m%d.txt", e, m); // Hard code here when changing the orbit
  output_file = fopen(buf, "w");
  if(output_file == NULL){
    printf("No file to write to");
  }

  if(input_file != NULL){
    while(fscanf(input_file, "%lf\t%lf\t%lf\t%lf", &t_p, &r_p, &phi_p, &ur) != EOF){

      xp.t = t_p;
      x.t = xp.t;
      xp.r = r_p;
      xp.theta = M_PI_2;
      xp.phi = phi_p;
       
      effsource_set_particle(&xp, E, L, ur);

      gsl_set_error_handler_off();

      double PhiS[2], dPhiS[8], ddPhiS[20], src[2];
      effsource_calc_m(m, &x, PhiS, dPhiS, ddPhiS, src);

      fprintf(output_file, "%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\n",
        x.t, x.r, x.theta, x.phi, src[0], src[1]);
    }
  }

  fclose(input_file);
  fclose(output_file);
  
  return 0;
}

  