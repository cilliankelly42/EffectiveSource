#include <stdlib.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <gsl/gsl_errno.h>
#include "effsource.h"
#include "./test/decompose.h"

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

  // Load basic parameters rmin,rmax,energy, angular_momentum, array_length, time_period, half_time_period
  char buffer[256];
  double rmin, rmax, E, L, Trad, halfTrad, rfield;
  int arrData_length;

  parameter_file = fopen("/home/cillian/.eff_source_J/notebooks/parameters.dat", "r");
  if(parameter_file == NULL){
    printf("Error finding the parameter file");
  }

  fgets(buffer, sizeof(buffer), parameter_file);
  sscanf(buffer, "%lf %lf %lf %lf %d %lf %lf %lf", &rmin, &rmax, &E, &L, &arrData_length, &Trad, &halfTrad, &rfield);

  fclose(parameter_file);

  /*
  double t_p[arrData_length];
  double r_p[arrData_length];
  double phi_p[arrData_length];

  char line[0x100];
  int i = 0;

  
  while (i < arrData_length && fgets(line, sizeof(line), input_file)) {
      double a, b, c;

      if (sscanf(line, "%lf %lf %lf", &a, &b, &c) == 3) {
          t_p[i] = a;
          r_p[i] = b;
          phi_p[i] = c;
          i++;
      }
      else {
          printf("error with line %d\n", i+1);
      }
  }
 

  fclose(input_file);
  */

  /* Mass and spin of the central black hole */
  const double a = 0.5;
  const double M = 1.0;
  effsource_init(M, a);

  double t_p, r_p, phi_p, ur; 
  struct coordinate xp;

  /* The point where we measure the singular field/effective source (ordered like {r,theta,phi,t})*/
  struct coordinate x = {10,  M_PI_2 + 0.00001, 0.0, 0.0}; // It seems this definition is unnecessary as the r and theta values are overwritten in the loop below and the phi & t values are not used
  
  int m = 1;
 
  /* Load the orbit data */
  input_file = fopen("/home/cillian/.eff_source_J/notebooks/trajectory_data_t.dat", "r");
  if(input_file == NULL){
    printf("Couldn't open input_file");
  }

  char buf[0x100];
  sprintf(buf, "./output_data/eccentric_08/output_t_m%d.txt", m); // Hard code here when changing the orbit
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

  /*
  for(int i = 0; i<arrData_length; i++){

    // Set the coordinate values of the particle
    xp.t = t_p[i]; //Not sure if this is used but could evolve it straightforwardly like this
    xp.r = r_p[i];
    xp.theta = M_PI_2;
    xp.phi = phi_p[i]; 
    x.t = xp.t;

    // The following loop determines if ur is positive or negative ASSUMING THE ORBIT BEGINS AT rmin 
    //for(int j = 0; j <= 18; j+=2)
    //{
     // if(xp.r == rmin || xp.r == rmax){
       // ur = 0;
      //}
      //else if((j * halfTrad< xp.t) && (xp.t < (j+1) * halfTrad)){
       // ur = sqrt(-1 + e*e - 2*(-(l-a*e)*(l-a*e)*M/(r_p[i]*r_p[i]*r_p[i]) + (l*l-a*a*(e*e-1))/(2.*r_p[i]*r_p[i]) - M/r_p[i]));
      //}
      //else if(((j+1) * halfTrad< xp.t) && (xp.t < (j+2) * halfTrad)){
       // ur  = -sqrt(-1 + e*e - 2*(-(l-a*e)*(l-a*e)*M/(r_p[i]*r_p[i]*r_p[i]) + (l*l-a*a*(e*e-1))/(2.*r_p[i]*r_p[i]) - M/r_p[i]));
      //}
    } 

    // initialising the particle
    effsource_set_particle(&xp, e, l, ur);

    // Disable the GSL error handler so that it doesn't abort due to roundoff errors //
    gsl_set_error_handler_off();

    char filename[100];
    sprintf(filename, "./data/datam%d_ecc03/loop_%d.txt", m, i);
    output_file = fopen(filename, "w");

    double PhiS[2], dPhiS[8], ddPhiS[20], src[2]; //src_num[2],
    effsource_calc_m(m, &x, PhiS, dPhiS, ddPhiS, src);
    //m_decompose(m, x, src_calc, src_num);

    fprintf(output_file, "%.15g\t%.15g\t%.15g\t%.15g\t%.15g\n",
    x.r, x.theta, xp.t, 
    src[0], src[1]);

    fclose(output_file);
    }
        for(double r=8.9; r<=11.1; r+=0.1)
    {
        x.r = r;
        for(double theta=M_PI_2-0.1; theta<=M_PI_2+0.1; theta+=0.011)
        {
        double PhiS[2], src_num[2], dPhiS[8], ddPhiS[20], src[2];
        x.theta     = theta;
        effsource_calc_m(m, &x, PhiS, dPhiS, ddPhiS, src);
        m_decompose(m, x, src_calc, src_num);

        fprintf(output_file, "%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\n",
            x.r, x.theta, xp.t,
            src[0], src[1], src_num[0], src_num[1]);
        }
    }
      
    fclose(output_file);
  }
*/
//}
