#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
/* 
Ball in air projectile motion with drag
   solutions depends on a single 'internal' parameter, 
   the weight-specific damping coef d = D/(mg) 
       where D = 1/2 C_d \rho A. (A is the 
       x-sectional area and rho is the density of air. 
 */ 
int
func (double t, const double y[], double f[],
      void *params)
{
  double d = *(double *)params; // is my 'd' parameter. 
  double v = pow(y[1]*y[1]+y[3]*y[3],0.5); 
  f[0] = y[1];
  f[1] = -d*v*y[1];
  f[2] = y[3]; 
  f[3] = -1.0-d*v*y[3]; 
  return GSL_SUCCESS;
}

int
jac (double t, const double y[], double *dfdy, 
     double dfdt[], void *params)
{
  double d = *(double *)params;
  double v = pow(y[1]*y[1]+y[3]*y[3],.5); 
  gsl_matrix_view dfdy_mat 
    = gsl_matrix_view_array (dfdy, 4, 4);
  gsl_matrix * m = &dfdy_mat.matrix; 
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 0, 2, 0.0);
  gsl_matrix_set (m, 0, 3, 0.0);
  gsl_matrix_set (m, 1, 0, 0.0);
  gsl_matrix_set (m, 1, 1, -d*(y[1]*y[1]/v+v));
  gsl_matrix_set (m, 1, 2, 0.0);
  gsl_matrix_set (m, 1, 3, -d*y[1]*y[3]/v);
  gsl_matrix_set (m, 2, 0, 0.0);
  gsl_matrix_set (m, 2, 3, 1.0);
  gsl_matrix_set (m, 2, 2, 0.0);
  gsl_matrix_set (m, 2, 1, 0.0);
  gsl_matrix_set (m, 3, 0, 0.0);
  gsl_matrix_set (m, 3, 3, -d*(y[3]*y[3]/v+v));
  gsl_matrix_set (m, 3, 2, 0.0);
  gsl_matrix_set (m, 3, 1, -d*y[1]*y[3]/v);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  dfdt[2] = 0.0; 
  dfdt[3] = 0.0; 
  return GSL_SUCCESS;
}

int
main (void)
{
  const gsl_odeiv_step_type * T 
    = gsl_odeiv_step_rk8pd;

  gsl_odeiv_step * s 
    = gsl_odeiv_step_alloc (T, 4);
  gsl_odeiv_control * c 
    = gsl_odeiv_control_y_new (1e-6, 0.0);
  gsl_odeiv_evolve * e 
    = gsl_odeiv_evolve_alloc (4);

  double d = .0085;// assuming a 1Kg ball mass, 2' diameter, etc. 
  //double d = 0.0; // for testing, this is no air drag 
  double g = 9.8; 
  gsl_odeiv_system sys = {func, jac, 4, &d};
  int ii=1, num = 100; // minimum number of points to have 
  double t = 0.0, t1 =1.0*g, ti;
  double h = 1e-6;
  double y[4] = { 0.0, 10.3, 0.0, 6.7};// initial consditions (x, v_0x, y_0, v_y0) 
  while (ii < num) 
    { 
      ti = ii*t1/num; 
      ii++;
  while (t < ti)
      {
      int status = gsl_odeiv_evolve_apply (e, c, s,
                                           &sys, 
                                           &t, ti,
                                           &h, y);

      if (status != GSL_SUCCESS)
          break;

      printf("%.5e %.5e %.5e  %.5e  %.5e\n", t, y[0]/g, y[1], y[2]/g, y[3]);
      }
    }
  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);
  return 0;
}
