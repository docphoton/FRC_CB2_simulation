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
Ball Launch
    Solving the equations of motion for the ball rolling-while-launching 
    in the CB2 Catapult. 
 */ 

// setup globals...they who say globals are bad programming practice don't understand economy in physics simulations!!! End of rant...

 double A, B, C, gamm, g, k, m, r, Iball,Iarm, l0, k, marmg_gam, g, det, theta0, pi; 

// now onto the differential equations, their jacobian and the main program...

int
func (double t, const double y[], double f[],
      void *params)
{
  double d = *(double *)params; // is my 'd' parameter. 
  double v1, v2; 
  v1 = -m*g*r*sin(y[2]+gamm); 
  v2 = -(marmg_gam+m*g*l0)*cos(y[2])-k*(y[2]-theta0)-m*g*r*y[0]*cos(y[2]+gamm)+m*g*r*sin(y[2]+gamm); 
  f[0] = y[1];
  f[1] = (B*v1-C*v2)/det;
  f[2] = y[3];
  f[3] = (-C*v1+A*v2)/det; 
  return GSL_SUCCESS;
}

int
jac (double t, const double y[], double *dfdy, 
     double dfdt[], void *params)
{
  double d = *(double *)params;
  double w1, w2, z1, z2; 
  w1 = -m*g*r*cos(y[2]+gamm); 
  w2 = -k+(m*g*l0+marmg_gam)*sin(y[2])+m*g*r*y[0]*sin(y[2]+gamm)+m*g*r*cos(y[2]+gamm); 
  z1 = C*cos(y[2]+gamm)*m*g*r;
  z2 =  -m*g*r*A*cos(y[2]+gamm);
  gsl_matrix_view dfdy_mat 
    = gsl_matrix_view_array (dfdy, 4, 4);
  gsl_matrix * m = &dfdy_mat.matrix; 
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 0, 2, 0.0);
  gsl_matrix_set (m, 0, 3, 0.0);
  gsl_matrix_set (m, 1, 0, z1);
  gsl_matrix_set (m, 1, 1, 0.0);
  gsl_matrix_set (m, 1, 2, (B*w1-C*w2)/det);
  gsl_matrix_set (m, 1, 3, 0.0);
  gsl_matrix_set (m, 2, 0, 0.0);
  gsl_matrix_set (m, 2, 3, 1.0);
  gsl_matrix_set (m, 2, 2, 0.0);
  gsl_matrix_set (m, 2, 1, 0.0);
  gsl_matrix_set (m, 3, 0, z2);
  gsl_matrix_set (m, 3, 3, 0.0);
  gsl_matrix_set (m, 3, 2, (-C*w1+A*w2)/det);
  gsl_matrix_set (m, 3, 1, 0.0);
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
  double d; 
  gsl_odeiv_system sys = {func, jac, 4, &d};
  int ii=1, num = 100; // minimum number of points to have 
  double t = 0.0, t1 =.140, ti;
  double h = 1e-6, velx, vely, vtot, vangle;
 
// initialization of all the shooter physical parameters. 
 pi = 4.0*atan(1.0); 
 m = 1.4;              //mass of the ball, Kg
 r = .3;               //radius of the ball, meters 
 Iball = 2.0*m/5.0*r*r;    // I of ball
 l0 = .45;            // distance between pivot and ball rest, m
 gamm = .5;         // angle of the wrist. radians
 Iarm = (1.75)*l0*l0/3.0; // moment of inertia of the catapult ar,
   k = 70;           // torsional spring constant of the catapult spring. j/rad/rad
 marmg_gam = 13.7;     // m_catapult_arm * g * cm distance to pivot. For the potential energy in the arm lift; 
 theta0 = 90*2*pi/360; // the zero of the force is at theta0. For this simulation, it is combined wiht k to get the starting torque and how fast that decreases with angle. 
 g = 9.8;            // gravity , m/s^2
 A =m*r*r+Iball;  
 B = Iball+Iarm+m*(l0*l0+r*r)-2*l0*r*sin(gamm);  // ignoring the non-linear terms here in the coef...not too systematic!!!
 C =  -Iball + m*r*(l0*sin(gamm)-r); 
 det = A*B-C*C; 
 // initial conditions in y = phi, phidot, theta, thetadot
  double y[4] = { 0.0, 0.0, -0.40, 0.0};// initial consditions (phi, phidot, theta, thetadot) 
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
      // have to use the release co-ordinates to determine the final velocity vector of the ball. 
      velx = -l0*sin(y[2])*y[3]+r*y[1]*cos(y[2]+gamm)-r*y[0]*y[3]*sin(y[2]+gamm)-r*y[3]*cos(y[2]+gamm); 
      vely = l0*y[3]*cos(y[2])+r*y[0]*y[3]*cos(y[2]+gamm)+r*y[1]*sin(y[2]+gamm)-r*y[3]*sin(y[2]+gamm); 
	//      printf("%.5e %.5e %.5e  %.5e  %.5e\n", t, y[0], y[1], y[2], y[3]);
      vtot = pow(velx*velx+vely*vely,.5); 
      vangle = -180*atan(vely/velx)/pi; 
      printf("%.5e %.5e %.5e  %.5e  %.5e  %.5e\n", t, y[0], 180*y[2]/pi, vtot, vangle, y[1]-y[3]);
      }
    }
  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);
  return 0;
}
