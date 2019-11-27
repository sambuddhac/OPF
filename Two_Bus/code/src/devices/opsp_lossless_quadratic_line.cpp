#include "opsp_lossless_quadratic_line.hpp"

void LosslessQuadraticLine::solve(const double &rho, const double &rho_old)
{ 
  terminals[0]->scale_dual_variables(rho_old/rho);
  terminals[1]->scale_dual_variables(rho_old/rho);
  
  double *__restrict__ p1 = terminals[0]->p();
  double *__restrict__ p2 = terminals[1]->p();
  const double *__restrict__ u1 = terminals[0]->u();
  const double *__restrict__ u2 = terminals[1]->u();
  double val1, val2, val3, val4, y;
  
  double *__restrict__ theta1 = terminals[0]->theta();
  const double *__restrict__ v1 = terminals[0]->v();
  double *__restrict__ theta2 = terminals[1]->theta();
  const double *__restrict__ v2 = terminals[1]->v();
  
  for(int i = 0; i < len; ++i) {
    //terminals[0]->set_u(terminals[0]->u(i)*rho_old/rho, i);
    //terminals[1]->set_u(terminals[1]->u(i)*rho_old/rho, i);
    
    val1 = (p1[i] - u1[i]);
    val2 = (p2[i] - u2[i]);
    val3 = (theta1[i] - v1[i]);
    val4 = (theta2[i] - v2[i]);
    
    y = (val1 - val2)/(2.0);
    y = (2.0*y + (val4 - val3))/3.0;

    p1[i] = y;//0.5*(val1 - val2);
    p2[i] = -y;//0.5*(val2 - val1);
    


    theta1[i] = (val3 + val4 - y)/2.0;//v1[i];
    theta2[i] = (val3 + val4 + y)/2.0;//v2[i];
    // theta1[i] -= v1[i];
    // theta2[i] -= v2[i];
  }
}