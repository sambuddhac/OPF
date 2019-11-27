#include "opsp_fixed_load.hpp"

FixedLoad::FixedLoad(double *__restrict__ l, int len) : Device(1, len)
{
  loss_profile = new double[len];
  for(int i = 0; i < len; i++) {
    loss_profile[i] = l[i];
  }  
}


void FixedLoad::solve(const double &rho, const double &rho_old)
{
  // minimize
  //   sum[t=0..T](1/2*square(p[t] - l[t]) + rho/2*square(p[t] - v[t]))
  // subject to
  //   p[t] >= l[t], t=0..T
  // end
  terminals[0]->scale_dual_variables(rho_old/rho);
  double *__restrict__ p = terminals[0]->p();
  double *__restrict__ u = terminals[0]->u();
  
  double *__restrict__ theta = terminals[0]->theta();
  const double *__restrict__ v = terminals[0]->v();
  double val;
  
  const double s = rho/(10.0 + rho);
  
  for(int j = 0; j < len; ++j) {
    //double val = terminals[0]->p(j) - terminals[0]->u(j);
    //val = (val >= loss_profile[j]) ? val : loss_profile[j];
    
    // fixed loads are now allowed to "oversaturate" and "shortfall"
    // minimze || p - v||_2^2
    // //s.t. p >= loss_profile
    
    // 0.5 || p - l || + 0.5 rho ||p - v||
    // (p - l) + rho(p-v) = 0
    // (1+rho)p -l - rhov = 0
    // p = (l + rho*v)/(1+rho)
    
    val = (p[j] - u[j])*s+ loss_profile[j]*(1.0-s);
    val = (val < 0) ? 0 : val;
        
		//val -= (0.1/rho);
    p[j] = loss_profile[j];//val;//((val <= loss_profile[j]) ? loss_profile[j] : val);
    
    
    //theta[j] -= v[j];
  }
}
