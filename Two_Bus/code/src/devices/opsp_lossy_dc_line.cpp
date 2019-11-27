#include "opsp_lossy_dc_line.hpp"

void LossyDCLine::solve(const double &rho, const double &rho_old)
{    
  terminals[0]->scale_dual_variables(rho_old/rho);
  terminals[1]->scale_dual_variables(rho_old/rho);
  
  const double s = rho/(rho+EPSILON);
  // not vectorizable, but hey, let's give the compiler as much help as possible
  
  double *__restrict__ p1 = terminals[0]->p();
  double *__restrict__ p2 = terminals[1]->p();
  const double *__restrict__ u1 = terminals[0]->u();
  const double *__restrict__ u2 = terminals[1]->u();
  double v1, v2, u, v;
  bool pos;

  for(int i = 0; i < len; ++i) {
    v1 = s*(p1[i] - u1[i]);
    v2 = s*(p2[i] - u2[i]);
    
    u = (v1 - v2)/(2*g);
    v = (v1 + v2)/(2*g);
    
    pos = (u >= 0);

    v = (v >= L) ? L : v;
    u = (u <= 0) ? -u : u;
    
    const double UB = BOUND + M*v/beta_square;
    const double Lu = L*u;
    
    double vv = beta_square*(v-1.0)*(v-1.0);
    double uu = u*u;
        
    if( Lu - u <= UB ) {
      u = M;
      v = L;
    }
    else if( vv + uu - beta_square <= 0 ) {
      // do nothing
    } else {
      // estimate upper bound for lambda
      // upper bound: fix v, project u to line v = L/M |u|
      double lambdaU = (Lu - M*v)/(M - Lu);
      // lower bound: find lambda that makes v = 0
      double lambdaL = MAX(0, -v);

      double lambda = (lambdaU + lambdaL)/2.0;

      const double beta4 = beta_square*beta_square;
      
      // no expensive divides in bisection loop
      while(lambdaU - lambdaL >= BISECTION_TOL)
      {
        lambda = (lambdaU + lambdaL)/2.0;
        //cout << lambdaU << " " << lambdaL << " "<< lambda << endl;

        double s1 = (lambda+beta_square)*(lambda+beta_square);
        double s2 = (1.0 + lambda)*(1.0 + lambda);  
        if( s1*vv + beta4*uu*s2 - beta_square*s1*s2 <= 0 )
          lambdaU = lambda;
        else
          lambdaL = lambda;
      }
      
      u = beta_square*u/(lambda + beta_square);
      v = (v + lambda)/(1.0+lambda);
    }
    
    u = (pos) ? u : -u;
    
    p1[i] = g*(v + u);
    p2[i] = g*(v - u);
  }
}
