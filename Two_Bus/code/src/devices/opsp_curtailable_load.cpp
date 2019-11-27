#include "opsp_curtailable_load.hpp"

void CurtailableLoad::solve(const double &rho, const double &rho_old)
{
  // minimize
  //   sum[t=0..T](gamma*pos(pdes - p[t]) + rho/2*square(p[t] - v[t]))
  // subject to
  //   p[t] >= 0, t=0..T
  // end

	// rho/2*square(p - val) + gamma*pos(pdes - p)
	// rho*(p - val) - gamma = 0
	// p =val + gamma/rho if (pdes > val)
	// orw, p = val if (pdes < val)
  
  terminals[0]->scale_dual_variables(rho_old/rho);
  double *__restrict__ p = terminals[0]->p();
  const double *__restrict__ u = terminals[0]->u();
  
  double *__restrict__ theta = terminals[0]->theta();
  const double *__restrict__ v = terminals[0]->v();
  
  double val;
  
  for(int j = 0; j < len; ++j) {
    val = p[j] - u[j];
    
    val = MAX(val,pdes) - MAX(pdes - gamma/rho - val, 0);
    val = (val <= 0) ? 0 : val;
    //printf("%f ", val);    
    p[j] = val;
    //theta[j] -= v[j];
  }
  //printf("\n");
}

const double CurtailableLoad::objective() const
{
  double obj = 0;
  const double *__restrict__ p = terminals[0]->p();
  for(int j = 0; j < len; ++j) {
    //std::cout << val << std::endl;
    obj += gamma*MAX(0, pdes - p[j]);
  }
  return obj;
}
