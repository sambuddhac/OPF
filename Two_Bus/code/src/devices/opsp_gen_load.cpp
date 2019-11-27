#include "opsp_gen_load.hpp"

GeneratingLoad::GeneratingLoad(double *__restrict__ l, double gamma, int len) : 
  Device(1, len), gamma(gamma)
{ 
	pdes = new double[len];
  for(int i = 0; i < len; i++) {
    pdes[i] = l[i];
  }
}

void GeneratingLoad::solve(const double &rho, const double &rho_old)
{
  // minimize
  //   sum[t=0..T](gamma*(p[t] - pdes) + rho/2*square(p[t] - v[t]))
  // subject to
  //   p[t] <= 0, t=0..T
	//	 p[t] >= pdes, t=0..T (pdes is < 0)
  // end

	// rho/2*square(p - val) + gamma*(p - pdes)
	// rho*(p - val) + gamma = 0
	// p = val - gamma/rho 

  
  terminals[0]->scale_dual_variables(rho_old/rho);
  double *__restrict__ p = terminals[0]->p();
  const double *__restrict__ u = terminals[0]->u();
  double *__restrict__ theta = terminals[0]->theta();
  const double *__restrict__ v = terminals[0]->v();
  double val;
  
  for(int j = 0; j < len; ++j) {
    val = p[j] - u[j];
    
		val = val - gamma/rho;
    val = (val >= 0) ? 0 : val;
		val = (val <= pdes[j]) ? pdes[j] : val;
    
    //printf("val: %f vars: %f\n", val, cur->vars.p[j][0]);
    p[j] = val;
    //theta[j] -= v[j];
  }
}

const double GeneratingLoad::objective() const
{
  double obj = 0;
  const double *__restrict__ p = terminals[0]->p();
  for(int j = 0; j < len; ++j) {
    //std::cout << val << std::endl;
    obj += gamma*MAX(0, p[j] - pdes[j]);
  }
  return obj;
}