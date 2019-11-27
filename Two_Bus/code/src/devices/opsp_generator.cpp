#include "opsp_generator.hpp"

Generator::Generator(double alpha, double beta, double pmax, double pmin, double slew, int len) : Device(1, len)
{
  this->gen = new CVX_Generator;
  
  gen->set_defaults();
  gen->setup_indexing();
  gen->settings.verbose = 0;
  gen->settings.eps = CVXGEN_EPS;
  gen->settings.resid_tol = CVXGEN_RESID_TOL;
  gen->settings.max_iters = CVXGEN_ITERS;
  gen->settings.refine_steps = CVXGEN_REFINE_STEPS;
  gen->settings.kkt_reg = CVXGEN_KKT_REG;
  
  gen->params.alpha[0] = alpha;  // from opf
  gen->params.beta[0] = beta;
  gen->params.pmax[0] = pmax;
  gen->params.pmin[0] = pmin;
  gen->params.S[0] = slew;  
}

void Generator::solve(const double &rho, const double &rho_old)
{    
  gen->params.rho[0] = rho;

  // TODO: really would just like to point over
  terminals[0]->copy_pvec(gen->params.v, rho_old/rho);

  gen->solve();
  if(!gen->work.converged) {
    printf("gen: did not converge!\n");
  }
  
  terminals[0]->set_pvec(gen->vars.p);
}

const double Generator::objective() const
{
  double obj = 0;
  const double *__restrict__ p = terminals[0]->p();
  for(int j = 0; j < len; ++j) {
    //std::cout << val << std::endl;
    obj += gen->params.alpha[0]*p[j]*p[j] - gen->params.beta[0]*p[j];
  }
  return obj;
}