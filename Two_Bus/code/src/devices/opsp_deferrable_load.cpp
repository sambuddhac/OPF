#include "opsp_deferrable_load.hpp"

DeferrableLoad::DeferrableLoad(double E, double pmin, double pmax, int Tstart, int Tend, int len) : Device(1, len)
{
  this->def = new CVX_DeferrableLoad;// *)solvers[i];

  def->set_defaults();
  def->setup_indexing();
  def->settings.verbose = 0;
  def->settings.eps = CVXGEN_EPS;
  def->settings.resid_tol = CVXGEN_RESID_TOL;
  def->settings.max_iters = CVXGEN_ITERS;
  def->settings.refine_steps = CVXGEN_REFINE_STEPS;
  def->settings.kkt_reg = CVXGEN_KKT_REG;
  
  def->params.E[0] = E;
  //def->params.pmin[0] = pmin;
  def->params.pmax[0] = pmax;

  for(int j = 0; j < len; ++j) {
   if(Tstart <= j && j <= Tend)
     def->params.c[j][0] = 1;
   else
     def->params.c[j][0] = 0;
  }
}

void DeferrableLoad::solve(const double &rho, const double &rho_old)
{    
  // TODO: really would just like to point over
  terminals[0]->copy_pvec(def->params.v, rho_old/rho);

  def->solve();

  terminals[0]->set_pvec(def->vars.p);
}
