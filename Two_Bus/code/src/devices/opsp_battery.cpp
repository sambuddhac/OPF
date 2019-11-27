#include "opsp_battery.hpp"

Battery::Battery(double q0, double qf, double P, double Q, int len) : Device(1,len)
{
  this->bat = new CVX_Battery;

  bat->set_defaults();
  bat->setup_indexing();
  bat->settings.verbose = 0;
  bat->settings.eps = CVXGEN_EPS;
  bat->settings.resid_tol = CVXGEN_RESID_TOL;
  bat->settings.max_iters = CVXGEN_ITERS;
  bat->settings.refine_steps = CVXGEN_REFINE_STEPS;
  bat->settings.kkt_reg = CVXGEN_KKT_REG;

  bat->params.q0[0] = q0;
  //bat->params.qf[0] = qf;
  bat->params.Q[0] = Q;
  bat->params.P[0] = P;
}

void Battery::solve(const double &rho, const double &rho_old)
{    
  terminals[0]->copy_pvec(bat->params.v, rho_old/rho);

  bat->solve();
    
  terminals[0]->set_pvec(bat->vars.p);
}

