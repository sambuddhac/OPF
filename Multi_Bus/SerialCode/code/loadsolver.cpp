//solver.cpp
/* Produced by CVXGEN, 2014-08-13 17:31:43 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: solver.c. */
/* Description: Main solver file. */
#include "loadsolver.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>
double Loadsolver::eval_gap(void) {
  int i;
  double gap;
  gap = 0;
  for (i = 0; i < 0; i++)
    gap += work.z[i]*work.s[i];
  return gap;
}
void Loadsolver::set_defaults(void) {
  settings.resid_tol = 1e-6;
  settings.eps = 1e-4;
  settings.max_iters = 25;
  settings.refine_steps = 1;
  settings.s_init = 1;
  settings.z_init = 1;
  settings.debug = 0;
  settings.verbose = 1;
  settings.verbose_refinement = 0;
  settings.better_start = 1;
  settings.kkt_reg = 1e-7;
}
void Loadsolver::setup_pointers(void) {
  work.y = work.x + 1;
  work.s = work.x + 1;
  work.z = work.x + 1;
  vars.Thetal = work.x + 0;
}
void Loadsolver::setup_indexing(void) {
  setup_pointers();
}
void Loadsolver::set_start(void) {
  int i;
  for (i = 0; i < 1; i++)
    work.x[i] = 0;
  for (i = 0; i < 0; i++)
    work.y[i] = 0;
  for (i = 0; i < 0; i++)
    work.s[i] = (work.h[i] > 0) ? work.h[i] : settings.s_init;
  for (i = 0; i < 0; i++)
    work.z[i] = settings.z_init;
}
double Loadsolver::eval_objv(void) {
  int i;
  double objv;
  /* Borrow space in work.rhs. */
  multbyP(work.rhs, work.x);
  objv = 0;
  for (i = 0; i < 1; i++)
    objv += work.x[i]*work.rhs[i];
  objv *= 0.5;
  for (i = 0; i < 1; i++)
    objv += work.q[i]*work.x[i];
  objv += work.quad_953246392320[0]-params.L[0];
  return objv;
}
void Loadsolver::fillrhs_aff(void) {
  int i;
  double *r1, *r2, *r3, *r4;
  r1 = work.rhs;
  r2 = work.rhs + 1;
  r3 = work.rhs + 1;
  r4 = work.rhs + 1;
  /* r1 = -A^Ty - G^Tz - Px - q. */
  multbymAT(r1, work.y);
  multbymGT(work.buffer, work.z);
  for (i = 0; i < 1; i++)
    r1[i] += work.buffer[i];
  multbyP(work.buffer, work.x);
  for (i = 0; i < 1; i++)
    r1[i] -= work.buffer[i] + work.q[i];
  /* r2 = -z. */
  for (i = 0; i < 0; i++)
    r2[i] = -work.z[i];
  /* r3 = -Gx - s + h. */
  multbymG(r3, work.x);
  for (i = 0; i < 0; i++)
    r3[i] += -work.s[i] + work.h[i];
  /* r4 = -Ax + b. */
  multbymA(r4, work.x);
  for (i = 0; i < 0; i++)
    r4[i] += work.b[i];
}
void Loadsolver::fillrhs_cc(void) {
  int i;
  double *r2;
  double *ds_aff, *dz_aff;
  double mu;
  double alpha;
  double sigma;
  double smu;
  double minval;
  r2 = work.rhs + 1;
  ds_aff = work.lhs_aff + 1;
  dz_aff = work.lhs_aff + 1;
  mu = 0;
  for (i = 0; i < 0; i++)
    mu += work.s[i]*work.z[i];
  /* Don't finish calculating mu quite yet. */
  /* Find min(min(ds./s), min(dz./z)). */
  minval = 0;
  for (i = 0; i < 0; i++)
    if (ds_aff[i] < minval*work.s[i])
      minval = ds_aff[i]/work.s[i];
  for (i = 0; i < 0; i++)
    if (dz_aff[i] < minval*work.z[i])
      minval = dz_aff[i]/work.z[i];
  /* Find alpha. */
  if (-1 < minval)
      alpha = 1;
  else
      alpha = -1/minval;
  sigma = 0;
  for (i = 0; i < 0; i++)
    sigma += (work.s[i] + alpha*ds_aff[i])*
      (work.z[i] + alpha*dz_aff[i]);
  sigma /= mu;
  sigma = sigma*sigma*sigma;
  /* Finish calculating mu now. */
  smu = sigma*mu;
  /* Fill-in the rhs. */
  for (i = 0; i < 1; i++)
    work.rhs[i] = 0;
  for (i = 1; i < 1; i++)
    work.rhs[i] = 0;
  for (i = 0; i < 0; i++)
    r2[i] = work.s_inv[i]*(smu - ds_aff[i]*dz_aff[i]);
}
void Loadsolver::refine(double *target, double *var) {
  int i, j;
  double *residual = work.buffer;
  double norm2;
  double *new_var = work.buffer2;
  for (j = 0; j < settings.refine_steps; j++) {
    norm2 = 0;
    matrix_multiply(residual, var);
    for (i = 0; i < 1; i++) {
      residual[i] = residual[i] - target[i];
      norm2 += residual[i]*residual[i];
    }
#ifndef ZERO_LIBRARY_MODE
    if (settings.verbose_refinement) {
      if (j == 0)
        printf("Initial residual before refinement has norm squared %.6g.\n", norm2);
      else
        printf("After refinement we get squared norm %.6g.\n", norm2);
    }
#endif
    /* Solve to find new_var = KKT \ (target - A*var). */
    ldl_solve(residual, new_var);
    /* Update var += new_var, or var += KKT \ (target - A*var). */
    for (i = 0; i < 1; i++) {
      var[i] -= new_var[i];
    }
  }
#ifndef ZERO_LIBRARY_MODE
  if (settings.verbose_refinement) {
    /* Check the residual once more, but only if we're reporting it, since */
    /* it's expensive. */
    norm2 = 0;
    matrix_multiply(residual, var);
    for (i = 0; i < 1; i++) {
      residual[i] = residual[i] - target[i];
      norm2 += residual[i]*residual[i];
    }
    if (j == 0)
      printf("Initial residual before refinement has norm squared %.6g.\n", norm2);
    else
      printf("After refinement we get squared norm %.6g.\n", norm2);
  }
#endif
}
double Loadsolver::calc_ineq_resid_squared(void) {
  /* Calculates the norm ||-Gx - s + h||. */
  double norm2_squared;
  int i;
  /* Find -Gx. */
  multbymG(work.buffer, work.x);
  /* Add -s + h. */
  for (i = 0; i < 0; i++)
    work.buffer[i] += -work.s[i] + work.h[i];
  /* Now find the squared norm. */
  norm2_squared = 0;
  for (i = 0; i < 0; i++)
    norm2_squared += work.buffer[i]*work.buffer[i];
  return norm2_squared;
}
double Loadsolver::calc_eq_resid_squared(void) {
  /* Calculates the norm ||-Ax + b||. */
  double norm2_squared;
  int i;
  /* Find -Ax. */
  multbymA(work.buffer, work.x);
  /* Add +b. */
  for (i = 0; i < 0; i++)
    work.buffer[i] += work.b[i];
  /* Now find the squared norm. */
  norm2_squared = 0;
  for (i = 0; i < 0; i++)
    norm2_squared += work.buffer[i]*work.buffer[i];
  return norm2_squared;
}
void Loadsolver::better_start(void) {
  /* Calculates a better starting point, using a similar approach to CVXOPT. */
  /* Not yet speed optimized. */
  int i;
  double *x, *s, *z, *y;
  double alpha;
  work.block_33[0] = -1;
  /* Make sure sinvz is 1 to make hijacked KKT system ok. */
  for (i = 0; i < 0; i++)
    work.s_inv_z[i] = 1;
  fill_KKT();
  ldl_factor();
  fillrhs_start();
  /* Borrow work.lhs_aff for the solution. */
  ldl_solve(work.rhs, work.lhs_aff);
  /* Don't do any refinement for now. Precision doesn't matter too much. */
  x = work.lhs_aff;
  s = work.lhs_aff + 1;
  z = work.lhs_aff + 1;
  y = work.lhs_aff + 1;
  /* Just set x and y as is. */
  for (i = 0; i < 1; i++)
    work.x[i] = x[i];
  for (i = 0; i < 0; i++)
    work.y[i] = y[i];
  /* Now complete the initialization. Start with s. */
  /* Must have alpha > max(z). */
  alpha = -1e99;
  for (i = 0; i < 0; i++)
    if (alpha < z[i])
      alpha = z[i];
  if (alpha < 0) {
    for (i = 0; i < 0; i++)
      work.s[i] = -z[i];
  } else {
    alpha += 1;
    for (i = 0; i < 0; i++)
      work.s[i] = -z[i] + alpha;
  }
  /* Now initialize z. */
  /* Now must have alpha > max(-z). */
  alpha = -1e99;
  for (i = 0; i < 0; i++)
    if (alpha < -z[i])
      alpha = -z[i];
  if (alpha < 0) {
    for (i = 0; i < 0; i++)
      work.z[i] = z[i];
  } else {
    alpha += 1;
    for (i = 0; i < 0; i++)
      work.z[i] = z[i] + alpha;
  }
}
void Loadsolver::fillrhs_start(void) {
  /* Fill rhs with (-q, 0, h, b). */
  int i;
  double *r1, *r2, *r3, *r4;
  r1 = work.rhs;
  r2 = work.rhs + 1;
  r3 = work.rhs + 1;
  r4 = work.rhs + 1;
  for (i = 0; i < 1; i++)
    r1[i] = -work.q[i];
  for (i = 0; i < 0; i++)
    r2[i] = 0;
  for (i = 0; i < 0; i++)
    r3[i] = work.h[i];
  for (i = 0; i < 0; i++)
    r4[i] = work.b[i];
}
long Loadsolver::solve(void) {
  int i;
  int iter;
  double *dx, *ds, *dy, *dz;
  double minval;
  double alpha;
  work.converged = 0;
  setup_pointers();
  pre_ops();
#ifndef ZERO_LIBRARY_MODE
  if (settings.verbose)
    printf("iter     objv        gap       |Ax-b|    |Gx+s-h|    step\n");
#endif
  fillq();
  fillh();
  fillb();
  if (settings.better_start)
    better_start();
  else
    set_start();
  for (iter = 0; iter < settings.max_iters; iter++) {
    for (i = 0; i < 0; i++) {
      work.s_inv[i] = 1.0 / work.s[i];
      work.s_inv_z[i] = work.s_inv[i]*work.z[i];
    }
    work.block_33[0] = 0;
    fill_KKT();
    ldl_factor();
    /* Affine scaling directions. */
    fillrhs_aff();
    ldl_solve(work.rhs, work.lhs_aff);
    refine(work.rhs, work.lhs_aff);
    /* Centering plus corrector directions. */
    fillrhs_cc();
    ldl_solve(work.rhs, work.lhs_cc);
    refine(work.rhs, work.lhs_cc);
    /* Add the two together and store in aff. */
    for (i = 0; i < 1; i++)
      work.lhs_aff[i] += work.lhs_cc[i];
    /* Rename aff to reflect its new meaning. */
    dx = work.lhs_aff;
    ds = work.lhs_aff + 1;
    dz = work.lhs_aff + 1;
    dy = work.lhs_aff + 1;
    /* Find min(min(ds./s), min(dz./z)). */
    minval = 0;
    for (i = 0; i < 0; i++)
      if (ds[i] < minval*work.s[i])
        minval = ds[i]/work.s[i];
    for (i = 0; i < 0; i++)
      if (dz[i] < minval*work.z[i])
        minval = dz[i]/work.z[i];
    /* Find alpha. */
    if (-0.99 < minval)
      alpha = 1;
    else
      alpha = -0.99/minval;
    /* Update the primal and dual variables. */
    for (i = 0; i < 1; i++)
      work.x[i] += alpha*dx[i];
    for (i = 0; i < 0; i++)
      work.s[i] += alpha*ds[i];
    for (i = 0; i < 0; i++)
      work.z[i] += alpha*dz[i];
    for (i = 0; i < 0; i++)
      work.y[i] += alpha*dy[i];
    work.gap = eval_gap();
    work.eq_resid_squared = calc_eq_resid_squared();
    work.ineq_resid_squared = calc_ineq_resid_squared();
#ifndef ZERO_LIBRARY_MODE
    if (settings.verbose) {
      work.optval = eval_objv();
      printf("%3d   %10.3e  %9.2e  %9.2e  %9.2e  % 6.4f\n",
          iter+1, work.optval, work.gap, sqrt(work.eq_resid_squared),
          sqrt(work.ineq_resid_squared), alpha);
    }
#endif
    /* Test termination conditions. Requires optimality, and satisfied */
    /* constraints. */
    if (   (work.gap < settings.eps)
        && (work.eq_resid_squared <= settings.resid_tol*settings.resid_tol)
        && (work.ineq_resid_squared <= settings.resid_tol*settings.resid_tol)
       ) {
      work.converged = 1;
      work.optval = eval_objv();
      return iter+1;
    }
  }
  return iter;
}

//util.cpp
/* Produced by CVXGEN, 2014-08-13 17:31:43 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: util.c. */
/* Description: Common utility file for all cvxgen code. */

void Loadsolver::tic(void) {
  tic_timestart = clock();
}
float Loadsolver::toc(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  printf("time: %8.2f.\n", (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC);
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}
float Loadsolver::tocq(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}
void Loadsolver::printmatrix(char *name, double *A, int m, int n, int sparse) {
  int i, j;
  printf("%s = [...\n", name);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
      if ((sparse == 1) && (A[i+j*m] == 0))
        printf("         0");
      else
        printf("  % 9.4f", A[i+j*m]);
    printf(",\n");
  }
  printf("];\n");
}
double Loadsolver::unif(double lower, double upper) {
  return lower + ((upper - lower)*rand())/RAND_MAX;
}
/* Next function is from numerical recipes in C. */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
float Loadsolver::ran1(long*idum, int reset) {
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (reset) {
    iy = 0;
  }
  if (*idum<=0||!iy) {
    if (-(*idum)<1)*idum=1;
    else *idum=-(*idum);
    for (j=NTAB+7; j>=0; j--) {
      k = (*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum<0)*idum+=IM;
      if (j<NTAB)iv[j]=*idum;
    }
    iy = iv[0];
  }
  k = (*idum)/IQ;
  *idum = IA*(*idum-k*IQ)-IR*k;
  if (*idum<0)*idum += IM;
  j = iy/NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy)> RNMX) return RNMX;
  else return temp;
}
/* Next function is from numerical recipes in C. */
float Loadsolver::randn_internal(long *idum, int reset) {
  static int iset=0;
  static float gset;
  float fac, rsq, v1, v2;
  if (reset) {
    iset = 0;
  }
  if (iset==0) {
    do {
      v1 = 2.0*ran1(idum, reset)-1.0;
      v2 = 2.0*ran1(idum, reset)-1.0;
      rsq = v1*v1+v2*v2;
    } while(rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }
}
double Loadsolver::randn(void) {
  return randn_internal(&global_seed, 0);
}
void Loadsolver::reset_rand(void) {
  srand(15);
  global_seed = 1;
  randn_internal(&global_seed, 1);
}

//matrix_support.cpp
/* Produced by CVXGEN, 2014-08-13 17:31:43 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: matrix_support.c. */
/* Description: Support functions for matrix multiplication and vector filling. */

void Loadsolver::multbymA(double *lhs, double *rhs) {
}
void Loadsolver::multbymAT(double *lhs, double *rhs) {
  lhs[0] = 0;
}
void Loadsolver::multbymG(double *lhs, double *rhs) {
}
void Loadsolver::multbymGT(double *lhs, double *rhs) {
  lhs[0] = 0;
}
void Loadsolver::multbyP(double *lhs, double *rhs) {
  /* TODO use the fact that P is symmetric? */
  /* TODO check doubling / half factor etc. */
  lhs[0] = rhs[0]*(2*work.frac_121674190848);
}
void Loadsolver::fillq(void) {
  work.q[0] = 2*(-params.Vl_N_avg[0]-params.Thetal_N_avg[0]+params.vl_N[0])*work.frac_121674190848;
}
void Loadsolver::fillh(void) {
}
void Loadsolver::fillb(void) {
}
void Loadsolver::pre_ops(void) {
  work.frac_121674190848 = params.rho[0];
  work.frac_121674190848 /= 2;
  work.quad_953246392320[0] = (-params.Vl_N_avg[0]-params.Thetal_N_avg[0]+params.vl_N[0])*work.frac_121674190848*(-params.Vl_N_avg[0]-params.Thetal_N_avg[0]+params.vl_N[0]);
}

//ldl.cpp
/* Produced by CVXGEN, 2014-08-13 17:31:43 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: ldl.c. */
/* Description: Basic test harness for solver.c. */

/* Be sure to place ldl_solve first, so storage schemes are defined by it. */
void Loadsolver::ldl_solve(double *target, double *var) {
  int i;
  /* Find var = (L*diag(work.d)*L') \ target, then unpermute. */
  /* Answer goes into var. */
  /* Forward substitution. */
  /* Include permutation as we retrieve from target. Use v so we can unpermute */
  /* later. */
  work.v[0] = target[0];
  /* Diagonal scaling. Assume correctness of work.d_inv. */
  for (i = 0; i < 1; i++)
    work.v[i] *= work.d_inv[i];
  /* Back substitution */
  /* Unpermute the result, from v to var. */
  var[0] = work.v[0];
#ifndef ZERO_LIBRARY_MODE
  if (settings.debug) {
    printf("Squared norm for solution is %.8g.\n", check_residual(target, var));
  }
#endif
}
void Loadsolver::ldl_factor(void) {
  work.d[0] = work.KKT[0];
  if (work.d[0] < 0)
    work.d[0] = settings.kkt_reg;
  else
    work.d[0] += settings.kkt_reg;
  work.d_inv[0] = 1/work.d[0];
#ifndef ZERO_LIBRARY_MODE
  if (settings.debug) {
    printf("Squared Frobenius for factorization is %.8g.\n", check_factorization());
  }
#endif
}
double Loadsolver::check_factorization(void) {
  /* Returns the squared Frobenius norm of A - L*D*L'. */
  double temp, residual;
  /* Only check the lower triangle. */
  residual = 0;
  temp = work.KKT[0]-1*work.d[0]*1;
  residual += temp*temp;
  return residual;
}
void Loadsolver::matrix_multiply(double *result, double *source) {
  /* Finds result = A*source. */
  result[0] = work.KKT[0]*source[0];
}
double Loadsolver::check_residual(double *target, double *multiplicand) {
  /* Returns the squared 2-norm of lhs - A*rhs. */
  /* Reuses v to find the residual. */
  int i;
  double residual;
  residual = 0;
  matrix_multiply(work.v, multiplicand);
  for (i = 0; i < 1; i++) {
    residual += (target[i] - work.v[i])*(target[i] - work.v[i]);
  }
  return residual;
}
void Loadsolver::fill_KKT(void) {
  work.KKT[0] = 2*work.frac_121674190848;
}

//testsolver.cpp
/* Produced by CVXGEN, 2014-08-13 17:31:43 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */

//Vars vars;
//Params params;
//Workspace work;
//Settings settings;
Loadsolver::Loadsolver(double Pload) {
  params.rho[0] = 1.0;
  params.Vl_N_avg[0] = 0.0;
  params.Thetal_N_avg[0] = 0.0;
  params.vl_N[0] = 0.0;
  params.L[0] = -Pload;
  global_seed = 1;
  Thiterate = NULL;
  set_defaults();
  setup_indexing();
}
Loadsolver::~Loadsolver() {
}
#define NUMTESTS 0
void Loadsolver::mainsolve(double lsRho, double Apriceavg, double Aavg, double Aprice) {
  int num_iters;
#if (NUMTESTS > 0)
  int i;
  double time;
  double time_per;
#endif
  load_default_data(lsRho, Apriceavg, Aavg, Aprice);
  /* Solve problem instance for the record. */
  settings.verbose = 1;
  num_iters = solve();
#ifndef ZERO_LIBRARY_MODE
#if (NUMTESTS > 0)
  /* Now solve multiple problem instances for timing purposes. */
  settings.verbose = 0;
  tic();
  for (i = 0; i < NUMTESTS; i++) {
    solve();
  }
  time = tocq();
  printf("Timed %d solves over %.3f seconds.\n", NUMTESTS, time);
  time_per = time / NUMTESTS;
  if (time_per > 1) {
    printf("Actual time taken per solve: %.3g s.\n", time_per);
  } else if (time_per > 1e-3) {
    printf("Actual time taken per solve: %.3g ms.\n", 1e3*time_per);
  } else {
    printf("Actual time taken per solve: %.3g us.\n", 1e6*time_per);
  }
#endif
#endif
}
void Loadsolver::load_default_data(double lsRho, double Angpriceavg, double Angavg, double Angprice) {
  params.rho[0] = lsRho;
  params.Vl_N_avg[0] = Angpriceavg;
  params.Thetal_N_avg[0] = Angavg;
  params.vl_N[0] = Angprice;
  //set_defaults();
  //setup_indexing();
}

double Loadsolver::getThetaSol(void) {
 Thiterate = vars.Thetal;
 return *Thiterate;
}
