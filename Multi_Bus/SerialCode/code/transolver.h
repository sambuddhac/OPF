/* Produced by CVXGEN, 2014-08-13 19:27:00 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: solver.h. */
/* Description: Header file with relevant definitions. */
#ifndef TRANSOLVER_H
#define TRANSOLVER_H
#include <time.h>
#include <stdlib.h>
//#include "mex.h"
/* Uncomment the next line to remove all library dependencies. */
/*#define ZERO_LIBRARY_MODE */
#ifdef MATLAB_MEX_FILE
/* Matlab functions. MATLAB_MEX_FILE will be defined by the mex compiler. */
/* If you are not using the mex compiler, this functionality will not intrude, */
/* as it will be completely disabled at compile-time. */
//#include "mex.h"
#else
#ifndef ZERO_LIBRARY_MODE
#include <stdio.h>
#endif
#endif
/* Space must be allocated somewhere (testsolver.c, csolve.c or your own */
/* program) for the global variables vars, params, work and settings. */
/* At the bottom of this file, they are externed. */
#ifndef ZERO_LIBRARY_MODE
#include <math.h>
#define pm(A, m, n) printmatrix(#A, A, m, n, 1)
#endif

class Transolver {
private:
typedef struct Params_t {
  double rho[1];
  double Pt_N_init[2];
  double Pt_N_avg[2];
  double ut_N[2];
  double Vt_N_avg[2];
  double Thetat_N_avg[2];
  double vt_N[2];
  double PtMin[1];
  double PtMax[1];
  double Xt[1];
} Params;
typedef struct Vars_t {
  double *Pt; /* 2 rows. */
  double *Thetat; /* 2 rows. */
} Vars;
typedef struct Workspace_t {
  double h[4];
  double s_inv[4];
  double s_inv_z[4];
  double b[2];
  double q[4];
  double rhs[14];
  double x[14];
  double *s;
  double *z;
  double *y;
  double lhs_aff[14];
  double lhs_cc[14];
  double buffer[14];
  double buffer2[14];
  double KKT[26];
  double L[15];
  double d[14];
  double v[14];
  double d_inv[14];
  double gap;
  double optval;
  double ineq_resid_squared;
  double eq_resid_squared;
  double block_33[1];
  /* Pre-op symbols. */
  double frac_121674190848;
  double quad_373534269440[1];
  double quad_166572511232[1];
  double frac_994642780160;
  int converged;
} Workspace;
typedef struct Settings_t {
  double resid_tol;
  double eps;
  int max_iters;
  int refine_steps;
  int better_start;
  /* Better start obviates the need for s_init and z_init. */
  double s_init;
  double z_init;
  int verbose;
  /* Show extra details of the iterative refinement steps. */
  int verbose_refinement;
  int debug;
  /* For regularization. Minimum value of abs(D_ii) in the kkt D factor. */
  double kkt_reg;
} Settings;
Vars vars;
Params params;
Workspace work;
Settings settings;
double *Piterate1;
double *Thiterate1;
double *Piterate2;
double *Thiterate2;
long global_seed;
clock_t tic_timestart;
/* Function definitions in ldl.c: */
public:
Transolver(double, double, double); // Transolver object constructor
~Transolver(); // Transolver object destructor
void ldl_solve(double *target, double *var);
void ldl_factor(void);
double check_factorization(void);
void matrix_multiply(double *result, double *source);
double check_residual(double *target, double *multiplicand);
void fill_KKT(void);

/* Function definitions in matrix_support.c: */
void multbymA(double *lhs, double *rhs);
void multbymAT(double *lhs, double *rhs);
void multbymG(double *lhs, double *rhs);
void multbymGT(double *lhs, double *rhs);
void multbyP(double *lhs, double *rhs);
void fillq(void);
void fillh(void);
void fillb(void);
void pre_ops(void);

/* Function definitions in solver.c: */
double eval_gap(void);
void set_defaults(void);
void setup_pointers(void);
void setup_indexing(void);
void set_start(void);
double eval_objv(void);
void fillrhs_aff(void);
void fillrhs_cc(void);
void refine(double *target, double *var);
double calc_ineq_resid_squared(void);
double calc_eq_resid_squared(void);
void better_start(void);
void fillrhs_start(void);
long solve(void);

/* Function definitions in testsolver.c: */
void mainsolve(double, double, double, double, double, double, double, double, double, double, double, double, double);
void load_default_data(double, double, double, double, double, double, double, double, double, double, double, double, double);
double getPSol1(void);
double getThetaSol1(void);
double getPSol2(void);
double getThetaSol2(void);

/* Function definitions in util.c: */
void tic(void);
float toc(void);
float tocq(void);
void printmatrix(char *name, double *A, int m, int n, int sparse);
double unif(double lower, double upper);
float ran1(long*idum, int reset);
float randn_internal(long *idum, int reset);
double randn(void);
void reset_rand(void);
};

#endif
