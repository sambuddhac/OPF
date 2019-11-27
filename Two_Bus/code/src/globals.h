#ifndef OPF_GLOBALS_H
#define OPF_GLOBALS_H

#define MAX_TIME 96 // time horizon
#define NUMTYPES 6  // number of bus types
#define ABSTOL 1e-3

// parameters for rho update
#define LAMBDA    0.0001//0.01
#define MU        0.0005//0.05

// quadratic term for penalty on lossless lines
// (for creating networks)
#define EPSILON   1e-2

#define CVXGEN_EPS          1e-8
#define CVXGEN_RESID_TOL    1e-8
#define CVXGEN_ITERS        20
#define CVXGEN_REFINE_STEPS 1
#define CVXGEN_KKT_REG      1e-12

#define BISECTION_TOL       1e-8

#endif
