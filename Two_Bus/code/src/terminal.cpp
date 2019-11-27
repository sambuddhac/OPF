#include "terminal.hpp"


//3.80241e-15 2.4398e-14 4.05635 -0.305277
//----
//1.82736e-11 1.11208e-12 1 15.4329
//8.62952e-13 5.55235e-14 1 14.4083

// void Terminal::update_p(const double &ratio)
// {
//   for(int i = 0; i < T; i++) {    
//     this->_u[i] *= ratio;
//     this->_p[i] -= this->_u[i];
//   }
// }

// define a cout operator
std::ostream& operator<<(std::ostream& out, const Terminal& term)
{
  out << "DC terms" << std::endl;
  for(int i = 0; i < term.T; i++)
    out << term._p[i] << " ";
  out << std::endl << "AC terms" << std::endl;
  for(int i = 0; i < term.T; i++)
    out << term._theta[i] << " ";
  return out;
}

void Terminal::scale_dual_variables(const double &ratio) {
  double *__restrict__ u = this->_u;
  double *__restrict__ v = this->_v;
  for(int i = 0; i < this->T; i++) {
    u[i] *= ratio;//(*(gen->vars.p[i]), i);
    v[i] *= ratio;
  }
}

void Terminal::set_pvec(double *__restrict__*__restrict__ from)
{
  // only sets p, but also does the update for theta!
  double *__restrict__ p = this->_p;
  
  for(int i = 0; i < this->T; i++) {
    #ifdef DEBUG
    assert(from[i]);
    #endif
    p[i] = *(from[i]);//(*(gen->vars.p[i]), i);
  }
}

void Terminal::copy_pvec(double*__restrict__*__restrict__ to, const double &ratio)
{
  // only copies p and also scales dual variables
  // also updates theta (since objects that call this don't have phase dependency)
  double *__restrict__ u = this->_u;
  double *__restrict__ v = this->_v;
  double *__restrict__ p = this->_p;
  double *__restrict__ theta = this->_theta;
  
  for(int i = 0; i < this->T; i++) {
    #ifdef DEBUG
    assert(to[i]);
    #endif
    u[i] *= ratio; // update the dual varaible
    v[i] *= ratio; // update the dual variable
    
    *(to[i]) = p[i] - u[i];// set the parameter properly
    //theta[i] -= v[i];
  }  
}


void Terminal::update(const double * __restrict__ pbar, const double * __restrict__ thetabar)
{
  //static const double alpha = 1.0;
  
  primal_residual = 0.0;
  dual_residual = 0.0;
  
  // allows vectorized code
  double *__restrict__ p = _p;
  double *__restrict__ u = _u;
  
  double *__restrict__ theta = _theta;
  double *__restrict__ v = _v;
  //double *__restrict__ d = _d;
  
  double *__restrict__ oldp = _oldp;
  double *__restrict__ oldtheta = _oldtheta;
  
  for(int i = 0; i < T; i++) {
    //assert(pbar[i]);
    
    // update p
    //if(i == 0) {
    dual_residual += (-oldp[i] - pbar[i] + p[i])*(-oldp[i] - pbar[i] + p[i]);
    primal_residual += pbar[i]*pbar[i];
    //}
    
    u[i] += pbar[i];
    //printf("(%f -> %f) ", _oldp[i], _p[i] - pbar[i]);
    
    oldp[i] = p[i] - pbar[i];  // p - pbar
    // for overrelaxation
    //d[i] = (p[i] - pt) + (1.0 - alpha)*d[i];
    
    p[i] -= pbar[i];//(-pt + (1.0 - alpha)*d[i]);// p - pbar + d
    
    // update thetas
    //if(i==0) {
    dual_residual += (thetabar[i] - oldtheta[i])*(thetabar[i] - oldtheta[i]);
    primal_residual += (theta[i] - thetabar[i])*(theta[i] - thetabar[i]);
  //}
    v[i] += (theta[i] - thetabar[i]);
    
    oldtheta[i] = theta[i];  
    theta[i] = thetabar[i];// - u[i];// - p[i]; /
    
    //printf("%d: %f %f\n", i+1, primal_residual, dual_residual);
  }
  //printf("%f\n", pbar[0]);
}
