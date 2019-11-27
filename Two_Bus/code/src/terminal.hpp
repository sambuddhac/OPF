#ifndef OPSP_TERMINAL_H
#define OPSP_TERMINAL_H

#include <iostream>
#include <assert.h>

// TODO: terminal = vector
// TODO: dual terminal = vector
class Terminal {
public:
  Terminal(int N) :
    _p(new double[N]()), _u(new double[N]()), 
    _theta(new double[N]()), _v(new double[N]()),
    /*_d(new double[N]()),*/ 
    _oldp(new double[N]()), _oldtheta(new double[N]()),
    primal_residual(0.0), dual_residual(0.0), T(N)
  { }
  
  inline const double &p(int i) const;
  inline const double &u(int i) const;
  inline const double &theta(int i) const;
  inline const double &v(int i) const;

  // these accessors basically expose everything
  inline double *&p() { return _p; }
  inline double *&u() { return _u; }
  inline double *&oldp() { return _oldp; }
  
  inline double *&theta() { return _theta; }
  inline double *&v() { return _v; }
  inline double *&oldtheta() { return _oldtheta; }
  
  //inline void set_oldp(const double &p, int i);
  //inline void set_p(const double &p, int i);
  //inline void set_u(const double &u, int i);
  void scale_dual_variables(const double &ratio);
  
  // define a cout operator
  friend std::ostream& operator<<(std::ostream& out, const Terminal& term);
  
  void set_pvec(double*__restrict__*__restrict__ from);
  void copy_pvec(double*__restrict__*__restrict__ to, const double &ratio);
  
  inline const double &r_square() const;
  inline const double &s_square() const;
  
  void update(const double *pbar, const double *thetabar);
  
  virtual ~Terminal()
  {
    if(_p) delete[] _p;
    if(_u) delete[] _u;
    if(_oldp) delete[] _oldp;
    if(_theta) delete[] _theta;
    if(_v) delete[] _v;
    if(_oldtheta) delete[] _oldtheta;
    //if(_d) delete[] _d;
  }
  
protected:
  double *_p;
  double *_u;
  double *_theta;
  double *_v;
  
  //double *_d; // for overrelaxtion
  
  double *_oldp;
  double *_oldtheta;
  
  double primal_residual;
  double dual_residual;
  
  const int T;
private:
  // disallow copy constructor
  Terminal(const Terminal&);
  Terminal& operator=(const Terminal&);
};


inline const double &Terminal::p(int i) const {
  #ifdef DEBUG
  assert(i < T && i >= 0);
  #endif
  return this->_p[i];
}

inline const double &Terminal::u(int i) const { 
  #ifdef DEBUG
  assert(i < T && i >= 0);
  #endif
  return this->_u[i];
}

inline const double &Terminal::theta(int i) const {
  #ifdef DEBUG
  assert(i < T && i >= 0);
  #endif
  return this->_theta[i];
}

inline const double &Terminal::v(int i) const {
  #ifdef DEBUG
  assert(i < T && i >= 0);
  #endif
  return this->_v[i];
}

// inline void Terminal::set_oldp(const double &p, int i) {
//   #ifdef DEBUG
//   assert(i < T && i >= 0);
//   #endif
//   this->_oldp[i] = p; 
// }
// 
// inline void Terminal::set_p(const double &p, int i) {
//   #ifdef DEBUG
//   assert(i < T && i >= 0);
//   #endif
//   this->_p[i] = p; 
// }
// 
// inline void Terminal::set_u(const double &u, int i) { 
//   #ifdef DEBUG
//   assert(i < T && i >= 0);
//   #endif
//   this->_u[i] = u; 
// }

inline const double &Terminal::r_square() const {
  return primal_residual;
}

inline const double &Terminal::s_square() const {  
  return dual_residual;
}
  
#endif // OPSP_TERMINAL_H
