#include "opsp_lossy_line.hpp"

#include <iostream>
using namespace std;
void LossyLine::solve(const double &rho, const double &rho_old)
{    
  terminals[0]->scale_dual_variables(rho_old/rho);
  terminals[1]->scale_dual_variables(rho_old/rho);
  
  const double s = 1.0;//rho/(rho+EPSILON);
  const double s2 = 1.0;//RATIO*rho/(RATIO*rho+EPSILON);
  
  // not vectorizable, but hey, let's give the compiler as much help as possible
  
  
  double *__restrict__ p1 = terminals[0]->p();
  double *__restrict__ p2 = terminals[1]->p();
  const double *__restrict__ u1 = terminals[0]->u();
  const double *__restrict__ u2 = terminals[1]->u();
  
  double *__restrict__ theta1 = terminals[0]->theta();
  double *__restrict__ theta2 = terminals[1]->theta();
  const double *__restrict__ v1 = terminals[0]->v();
  const double *__restrict__ v2 = terminals[1]->v();
  
  double val1, val2, val3, val4, x, y;
  bool pos;



  /*
   * minimize ||p1 - val1||_2^2 + ||p2 - val2||_2^2 + ||theta1 - val3||_2^2 + ||theta2 - val4||_2^2
   * s.t.
   *     abs(p1 - p2)/(2*b) <= (Cmax/b)
   *     (p1 + p2)/2g = 1 - cos(theta2 - theta1)
   *         equiv to
   *     (p1+p2)/2g = 1/2*((p1+p2)/2g)^2 + (1/2)*((p1-p2)/2b)^2
   *     (p1 - p2)/2b = theta2 - theta1
   *     
   *     x = p1+p2/2g, y = p1-p2/2b
   *     constraints are
   *     abs(y) <= Cmax/b
   *     x = 1/2x^2 + 1/2*y^2   -> (x-1)^2 - 1 + y^2 = 0
   *     y = theta2 - theta1
   */  
  
  //printf("%f %f %f %f\n", p1[0], p2[0], theta1[0], theta2[0]);
  for(int i = 0; i < len; ++i) {
    /*
    val1 = p1[i] - u1[i];
    val2 = p2[i] - u2[i];
    val3 = theta1[i] - v1[i];
    val4 = theta2[i] - v2[i];
    
    y = (val1 - val2)/(2.0*b);
    y = (2.0*y + RATIO*(val4 - val3))/(2.0+RATIO);

    p1[i] = b*y;//0.5*(val1 - val2);
    p2[i] = -b*y;//0.5*(val2 - val1);
    


    theta1[i] = (val3 + val4 - y)/2.0;//v1[i];
    theta2[i] = (val3 + val4 + y)/2.0;//v2[i];*/
    
    val1 = s*(p1[i] - u1[i]);
    val2 = s*(p2[i] - u2[i]);
    //theta1[i] *= M_PI/180.0;
    //theta2[i] *= M_PI/180.0;
    
    // assumes \bar v = 0, so \tilde v = v
    // this drives consensus... but we need a different rho weighting
    // *or* we need to put AC terminals on *every* device.
    //
    // i'm going to do the latter
    val3 = s2*(theta1[i] - v1[i]); val4 = s2*(theta2[i] - v2[i]);
    
    x = (val1 + val2)/(2.0*g);
    y = (val1 - val2)/(2.0*b);
    y = (2.0*y + RATIO*(val4 - val3))/(2.0+RATIO);    
    
    pos = (y >= 0);         // store the sign

    x = (x >= L) ? L : x;
    y = (y <= 0) ? -y : y;  // work with absolute value
    
    // DC CASE is commented out
    //const double xx = (x-1.0)*(x-1.0);
    const double yy = y*y;
    // 
    // // just normalize (when no phase)
    // const double d = sqrt((x-1.0)*(x-1.0) + y*y);
    // 
    // x = 1.0+(x-1.0)/d;
    // y = y/d;
    
    
    //const double yy = y*y;
    
    // this is the convex case
    // if we use the nonconvex version, projection is exact
    
    if( p*y - L*y  + M*x - p*M >= 0) {
      y = M;
      x = L;
    }
    else if( (x-1.0)*(x-1.0) + yy  - 1.0 <= 0 ) {
      //cout << "i'm inside.... " << endl;
      // do nothing
    } else {
      // XXX: use bisection for better accuracy... inaccuracies may be killing us
      // const double xx = (x - p)*(x - p);
      // const double d = sqrt(xx + yy);
      // 
      // const double u1 = y/d;
      // const double u2 = (x - p)/d;
      // 
      // // p = 1.5 - 0.5 L
      // // b = 2*(p-1) = 2*(0.5 - 0.5L) = 1 - L
      // // c = (p-1)^2 - 1 = (0.5 - 0.5L)^2 - 1 # not worth simplifying
      // // (1-L)*(1-L) = 1-M*M
      // 
      // // t = (-b*u2 + sqrt( (b*u2)^2 - 4*c ))/2
      // // t is distance from the point (0, p) in direction u = (u1,u2)
      // const double t = 0.5*(-(1.0-L)*u2 + sqrt( (1.0-M*M)*u2*u2 - 4.0*(p-1)*(p-1) + 4.0 ));
      // 
      // // (approximate) solution
      // y = t*u1;
      // x = t*u2 + p;

            
        // estimate upper bound for lambda
        // upper bound: fix v, project u to line x = L/M |y|
        double lambdaU = (M*x - L*y)/(L*y - M);
        double newx = x; 
        double newy = y;
        // lower bound: find lambda that makes v = 0
        double lambdaL = 0;//MAX(0,-x);//MAX(0, -x);
            
        double lambda = (lambdaU + lambdaL)/2.0;
        //printf("projecting %f %f\n", x, y);
        
        // no expensive divides in bisection loop
        //cout << "begin!" << endl;
        while(lambdaU - lambdaL >= BISECTION_TOL)
        {
          lambda = (lambdaU + lambdaL)/2.0;
          //cout << lambdaU << " " << lambdaL << " "<< lambda << endl;
          
          newx = (x + lambda)/(1.0 + lambda);
          newy = (RATIO+2)*y/((RATIO+2) + 2.0*lambda);
          
          //printf("new: %f,%f\n", newx, newy);
          // double s1 = (x - 1.0);
          // double s2 = y;  
          if( (newx - 1.0)*(newx - 1.0) + newy*newy < 1.0-BISECTION_TOL ) {
            lambdaU = lambda;
          }
          else
            lambdaL = lambda;
        }
      
        y = newy;
        x = newx;
        //printf("%f %f, with lambda=%f\n", x, y,lambda);
        
    }
    
    
    //cout << "x: " << x << "x1: " << x1 << endl;
    //cout << "y: " << y << "y1: " << y1 << endl;
    
    y = (pos) ? y : -y;
    
    //cout << sqrt((x-1.0)*(x-1.0) + y*y) << endl;

// XXX: it's somehow periodic... i wonder if this is because of loops?
    double a1 = (g*x + b*y);
    double a2 = (g*x - b*y);
    // eliminates edge cases
    if (a1 + a2 <= 0) {
      a1 = -a2;
    }
    p1[i] = a1;
    p2[i] = a2;
    theta1[i] = (val3 + val4 - y)/2.0;
    theta2[i] = (val3 + val4 + y)/2.0; 
    //theta1[i] *= 180.0/M_PI;
    //theta2[i] *= 180.0/M_PI;
   // if(0.5*(p1[i]+p2[i]) > 0.001) { 
   /* cout << x << " =? " << 0.5*x*x + 0.5*y*y << endl;
    cout << 0.5*(p1[i] - p2[i]) << " =? " << b*(theta2[i] - theta1[i]) << endl;
    cout << "flow: " << 0.5*(p1[i] - p2[i]) << ", loss: " << 0.5*(p1[i]+p2[i]) << endl; 
    cout << 0.5*(p1[i]+p2[i]) << " =? " << g*(1.0 - cos(theta2[i] - theta1[i])) << endl;
    cout << "theta1: " << theta1[i] << " theta2: " << theta2[i] << endl;
    cout << endl; */
    //}
  }
  
  //printf("%f %f %f %f\n\n", p1[0], p2[0], theta1[0], theta2[0]);
  
}
