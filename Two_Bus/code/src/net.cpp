#include "net.hpp"
#include "math.h"

bool Net::address_exists(const Terminal *t, int low, int high)
{
  if(high < low) return false;
  int mid = (low + high)/2;
  if(terminals[mid] > t)
    return address_exists(t, low, mid-1);
  else if(terminals[mid] < t)
    return address_exists(t, mid+1, high);
  else
    return true;
}


void Net::add_terminal(const Terminal *t)
{
  // does nothing if you pass in a null pointer
  if(t != NULL) {
    //printf("terminal array!\n");
    int i = N;
    bool unique = true;
    if(N == 0) {
        terminals = (Terminal **) malloc(sizeof(Terminal *));
    } else {
      // O(log N) search to maintain uniqueness
      unique = !address_exists(t, 0, N);
    
      terminals = (Terminal **) realloc(terminals, (N+1)*sizeof(Terminal *));
      // O(N) insertion to maintain list of sorted addresses
      while(terminals[i-1] > t && i > 0) {
        terminals[i] = terminals[i-1];
        i--;
      }
    }
    if(unique) {
      terminals[i] = (Terminal *) t;  // violates const-correctness!
      N++;
    }
  }
}

void Net::reduce(double &primal_residual, double &dual_residual)
{
  // TODO: terminals connected to net must all have the same size
  // TODO: get length of terminal here
  double *__restrict__ pbar = new double[MAX_TIME]();
  double *__restrict__ thetabar = new double[MAX_TIME]();
  const double *__restrict__ p;
  const double *__restrict__ theta;
  
  
  /* BEGIN DC TERMINALS */
  // printf("nAC: %d\n", nAC);
  // printf("DC: ");
  //printf("\n(p, theta, u, v)\n");
  for(int i = 0; i < N; i++)
  {
    p = terminals[i]->p();
    theta = terminals[i]->theta();
    //printf("%f ", p[0]);
    for(int j = 0; j < MAX_TIME; j++) {
      pbar[j] += p[j]/N;
      thetabar[j] += theta[j]/N;
    }
  }
  //printf("BAR (%f, %f)\n", pbar[29], thetabar[29]);
  //printf(" = %f\n", pbar[0]);
  
  // printf("******\n");
  // printf("%d : ", N);
  // for(int i = 0; i < N; i++)
  // {
  //   for(int j = 0; j < MAX_TIME; j++)
  //     printf("%f ", terminals[i]->p(j));
  //   printf("\n");
  // }
  // printf("--\n");
  // for(int j = 0; j < MAX_TIME; j++) 
  //   printf("%f ", pbar[j]);
  // printf("******\n");
  //printf("before DC %f %f\n", sqrt(primal_residual), sqrt(dual_residual));
  
  // can't vectorize what follows
  for(int i = 0; i < N; i++)
  {
    //std::cout << pbar[i] << std::endl;
    terminals[i]->update(pbar, thetabar);
    //printf("(%f, %f %f, %f)\n", terminals[i]->p()[29], terminals[i]->theta()[29], terminals[i]->u()[29], terminals[i]->v()[29]);
    //printf("r:%f, s:%f\n", terminals[i]->r_square(), terminals[i]->s_square());
    primal_residual += terminals[i]->r_square();
    dual_residual += terminals[i]->s_square();
  }
  //printf("%f %f\n", primal_residual, dual_residual);
  
  delete[] pbar;
  delete[] thetabar;
}

void Net::print()
{
  std::cout << "Net with " << N << " terminals." << std::endl;
  for(int i = 0; i < N; i++)
  {
    std::cout << *terminals[i] << std::endl;
  }
}
