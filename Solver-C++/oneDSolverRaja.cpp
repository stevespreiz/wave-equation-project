/* n-Dimensional Wave Equation Solver
*
*  Currently being built to support 1D then 2D
*
*/

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>

#include "RAJA/RAJA.hpp"

using namespace std;

// declaring LAPACK linear system solver (compile with  -llapack  flag)
// dgesv_ is a symbol in the LAPACK library files
extern "C" {
extern int dgesv_(int*,int*,double*,int*,int*,double*,int*,int*);
}

// Initial Position
class F{
public:
  double operator() (double x) { return sin(3*M_PI*x/2.); }
};

// Initial Velocity
class G{
public:
  double operator() (double x) { return 0; }
};

// Left side BC
class L{
public:
  double operator() (double x) { return 0; }
};

// Right side BC
class R{
public:
  double operator() (double x) { return 0; }
};

class Y{
public:
  double operator() (double x, double t, F f, int c = 1) {return 0.5*(f(x-c*t)+f(x+c*t));}
};


struct Definition {
  double a;
  double b;
  double c;
  int N;
  F f;
  G g;
  L l;
  R r;
};

void firstStep(Definition* def, double sigma, double* x, double dt, const int ja, const int jb,
        double* unm1, double* un, int nD, int oacc){
  if(nD ==  1){
    RAJA::forall<RAJA::loop_exec>(RAJA::RangeSegment(ja,jb+1),[=](int i){
      un[i] = unm1[i]
            + dt*def->g(x[i])
            + pow(sigma,2)/2*(unm1[i+1]-2*unm1[i]+unm1[i-1]);
    });
  }
}

void timeStep(double sigma, const int ja, const int jb, double* unm1, double* un, double* unp1,
        int nD, int oacc){
  if(nD == 1){
    RAJA::forall<RAJA::loop_exec>(RAJA::RangeSegment(ja,jb+1),[=] (int i){
      unp1[i] = 2*un[i] - unm1[i] + pow(sigma,2)*(un[i+1]-2*un[i]+un[i-1]);
    });
  }
}

void BC(Definition* def, double sigma, double* x, int n, double dt, const int ja, const int jb,
        double* unp1, int nD, int iCase, int oacc){
  double dx = x[1]-x[0];
  if(nD == 1){
    // Dirchlet left (assume l_tt = 0), Neumann right
    if(iCase == 1){
      if(oacc == 2){
        unp1[jb+1] = unp1[jb-1] + 2*dx*def->r(n*dt);
        unp1[ja-1] = 2*unp1[ja]-unp1[ja+1];
      }
    }
    // Dirchlet both sides (assume l_tt = r_tt = 0)
    if(iCase == 2){
      if(oacc == 2){
        unp1[ja-1] = 2*unp1[ja]-unp1[ja+1];
        unp1[jb+1] = 2*unp1[jb]-unp1[jb-1];
      }
    }
  }
}



int main(int argc, char* argv[]){
  /*  Inputs:
  *   defintion - struct - containing a,b,c,N,f,g,and BC functions
  *   sigma     - double - CFL parameter
  *   tf        - double - final time
  *   nD        - int    - number of dimensions
  *   icase     - int    - case to solve
  *   oacc      - int    - order of accuracy
  */

  ofstream fout("out.txt");

  // Problem Initialization - eventually move these to command line arguments
  Definition* def = new Definition;
  F f;
  G g;
  L l;
  R r;
  Y y;

  def->a = 0;
  def->b = 1;
  def->c = 1;
  def->N = atoi(argv[6]); // just to make testing easier
  def->f = f;
  def->g = g;
  def->l = l;
  def->r = r;

  double sigma = stod(argv[1]);
  double tf    = stod(argv[2]);
  int nD       = atoi(argv[3]);
  int icase    = atoi(argv[4]);
  int oacc     = atoi(argv[5]);

  /////////////////////////////////////////////////////////////////////////////
  //  Setup

  //  Set indexing variabls
  const int ja = oacc/2;
  const int jb = def->N + oacc/2;

  //  Steps in space
  double dx = (def->b - def->a)/def->N;
  int arrSize = def->N + 1 + oacc;
  double* x = new double[arrSize];
  for(int i = 0; i < arrSize; i++){
    x[i] = (i - oacc/2)*dx;
    fout << x[i] << "\t";
  }
  fout << endl;

  //  Step in time
  double dttilde = sigma*dx/def->c;
  int nt         = ceil(tf/dttilde);
  double dt      = tf/nt;

  //  Update sigma
  sigma  = def->c * dt/dx;

  //  Initialize grid
  double* unm1 = new double[arrSize];
  double* un   = new double[arrSize];
  double* unp1 = new double[arrSize];

  /////////////////////////////////////////////////////////////////////////////
  //  Initial Condition
  for(int i = 0; i < arrSize; i++){
    unm1[i] = def->f(x[i]);
    fout << unm1[i] << "\t";
  }
  fout << endl;

  /////////////////////////////////////////////////////////////////////////////
  //  First Time Step

  firstStep(def, sigma, x, dt, ja, jb, unm1, un, nD, oacc);
  BC(def,sigma,x,1,dt,ja,jb,un,nD,icase,oacc);

  for(int i = 0; i < arrSize ; i++)
    fout<<un[i] << "\t";
  fout << endl;
  /////////////////////////////////////////////////////////////////////////////
  //  Rest of the time steps
  int n = 2;
  while(n*dt <= tf){
    // Update middle values
    timeStep(sigma,ja,jb,unm1,un,unp1,nD,oacc);

    // Update boundary conditions
    BC(def,sigma,x,n,dt,ja,jb,unp1,nD,icase,oacc);

    // Update arrays
    for(int i = 0; i<arrSize; i++){
      unm1[i] = un[i];
      un[i] = unp1[i];
    }


    for(int i = 0; i< arrSize; i++){
      fout << un[i] << "\t";
    }
    fout << endl;

    n++;
  }

  /////////////////////////////////////////////////////////////////////////////
  //  Error
  double e = 0;
  for(int i = ja; i <= jb; i++){
    e = max(e,unp1[i]-y(x[i],tf,f));
  }
  cout << "inf norm err: " << e << endl;

  /////////////////////////////////////////////////////////////////////////////
  //  Output
  // for(int i = ja; i <= jb; i++){
  //   cout << "x: " << x[i] << "\tu = " << un[i] << endl;
  // }

  fout.close();
  /////////////////////////////////////////////////////////////////////////////
  //  Memory Cleanup
  delete[] unm1;
  delete[] un;
  delete[] unp1;
  delete[] x;
  delete def;

  return 0;
}
