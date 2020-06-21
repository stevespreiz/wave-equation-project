/* n-Dimensional Wave Equation Solver
*
*  Currently being built to support 1D then 2
*
*/

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

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

void timeStep(double sigma, int ja, int jb, double* unm1, double* un, double* unp1,
        int nD, int oacc){
  if(nD == 1){
    for(int i = ja; i <= jb; i++){
      unp1[i] = 2*un[i] - unm1[i] + pow(sigma,2)*(un[i+1]-2*un[i]+un[i-1]);
      if(oacc > 2){
        unp1[i] += (pow(sigma,2)-pow(sigma,4))/12*(un[i+2]-4*un[i+1]+6*un[i]-4*un[i-1]+un[i-2]);
      }
    }
  }
}

void BC(Definition* def, double sigma, double* x, int n, double dt, int ja, int jb,
        double* unp1, int nD, int iCase, int oacc){
  if(nD == 1){
    // Dirchlet left, Neumann right
    if(iCase == 1){
      if(oacc == 2){
        unp1[jb+1] = unp1[jb-1];
        unp1[ja-1] = 2*unp1[ja]-unp1[ja+1];
      }
    }
    // Dirchlet both sides
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

  def->a = 0;
  def->b = 1;
  def->c = 1;
  def->N = 100;
  def->f = f;
  def->g = g;
  def->l = l;
  def->r = r;

  double sigma = 0.9;
  double tf    = 3;
  int nD       = 1;
  int icase    = 1;
  int oacc     = 2;

  /////////////////////////////////////////////////////////////////////////////
  //  Setup

  //  Set indexing variabls
  int ja = oacc/2;
  int jb = def->N + oacc/2;

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
  double dttilde = sigma*dx*def->c;
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
  if(nD == 1){
    if(oacc == 2){
      for(int i = ja; i <= jb; i++){
        un[i] = (1-pow(sigma,2))*unm1[i] + dt*def->g(x[i])+pow(sigma,2)/2*(unm1[i-1]+unm1[i+1]);

      }
    }
    else if(oacc == 4){
      for(int i = ja; i<= jb; i++){

      }
    }
  }


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

    // Update array pointers
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
  //  Output
  for(int i = ja; i <= jb; i++){
    cout << "x: " << x[i] << "\tu = " << un[i] << endl;
  }

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
