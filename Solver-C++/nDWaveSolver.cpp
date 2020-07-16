/* n-Dimensional Wave Equation Solver
*
*  Currently being built to support 1D then 2D
*
*/

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
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

void firstStep(Definition* def, double sigma, double* x, double dt, int ja, int jb,
        double* unm1, double* un, int nD, int oacc){
  if(nD ==  1){
    for(int i = ja; i<=jb; i++){
      un[i] = unm1[i]
            + dt*def->g(x[i])
            + pow(sigma,2)/2*(unm1[i+1]-2*unm1[i]+unm1[i-1]);
      if(oacc > 2){
        un[i] += -1*pow(sigma,2)/2*(unm1[i+2]-4*unm1[i+1]+6*unm1[i]-4*unm1[i-1]+unm1[i-2])/12
               + dt*pow(sigma,2)/6*( (def->g(x[i+1])-2*def->g(x[i])+def->g(x[i-1])) - (def->g(x[i+2])-4*def->g(x[i+1])+6*def->g(x[i])-4*def->g(x[i-1])+def->g(x[i-2]))/12)
               + pow(sigma,4)/24*(unm1[i+2]-4*unm1[i+1]+6*unm1[i]-4*unm1[i-1]+unm1[i-2]);
        if(oacc > 4){
          un[i] += pow(sigma,2)/2*(unm1[i+3]-6*unm1[i+2]+15*unm1[i+1]-20*unm1[i]+15*unm1[i-1]-6*unm1[i-2]+unm1[i-3])/90+
                   dt*pow(sigma,2)/6*(def->g(x[i+3])-6*def->g(x[i+2])+15*def->g(x[i+1])-20*def->g(x[i])+15*def->g(x[i-1])-6*def->g(x[i-2])+def->g(x[i-3]))/90+
                   -pow(sigma,4)/24*(unm1[i+3]-6*unm1[i+2]+15*unm1[i+1]-20*unm1[i]+15*unm1[i-1]-6*unm1[i-2]+unm1[i-3])/72+
                   -dt*pow(sigma,4)/120*(def->g(x[i+3])-6*def->g(x[i+2])+15*def->g(x[i+1])-20*def->g(x[i])+15*def->g(x[i-1])-6*def->g(x[i-2])+def->g(x[i-3]))+
                   pow(sigma,6)/720*(unm1[i+3]-6*unm1[i+2]+15*unm1[i+1]-20*unm1[i]+15*unm1[i-1]-6*unm1[i-2]+unm1[i-3]);
        }
      }
    }
  }
}

void timeStep(double sigma, int ja, int jb, double* unm1, double* un, double* unp1,
        int nD, int oacc){
  if(nD == 1){
    for(int i = ja; i <= jb; i++){
      unp1[i] = 2*un[i] - unm1[i] + pow(sigma,2)*(un[i+1]-2*un[i]+un[i-1]);
      if(oacc > 2){
        unp1[i] -= (pow(sigma,2)-pow(sigma,4))/12*(un[i+2]-4*un[i+1]+6*un[i]-4*un[i-1]+un[i-2]);
        if(oacc > 4){
          unp1[i] += (pow(sigma,2)/90-pow(sigma,4)/72+pow(sigma,6)/360)*(un[i+3]-6*un[i+2]+15*un[i+1]-20*un[i]+15*un[i-1]-6*un[i-2]+un[i-3]);
        }
      }
    }
  }
}

void BC(Definition* def, double sigma, double* x, int n, double dt, int ja, int jb,
        double* unp1, int nD, int iCase, int oacc){
  double dx = x[1]-x[0];
  if(nD == 1){
    // Dirchlet left (assume l_tt = 0), Neumann right
    if(iCase == 1){
      if(oacc == 2){
        unp1[jb+1] = unp1[jb-1] + 2*dx*def->r(n*dt);
        unp1[ja-1] = 2*unp1[ja]-unp1[ja+1];
      }
      else if(oacc == 4){
        // Left hand Dirchlet DDFA
        double* f0 = new double[2];
        double* f1 = new double[2];
        double* f2 = new double[2];

        f0[0] = pow(def->c,2)/pow(dx,2)*(-2*unp1[ja]+unp1[ja+1] - (6*unp1[ja]-4*unp1[ja+1]+unp1[ja+2])/12);
        f0[1] = pow(def->c,2)/pow(dx,4)*(6*unp1[ja]-4*unp1[ja+1]+unp1[ja+2]);

        f1[0] = pow(def->c,2)/pow(dx,2)*(1-2*unp1[ja]+unp1[ja+1] - (-4+6*unp1[ja]-4*unp1[ja+1]+unp1[ja+2])/12);
        f1[1] = pow(def->c,2)/pow(dx,4)*(-4+6*unp1[ja]-4*unp1[ja+1]+unp1[ja+2]);

        f2[0] = pow(def->c,2)/pow(dx,2)*(-2*unp1[ja]+unp1[ja+1] - (1+6*unp1[ja]-4*unp1[ja+1]+unp1[ja+2])/12);
        f2[1] = pow(def->c,2)/pow(dx,4)*(1+6*unp1[ja]-4*unp1[ja+1]+unp1[ja+2]);

        int n = 2;
        int nrhs = 1;
        double*A = new double[n*n];
        int lda = 2;
        double* b = new double[n];
        int ldb = 2;
        int* ipiv = new int [n];
        int info;

        A[0] = f1[0]-f0[0];
        A[1] = f1[1]-f0[1];
        A[2] = f2[0]-f0[0];
        A[3] = f2[1]-f0[1];

        b[0] = -1*f0[0];
        b[1] = -1*f0[1];

        dgesv_(&n,&nrhs,A,&lda,ipiv,b,&ldb,&info);

        unp1[ja-1] = b[0];
        unp1[ja-2] = b[1];

        // Right hand Neumann
        f0[0] = 0 - unp1[jb-1] - 1/3*(0 - 2*0 + 2*unp1[jb-1]-unp1[jb-2]);
        f0[1] =                       0 - 2*0 + 2*unp1[jb-1]-unp1[jb-2];

        f1[0] = 1-unp1[jb-1] - 1/3*(0 - 2*1 +2 *unp1[jb-1] - unp1[jb-2]);
        f1[1] =                     0 - 2*1 + 2*unp1[jb-1] - unp1[jb-2];

        f2[0] = 0 - unp1[jb-1] - 1/3*(1- 2*0 + 2*unp1[jb-1] - unp1[jb-2]);
        f2[1] =                      1 - 2*0 + 2*unp1[jb-1] - unp1[jb-2];

        A[0] = f1[0]-f0[0];
        A[1] = f1[1]-f0[1];
        A[2] = f2[0]-f0[0];
        A[3] = f2[1]-f0[1];

        b[0] = -1*f0[0];
        b[1] = -1*f0[1];

        dgesv_(&n,&nrhs,A,&lda,ipiv,b,&ldb,&info);

        unp1[jb+1] = b[0];
        unp1[jb+2] = b[1];

        delete[] A;
        delete[] b;
        delete[] ipiv;
        delete[] f0;
        delete[] f1;
        delete[] f2;
      }
      else if(oacc == 6){
        // Left hand Dirchlet
        double* f0 = new double[3];
        double* f1 = new double[3];
        double* f2 = new double[3];
        double* f3 = new double[3];

        f0[0] = 1/180.*(2*0 - 27*0 + 270*0 - 490*unp1[ja]  + 270*unp1[ja+1] - 27*unp1[ja+2] + 2*unp1[ja+3]);
        f0[1] = 1/6.*  ( -0 + 12*0 - 39*0  + 56* unp1[ja]  - 39* unp1[ja+1] + 12*unp1[ja+2] -  unp1[ja+3]);
        f0[2] =           0 - 6*0  + 15*0  - 20* unp1[ja]  + 15* unp1[ja+1] - 6* unp1[ja+2] +  unp1[ja+3];

        f1[0] = 1/180.*(2*0 - 27*0 + 270*1 - 490*unp1[ja]  + 270*unp1[ja+1] - 27*unp1[ja+2] + 2*unp1[ja+3]);
        f1[1] = 1/6.*  ( -0 + 12*0 - 39*1  + 56* unp1[ja]  - 39* unp1[ja+1] + 12*unp1[ja+2] -   unp1[ja+3]);
        f1[2] =           0 - 6*0  + 15*1  - 20* unp1[ja]  + 15* unp1[ja+1] - 6* unp1[ja+2] +   unp1[ja+3];

        f2[0] = 1/180.*(2*0 - 27*1 + 270*0 - 490*unp1[ja]  + 270*unp1[ja+1] - 27*unp1[ja+2] + 2*unp1[ja+3]);
        f2[1] = 1/6.*  ( -0 + 12*1 - 39*0  + 56* unp1[ja]  - 39* unp1[ja+1] + 12*unp1[ja+2] -   unp1[ja+3]);
        f2[2] =           0 - 6* 1 + 15*0  - 20* unp1[ja]  + 15* unp1[ja+1] - 6* unp1[ja+2] +   unp1[ja+3];

        f3[0] = 1/180.*(2*1 - 27*0 + 270*0 - 490*unp1[ja]  + 270*unp1[ja+1] - 27*unp1[ja+2] + 2*unp1[ja+3]);
        f3[1] = 1/6.*  ( -1 + 12*0 - 39*0  + 56* unp1[ja]  - 39* unp1[ja+1] + 12*unp1[ja+2] -   unp1[ja+3]);
        f3[2] =           1 - 6*0  + 15*0  - 20* unp1[ja]  + 15* unp1[ja+1] - 6* unp1[ja+2] +   unp1[ja+3];


        int n = 3;
        int nrhs = 1;
        double*A = new double[n*n];
        int lda = n;
        double* b = new double[n];
        int ldb = n;
        int* ipiv = new int [n];
        int info;

        A[0] = f1[0]-f0[0];
        A[1] = f1[1]-f0[1];
        A[2] = f1[2]-f0[2];
        A[3] = f2[0]-f0[0];
        A[4] = f2[1]-f0[1];
        A[5] = f2[2]-f0[2];
        A[6] = f3[0]-f0[0];
        A[7] = f3[1]-f0[1];
        A[8] = f3[2]-f0[2];

        b[0] = -1*f0[0];
        b[1] = -1*f0[1];
        b[2] = -1*f0[2];

        dgesv_(&n,&nrhs,A,&lda,ipiv,b,&ldb,&info);

        unp1[ja-1] = b[0];
        unp1[ja-2] = b[1];
        unp1[ja-3] = b[2];


        // Right hand Neumann
        f0[0] = 1/60.*(-1*unp1[jb-3] + 9*unp1[jb-2] - 45*unp1[jb-1] + 45*0 - 9*0 + 1*0);
        f0[1] = 1/8.* ( 1*unp1[jb-3] - 8*unp1[jb-2] + 13*unp1[jb-1] - 13*0 + 8*0 - 1*0);
        f0[2] = 1/2.* (-1*unp1[jb-3] + 4*unp1[jb-2] - 5* unp1[jb-1] + 5*0  - 4*0 + 1*0);

        f1[0] = 1/60.*(-1*unp1[jb-3] + 9*unp1[jb-2] - 45*unp1[jb-1] + 45*1 - 9*0 + 1*0);
        f1[1] = 1/8.* ( 1*unp1[jb-3] - 8*unp1[jb-2] + 13*unp1[jb-1] - 13*1 + 8*0 - 1*0);
        f1[2] = 1/2.* (-1*unp1[jb-3] + 4*unp1[jb-2] - 5* unp1[jb-1] + 5* 1 - 4*0 + 1*0);

        f2[0] = 1/60.*(-1*unp1[jb-3] + 9*unp1[jb-2] - 45*unp1[jb-1] + 45*0 - 9*1 + 1*0);
        f2[1] = 1/8.* ( 1*unp1[jb-3] - 8*unp1[jb-2] + 13*unp1[jb-1] - 13*0 + 8*1 - 1*0);
        f2[2] = 1/2.* (-1*unp1[jb-3] + 4*unp1[jb-2] - 5* unp1[jb-1] + 5* 0 - 4*1 + 1*0);

        f3[0] = 1/60.*(-1*unp1[jb-3] + 9*unp1[jb-2] - 45*unp1[jb-1] + 45*0 - 9*0 + 1*1);
        f3[1] = 1/8.* ( 1*unp1[jb-3] - 8*unp1[jb-2] + 13*unp1[jb-1] - 13*0 + 8*0 - 1*1);
        f3[2] = 1/2.* (-1*unp1[jb-3] + 4*unp1[jb-2] - 5* unp1[jb-1] + 5* 0 - 4*0 + 1*1);


        A[0] = f1[0]-f0[0];
        A[1] = f1[1]-f0[1];
        A[2] = f1[2]-f0[2];
        A[3] = f2[0]-f0[0];
        A[4] = f2[1]-f0[1];
        A[5] = f2[2]-f0[2];
        A[6] = f3[0]-f0[0];
        A[7] = f3[1]-f0[1];
        A[8] = f3[2]-f0[2];

        b[0] = -1*f0[0];
        b[1] = -1*f0[1];
        b[2] = -1*f0[2];

        dgesv_(&n,&nrhs,A,&lda,ipiv,b,&ldb,&info);

        unp1[jb+1] = b[0];
        unp1[jb+2] = b[1];
        unp1[jb+3] = b[2];

        delete[] A;
        delete[] b;
        delete[] ipiv;
        delete[] f0;
        delete[] f1;
        delete[] f2;
        delete[] f3;
      }
    }
    // Dirchlet both sides (assume l_tt = r_tt = 0)
    if(iCase == 2){
      if(oacc == 2){
        unp1[ja-1] = 2*unp1[ja]-unp1[ja+1];
        unp1[jb+1] = 2*unp1[jb]-unp1[jb-1];
      }
      else if(oacc == 4){
        // Left hand Dirchlet DDFA
        double* f0 = new double[2];
        double* f1 = new double[2];
        double* f2 = new double[2];

        f0[0] = pow(def->c,2)/pow(dx,2)*(-2*unp1[ja]+unp1[ja+1] - (6*unp1[ja]-4*unp1[ja+1]+unp1[ja+2])/12);
        f0[1] = pow(def->c,2)/pow(dx,4)*(6*unp1[ja]-4*unp1[ja+1]+unp1[ja+2]);

        f1[0] = pow(def->c,2)/pow(dx,2)*(1-2*unp1[ja]+unp1[ja+1] - (-4+6*unp1[ja]-4*unp1[ja+1]+unp1[ja+2])/12);
        f1[1] = pow(def->c,2)/pow(dx,4)*(-4+6*unp1[ja]-4*unp1[ja+1]+unp1[ja+2]);

        f2[0] = pow(def->c,2)/pow(dx,2)*(-2*unp1[ja]+unp1[ja+1] - (1+6*unp1[ja]-4*unp1[ja+1]+unp1[ja+2])/12);
        f2[1] = pow(def->c,2)/pow(dx,4)*(1+6*unp1[ja]-4*unp1[ja+1]+unp1[ja+2]);

        int n = 2;
        int nrhs = 1;
        double*A = new double[n*n];
        int lda = 2;
        double* b = new double[n];
        int ldb = 2;
        int* ipiv = new int [n];
        int info;

        A[0] = f1[0]-f0[0];
        A[1] = f1[1]-f0[1];
        A[2] = f2[0]-f0[0];
        A[3] = f2[1]-f0[1];

        b[0] = -1*f0[0];
        b[1] = -1*f0[1];

        dgesv_(&n,&nrhs,A,&lda,ipiv,b,&ldb,&info);

        unp1[ja-1] = b[0];
        unp1[ja-2] = b[1];

        // Right hand Dirchlet
      }
      else if(oacc == 6){
        // Left hand Dirchlet

        // Right hand Dirchlet

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
  double dttilde = sigma*dx/def->c;
  int nt         = ceil(tf/dttilde);
  double dt      = tf/nt;

  //  Update sigma
  sigma = def->c * dt/dx;

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
