/* n-Dimensional Wave Equation Solver
*
*  Currently being built to support 1D then 2D
*
*/

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>

// #include <RAJA/RAJA.hpp>
using namespace std;

// Indexing macro
#define ind2D(r,c,nx) r*nx + c

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
  double operator() (double t) { return 0; }
};

// Right side BC
class R{
public:
  double operator() (double t) { return 0; }
};

// Initial Position
class F2{
public:
  double operator() (double x,double y) { return sin(2*M_PI*x)*sin(2*M_PI*y); }
};

// Initial Velocity
class G2{
public:
  double operator() (double x,double y) { return 0; }
};

// Left side BC
class L2{
public:
  double operator() (double t) { return 0; }
};

// Right side BC
class R2{
public:
  double operator() (double t) { return 0; }
};

class T2{
public:
  double operator() (double t) { return 0; }
};

// Right side BC
class B2{
public:
  double operator() (double t) { return 0; }
};


double y1(double x, double t, F f, int c = 1) {
  return 0.5*(f(x-c*t)+f(x+c*t));
}

double y2(double x, double y, double t, F2 f, int c = sqrt(1/2)){
  return cos(2*M_PI*t)*sin(2*M_PI*x)*sin(2*M_PI*y);
}


struct Definition {
  double a_x,a_y,a_z;
  double b_x,b_y,b_z;
  double c;
  int N_x,N_y,N_z;
  F f;
  G g;
  L l;
  R r;

  F2 f2;
  G2 g2;
  L2 l2;
  R2 r2;
  T2 t2;
  B2 b2;
};

struct Indices{
  int ja_x,ja_y,ja_z;
  int jb_x,jb_y,jb_z;
};

struct Sigma{
  double sigma_x,sigma_y,sigma_z;
};

void firstStep(Definition* def, Sigma* sig, double* x, double* y,double dt, Indices* ind,
        double* unm1, double* un, int nD, int oacc){
  if(nD ==  1){
    double sigma = sig->sigma_x;
    int ja = ind->ja_x;
    int jb = ind->jb_x;

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
  else if(nD == 2){
    int nx = def->N_x+1+oacc;
    for(int r = ind->ja_x; r <= ind->jb_x; r++){
      for(int c = ind->ja_y; c <= ind->jb_y; c++){
        un[ind2D(r,c,nx)] = unm1[ind2D(r,c,nx)] +
                dt*def->g2(x[r],y[c]) +
                pow(sig->sigma_x,2)/2.*(unm1[ind2D(r-1,c,nx)]-2*unm1[ind2D(r,c,nx)]+unm1[ind2D(r+1,c,nx)]) +
                pow(sig->sigma_y,2)/2.*(unm1[ind2D(r,c-1,nx)]-2*unm1[ind2D(r,c,nx)]+unm1[ind2D(r,c+1,nx)]);
      }
    }
  }
}

void timeStep(Definition* def, Sigma* sig, Indices* ind, double* unm1, double* un, double* unp1,
        int nD, int oacc){
  if(nD == 1){
    double sigma = sig->sigma_x;
    int ja = ind->ja_x;
    int jb = ind->jb_x;

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
  else if(nD == 2){
    int nx = def->N_x+1+oacc;
    for(int r = ind->ja_x; r <= ind->jb_x; r++){
      for(int c = ind->ja_y; c <= ind->jb_y; c++){
        unp1[ind2D(r,c,nx)] = 2*un[ind2D(r,c,nx)] - unm1[ind2D(r,c,nx)]
                  + pow(sig->sigma_x,2)*(unm1[ind2D(r-1,c,nx)]-2*unm1[ind2D(r,c,nx)]+unm1[ind2D(r+1,c,nx)])
                  + pow(sig->sigma_y,2)*(unm1[ind2D(r,c-1,nx)]-2*unm1[ind2D(r,c,nx)]+unm1[ind2D(r,c+1,nx)]);
      }
    }
  }
}

void BC(Definition* def, Sigma* sig, double* x, int n, double dt, Indices* ind,
        double* unp1, int nD, int iCase, int oacc){

  if(nD == 1){
    double dx = x[1]-x[0];
    //double sigma = sig->sigma_x;
    int ja = ind->ja_x;
    int jb = ind->jb_x;
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
  else if(nD == 2){
    if(iCase == 1){

    }
    else if(iCase == 2){
      int ja_x = ind->ja_x;
      int ja_y = ind->ja_y;
      int jb_x = ind->jb_x;
      int jb_y = ind->jb_y;
      int nx   = def->N_x+1+oacc;

      if(oacc == 2){
        for(int i = ja_x; i <= jb_x; i++){
          unp1[ind2D(i,ja_y-1,nx)] = 2*unp1[ind2D(i,ja_y,nx)] - unp1[ind2D(i, ja_y+1, nx)];
          unp1[ind2D(i,jb_y+1,nx)] = 2*unp1[ind2D(i,jb_y,nx)] - unp1[ind2D(i, jb_y-1, nx)];
        }
        for(int i = ja_y; i <= jb_y; i++){
          unp1[ind2D(ja_x-1,i,nx)] = 2*unp1[ind2D(ja_x,i,nx)] - unp1[ind2D(ja_x+1,i, nx)];
          unp1[ind2D(jb_x+1,i,nx)] = 2*unp1[ind2D(jb_x,i,nx)] - unp1[ind2D(jb_x-1,i, nx)];
        }
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
  double sigma = stod(argv[1]);
  double tf    = stod(argv[2]);
  int nD       = atoi(argv[3]);
  int icase    = atoi(argv[4]);
  int oacc     = atoi(argv[5]);

  Definition* def = new Definition;
  F f;
  G g;
  L l;
  R r;
  F2 f2;
  G2 g2;
  L2 l2;
  R2 r2;
  T2 t2;
  B2 b2;

  // Problem Initialization
  if(nD == 1){
    def->a_x = 0;
    def->b_x = 1;
    def->c = 1;
    def->N_x = atoi(argv[6]);
    def->f = f;
    def->g = g;
    def->l = l;
    def->r = r;
  }
  else if(nD == 2){
    def->a_x = 0;
    def->b_x = 1;
    def->a_y = 0;
    def->b_y = 1;
    def->c = sqrt(1/2.);
    def->N_x = atoi(argv[6]);
    def->N_y = def->N_x; // same for now
    def->f2 = f2;
    def->g2 = g2;
    def->l2 = l2;
    def->r2 = r2;
    def->t2 = t2;
    def->b2 = b2;
  }


  /////////////////////////////////////////////////////////////////////////////
  // Setup
  /////////////////////////////////////////////////////////////////////////////
  int arrSize,matSize;
  double dt;
  Indices* ind = new Indices;
  double* x;
  double* y;
  // double* z;
  Sigma* sig = new Sigma;
  int nx;

  if(nD == 1){
    //  Set indexing variabls
    ind->ja_x = oacc/2;
    ind->jb_x = def->N_x+ oacc/2;

    //  Steps in space
    double dx = (def->b_x - def->a_x)/(double)def->N_x;
    arrSize = def->N_x+ 1 + oacc;
    x = new double[arrSize];
    for(int i = 0; i < arrSize; i++){
      x[i] = (i - oacc/2)*dx;
    //  fout << x[i] << "\t";
    }
    //fout << endl;

    //  Step in time
    double dttilde = sigma*dx/def->c;
    int nt         = ceil(tf/dttilde);
    dt      = tf/nt;

    nx = arrSize;
    matSize = arrSize;

    //  Update sigma
    sig->sigma_x = def->c * dt/dx;
  }
  else if(nD == 2){
    //  Set indexing variabls
    ind->ja_x = oacc/2;
    ind->jb_x = def->N_x + oacc/2;
    ind->ja_y = oacc/2;
    ind->jb_y = def->N_y + oacc/2;

    // Steps in space
    double dx = (def->b_x - def->a_x)/(double)def->N_x;
    arrSize = def->N_x+ 1 + oacc;
    x = new double[arrSize];
    for(int i = 0; i < arrSize; i++){
      x[i] = (i - oacc/2)*dx;
    }

    double dy = (def->b_y - def->a_y)/(double)def->N_y;
    arrSize = def->N_y+ 1 + oacc;
    y = new double[arrSize];
    for(int i = 0; i < arrSize; i++){
      y[i] = (i - oacc/2)*dy;
    }

    //  Step in time
    double dttilde = sigma*dx*dy/def->c*sqrt(1/dx/dx+1/dy/dy);
    int nt         = ceil(tf/dttilde);
    dt             = tf/(double)nt;

    //  Update sigma
    sig->sigma_x = def->c * dt/dx;
    sig->sigma_y = def->c * dt/dy;

    nx = arrSize;

    matSize = arrSize*arrSize;
  }

  //  Initialize grid - stored in row major order
  double* unm1 = new double[matSize];
  double* un   = new double[matSize];
  double* unp1 = new double[matSize];


  /////////////////////////////////////////////////////////////////////////////
  // Initial Condition
  /////////////////////////////////////////////////////////////////////////////
  if(nD == 1){
    for(int i = 0; i < arrSize; i++){
      unm1[i] = def->f(x[i]);
      // fout << unm1[i] << "\t";
    }
    // fout << endl;
  }
  else if(nD == 2){
    for(int i = 0; i < nx; i++){
      for(int j = 0; j < nx; j++){
        unm1[ind2D(i,j,nx)] = def->f2(x[i],y[j]);
        // fout << unm1[i*nx+j] << "\t";
      }
      //fout << endl;
    }
  }
  cout << "IC" << endl;
  // fout << endl << endl << endl;


  /////////////////////////////////////////////////////////////////////////////
  // First Time Step
  /////////////////////////////////////////////////////////////////////////////

  // firstStep(def, sig, x, y, dt, ind, unm1, un, nD, oacc);
  // BC(def,sig,x,1,dt,ind,un,nD,icase,oacc);

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < nx; j++){
      un[ind2D(i,j,nx)] = y2(x[i],y[j],dt,def->f2);
    }
  }


  cout << "First Step" << endl;
  // for(int i = 0; i < matSize ; i++)
  //    fout<<un[i] << "\t";
  // fout << endl;

  // for(int i = 0; i < nx; i++){
  //   for(int j = 0; j < nx; j++){
  //     fout << un[i*nx+j] << "\t";
  //   }
  //   fout << endl;
  // }
  // fout << endl;

  // double err = 0;
  // if(nD == 1){
  //   for(int i = ind->ja_x; i <= ind->jb_x; i++){
  //     err = max(err,unp1[i]-y1(x[i],0,def->f));
  //   }
  //
  //   // RAJA::ReduceMax<RAJA::seq_reduce, double> err(-1.0);
  //   // RAJA::forall<RAJA::seq_exec>(RAJA::RangeSegment(ja,jb+1),[=] (int i) {
  //   //   double myErr = abs(unp1[i]-y(x[i],tf,f,def->c));
  //   //   err.max(myErr);
  //   // });
  // }
  // else if(nD == 2){
  //   for(int r = ind->ja_x; r <= ind->jb_x; r++){
  //     for(int c = ind->ja_y; c <= ind-> jb_y; c++){
  //       err = max(err, abs( un[ind2D(r,c,nx)] - y2(x[r],y[c],dt,def->f2) ) );
  //     }
  //   }
  // }





  /////////////////////////////////////////////////////////////////////////////
  //  Rest of the time steps
  /////////////////////////////////////////////////////////////////////////////
  int n = 2;
  while(n*dt <= tf){
    // Update middle values
    timeStep(def,sig,ind,unm1,un,unp1,nD,oacc);
    // Update boundary conditions
    BC(def,sig,x,n,dt,ind,unp1,nD,icase,oacc);
    // Update arrays
    for(int i = 0; i<arrSize; i++){
      unm1[i] = un[i];
      un[i] = unp1[i];
    }


    // for(int i = 0; i< arrSize; i++){
    //   fout << un[i] << "\t";
    // }
    // fout << endl;
    cout << n << endl;
    n++;
  }
  cout << "Final time reached" << endl;
  /////////////////////////////////////////////////////////////////////////////
  //  Error
  /////////////////////////////////////////////////////////////////////////////
  double err = 0;
  if(nD == 1){
    for(int i = ind->ja_x; i <= ind->jb_x; i++){
      err = max(err,unp1[i]-y1(x[i],tf,def->f));
    }

    // RAJA::ReduceMax<RAJA::seq_reduce, double> err(-1.0);
    // RAJA::forall<RAJA::seq_exec>(RAJA::RangeSegment(ja,jb+1),[=] (int i) {
    //   double myErr = abs(unp1[i]-y(x[i],tf,f,def->c));
    //   err.max(myErr);
    // });
  }
  else if(nD == 2){
    for(int r = ind->ja_x; r <= ind->jb_x; r++){
      for(int c = ind->ja_y; c <= ind-> jb_y; c++){
        err = max(err, abs( unp1[ind2D(r,c,nx)] - y2(x[r],y[c],tf,def->f2) ) );
      }
    }
  }
  cout << "inf norm err: " << err << endl;

  /////////////////////////////////////////////////////////////////////////////
  //  Output
  /////////////////////////////////////////////////////////////////////////////
  // for(int i = ja; i <= jb; i++){
  //   cout << "x: " << x[i] << "\tu = " << un[i] << endl;
  // }
  //
  fout.close();
  /////////////////////////////////////////////////////////////////////////////
  //  Memory Cleanup
  /////////////////////////////////////////////////////////////////////////////

  delete[] unm1;
  delete[] un;
  delete[] unp1;
  delete[] x;
  delete[] y;

  cout << "Arr Cleaned" << endl;
  // delete[] z;
  delete def;
  delete ind;
  delete sig;
  cout << "Done" << endl;

  return 0;
}
