#include <iostream>
#include <fstream>


using namespace std;

// dgesv_ is a symbol in the LAPACK library files
extern "C" {
extern int dgesv_(int*,int*,double*,int*,int*,double*,int*,int*);
}

int main(int argc,char* argv[]){
  int n = 2;
  int nrhs = 1;
  double* A = new double[n*n];
  int lda = 2;
  double* b = new double[n];
  int ldb = 2;
  int* ipiv = new int[n];
  int info;

  A[0] = 2;
  A[1] = 3;
  A[2] = -5;
  A[3] = 1;

  b[0] = 15;
  b[1] = 31;

  dgesv_(&n,&nrhs,A,&lda,ipiv,b,&ldb,&info);

  cout << b[0] << endl << b[1] << endl;

  return 0;
}
