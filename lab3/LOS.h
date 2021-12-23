#pragma once
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

typedef double real;

class LOS {
   int N;
   int maxiter;
   double eps, nev, nev_next = 0;
   real* ggl, * ggu, * di, * pr, * L, * U, * D;
   int* ig, * jg;
   real* r, * z, * p, * x;
   real** matrixG;
   real alpha, beta;
public:
   void input(ifstream& kuslau, ifstream& ig, ifstream& jg, ifstream& ggl, ifstream& ggu, ifstream& di, ifstream& pr);
   void output(ofstream& output);
   void los(ofstream& iteration);
   void los_d(ofstream& iteration);
   void los_LU(ofstream& iteration);
   void Gilbert(int N);
private:
   void Forward(real* vec, real* res);
   void Backward(real* vec, real* res);
   void LU(real* L, real* U, real* D);
   void r0();
   void xk();
   void rk();
   void zk(real* buf);
   void pk(real* Ak);
   void calcAlpha();
   void calcBeta(real* Ar);
   void copyVec(real* a, real* b, int size);
   void matMul(real* vec, real* res);
   void vecDivD(real* res);
   void prof();
   real norm(real* vec);
   real scalar(real* a, real* b, int size);
};

